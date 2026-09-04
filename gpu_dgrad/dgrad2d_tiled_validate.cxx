// Stage A0 gate: the tiled 2D fused gradient (tile spine + pinned/streams
// ring) must be byte-identical to the in-core fused path and to the CPU
// reference for every tile size and residency, and must run an image larger
// than device memory at a footprint independent of the image size.
//
// Gates (exit code 0 iff all pass):
//   1. tiled vs in-core ComputeDiscreteGradient2DFused: memcmp = 0 for
//      T in {64, 256, 512, 1000} x (slots, tilesPerLaunch) in {(1,1),(3,1),(3,4)}
//      on the dgrad2d_validate case set + 4096^2.
//   2. tiled vs CPU MyRobinsNoalloc::ComputePairing: memcmp = 0 (same cases).
//   3. procedural 8192^2 and 32768^2 (no host image, no host gradient): every
//      byte assigned, Euler characteristic == 1, device free memory sampled per
//      batch is flat and independent of image size, no device leak.
//
// Host-only TU (MSVC): may include the gi_* stack to build the CPU reference.

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#endif

#include "gi_discrete_gradient_computer.h"
#include "dgrad_gpu_api.h"
#include "dgrad2d_tables.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <random>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace GInt;
typedef Accurate2D::GridType GridType;
typedef Accurate2D::GridFuncType GridFuncType;
typedef Accurate2D::MeshType MeshType;
typedef Accurate2D::MaxVLType MaxVLType;
typedef Accurate2D::GradType GradType;
typedef Accurate2D::RobinsType RobinsType;

static double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0)
        .count();
}

static double host_peak_mb() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return pmc.PeakWorkingSetSize / 1e6;
#endif
    return -1.0;
}
static double host_now_mb() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        return pmc.WorkingSetSize / 1e6;
#endif
    return -1.0;
}

struct Config {
    int T, slots, tpl;
};
static const Config kConfigs[] = {
    {64, 1, 1}, {64, 3, 1}, {64, 3, 4},
    {256, 1, 1}, {256, 3, 1}, {256, 3, 4},
    {512, 1, 1}, {512, 3, 1}, {512, 3, 4},
    {1000, 1, 1}, {1000, 3, 1}, {1000, 3, 4},
};

// Gates 1 + 2 on one image. Returns the number of mismatched bytes summed
// over all configurations (0 = pass).
static long long RunCase(const char* name, int X, int Y,
                         const std::function<float(int, int)>& f,
                         const gpu::Dgrad2DTables& tables, bool timing) {
    printf("\n=== case %s (%dx%d) ===\n", name, X, Y);
    std::vector<float> values((size_t)X * Y);
    for (int gy = 0; gy < Y; gy++)
        for (int gx = 0; gx < X; gx++) values[(size_t)gy * X + gx] = f(gx, gy);
    const size_t n_cells = (size_t)(2 * X - 1) * (size_t)(2 * Y - 1);

    // in-core reference
    std::vector<uint8_t> ref(n_cells, 0);
    gpu::Dgrad2DTimings tm_ref{};
    if (!gpu::ComputeDiscreteGradient2DFused(values.data(), X, Y, tables, ref.data(), &tm_ref) ||
        (timing && !gpu::ComputeDiscreteGradient2DFused(values.data(), X, Y, tables,
                                                        ref.data(), &tm_ref))) {
        printf("GPU FAIL: in-core fused returned false\n");
        return 1;
    }

    // CPU reference
    GridType grid(Vec2i(X, Y), Vec2b(false, false));
    GridFuncType func(&grid, values.data());
    MeshType mesh(&grid);
    MaxVLType maxv(&mesh, &func);
    maxv.ComputeOutput();
    GradType grad(&mesh);
    grad.ClearAllGradient();
    RobinsType robins(&mesh, &maxv, &grad);
    robins.ComputePairing();
    std::vector<uint8_t> cpu(n_cells);
    for (size_t c = 0; c < n_cells; c++) cpu[c] = (uint8_t)grad.getAsChar((INDEX_TYPE)c);

    long long bad = 0;
    const size_t nCfg = sizeof(kConfigs) / sizeof(kConfigs[0]);
    for (size_t k = 0; k < nCfg; k++) {
        const Config& c = kConfigs[k];
        const gpu::TileGrid2D tg = gpu::MakeTileGrid2D(X, Y, c.T);
        gpu::TiledOptions2D opt;
        opt.slots = c.slots;
        opt.tilesPerLaunch = c.tpl;
        std::vector<uint8_t> out(n_cells, 0xEE);  // poison: unwritten bytes show
        gpu::TileStats2D st;
        const bool ok = gpu::ComputeDiscreteGradient2DTiled(
            gpu::HostImageSource2D(values.data(), X), tg, tables,
            gpu::HostGradSink2D(out.data(), X), opt, &st);
        long long mmIn = 0, mmCpu = 0;
        size_t first = n_cells;
        if (!ok) {
            mmIn = (long long)n_cells;
        } else {
            for (size_t i = 0; i < n_cells; i++) {
                if (out[i] != ref[i]) { if (first == n_cells) first = i; mmIn++; }
                if (out[i] != cpu[i]) mmCpu++;
            }
        }
        bad += mmIn + mmCpu;
        printf("  T=%4d slots=%d tpl=%d tiles=%5lld launches=%5lld : vs in-core %lld, vs CPU %lld %s",
               c.T, c.slots, c.tpl, (long long)tg.numTiles, (long long)st.launches, mmIn,
               mmCpu, (mmIn == 0 && mmCpu == 0 && ok) ? "PASS" : "FAIL");
        if (first != n_cells)
            printf("  (first mismatch cell %zu at (%lld,%lld): tiled=0x%02x ref=0x%02x)",
                   first, (long long)(first % (size_t)(2 * X - 1)),
                   (long long)(first / (size_t)(2 * X - 1)), out[first], ref[first]);
        printf("\n");
        if (timing && c.slots == 3) {
            printf("      h2d/kern/d2h sum %.2f / %.2f / %.2f ms, wall %.2f ms, "
                   "H2D %.1f GB/s D2H %.1f GB/s, pinned %.1f MB dev %.1f MB\n",
                   st.h2d_ms, st.kernel_ms, st.d2h_ms, st.wall_ms,
                   st.bytesH2D / (st.h2d_ms * 1e6), st.bytesD2H / (st.d2h_ms * 1e6),
                   st.hostPinnedBytes / 1e6, st.deviceBytesAllocated / 1e6);
        }
    }
    // in-core vs CPU (the pre-existing gate, re-checked since the kernel body changed)
    long long mmRef = 0;
    for (size_t i = 0; i < n_cells; i++) if (ref[i] != cpu[i]) mmRef++;
    bad += mmRef;
    printf("  in-core fused vs CPU : %lld %s\n", mmRef, mmRef == 0 ? "PASS" : "FAIL");
    if (timing)
        printf("  in-core fused h2d/kern/d2h : %.2f / %.2f / %.2f ms\n", tm_ref.h2d_ms,
               tm_ref.kernel_ms, tm_ref.d2h_ms);
    return bad;
}

// Deterministic procedural field: smooth bumps + hashed noise, no storage.
static inline float proc_value(int64_t x, int64_t y) {
    uint32_t h = (uint32_t)(x * 73856093ull ^ y * 19349663ull ^ 0x9E3779B9u);
    h ^= h >> 16; h *= 0x7feb352dU; h ^= h >> 15; h *= 0x846ca68bU; h ^= h >> 16;
    const float noise = (h & 0xFFFFFF) / 16777216.0f;
    const float fx = (float)(x % 1024) / 1024.f, fy = (float)(y % 1024) / 1024.f;
    const float bump = std::sin(6.2831853f * fx * 3.f) * std::cos(6.2831853f * fy * 2.f);
    return bump + 0.25f * noise;
}

struct FootprintResult {
    bool ok = false;
    long long crit[3] = {0, 0, 0};
    long long unassigned = 0;
    long long cells = 0;
    gpu::TileStats2D st;
    int64_t devFreeBefore = 0, devFreeAfter = 0;
    double hostPeakBefore = 0, hostPeakAfter = 0, hostNowMax = 0;
};

static FootprintResult RunProcedural(int64_t N, int T, const gpu::Dgrad2DTables& tables) {
    FootprintResult r;
    printf("\n=== procedural %lldx%lld, T=%d, slots=3, tpl=4 (no host image, no host gradient) ===\n",
           (long long)N, (long long)N, T);
    gpu::QueryDeviceMemory(&r.devFreeBefore, nullptr);
    r.hostPeakBefore = host_peak_mb();

    gpu::TileSource2D src = [](int64_t x0, int64_t y0, int w, int h, float* dst, int pitch) {
#pragma omp parallel for schedule(static)
        for (int rr = 0; rr < h; rr++) {
            float* d = dst + (int64_t)rr * pitch;
            for (int c = 0; c < w; c++) d[c] = proc_value(x0 + c, y0 + rr);
        }
        return true;
    };
    const int64_t W = 2 * N - 1;
    gpu::ByteSink2D sink = [&](int64_t cx0, int64_t cy0, int w, int h, const uint8_t* s, int pitch) {
        long long c0 = 0, c1 = 0, c2 = 0, un = 0;
#pragma omp parallel for schedule(static) reduction(+ : c0, c1, c2, un)
        for (int rr = 0; rr < h; rr++) {
            const uint8_t* row = s + (int64_t)rr * pitch;
            const int py = (int)((cy0 + rr) & 1);
            for (int c = 0; c < w; c++) {
                const uint8_t b = row[c];
                if (!(b & 1)) un++;
                if ((b >> 2 & 7) == 7) {
                    const int d = py + (int)((cx0 + c) & 1);
                    if (d == 0) c0++; else if (d == 1) c1++; else c2++;
                }
            }
        }
        r.crit[0] += c0; r.crit[1] += c1; r.crit[2] += c2; r.unassigned += un;
        r.cells += (long long)w * h;
        const double hn = host_now_mb();
        if (hn > r.hostNowMax) r.hostNowMax = hn;
        (void)W;
        return true;
    };
    const gpu::TileGrid2D tg = gpu::MakeTileGrid2D(N, N, T);
    gpu::TiledOptions2D opt;
    opt.slots = 3;
    opt.tilesPerLaunch = 4;
    r.ok = gpu::ComputeDiscreteGradient2DTiled(src, tg, tables, sink, opt, &r.st);
    gpu::QueryDeviceMemory(&r.devFreeAfter, nullptr);
    r.hostPeakAfter = host_peak_mb();

    const long long euler = r.crit[0] - r.crit[1] + r.crit[2];
    const bool coverage = r.cells == W * W;
    printf("  ok=%d cells=%lld (expected %lld) unassigned=%lld crit=%lld/%lld/%lld euler=%lld\n",
           (int)r.ok, r.cells, (long long)(W * W), r.unassigned, r.crit[0], r.crit[1],
           r.crit[2], euler);
    printf("  tiles=%lld launches=%lld bytes H2D %.2f GB D2H %.2f GB, wall %.1f ms\n",
           (long long)r.st.tiles, (long long)r.st.launches, r.st.bytesH2D / 1e9,
           r.st.bytesD2H / 1e9, r.st.wall_ms);
    printf("  h2d/kern/d2h sums %.1f / %.1f / %.1f ms -> H2D %.1f GB/s, D2H %.1f GB/s\n",
           r.st.h2d_ms, r.st.kernel_ms, r.st.d2h_ms, r.st.bytesH2D / (r.st.h2d_ms * 1e6),
           r.st.bytesD2H / (r.st.d2h_ms * 1e6));
    printf("  device: allocated %.1f MB, free min/max during run %.1f / %.1f MB, "
           "free before/after %.1f / %.1f MB (delta %lld B)\n",
           r.st.deviceBytesAllocated / 1e6, r.st.deviceFreeMin / 1e6, r.st.deviceFreeMax / 1e6,
           r.devFreeBefore / 1e6, r.devFreeAfter / 1e6,
           (long long)(r.devFreeAfter - r.devFreeBefore));
    printf("  host: pinned %.1f MB, working set max %.1f MB, peak before/after %.1f / %.1f MB\n",
           r.st.hostPinnedBytes / 1e6, r.hostNowMax, r.hostPeakBefore, r.hostPeakAfter);
    r.ok = r.ok && coverage && r.unassigned == 0 && euler == 1;
    printf("  -> %s\n", r.ok ? "PASS" : "FAIL");
    return r;
}

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
    gpu::Dgrad2DTables tables;
    {
        GridType grid(Vec2i(8, 8), Vec2b(false, false));
        MeshType mesh(&grid);
        if (!gpu::BuildDgrad2DTablesFromMesh(mesh, tables)) {
            printf("TABLE FAIL\n");
            return 2;
        }
    }

    std::mt19937 rng(12345);
    std::uniform_real_distribution<float> uni(0.f, 1.f);
    std::uniform_int_distribution<int> quant(0, 4);
    auto noise = [&](int, int) { return uni(rng); };
    auto qnoise = [&](int, int) { return (float)quant(rng); };
    auto constant = [](int, int) { return 42.0f; };
    auto xramp = [](int gx, int) { return (float)gx; };
    auto bumps = [](int gx, int gy) {
        auto g = [&](float cx, float cy, float s) {
            const float dx = gx - cx, dy = gy - cy;
            return std::exp(-(dx * dx + dy * dy) / (2 * s * s));
        };
        return g(50, 60, 25) + 0.8f * g(300, 150, 60) + 0.5f * g(120, 300, 40);
    };

    long long total = 0;
    total += RunCase("constant", 64, 64, constant, tables, false);
    total += RunCase("xramp", 128, 96, xramp, tables, false);
    total += RunCase("noise", 256, 256, noise, tables, false);
    total += RunCase("qnoise-odd", 251, 247, qnoise, tables, false);
    total += RunCase("bumps", 512, 384, bumps, tables, false);
    total += RunCase("tiny-3x3", 3, 3, noise, tables, false);
    total += RunCase("tiny-5x3", 5, 3, qnoise, tables, false);
    total += RunCase("odd-1001x777", 1001, 777, noise, tables, false);
    if (!quick) {
        total += RunCase("noise-1024", 1024, 1024, noise, tables, true);
        total += RunCase("qnoise-4096", 4096, 4096, qnoise, tables, true);
    }

    // Gate 3: footprint independence. The 8K^2 image would fit; the 32K^2
    // one (1.07 Gpix: 4.3 GB values + 4.3 GB gradient) does not fit the 8 GB
    // device and is never materialized on the host either.
    bool footprint_ok = true;
    if (!quick) {
        FootprintResult a = RunProcedural(8192, 512, tables);
        FootprintResult b = RunProcedural(32768, 512, tables);
        footprint_ok = a.ok && b.ok;
        const int64_t freeDiff = a.st.deviceFreeMin - b.st.deviceFreeMin;
        const int64_t leakA = a.devFreeAfter - a.devFreeBefore;
        const int64_t leakB = b.devFreeAfter - b.devFreeBefore;
        printf("\n  footprint: device free-min 8K vs 32K differ by %lld B (limit 64 MB); "
               "leak 8K %lld B, 32K %lld B; host peak after 32K %.1f MB\n",
               (long long)freeDiff, (long long)leakA, (long long)leakB, b.hostPeakAfter);
        if (freeDiff > (64ll << 20) || freeDiff < -(64ll << 20)) footprint_ok = false;
        if (leakA != 0 || leakB != 0) footprint_ok = false;
        printf("  footprint gate -> %s\n", footprint_ok ? "PASS" : "FAIL");
    }

    const bool pass = total == 0 && footprint_ok;
    printf("\n=== TOTAL mismatches: %lld, footprint %s -> %s ===\n", total,
           quick ? "skipped" : (footprint_ok ? "ok" : "FAIL"), pass ? "ALL PASS" : "FAIL");
    return pass ? 0 : 1;
}
