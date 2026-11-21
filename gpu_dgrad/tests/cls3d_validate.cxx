// Gate for the 3D GPU lower-star gradient (cls3d_kernel.cu / cls3d_tiled.cu).
//
// Gates (exit code 0 iff all pass):
//   1. tables: generated from Explicit3x3x3SmallRegularGrid iterators and
//      identical to what TopologicalRegularGrid3D's iterators give.
//   2. GPU in-core == CPU CLSSlidingGradient::ComputePairing_sliding,
//      memcmp over the full cell array (domain faces, edges and corners
//      included), on random / constant / ramp / quantized fields, odd and
//      degenerate sizes and 256^3. On the smaller cases the CPU sliding result
//      is also checked against CLSGradient<...,4,6>::ComputePairing.
//   3. tiled == in-core for several tile shapes (cubes, anisotropic, slabs),
//      residencies and tiles per launch.
//   4. constrained: GPU (in-core and tiled) == CPU 4-argument
//      CLSGradient<...,4,6> with restriction labels, including labels whose
//      key class*3+bv collides across classes and is clamped at 23; all-zero
//      labels == the unconstrained output.
//   5. fixed footprint: a procedural volume streamed tile by tile (no host
//      volume, no host gradient): every cell assigned, Euler characteristic 1,
//      an order-independent hash equal to the in-core result (512^3), device
//      free memory flat across volume sizes, 0-byte device leak.
// Then a timing table against the CPU (16 OpenMP threads).
//
// Host-only translation unit (MSVC): it may include the gi_* stack.
// Usage: cls3d_validate [--quick]   (--quick skips 256^3, timings and gate 5)

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#endif

#include "gi_basic_types.h"
#include "gi_timing.h"
#include "gi_regular_grid.h"
#include "gi_topological_regular_grid.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_labeling.h"
#include "gi_discrete_gradient_labeling.h"
#include "gi_max_vertex_labeling.h"
#include "gradient/gi_cls_gradient.h"
#include "gradient/gi_cls_sliding_gradient.h"

#include "dgrad_gpu_api.h"
#include "cls/cls3d_tables.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <random>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace GInt;
typedef RegularGrid3D GridType;
typedef RegularGridTrilinearFunction GridFuncType;
typedef TopologicalRegularGrid3D MeshType;
typedef DiscreteGradientLabeling<MeshType> GradType;
typedef RegularGridMaxMinVertexLabeling3D<MeshType, GridFuncType> MaxVLType;
typedef CLSGradient<MeshType, MaxVLType, GradType, 4, 6> CLSType;   // ondemandaccurate
typedef CLSSlidingGradient<GridType, GridFuncType, MeshType, MaxVLType, GradType> SlidingType;
static const int kKeyClamp = 4 * 6 - 1;   // LMAXDIM * LMAX_LABEL - 1 of CLSType

static double MsSince(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

struct Vol {
    int64_t X, Y, Z;
    std::vector<float> v;
    int64_t NC() const { return (2 * X - 1) * (2 * Y - 1) * (2 * Z - 1); }
};

// ---------------------------------------------------------------------------
// CPU references
// ---------------------------------------------------------------------------

struct CpuTiming {
    double label_ms = 0, pair_ms = 0;
};

// The ondemandaccurate first pass: max/min labeling + ComputePairing_sliding.
static std::vector<uint8_t> CpuSliding(Vol& vol, CpuTiming* tm) {
    GridType grid(Vec3i((int)vol.X, (int)vol.Y, (int)vol.Z), Vec3b(false, false, false));
    GridFuncType func(&grid, vol.v.data());
    MeshType mesh(&grid);
    auto t0 = std::chrono::steady_clock::now();
    MaxVLType maxv(&mesh, &func);
    maxv.ComputeOutput();
    const double lab = MsSince(t0);
    GradType grad(&mesh);
    grad.ClearAllGradient();
    t0 = std::chrono::steady_clock::now();
    {
        SlidingType robins(&grid, &func, &mesh, &maxv, &grad);
        robins.ComputePairing_sliding();
    }
    const double pr = MsSince(t0);
    if (tm) { tm->label_ms = lab; tm->pair_ms = pr; }
    std::vector<uint8_t> out((size_t)vol.NC());
    for (int64_t c = 0; c < vol.NC(); c++) out[(size_t)c] = (uint8_t)grad.getAsChar((INDEX_TYPE)c);
    return out;
}

// CLSGradient<...,4,6>::ComputePairing, optionally constrained.
static std::vector<uint8_t> CpuRobins(Vol& vol, const std::vector<uint8_t>* classes) {
    GridType grid(Vec3i((int)vol.X, (int)vol.Y, (int)vol.Z), Vec3b(false, false, false));
    GridFuncType func(&grid, vol.v.data());
    MeshType mesh(&grid);
    MaxVLType maxv(&mesh, &func);
    maxv.ComputeOutput();
    GradType grad(&mesh);
    grad.ClearAllGradient();
    DenseLabeling<char>* lab = nullptr;
    if (classes) {
        lab = new DenseLabeling<char>(mesh.numCells());
        for (int64_t c = 0; c < vol.NC(); c++) lab->SetLabel((INDEX_TYPE)c, (char)(*classes)[(size_t)c]);
    }
    {
        CLSType robins(&mesh, &maxv, lab, &grad);
        robins.ComputePairing();
    }
    delete lab;
    std::vector<uint8_t> out((size_t)vol.NC());
    for (int64_t c = 0; c < vol.NC(); c++) out[(size_t)c] = (uint8_t)grad.getAsChar((INDEX_TYPE)c);
    return out;
}

// ---------------------------------------------------------------------------
// compare helpers
// ---------------------------------------------------------------------------

static int64_t Diff(const char* what, const Vol& vol, const std::vector<uint8_t>& a,
                    const std::vector<uint8_t>& b, bool print = true) {
    if (a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size()) == 0) return 0;
    int64_t n = 0, first = -1;
    for (size_t i = 0; i < a.size(); i++)
        if (a[i] != b[i]) { if (first < 0) first = (int64_t)i; n++; }
    if (print && first >= 0) {
        const int64_t W = 2 * vol.X - 1, H = 2 * vol.Y - 1;
        printf("    %s: %lld cells differ, first cell (%lld,%lld,%lld) a=0x%02x b=0x%02x\n", what,
               (long long)n, (long long)(first % W), (long long)((first / W) % H),
               (long long)(first / (W * H)), a[(size_t)first], b[(size_t)first]);
    }
    return n == 0 ? 1 : n;
}

static bool TiledRun(const Vol& vol, const std::vector<uint8_t>* classes, int Tx, int Ty, int Tz,
                     int slots, int tpl, int shape, const gpu::Dgrad3DTables& tables,
                     std::vector<uint8_t>& out, gpu::TileStats3D* st) {
    out.assign((size_t)vol.NC(), 0xEE);   // poison: an unwritten byte shows
    const gpu::TileGrid3D tg = gpu::MakeTileGrid3D(vol.X, vol.Y, vol.Z, Tx, Ty, Tz);
    gpu::TiledOptions3D opt;
    opt.slots = slots;
    opt.tilesPerLaunch = tpl;
    opt.block_shape = shape;
    opt.packOutput = ((Tx + Ty + Tz + slots + tpl) & 1) != 0;   // cover both download forms
    const auto src = gpu::HostVolumeSource3D(vol.v.data(), vol.X, vol.Y);
    const auto sink = gpu::HostGradSink3D(out.data(), vol.X, vol.Y);
    if (classes)
        return gpu::ComputeCLSGradient3DConstrainedTiled(
            src, gpu::HostCellSource3D(classes->data(), vol.X, vol.Y), tg, tables, kKeyClamp,
            sink, opt, st);
    return gpu::ComputeCLSGradient3DTiled(src, tg, tables, sink, opt, st);
}

struct TileCfg {
    int Tx, Ty, Tz, slots, tpl, shape;
};

static std::vector<TileCfg> ConfigsFor(const Vol& vol) {
    std::vector<TileCfg> c;
    const int64_t nv = vol.X * vol.Y * vol.Z;
    const int X = (int)vol.X, Y = (int)vol.Y;
    if (nv <= 4096) c.push_back({1, 1, 1, 3, 7, 0});
    if (nv <= 200000) {
        c.push_back({2, 2, 2, 3, 5, 1});
        c.push_back({3, 3, 3, 1, 1, 0});
        c.push_back({5, 3, 2, 3, 4, 2});
        c.push_back({7, 7, 7, 3, 1, 0});
        c.push_back({16, 16, 16, 3, 2, 1});
        c.push_back({X, Y, 1, 3, 3, 0});     // one-slice slabs
        c.push_back({X, 4, 3, 2, 2, 2});
    } else {
        c.push_back({64, 64, 64, 1, 1, 0});
        c.push_back({64, 64, 64, 3, 1, 0});
        c.push_back({64, 64, 64, 3, 4, 1});
        c.push_back({100, 37, 61, 3, 2, 2});
        c.push_back({128, 128, 32, 3, 2, 0});
        c.push_back({X, Y, 16, 3, 1, 0});
    }
    return c;
}

// ---------------------------------------------------------------------------
// gates 2 + 3 (+ 4) on one volume
// ---------------------------------------------------------------------------

static int64_t RunCase(const std::string& name, Vol& vol, const gpu::Dgrad3DTables& tables,
                       bool constrained) {
    printf("\n=== %s (%lldx%lldx%lld) ===\n", name.c_str(), (long long)vol.X, (long long)vol.Y,
           (long long)vol.Z);
    int64_t bad = 0;
    CpuTiming ct;
    const std::vector<uint8_t> cpu = CpuSliding(vol, &ct);

    const int64_t nv = vol.X * vol.Y * vol.Z;
    if (nv <= 300000) {
        const std::vector<uint8_t> rob = CpuRobins(vol, nullptr);
        const int64_t d = Diff("CPU sliding vs CLSGradient", vol, cpu, rob);
        printf("  CPU sliding vs CLSGradient::ComputePairing : %s\n", d ? "DIFFER" : "equal");
        bad += d;
    }

    std::vector<uint8_t> gin((size_t)vol.NC(), 0xEE);
    gpu::Dgrad3DTimings tm{};
    if (!gpu::ComputeCLSGradient3D(vol.v.data(), vol.X, vol.Y, vol.Z, tables, gin.data(), &tm)) {
        printf("  GPU in-core FAILED to run\n");
        return bad + 1;
    }
    const int64_t d0 = Diff("GPU in-core vs CPU", vol, cpu, gin);
    printf("  GPU in-core vs CPU sliding : %s   (CPU label %.1f + pair %.1f ms; GPU wall %.1f, kernel %.2f ms)\n",
           d0 ? "FAIL" : "PASS", ct.label_ms, ct.pair_ms, tm.wall_ms, tm.kernel_ms);
    bad += d0;

    // in-core host modes: direct z-slab DMA (arrays page-locked for the call)
    // and z-slabs through the pinned ring, several slab depths
    {
        const int modes[4] = {2, 1, 2, 1};
        const int depths[4] = {1, 2, 0, 3};
        for (int k = 0; k < 4; k++) {
            gpu::InCoreOptions3D o;
            o.hostMode = modes[k];
            o.slabDepth = depths[k];
            o.packOutput = k == 3;   // packed download, expanded into the caller's array
            std::vector<uint8_t> out((size_t)vol.NC(), 0xEE);
            const bool ok = gpu::ComputeCLSGradient3DEx(vol.v.data(), nullptr, 0, vol.X, vol.Y, vol.Z,
                                                             tables, out.data(), o, nullptr);
            const int64_t d = ok ? Diff("in-core host mode vs in-core", vol, gin, out) : 1;
            printf("  in-core %s slabDepth=%d%s : %s\n", modes[k] == 2 ? "direct" : "copy  ", depths[k],
                   o.packOutput ? " packed" : "", !ok ? "RUN FAILED" : (d ? "FAIL" : "PASS"));
            bad += d;
        }
    }

    for (const TileCfg& c : ConfigsFor(vol)) {
        std::vector<uint8_t> out;
        gpu::TileStats3D st;
        const bool ok = TiledRun(vol, nullptr, c.Tx, c.Ty, c.Tz, c.slots, c.tpl, c.shape, tables, out, &st);
        const int64_t d = ok ? Diff("tiled vs in-core", vol, gin, out) : 1;
        printf("  tiled T=(%d,%d,%d) slots=%d tpl=%d shape=%d tiles=%lld launches=%lld : %s\n", c.Tx,
               c.Ty, c.Tz, c.slots, c.tpl, c.shape, (long long)st.tiles, (long long)st.launches,
               !ok ? "RUN FAILED" : (d ? "FAIL" : "PASS"));
        bad += d;
    }

    if (constrained) {
        // labels: spatial blocks of classes 0..3 plus a sprinkle of 4..7 (keys
        // above the clamp), so subsets split across classes, collide (class 1,
        // bv 0 vs class 0, bv 3) and clamp
        std::vector<uint8_t> cls((size_t)vol.NC());
        std::mt19937 rng(777);
        std::uniform_int_distribution<int> sprinkle(0, 99);
        const int64_t W = 2 * vol.X - 1, H = 2 * vol.Y - 1;
        for (int64_t c = 0; c < vol.NC(); c++) {
            const int64_t x = c % W, y = (c / W) % H, z = c / (W * H);
            int k = (int)(((x / 5) + 2 * (y / 7) + 3 * (z / 3)) % 4);
            const int r = sprinkle(rng);
            if (r < 6) k = 4 + r % 4;
            cls[(size_t)c] = (uint8_t)k;
        }
        const std::vector<uint8_t> ccpu = CpuRobins(vol, &cls);
        std::vector<uint8_t> cg((size_t)vol.NC(), 0xEE);
        bool ok = gpu::ComputeCLSGradient3DConstrained(vol.v.data(), cls.data(), vol.X, vol.Y,
                                                            vol.Z, tables, kKeyClamp, cg.data(), nullptr);
        int64_t d = ok ? Diff("constrained GPU vs CPU", vol, ccpu, cg) : 1;
        printf("  constrained GPU in-core vs CPU CLSGradient(labels) : %s\n", d ? "FAIL" : "PASS");
        bad += d;
        const TileCfg tc = ConfigsFor(vol).back();
        std::vector<uint8_t> ct2;
        ok = TiledRun(vol, &cls, tc.Tx, tc.Ty, tc.Tz, 3, 2, 0, tables, ct2, nullptr);
        d = ok ? Diff("constrained tiled vs in-core", vol, cg, ct2) : 1;
        printf("  constrained tiled T=(%d,%d,%d) vs in-core : %s\n", tc.Tx, tc.Ty, tc.Tz, d ? "FAIL" : "PASS");
        bad += d;
        const TileCfg tc2 = ConfigsFor(vol)[1];
        ok = TiledRun(vol, &cls, tc2.Tx, tc2.Ty, tc2.Tz, 3, 4, 1, tables, ct2, nullptr);
        d = ok ? Diff("constrained tiled vs in-core", vol, cg, ct2) : 1;
        printf("  constrained tiled T=(%d,%d,%d) vs in-core : %s\n", tc2.Tx, tc2.Ty, tc2.Tz, d ? "FAIL" : "PASS");
        bad += d;
        {
            gpu::InCoreOptions3D o;
            o.hostMode = 2;
            o.slabDepth = 2;
            std::vector<uint8_t> out((size_t)vol.NC(), 0xEE);
            ok = gpu::ComputeCLSGradient3DEx(vol.v.data(), cls.data(), kKeyClamp, vol.X, vol.Y, vol.Z,
                                                  tables, out.data(), o, nullptr);
            d = ok ? Diff("constrained direct vs in-core", vol, cg, out) : 1;
            printf("  constrained in-core direct slabDepth=2 vs in-core : %s\n", d ? "FAIL" : "PASS");
            bad += d;
        }
        // all-zero classes reproduce the unconstrained gradient
        std::vector<uint8_t> zero((size_t)vol.NC(), 0);
        ok = gpu::ComputeCLSGradient3DConstrained(vol.v.data(), zero.data(), vol.X, vol.Y, vol.Z,
                                                       tables, kKeyClamp, cg.data(), nullptr);
        d = ok ? Diff("zero-class constrained vs unconstrained", vol, gin, cg) : 1;
        printf("  constrained (all-zero labels) vs unconstrained : %s\n", d ? "FAIL" : "PASS");
        bad += d;
    }
    return bad;
}

// ---------------------------------------------------------------------------
// gate 5: fixed device-memory footprint
// ---------------------------------------------------------------------------

static inline float ProcValue(int64_t x, int64_t y, int64_t z) {
    uint32_t h = (uint32_t)(x * 73856093ull ^ y * 19349663ull ^ z * 83492791ull ^ 0x9E3779B9u);
    h ^= h >> 16; h *= 0x7feb352dU; h ^= h >> 15; h *= 0x846ca68bU; h ^= h >> 16;
    const float noise = (float)(h & 0xFF) / 256.0f;   // coarse: plenty of exact ties
    const float fx = (float)(x % 256) / 256.f, fy = (float)(y % 256) / 256.f,
                fz = (float)(z % 256) / 256.f;
    const float bump = std::sin(6.2831853f * fx * 2.f) * std::cos(6.2831853f * fy * 3.f) *
                       std::sin(6.2831853f * fz * 1.f + 0.3f);
    return bump + 0.05f * noise;
}

static inline uint64_t Mix(uint64_t i, uint8_t b) {
    uint64_t h = (i * 0x9E3779B97F4A7C15ull) ^ ((uint64_t)b * 0xC2B2AE3D27D4EB4Full);
    h ^= h >> 31; h *= 0xBF58476D1CE4E5B9ull; h ^= h >> 29;
    return h;
}

struct ProcResult {
    bool ok = false;
    int64_t crit[4] = {0, 0, 0, 0};
    int64_t unassigned = 0, cells = 0;
    uint64_t hash = 0;
    gpu::TileStats3D st;
    int64_t freeBefore = 0, freeAfter = 0;
};

static ProcResult RunProcedural(int64_t N, int Tx, int Ty, int Tz, int slots, int tpl,
                                const gpu::Dgrad3DTables& tables) {
    ProcResult r;
    printf("\n=== procedural %lld^3, T=(%d,%d,%d), slots=%d, tpl=%d (no host volume, no host gradient) ===\n",
           (long long)N, Tx, Ty, Tz, slots, tpl);
    gpu::QueryDeviceMemory(&r.freeBefore, nullptr);
    gpu::TileSource3D src = [](const gpu::Box3D& b, float* dst, int64_t pitch, int64_t spitch) {
#pragma omp parallel for schedule(static)
        for (int z = 0; z < b.d; z++)
            for (int y = 0; y < b.h; y++) {
                float* d = dst + z * spitch + y * pitch;
                for (int x = 0; x < b.w; x++) d[x] = ProcValue(b.x0 + x, b.y0 + y, b.z0 + z);
            }
        return true;
    };
    const int64_t W = 2 * N - 1;
    gpu::ByteSink3D sink = [&](const gpu::Box3D& b, const uint8_t* s, int64_t pitch, int64_t spitch) {
        long long c0 = 0, c1 = 0, c2 = 0, c3 = 0, un = 0;
        uint64_t hs = 0;
#pragma omp parallel for schedule(static) reduction(+ : c0, c1, c2, c3, un, hs)
        for (int z = 0; z < b.d; z++)
            for (int y = 0; y < b.h; y++) {
                const uint8_t* row = s + z * spitch + y * pitch;
                const int64_t gz = b.z0 + z, gy = b.y0 + y;
                for (int x = 0; x < b.w; x++) {
                    const uint8_t v = row[x];
                    const int64_t gx = b.x0 + x;
                    hs += Mix((uint64_t)((gz * W + gy) * W + gx), v);
                    if (!(v & 1)) un++;
                    if ((v >> 2 & 7) == 7) {
                        const int dim = (int)(gx & 1) + (int)(gy & 1) + (int)(gz & 1);
                        if (dim == 0) c0++; else if (dim == 1) c1++; else if (dim == 2) c2++; else c3++;
                    }
                }
            }
        r.crit[0] += c0; r.crit[1] += c1; r.crit[2] += c2; r.crit[3] += c3;
        r.unassigned += un;
        r.cells += b.Count();
        r.hash += hs;
        return true;
    };
    const gpu::TileGrid3D tg = gpu::MakeTileGrid3D(N, N, N, Tx, Ty, Tz);
    gpu::TiledOptions3D opt;
    opt.slots = slots;
    opt.tilesPerLaunch = tpl;
    r.ok = gpu::ComputeCLSGradient3DTiled(src, tg, tables, sink, opt, &r.st);
    gpu::QueryDeviceMemory(&r.freeAfter, nullptr);
    const int64_t euler = r.crit[0] - r.crit[1] + r.crit[2] - r.crit[3];
    const bool coverage = r.cells == W * W * W;
    printf("  ok=%d cells=%lld (expected %lld) unassigned=%lld crit=%lld/%lld/%lld/%lld euler=%lld\n",
           (int)r.ok, (long long)r.cells, (long long)(W * W * W), (long long)r.unassigned,
           (long long)r.crit[0], (long long)r.crit[1], (long long)r.crit[2], (long long)r.crit[3],
           (long long)euler);
    printf("  tiles=%lld launches=%lld H2D %.2f GB D2H %.2f GB, wall %.1f ms (source %.1f, sink %.1f host ms)\n",
           (long long)r.st.tiles, (long long)r.st.launches, r.st.bytesH2D / 1e9, r.st.bytesD2H / 1e9,
           r.st.wall_ms, r.st.source_ms, r.st.sink_ms);
    printf("  device: allocated %.1f MB, free min/max during run %.1f / %.1f MB, free before/after %.1f / %.1f MB (delta %lld B)\n",
           r.st.deviceBytesAllocated / 1e6, r.st.deviceFreeMin / 1e6, r.st.deviceFreeMax / 1e6,
           r.freeBefore / 1e6, r.freeAfter / 1e6, (long long)(r.freeAfter - r.freeBefore));
    r.ok = r.ok && coverage && r.unassigned == 0 && euler == 1;
    printf("  -> %s\n", r.ok ? "PASS" : "FAIL");
    return r;
}

// ---------------------------------------------------------------------------
// timings
// ---------------------------------------------------------------------------

static void Timings(Vol& vol, const gpu::Dgrad3DTables& tables) {
    printf("\n=== timing %lld^3 random, CPU with 16 OpenMP threads ===\n", (long long)vol.X);
#ifdef _OPENMP
    const int prev = omp_get_max_threads();
    omp_set_num_threads(16);
#endif
    CpuTiming best{1e30, 1e30};
    std::vector<uint8_t> cpu;
    for (int r = 0; r < 3; r++) {
        CpuTiming t;
        cpu = CpuSliding(vol, &t);
        best.label_ms = std::min(best.label_ms, t.label_ms);
        best.pair_ms = std::min(best.pair_ms, t.pair_ms);
    }
#ifdef _OPENMP
    omp_set_num_threads(prev);
#endif
    struct Row {
        std::string label;
        double wall, kern, h2d, d2h;
        bool same;
    };
    std::vector<Row> rows;
    std::vector<uint8_t> out((size_t)vol.NC());
    {
        Row best_row{"GPU in-core (pageable, z-slabs)", 1e30, 0, 0, 0, true};
        for (int r = 0; r < 4; r++) {
            gpu::Dgrad3DTimings tm{};
            gpu::ComputeCLSGradient3D(vol.v.data(), vol.X, vol.Y, vol.Z, tables, out.data(), &tm);
            if (r > 0 && tm.wall_ms < best_row.wall)
                best_row = Row{best_row.label, tm.wall_ms, tm.kernel_ms, tm.h2d_ms, tm.d2h_ms, true};
        }
        best_row.same = out == cpu;
        rows.push_back(best_row);
    }
    const TileCfg cfgs[] = {{64, 64, 64, 3, 4, 0}, {128, 128, 32, 3, 2, 0}, {(int)vol.X, (int)vol.Y, 16, 3, 1, 0},
                            {(int)vol.X, (int)vol.Y, 32, 3, 1, 0}, {(int)vol.X, (int)vol.Y, 32, 4, 1, 1},
                            {(int)vol.X, (int)vol.Y, 32, 4, 1, 2}};
    for (const TileCfg& c : cfgs) {
        char lab[128];
        std::snprintf(lab, sizeof(lab), "GPU tiled T=(%d,%d,%d) slots=%d tpl=%d shape=%d", c.Tx, c.Ty, c.Tz,
                      c.slots, c.tpl, c.shape);
        Row br{lab, 1e30, 0, 0, 0, true};
        for (int r = 0; r < 4; r++) {
            gpu::TileStats3D st;
            TiledRun(vol, nullptr, c.Tx, c.Ty, c.Tz, c.slots, c.tpl, c.shape, tables, out, &st);
            if (r > 0 && st.wall_ms < br.wall) br = Row{lab, st.wall_ms, st.kernel_ms, st.h2d_ms, st.d2h_ms, true};
        }
        br.same = out == cpu;
        rows.push_back(br);
    }
    printf("\n  | path | ms | vs CPU pairing | vs CPU label+pairing | kernel ms (sum) | h2d ms (sum) | d2h ms (sum) | == CPU |\n");
    printf("  |---|---:|---:|---:|---:|---:|---:|:---:|\n");
    printf("  | CPU max/min labeling (16 thr) | %.1f | | | | | | |\n", best.label_ms);
    printf("  | CPU ComputePairing_sliding (16 thr) | %.1f | 1.0x | | | | | |\n", best.pair_ms);
    printf("  | CPU labeling + pairing | %.1f | | 1.0x | | | | |\n", best.label_ms + best.pair_ms);
    for (const Row& r : rows)
        printf("  | %s | %.1f | %.1fx | %.1fx | %.2f | %.2f | %.2f | %s |\n", r.label.c_str(), r.wall,
               best.pair_ms / r.wall, (best.label_ms + best.pair_ms) / r.wall, r.kern, r.h2d, r.d2h,
               r.same ? "yes" : "NO");
}

// ---------------------------------------------------------------------------

typedef std::function<float(int64_t, int64_t, int64_t)> Field;

static Vol MakeVol(int64_t X, int64_t Y, int64_t Z, const Field& f) {
    Vol v{X, Y, Z, std::vector<float>((size_t)(X * Y * Z))};
    for (int64_t z = 0; z < Z; z++)
        for (int64_t y = 0; y < Y; y++)
            for (int64_t x = 0; x < X; x++) v.v[(size_t)((z * Y + y) * X + x)] = f(x, y, z);
    return v;
}

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
    gpu::DeviceInfo di;
    if (!gpu::QueryDevice(di)) {
        printf("no CUDA device\n");
        return 2;
    }
    printf("device: %s, %d SMs, %.1f GB\n", di.name, di.sm_count, di.total_mem_bytes / 1e9);

    // gate 1: tables
    gpu::Dgrad3DTables tables;
    {
        Explicit3x3x3SmallRegularGrid sgrid;
        if (!gpu::BuildDgrad3DTablesFromSmallMesh<Explicit3x3x3SmallRegularGrid, Vec3l>(sgrid, tables)) {
            printf("TABLE FAIL: small-grid iterators break a kernel convention\n");
            return 1;
        }
        GridType g5(Vec3i(5, 5, 5), Vec3b(false, false, false));
        MeshType m5(&g5);
        if (!gpu::CheckDgrad3DTablesAgainstMesh<MeshType, Vec3l>(m5, tables)) {
            printf("TABLE FAIL: TopologicalRegularGrid3D iterators differ from the small grid\n");
            return 1;
        }
        printf("tables: small-grid == TopologicalRegularGrid3D iterator orders: PASS\n");
    }

    std::mt19937 rng(12345);
    std::uniform_real_distribution<float> uni(0.f, 1.f);
    std::uniform_int_distribution<int> q5(0, 4), q2(0, 1);
    const Field noise = [&](int64_t, int64_t, int64_t) { return uni(rng); };
    const Field qnoise = [&](int64_t, int64_t, int64_t) { return (float)q5(rng); };
    const Field binary = [&](int64_t, int64_t, int64_t) { return (float)q2(rng); };
    const Field constant = [](int64_t, int64_t, int64_t) { return 42.0f; };
    const Field xramp = [](int64_t x, int64_t, int64_t) { return (float)x; };
    const Field dramp = [](int64_t x, int64_t y, int64_t z) { return (float)(x + y + z); };
    const Field zdown = [](int64_t, int64_t, int64_t z) { return -(float)z; };
    const Field terrace = [](int64_t x, int64_t y, int64_t z) { return (float)((x + 2 * y + 3 * z) / 4); };
    const Field signedzero = [&](int64_t, int64_t, int64_t) { return q2(rng) ? 0.0f : -0.0f; };
    const Field bumps = [](int64_t x, int64_t y, int64_t z) {
        auto g = [&](float cx, float cy, float cz, float s) {
            const float dx = x - cx, dy = y - cy, dz = z - cz;
            return std::exp(-(dx * dx + dy * dy + dz * dz) / (2 * s * s));
        };
        return g(10, 12, 9, 6) + 0.8f * g(30, 20, 14, 9) + 0.5f * g(20, 33, 5, 7);
    };

    int64_t bad = 0;
    struct C { const char* name; int64_t X, Y, Z; const Field* f; bool constrained; };
    const C cases[] = {
        {"noise", 24, 20, 18, &noise, true},
        {"qnoise-odd", 37, 29, 23, &qnoise, true},
        {"binary", 17, 19, 21, &binary, true},
        {"constant", 16, 16, 16, &constant, true},
        {"xramp", 21, 9, 7, &xramp, false},
        {"diag-ramp", 13, 15, 11, &dramp, true},
        {"z-down", 9, 10, 11, &zdown, false},
        {"terrace", 25, 18, 14, &terrace, true},
        {"signed-zero", 11, 12, 13, &signedzero, false},
        {"bumps", 40, 44, 20, &bumps, true},
        {"tiny-1x1x1", 1, 1, 1, &noise, false},
        {"tiny-2x2x2", 2, 2, 2, &noise, true},
        {"tiny-3x3x3", 3, 3, 3, &qnoise, true},
        {"line-1x1x9", 1, 1, 9, &noise, false},
        {"line-7x1x1", 7, 1, 1, &qnoise, false},
        {"sheet-1x6x5", 1, 6, 5, &noise, true},
        {"sheet-4x1x3", 4, 1, 3, &qnoise, false},
        {"sheet-9x8x1", 9, 8, 1, &noise, true},
        {"thin-2x9x2", 2, 9, 2, &binary, true},
        {"thin-64x3x2", 64, 3, 2, &noise, false},
        {"qnoise-65x66x67", 65, 66, 67, &qnoise, true},
    };
    for (const C& c : cases) {
        Vol v = MakeVol(c.X, c.Y, c.Z, *c.f);
        bad += RunCase(c.name, v, tables, c.constrained);
    }
    if (!quick) {
        Vol big = MakeVol(256, 256, 256, noise);
        bad += RunCase("noise-256^3", big, tables, false);
        Vol bigq = MakeVol(256, 256, 256, qnoise);
        bad += RunCase("qnoise-256^3", bigq, tables, true);
        Vol odd = MakeVol(301, 157, 203, binary);
        bad += RunCase("binary-301x157x203", odd, tables, false);
    }

    bool footprint_ok = true;
    if (!quick) {
        // in-core reference hash at 512^3
        const int64_t N0 = 512;
        uint64_t refHash = 0;
        {
            Vol p{N0, N0, N0, std::vector<float>((size_t)(N0 * N0 * N0))};
#pragma omp parallel for schedule(static)
            for (int64_t z = 0; z < N0; z++)
                for (int64_t y = 0; y < N0; y++)
                    for (int64_t x = 0; x < N0; x++) p.v[(size_t)((z * N0 + y) * N0 + x)] = ProcValue(x, y, z);
            std::vector<uint8_t> g((size_t)p.NC());
            gpu::Dgrad3DTimings tm{};
            if (!gpu::ComputeCLSGradient3D(p.v.data(), N0, N0, N0, tables, g.data(), &tm)) {
                printf("in-core 512^3 FAILED\n");
                footprint_ok = false;
            }
            const int64_t W = 2 * N0 - 1;
            uint64_t hs = 0;
#pragma omp parallel for schedule(static) reduction(+ : hs)
            for (int64_t i = 0; i < p.NC(); i++) hs += Mix((uint64_t)i, g[(size_t)i]);
            refHash = hs;
            (void)W;
            printf("\nin-core 512^3 procedural: wall %.1f ms, kernel %.1f ms\n", tm.wall_ms, tm.kernel_ms);
        }
        ProcResult a = RunProcedural(512, 128, 128, 64, 3, 2, tables);
        ProcResult b = RunProcedural(1024, 128, 128, 64, 3, 2, tables);
        const bool hashOk = a.hash == refHash;
        printf("\n  tiled 512^3 hash == in-core hash: %s\n", hashOk ? "PASS" : "FAIL");
        const int64_t fd1 = a.st.deviceFreeMin - b.st.deviceFreeMin;
        const int64_t la = a.freeAfter - a.freeBefore, lb = b.freeAfter - b.freeBefore;
        printf("  footprint: device free-min 512^3 vs 1024^3 differ by %lld B (limit 64 MB); "
               "leak %lld / %lld B\n",
               (long long)fd1, (long long)la, (long long)lb);
        footprint_ok = footprint_ok && a.ok && b.ok && hashOk;
        if (std::llabs(fd1) > (64ll << 20)) footprint_ok = false;
        if (la != 0 || lb != 0) footprint_ok = false;
        printf("  footprint gate -> %s\n", footprint_ok ? "PASS" : "FAIL");

        Vol tv = MakeVol(256, 256, 256, noise);
        Timings(tv, tables);
    }

    const bool pass = bad == 0 && footprint_ok;
    printf("\n=== TOTAL mismatching cases/cells: %lld, footprint %s -> %s ===\n", (long long)bad,
           quick ? "skipped" : (footprint_ok ? "ok" : "FAIL"), pass ? "ALL PASS" : "FAIL");
    return pass ? 0 : 1;
}
