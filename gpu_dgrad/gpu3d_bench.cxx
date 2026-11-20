// GPU-only benchmark of the tiled 3D gradient paths on an N^3 volume
// (default N = 1024: 1.07 G vertices, 8.6 G cells - larger than the 8 GB
// device, so only the tiled drivers can run it):
//   Robins lower star (ComputeDiscreteGradient3DTiled, float values),
//   DMT PLS / PLS-BR (ComputeLowerStarGradient3DTiled, uint32 ranks),
//   DMT DMF first pass only and with the second pass (ComputeDMFGradient3DTiled).
// The volume and the gradient live in host memory (values and ranks generated
// once; the gradient buffer is reused), copied tile by tile on OpenMP threads.
// No CPU gradient is run.
//
// Usage: gpu3d_bench [N]

#include "gi_basic_types.h"
#include "gi_regular_grid.h"
#include "gi_topological_regular_grid.h"
#include "dgrad_gpu_api.h"
#include "dgrad3d_tables.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace GInt;
namespace gpu = GInt::gpu;

static inline uint32_t Fmix32(uint32_t h) {
    h ^= h >> 16; h *= 0x85ebca6bu; h ^= h >> 13; h *= 0xc2b2ae35u; h ^= h >> 16;
    return h;
}

static inline float ProcValue(int64_t x, int64_t y, int64_t z) {
    const uint32_t h = Fmix32((uint32_t)(x * 73856093ull ^ y * 19349663ull ^ z * 83492791ull));
    const float noise = (float)(h & 0xFF) / 256.0f;
    const float fx = (float)(x % 256) / 256.f, fy = (float)(y % 256) / 256.f, fz = (float)(z % 256) / 256.f;
    return std::sin(6.2831853f * fx * 2.f) * std::cos(6.2831853f * fy * 3.f) * std::sin(6.2831853f * fz + 0.3f) +
           0.05f * noise;
}

struct Row {
    std::string method, tile;
    bool ok;
    gpu::TileStats3D st;
    int64_t crit[4];
    int64_t freeDelta;
};

int main(int argc, char** argv) {
    const int64_t N = argc > 1 ? std::atoll(argv[1]) : 1024;
    gpu::DeviceInfo di;
    if (!gpu::QueryDevice(di)) { printf("no CUDA device\n"); return 2; }
    printf("device: %s, %d SMs, %.1f GB; volume %lld^3\n", di.name, di.sm_count, di.total_mem_bytes / 1e9, (long long)N);

    gpu::Dgrad3DTables tables;
    {
        Explicit3x3x3SmallRegularGrid sgrid;
        if (!gpu::BuildDgrad3DTablesFromSmallMesh<Explicit3x3x3SmallRegularGrid, Vec3l>(sgrid, tables)) return 1;
    }

    const int64_t nv = N * N * N, W = 2 * N - 1, nc = W * W * W;
    auto t0 = std::chrono::steady_clock::now();
    std::vector<float> vals((size_t)nv);
    std::vector<uint32_t> ranks((size_t)nv);
    std::vector<uint8_t> grad((size_t)nc);
#pragma omp parallel for schedule(static)
    for (long long z = 0; z < (long long)N; z++)
        for (int64_t y = 0; y < N; y++)
            for (int64_t x = 0; x < N; x++) {
                const int64_t i = (z * N + y) * N + x;
                vals[(size_t)i] = ProcValue(x, y, z);
                ranks[(size_t)i] = Fmix32((uint32_t)i);   // a bijection: pairwise distinct
            }
#pragma omp parallel for schedule(static)
    for (long long i = 0; i < (long long)nc; i++) grad[(size_t)i] = 0;
    printf("host volume + gradient buffers ready (%.1f GB) in %.0f ms\n",
           (nv * 8.0 + nc) / 1e9, std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());

    // parallel host adapters
    gpu::TileSource3D vsrc = [&](const gpu::Box3D& b, float* dst, int64_t pitch, int64_t spitch) {
#pragma omp parallel for schedule(static)
        for (int z = 0; z < b.d; z++)
            for (int y = 0; y < b.h; y++)
                std::memcpy(dst + z * spitch + y * pitch, &vals[(size_t)(((b.z0 + z) * N + (b.y0 + y)) * N + b.x0)],
                            sizeof(float) * (size_t)b.w);
        return true;
    };
    gpu::RankSource3D rsrc = [&](const gpu::Box3D& b, uint32_t* dst, int64_t pitch, int64_t spitch) {
#pragma omp parallel for schedule(static)
        for (int z = 0; z < b.d; z++)
            for (int y = 0; y < b.h; y++)
                std::memcpy(dst + z * spitch + y * pitch, &ranks[(size_t)(((b.z0 + z) * N + (b.y0 + y)) * N + b.x0)],
                            sizeof(uint32_t) * (size_t)b.w);
        return true;
    };
    gpu::ByteSink3D sink = [&](const gpu::Box3D& b, const uint8_t* s, int64_t pitch, int64_t spitch) {
#pragma omp parallel for schedule(static)
        for (int z = 0; z < b.d; z++)
            for (int y = 0; y < b.h; y++)
                std::memcpy(&grad[(size_t)(((b.z0 + z) * W + (b.y0 + y)) * W + b.x0)], s + z * spitch + y * pitch,
                            (size_t)b.w);
        return true;
    };
    auto census = [&](int64_t crit[4]) {
        long long c0 = 0, c1 = 0, c2 = 0, c3 = 0;
#pragma omp parallel for schedule(static) reduction(+ : c0, c1, c2, c3)
        for (long long z = 0; z < (long long)W; z++)
            for (int64_t y = 0; y < W; y++) {
                const uint8_t* row = &grad[(size_t)((z * W + y) * W)];
                const int pyz = (int)(y & 1) + (int)(z & 1);
                for (int64_t x = 0; x < W; x++)
                    if ((row[x] >> 2 & 7) == 7) {
                        const int d = pyz + (int)(x & 1);
                        if (d == 0) c0++; else if (d == 1) c1++; else if (d == 2) c2++; else c3++;
                    }
            }
        crit[0] = c0; crit[1] = c1; crit[2] = c2; crit[3] = c3;
    };

    // warm-up: the first launch of each kernel loads its module (lazy loading),
    // which would show up as device memory "lost" by the first timed run
    {
        const gpu::TileGrid3D wg = gpu::MakeTileGrid3D(64, 64, 64, 32, 32, 32);
        gpu::TiledOptions3D o;
        gpu::DmtTiledOptions3D d;
        gpu::DmtLowerStarStats ls;
        gpu::DmtDMFStats ds;
        bool w = gpu::ComputeDiscreteGradient3DTiled(vsrc, wg, tables, sink, o, nullptr);
        w = w && gpu::ComputeLowerStarGradient3DTiled(rsrc, wg, gpu::DmtLowerStarRule::PLS, sink, d, &ls, nullptr);
        w = w && gpu::ComputeLowerStarGradient3DTiled(rsrc, wg, gpu::DmtLowerStarRule::PLS_BR, sink, d, &ls, nullptr);
        w = w && gpu::ComputeDMFGradient3DTiled(rsrc, wg, true, sink, nullptr, d, &ds, nullptr);
        if (!w) { printf("warm-up failed\n"); return 1; }
    }

    struct T3 { int x, y, z; };
    const T3 tiles[] = {{128, 128, 64}, {256, 256, 32}, {(int)N, (int)N, 8}};
    std::vector<Row> rows;
    enum M { ROBINS, PLS, PLSBR, DMF1, DMF2 };
    const char* names[] = {"Robins (MSCEER lower star)", "DMT PLS", "DMT PLS-BR", "DMT DMF, alg 1 only", "DMT DMF, alg 1 + 2"};
    for (int m = ROBINS; m <= DMF2; m++) {
        for (const T3& t : tiles) {
            if (m >= DMF1 && t.z == 8) continue;   // DMF slab: local complex too large for 32-bit at this N
            const gpu::TileGrid3D tg = gpu::MakeTileGrid3D(N, N, N, t.x, t.y, t.z);
            Row r;
            r.method = names[m];
            char buf[64];
            std::snprintf(buf, sizeof(buf), "(%d,%d,%d)", t.x, t.y, t.z);
            r.tile = buf;
            int64_t f0 = 0, f1 = 0;
            gpu::QueryDeviceMemory(&f0, nullptr);
            if (m == ROBINS) {
                gpu::TiledOptions3D o;
                o.slots = 3;
                o.tilesPerLaunch = 1;
                r.ok = gpu::ComputeDiscreteGradient3DTiled(vsrc, tg, tables, sink, o, &r.st);
            } else if (m == PLS || m == PLSBR) {
                gpu::DmtTiledOptions3D o;
                gpu::DmtLowerStarStats s;
                r.ok = gpu::ComputeLowerStarGradient3DTiled(rsrc, tg, m == PLS ? gpu::DmtLowerStarRule::PLS : gpu::DmtLowerStarRule::PLS_BR,
                                                            sink, o, &s, &r.st);
                r.ok = r.ok && s.unassigned == 0;
            } else {
                gpu::DmtTiledOptions3D o;
                gpu::DmtDMFStats s;
                r.ok = gpu::ComputeDMFGradient3DTiled(rsrc, tg, m == DMF2, sink, nullptr, o, &s, &r.st);
            }
            gpu::QueryDeviceMemory(&f1, nullptr);
            r.freeDelta = f1 - f0;
            census(r.crit);
            const int64_t euler = r.crit[0] - r.crit[1] + r.crit[2] - r.crit[3];
            printf("%-28s T=%-14s ok=%d wall %8.1f ms | kernel %8.1f h2d %7.1f d2h %7.1f (summed) | source %7.1f sink %7.1f | dev %6.1f MB | euler %lld\n",
                   r.method.c_str(), r.tile.c_str(), (int)r.ok, r.st.wall_ms, r.st.kernel_ms, r.st.h2d_ms, r.st.d2h_ms,
                   r.st.source_ms, r.st.sink_ms, r.st.deviceBytesAllocated / 1e6, (long long)euler);
            r.ok = r.ok && euler == 1 && r.freeDelta == 0;
            rows.push_back(r);
        }
    }

    printf("\n| Method | Tile | Wall ms | G vertices/s | Kernel ms (sum) | H2D GB | D2H GB | Host copies ms (src+sink) | Device MB | Check |\n");
    printf("|---|---|---:|---:|---:|---:|---:|---:|---:|:---:|\n");
    bool all = true;
    for (const Row& r : rows) {
        all = all && r.ok;
        printf("| %s | %s | %.0f | %.2f | %.0f | %.2f | %.2f | %.0f | %.0f | %s |\n", r.method.c_str(), r.tile.c_str(),
               r.st.wall_ms, nv / (r.st.wall_ms * 1e6), r.st.kernel_ms, r.st.bytesH2D / 1e9, r.st.bytesD2H / 1e9,
               r.st.source_ms + r.st.sink_ms, r.st.deviceBytesAllocated / 1e6, r.ok ? "ok" : "FAIL");
    }
    printf("\n(Check: run succeeded, Euler characteristic of the critical cells == 1, device free memory unchanged.)\n");
    return all ? 0 : 1;
}
