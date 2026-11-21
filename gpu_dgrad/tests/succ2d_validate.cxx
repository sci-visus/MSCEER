// Stage A1 gate: the tiled SuccessorResolver (Label2DFlowTerminalsTiled) vs
// the in-core Label2DFlowTerminals.
//
//   1. lattice-index mode: memcmp = 0 vs in-core, both directions, for
//      T in {64, 256, 512} x (slots, tilesPerLaunch) in {(1,1), (3,4)}, on the
//      cls2d_validate case set + 4096^2 + a "staircase" ramp whose V-paths
//      cross dozens of tiles.
//   2. compact-id mode: the partition equals BuildBaseRegions2D's up to a
//      bijection of ids, and terminal count == region count.
//   3. frontier slots per tile <= 4T+4.
//   4. determinism: run twice, memcmp = 0.
//
// Host-only TU (MSVC); includes the gi_* stack for BuildBaseRegions2D.

#include "gi_discrete_gradient_computer.h"
#include "gi_extremum_arc_extract_2d.h"
#include "dgrad_gpu_api.h"
#include "cls/cls2d_tables.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <random>
#include <vector>

using namespace GInt;
typedef Accurate2D::GridType GridType;
typedef Accurate2D::MeshType MeshType;

struct Config { int T, slots, tpl; };
static const Config kConfigs[] = {{64, 1, 1}, {64, 3, 4}, {256, 1, 1}, {256, 3, 4},
                                  {512, 1, 1}, {512, 3, 4}};

static double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// compact partition vs BuildBaseRegions2D: consistent bijection check
static long long CheckCompact(const std::vector<int32_t>& compact, const std::vector<int32_t>& ref,
                              int X, int Y, bool ascending, int64_t terminals) {
    ExtNet2D::BaseRegions2D regions;
    ExtNet2D::BuildBaseRegions2D(ref.data(), X, Y, ascending, regions);
    const int32_t nr = regions.count();
    std::vector<int32_t> c2n((size_t)(terminals > 0 ? terminals : 0), -1);
    std::vector<int32_t> n2c((size_t)nr, -1);
    long long bad = 0;
    if ((int64_t)nr != terminals) bad++;
    for (size_t i = 0; i < compact.size(); i++) {
        const int32_t c = compact[i];
        const int32_t n = regions.latToNode[(size_t)ref[i]];
        if (c < 0 || c >= terminals || n < 0) { bad++; continue; }
        if (c2n[(size_t)c] < 0) c2n[(size_t)c] = n; else if (c2n[(size_t)c] != n) bad++;
        if (n2c[(size_t)n] < 0) n2c[(size_t)n] = c; else if (n2c[(size_t)n] != c) bad++;
    }
    return bad;
}

static long long RunCase(const char* name, int X, int Y,
                         const std::function<float(int, int)>& f,
                         const gpu::Dgrad2DTables& tables, bool timing) {
    printf("\n=== case %s (%dx%d) ===\n", name, X, Y);
    std::vector<float> values((size_t)X * Y);
    for (int gy = 0; gy < Y; gy++)
        for (int gx = 0; gx < X; gx++) values[(size_t)gy * X + gx] = f(gx, gy);
    const size_t n_cells = (size_t)(2 * X - 1) * (size_t)(2 * Y - 1);
    std::vector<uint8_t> grad(n_cells);
    if (!gpu::ComputeCLSGradient2DFused(values.data(), X, Y, tables, grad.data(), NULL)) {
        printf("GPU FAIL: gradient\n");
        return 1;
    }
    const size_t nv = (size_t)X * Y, nq = (size_t)(X - 1) * (Y - 1);
    std::vector<int32_t> refV(nv), refQ(nq);
    gpu::Dgrad2DTimings tm{};
    if (!gpu::Label2DFlowTerminals(grad.data(), X, Y, tables, refV.data(), refQ.data(), &tm) ||
        (timing && !gpu::Label2DFlowTerminals(grad.data(), X, Y, tables, refV.data(), refQ.data(), &tm))) {
        printf("GPU FAIL: in-core labels\n");
        return 1;
    }

    long long bad = 0;
    for (int dir = 0; dir < 2; dir++) {
        const bool asc = dir == 0;
        const std::vector<int32_t>& ref = asc ? refV : refQ;
        const int64_t NX = asc ? X : X - 1;
        for (const Config& c : kConfigs) {
            const gpu::TileGrid2D tg = gpu::MakeTileGrid2D(X, Y, c.T);
            gpu::ResolveOptions2D opt;
            opt.slots = c.slots;
            opt.tilesPerLaunch = c.tpl;
            // lattice mode, twice (determinism)
            std::vector<int32_t> out(ref.size(), -7), out2(ref.size(), -7);
            gpu::ResolveStats2D st, st2;
            const bool ok = gpu::Label2DFlowTerminalsTiled(gpu::HostGradSource2D(grad.data(), X), tg,
                                                           tables, asc, gpu::HostWordSink2D(out.data(), NX),
                                                           opt, &st) &&
                            gpu::Label2DFlowTerminalsTiled(gpu::HostGradSource2D(grad.data(), X), tg,
                                                           tables, asc, gpu::HostWordSink2D(out2.data(), NX),
                                                           opt, &st2);
            long long mm = 0, mm2 = 0;
            size_t first = ref.size();
            if (!ok) mm = (long long)ref.size();
            else
                for (size_t i = 0; i < ref.size(); i++) {
                    if (out[i] != ref[i]) { if (first == ref.size()) first = i; mm++; }
                    if (out[i] != out2[i]) mm2++;
                }
            // compact mode
            opt.compactIds = true;
            std::vector<int32_t> cmp(ref.size(), -7);
            gpu::ResolveStats2D stc;
            const bool okc = gpu::Label2DFlowTerminalsTiled(gpu::HostGradSource2D(grad.data(), X), tg,
                                                            tables, asc, gpu::HostWordSink2D(cmp.data(), NX),
                                                            opt, &stc);
            const long long badc = okc ? CheckCompact(cmp, ref, X, Y, asc, stc.terminals) : (long long)ref.size();
            const bool slotsOk = st.frontierSlotsMaxTile <= 4LL * c.T + 4;
            bad += mm + mm2 + badc + (slotsOk ? 0 : 1);
            printf("  %s T=%3d slots=%d tpl=%d tiles=%5lld: vs in-core %lld, run2 %lld, compact %lld, "
                   "slots max/tot %lld/%lld, passes %d, stitch hops %lld %s",
                   asc ? "vert" : "quad", c.T, c.slots, c.tpl, (long long)tg.numTiles, mm, mm2, badc,
                   (long long)st.frontierSlotsMaxTile, (long long)st.frontierSlots, st.passes,
                   (long long)st.stitchSteps,
                   (mm == 0 && mm2 == 0 && badc == 0 && slotsOk) ? "PASS" : "FAIL");
            if (first != ref.size())
                printf(" (first mismatch lat %zu: tiled=%d ref=%d)", first, out[first], ref[first]);
            printf("\n");
            if (timing && c.slots == 3)
                printf("      phase1 %.2f + stitch %.2f + phase2 %.2f = wall %.2f ms; dev %.1f MB pinned %.1f MB host-stitch %.2f MB\n",
                       st.phase1_ms, st.stitch_ms, st.phase2_ms, st.wall_ms,
                       st.deviceBytesAllocated / 1e6, st.hostPinnedBytes / 1e6, st.hostStitchBytes / 1e6);
        }
    }
    if (timing)
        printf("  in-core Label2DFlowTerminals (both dirs) h2d/kern/d2h %.2f / %.2f / %.2f ms\n",
               tm.h2d_ms, tm.kernel_ms, tm.d2h_ms);
    return bad;
}

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
    gpu::Dgrad2DTables tables;
    {
        GridType grid(Vec2i(8, 8), Vec2b(false, false));
        MeshType mesh(&grid);
        if (!gpu::BuildDgrad2DTablesFromMesh(mesh, tables)) { printf("TABLE FAIL\n"); return 2; }
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
    // staircase: a ramp with a slow diagonal so every V-path runs to the far
    // corner across all tiles (chains of length ~X+Y)
    auto stair = [](int gx, int gy) { return (float)gx + 0.001f * (float)gy; };

    long long total = 0;
    total += RunCase("constant", 64, 64, constant, tables, false);
    total += RunCase("xramp", 128, 96, xramp, tables, false);
    total += RunCase("noise", 256, 256, noise, tables, false);
    total += RunCase("qnoise-odd", 251, 247, qnoise, tables, false);
    total += RunCase("bumps", 512, 384, bumps, tables, false);
    total += RunCase("tiny-3x3", 3, 3, noise, tables, false);
    total += RunCase("tiny-5x3", 5, 3, qnoise, tables, false);
    total += RunCase("odd-1001x777", 1001, 777, noise, tables, false);
    total += RunCase("staircase-2048x512", 2048, 512, stair, tables, false);
    if (!quick) {
        total += RunCase("noise-1024", 1024, 1024, noise, tables, true);
        total += RunCase("qnoise-4096", 4096, 4096, qnoise, tables, true);
    }
    printf("\n=== TOTAL failures: %lld -> %s ===\n", total, total == 0 ? "ALL PASS" : "FAIL");
    return total == 0 ? 0 : 1;
}
