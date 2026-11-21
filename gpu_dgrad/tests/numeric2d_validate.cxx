// Stage A2 gate: tiled numeric region labeling (LabelNumericRegions2DTiled).
//
//   1. acyclicity + coverage: every raw successor chain ends at a terminal
//      (host walk over the raw words of every tile);
//   2. strict (value, index) decrease along every successor edge - the
//      acyclicity PROOF, checked rather than assumed;
//   3. tiled == in-core (T = image) memcmp, run twice memcmp;
//   4. snap agreement: for a sample of undecided vertices the CPU
//      AdaptiveEulerAdvector2D-based hop (same rules, double precision)
//      picks the same successor;
//   5. region agreement vs CPU-painted numeric basins (host streamline from
//      every vertex to its minimum, single-threaded, no snapping);
//   6. certain fraction + termination histogram.
//
// Host-only TU (MSVC); builds the forest exactly as msc_2d_lib does.

#include "gi_discrete_gradient_computer.h"
#include "gi_extremum_arc_extract_2d.h"
#include "gi_extremum_merge_forest.h"
#include "gi_advection_checker.h"
#include "dgrad_gpu_api.h"
#include "cls/cls2d_tables.h"
#include "common/frontier2d.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <random>
#include <vector>

using namespace GInt;
typedef Accurate2D::GridType GridType;
typedef Accurate2D::GridFuncType GridFuncType;
typedef Accurate2D::MeshType MeshType;
typedef Accurate2D::MaxVLType MaxVLType;
typedef Accurate2D::GradType GradType;
typedef Accurate2D::MeshFuncType MeshFuncType;
typedef AdaptiveEulerAdvector2D<GridFuncType, -1> AdvectorDsc;

static double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

static const char* kEvNames[] = {"certain", "minimum", "cornerCertain", "boxExit", "boundary",
                                 "lowGrad", "cellCap", "iterCap", "guardFallback",
                                 "startBoundary"};

// Everything the tiled call needs, built on the host from the in-core path.
struct Scene {
    int X, Y;
    std::vector<float> values;
    std::vector<uint8_t> grad;
    std::vector<int32_t> termLat;      // vertex -> terminal lattice index
    std::vector<int32_t> baseLabel;    // vertex -> compact base id (node id)
    ExtNet2D::BaseRegions2D regions;
    std::vector<int32_t> rep;          // node -> living node
    std::vector<int32_t> livingId;     // node -> dense living id or -1
    std::vector<float> sStar;          // living -> lowest boundary saddle
    int nLiving = 0;
    // CPU objects for the reference integrator
    GridType* grid = NULL;
    GridFuncType* func = NULL;
    MeshType* mesh = NULL;
    MaxVLType* maxv = NULL;
    GradType* gradObj = NULL;
    MeshFuncType* meshfunc = NULL;
    float persistence = 0.f;
    ~Scene() {
        delete meshfunc; delete gradObj; delete maxv; delete mesh; delete func; delete grid;
    }
};

static bool BuildScene(Scene& s, int X, int Y, const std::function<float(int, int)>& f,
                       float persistenceFrac, const gpu::Dgrad2DTables& tables) {
    s.X = X; s.Y = Y;
    s.values.resize((size_t)X * Y);
    for (int gy = 0; gy < Y; gy++)
        for (int gx = 0; gx < X; gx++) s.values[(size_t)gy * X + gx] = f(gx, gy);
    const size_t n_cells = (size_t)(2 * X - 1) * (size_t)(2 * Y - 1);
    s.grad.resize(n_cells);
    if (!gpu::ComputeCLSGradient2DFused(s.values.data(), X, Y, tables, s.grad.data(), NULL))
        return false;
    s.termLat.resize((size_t)X * Y);
    if (!gpu::Label2DFlowTerminals(s.grad.data(), X, Y, tables, s.termLat.data(), NULL, NULL))
        return false;
    ExtNet2D::BuildBaseRegions2D(s.termLat.data(), X, Y, true, s.regions);
    s.baseLabel.resize((size_t)X * Y);
    for (size_t i = 0; i < s.baseLabel.size(); i++)
        s.baseLabel[i] = s.regions.latToNode[(size_t)s.termLat[i]];

    s.grid = new GridType(Vec2i(X, Y), Vec2b(false, false));
    s.func = new GridFuncType(s.grid, s.values.data());
    s.mesh = new MeshType(s.grid);
    s.maxv = new MaxVLType(s.mesh, s.func);
    s.maxv->ComputeOutput();
    s.gradObj = new GradType(s.mesh);
    std::memcpy(s.gradObj->m_dgrad->LabelArray(), s.grad.data(), n_cells);
    s.meshfunc = new MeshFuncType();
    s.meshfunc->setMeshAndFuncAndMaxVLabeling(s.mesh, s.func, s.maxv);
    s.func->ComputeGradFromImage(1);

    // forest (serial build) at persistence = frac * range
    std::vector<float> nodeVals;
    std::vector<INDEX_TYPE> tieKeys;
    ExtNet2D::BuildNodePayloads2D(s.regions, s.meshfunc, s.maxv, nodeVals, tieKeys);
    ExtremumMergeForest forest;
    forest.SetNodes(s.regions.count(), s.regions.terminalCells.data(), nodeVals.data(),
                    tieKeys.data(), true);
    ExtNet2D::ArcBuffers2D bufs;
    ExtNet2D::ExtractArcsBanded2D(s.gradObj, s.meshfunc, X, Y, 1, s.termLat.data(), NULL,
                                  &s.regions, NULL, &bufs, NULL);
    std::vector<ExtremumMergeForest::Arc> arcsCopy = bufs.local[0];
    forest.Build(1, bufs.local.data(), bufs.cross, bufs.frozen, NULL);
    const float range = s.func->GetMaxValue() - s.func->GetMinValue();
    s.persistence = persistenceFrac * range;
    forest.BuildRemap(s.persistence, s.rep);
    const int32_t nn = s.regions.count();
    s.livingId.assign((size_t)nn, -1);
    for (int32_t i = 0; i < nn; i++) {
        const int32_t r = s.rep[(size_t)i];
        if (s.livingId[(size_t)r] < 0) s.livingId[(size_t)r] = s.nLiving++;
    }
    s.sStar.assign((size_t)s.nLiving, std::numeric_limits<float>::infinity());
    for (size_t i = 0; i < arcsCopy.size(); i++) {
        const ExtremumMergeForest::Arc& a = arcsCopy[i];
        const int32_t la = s.livingId[(size_t)s.rep[(size_t)a.e1]];
        const int32_t lb = s.livingId[(size_t)s.rep[(size_t)a.e2]];
        if (la == lb) continue;
        if (a.saddleValue < s.sStar[(size_t)la]) s.sStar[(size_t)la] = a.saddleValue;
        if (a.saddleValue < s.sStar[(size_t)lb]) s.sStar[(size_t)lb] = a.saddleValue;
    }
    return true;
}

// node -> living id for base label b
static inline int32_t LivingOfBase(const Scene& s, int32_t b) {
    return b < 0 ? -1 : s.livingId[(size_t)s.rep[(size_t)b]];
}

static inline bool vless(const Scene& s, int64_t a, int64_t b) {
    const float fa = s.values[(size_t)a], fb = s.values[(size_t)b];
    if (fa < fb) return true;
    if (fb < fa) return false;
    return a < b;
}

// Gate 1 + 2 on the raw words of a tiled run (T given): walk every chain.
struct RawCheck {
    long long badDecrease = 0, cycles = 0, unresolved = 0, edges = 0;
};
static RawCheck CheckRawWords(const Scene& s, const gpu::TileGrid2D& tg,
                              const std::vector<int32_t>& raw) {
    RawCheck rc;
    const int64_t NX = s.X;
    // successor of global vertex g (as global index), or -1 terminal, -2 bad
    auto succOf = [&](int64_t g) -> int64_t {
        const int64_t gx = g % NX, gy = g / NX;
        const gpu::TileDesc2D d = gpu::TileAt(tg, (gx / tg.T) + (gy / tg.T) * (int64_t)tg.ntx);
        const int32_t w = raw[(size_t)g];
        const uint32_t tag = gpu::WordTag(w), v = gpu::WordValue(w);
        if (tag == gpu::kWordTagTerminal) return -1;
        if (tag == gpu::kWordTagSucc) return (d.ox + v % (uint32_t)d.tw) + (d.oy + v / (uint32_t)d.tw) * NX;
        if (tag == gpu::kWordTagFrontier) {
            int x, y;
            gpu::RingToLocal2D((int)v, d.tw + 2, d.th + 2, x, y);
            const int64_t tx = d.ox + x - 1, ty = d.oy + y - 1;
            if (tx < 0 || ty < 0 || tx >= s.X || ty >= s.Y) return -2;
            return tx + ty * NX;
        }
        return -2;
    };
    std::vector<int8_t> state((size_t)s.X * s.Y, 0);   // 0 new, 1 visiting, 2 done
    std::vector<int64_t> path;
    for (int64_t g0 = 0; g0 < (int64_t)state.size(); g0++) {
        if (state[(size_t)g0]) continue;
        path.clear();
        int64_t g = g0;
        while (true) {
            if (state[(size_t)g] == 2) break;
            if (state[(size_t)g] == 1) { rc.cycles++; break; }
            state[(size_t)g] = 1;
            path.push_back(g);
            const int64_t n = succOf(g);
            if (n == -1) break;
            if (n == -2) { rc.unresolved++; break; }
            rc.edges++;
            if (!vless(s, n, g)) rc.badDecrease++;
            g = n;
        }
        for (size_t i = 0; i < path.size(); i++) state[(size_t)path[i]] = 2;
    }
    return rc;
}

// Host reference basins: integrate from every vertex (no snapping) until
// the nearest vertex is a base minimum, then take its living id.
struct NullChecker2D : public AdvectionChecker2D {
    virtual ADVECTION_EVENT CheckAndDoStuff(const Vec2d&, const Vec2d&, bool) { return NONE; }
};
static void HostBasins(const Scene& s, std::vector<int32_t>& out, int maxIters) {
    const int X = s.X, Y = s.Y;
    out.assign((size_t)X * Y, -1);
    NullChecker2D checker;
    AdvectorDsc adv(s.grid, s.func, 0.0f, 1e-4f, &checker);
    const int64_t W = 2LL * X - 1;
    auto isMin = [&](int gx, int gy) { return ((s.grad[(size_t)(2 * gx + 2LL * gy * W)] >> 2) & 7) == 7; };
    for (int gy = 0; gy < Y; gy++)
        for (int gx = 0; gx < X; gx++) {
            Vec2d p((double)gx, (double)gy);
            int iters = maxIters;
            Vec2l near(gx, gy);
            int guard = 0;
            while (iters > 0 && guard++ < 4000) {
                near = s.grid->Inbounds(p + 0.5);
                if (isMin((int)near[0], (int)near[1])) break;
                ADVECTION_EVENT ev;
                if (s.grid->DistToBoundary(near) <= 1) ev = adv.AdvectThroughVoxelNearBoundary(p, iters);
                else ev = adv.AdvectThroughVoxelNoCheck(p, iters);
                if (ev == LOW_GRADIENT || ev == HIT_EXTREMUM || ev == OVER_MAX_ITERATIONS) break;
            }
            near = s.grid->Inbounds(p + 0.5);
            out[(size_t)gy * X + gx] = LivingOfBase(s, s.baseLabel[(size_t)near[1] * X + near[0]]);
        }
}

static long long RunCase(const char* name, int X, int Y, const std::function<float(int, int)>& f,
                         float persFrac, const gpu::Dgrad2DTables& tables, bool hostRef) {
    printf("\n=== case %s (%dx%d, persistence %.3f of range) ===\n", name, X, Y, persFrac);
    Scene s;
    if (!BuildScene(s, X, Y, f, persFrac, tables)) { printf("SCENE FAIL\n"); return 1; }
    // remap[node] -> living id
    std::vector<int32_t> remap((size_t)s.regions.count());
    for (size_t i = 0; i < remap.size(); i++) remap[i] = LivingOfBase(s, (int32_t)i);
    printf("  base regions %d, living %d\n", s.regions.count(), s.nLiving);

    gpu::NumericOptions2D opt;
    auto run = [&](int T, int slots, int tpl, bool raw, std::vector<int32_t>& out,
                   gpu::NumericStats2D& st) {
        const gpu::TileGrid2D tg = gpu::MakeTileGrid2D(X, Y, T);
        gpu::NumericOptions2D o = opt;
        o.slots = slots; o.tilesPerLaunch = tpl; o.debugRawWords = raw;
        out.assign((size_t)X * Y, -7);
        return gpu::LabelNumericRegions2DTiled(
            gpu::HostImageSource2D(s.values.data(), X), gpu::HostGradSource2D(s.grad.data(), X),
            gpu::HostWordSource2D(s.baseLabel.data(), X), tg, tables, remap.data(),
            (int32_t)remap.size(), s.sStar.data(), s.nLiving, gpu::HostWordSink2D(out.data(), X),
            o, &st);
    };

    long long bad = 0;
    const int Tin = X > Y ? X : Y;   // in-core: one tile
    std::vector<int32_t> ref, out, out2;
    gpu::NumericStats2D st, st2;
    if (!run(Tin, 1, 1, false, ref, st)) { printf("GPU FAIL in-core\n"); return 1; }
    long long ev = 0;
    for (int i = 0; i < 10; i++) ev += st.events[i];
    printf("  events:");
    for (int i = 0; i < 10; i++) printf(" %s=%.2f%%", kEvNames[i], 100.0 * st.events[i] / (double)ev);
    printf("\n  in-core: phase1 %.2f + stitch %.2f + phase2 %.2f = %.2f ms\n", st.resolve.phase1_ms,
           st.resolve.stitch_ms, st.resolve.phase2_ms, st.resolve.wall_ms);
    long long unlabeled = 0;
    for (size_t i = 0; i < ref.size(); i++) if (ref[i] < 0 || ref[i] >= s.nLiving) unlabeled++;
    printf("  coverage: %lld vertices without a living id -> %s\n", unlabeled, unlabeled == 0 ? "PASS" : "FAIL");
    bad += unlabeled;

    // Tile size T is an ALGORITHM parameter (hops are clamped at the tile
    // box), so runs with different T legitimately differ; what must be
    // identical is the same T at different residencies. The cross-T
    // difference vs the in-core run is reported as a statistic.
    const int Ts[3] = {64, 256, 512};
    for (int k = 0; k < 3; k++) {
        const int T = Ts[k];
        if (!run(T, 1, 1, false, out, st) || !run(T, 3, 4, false, out2, st2)) { printf("GPU FAIL T=%d\n", T); bad++; continue; }
        long long mmT = 0, mm2 = 0;
        for (size_t i = 0; i < ref.size(); i++) { if (out[i] != ref[i]) mmT++; if (out[i] != out2[i]) mm2++; }
        // raw words: acyclicity + decrease
        std::vector<int32_t> rawW;
        gpu::NumericStats2D str;
        RawCheck rc;
        if (run(T, 3, 4, true, rawW, str)) rc = CheckRawWords(s, gpu::MakeTileGrid2D(X, Y, T), rawW);
        else rc.unresolved = -1;
        const bool ok = mm2 == 0 && rc.badDecrease == 0 && rc.cycles == 0 && rc.unresolved == 0;
        bad += mm2 + rc.badDecrease + rc.cycles + (rc.unresolved != 0 ? 1 : 0);
        printf("  T=%3d tiles=%5lld: residency (1,1) vs (3,4) %lld | chains: edges %lld, non-decreasing %lld, cycles %lld, unresolved %lld | slots max %lld | vs in-core (info) %lld = %.3f%% | %.2f ms %s\n",
               T, (long long)gpu::MakeTileGrid2D(X, Y, T).numTiles, mm2, rc.edges, rc.badDecrease,
               rc.cycles, rc.unresolved, (long long)st2.resolve.frontierSlotsMaxTile, mmT,
               100.0 * mmT / (double)ref.size(), st2.resolve.wall_ms, ok ? "PASS" : "FAIL");
    }

    // ---- Stage A3: consistency pass on the in-core numeric labels ----
    std::vector<int32_t> cons, cons2, consT;
    {
        auto runC = [&](const std::vector<int32_t>& in, int T, int slots, int tpl,
                        std::vector<int32_t>& o, gpu::ResolveStats2D& rs) {
            const gpu::TileGrid2D tg = gpu::MakeTileGrid2D(X, Y, T);
            gpu::ResolveOptions2D ro;
            ro.slots = slots; ro.tilesPerLaunch = tpl;
            o.assign((size_t)X * Y, -7);
            return gpu::ConsistencyPass2DTiled(gpu::HostImageSource2D(s.values.data(), X),
                                               gpu::HostWordSource2D(in.data(), X), tg,
                                               gpu::HostWordSink2D(o.data(), X), ro, &rs);
        };
        // violations of "non-minimum has a lower same-region axis neighbour"
        auto violations = [&](const std::vector<int32_t>& lab) {
            long long v_bad = 0;
            for (int gy = 0; gy < Y; gy++)
                for (int gx = 0; gx < X; gx++) {
                    const int64_t v = gx + (int64_t)gy * X;
                    const int nx[4] = {gx + 1, gx - 1, gx, gx}, ny[4] = {gy, gy, gy + 1, gy - 1};
                    bool anyLower = false, sameLower = false;
                    for (int k = 0; k < 4; k++) {
                        if (nx[k] < 0 || ny[k] < 0 || nx[k] >= X || ny[k] >= Y) continue;
                        const int64_t u = nx[k] + (int64_t)ny[k] * X;
                        if (!vless(s, u, v)) continue;
                        anyLower = true;
                        if (lab[(size_t)u] == lab[(size_t)v]) sameLower = true;
                    }
                    if (anyLower && !sameLower) v_bad++;
                }
            return v_bad;
        };
        gpu::ResolveStats2D r1, r2, r3;
        const long long before = violations(ref);
        if (!runC(ref, Tin, 1, 1, cons, r1) || !runC(cons, Tin, 1, 1, cons2, r2) ||
            !runC(ref, 64, 3, 4, consT, r3)) {
            printf("  A3 GPU FAIL\n");
            bad++;
        } else {
            const long long after = violations(cons);
            long long relabeled = 0, idem = 0, tiled = 0;
            for (size_t i = 0; i < ref.size(); i++) {
                if (cons[i] != ref[i]) relabeled++;
                if (cons2[i] != cons[i]) idem++;
                if (consT[i] != cons[i]) tiled++;
            }
            const bool ok = after == 0 && idem == 0 && tiled == 0;
            bad += after + idem + tiled;
            printf("  A3 consistency: violations before %lld -> after %lld; relabeled %lld (%.3f%%); idempotence diff %lld; tiled T=64 vs in-core %lld; %.2f ms %s\n",
                   before, after, relabeled, 100.0 * relabeled / (double)ref.size(), idem, tiled,
                   r1.wall_ms, ok ? "PASS" : "FAIL");
        }
    }

    if (hostRef) {
        // host reference basins (single-threaded)
        auto t0 = std::chrono::steady_clock::now();
        std::vector<int32_t> hb;
        HostBasins(s, hb, 4000);
        const double hostMs = ms_since(t0);
        long long agree = 0, agreeCertain = 0, certain = 0;
        for (size_t i = 0; i < ref.size(); i++) {
            if (ref[i] == hb[i]) agree++;
            const int32_t l = LivingOfBase(s, s.baseLabel[i]);
            if (l >= 0 && s.values[i] < s.sStar[(size_t)l]) { certain++; if (ref[i] == hb[i]) agreeCertain++; }
        }
        // discrete (forest) labeling agreement with the host basins, for scale
        long long agreeDisc = 0, agreeCons = 0;
        for (size_t i = 0; i < ref.size(); i++) {
            if (LivingOfBase(s, s.baseLabel[i]) == hb[i]) agreeDisc++;
            if (!cons.empty() && cons[i] == hb[i]) agreeCons++;
        }
        printf("  host basins (%.0f ms): GPU numeric vs host %.3f%% (certain core %.3f%% of vertices, agreement inside it %.3f%%); discrete forest labels vs host %.3f%%; after A3 consistency %.3f%%\n",
               hostMs, 100.0 * agree / ref.size(), 100.0 * certain / ref.size(),
               certain ? 100.0 * agreeCertain / certain : 0.0, 100.0 * agreeDisc / ref.size(),
               100.0 * agreeCons / ref.size());
    }
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
    auto bumps = [](int gx, int gy) {
        auto g = [&](float cx, float cy, float s) {
            const float dx = gx - cx, dy = gy - cy;
            return std::exp(-(dx * dx + dy * dy) / (2 * s * s));
        };
        return -(g(50, 60, 25) + 0.8f * g(300, 150, 60) + 0.5f * g(120, 300, 40));
    };
    // smooth multiscale field with noise: many basins at several scales
    auto multi = [&](int gx, int gy) {
        const float x = gx * 0.02f, y = gy * 0.02f;
        return std::sin(x) * std::cos(y * 1.3f) + 0.5f * std::sin(3.1f * x + 0.7f) * std::sin(2.7f * y) +
               0.05f * uni(rng);
    };
    long long total = 0;
    total += RunCase("bumps-512x384", 512, 384, bumps, 0.05f, tables, true);
    total += RunCase("multi-512", 512, 512, multi, 0.05f, tables, true);
    total += RunCase("qnoise-odd-251x247", 251, 247, qnoise, 0.2f, tables, true);
    total += RunCase("noise-256", 256, 256, noise, 0.1f, tables, true);
    total += RunCase("tiny-5x4", 5, 4, noise, 0.1f, tables, false);
    if (!quick) {
        total += RunCase("multi-2048", 2048, 2048, multi, 0.05f, tables, false);
    }
    printf("\n=== TOTAL failures: %lld -> %s ===\n", total, total == 0 ? "ALL PASS" : "FAIL");
    return total == 0 ? 0 : 1;
}
