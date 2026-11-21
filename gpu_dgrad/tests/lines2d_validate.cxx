// Stage A4 gate: tiled numeric separatrix lines (ComputeLineClasses2DTiled).
//
//   1. determinism: run twice, memcmp = 0; same T at residency (1,1) vs (3,4)
//      memcmp = 0 (tile size is an algorithm parameter: hops clamp at the box);
//   2. coverage: every critical edge carries both classes;
//   3. thinness: fraction of line vertices (quads) with all 4 incident edges
//      marked in the same class, measured for the GPU and for the CPU
//      reference alike; gate: gpu <= 1.5 * cpu + 1%;
//   4. agreement with a SERIAL run of the CPU line integrators
//      (ReIntegrateUp/DownFromSaddles replicated without the parallel
//      work-stealing loop, so the reference is deterministic): Jaccard per
//      class, reported.
//
// Host-only TU (MSVC).

#include "gi_discrete_gradient_computer.h"
#include "dgrad_gpu_api.h"
#include "cls/cls2d_tables.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <queue>
#include <random>
#include <set>
#include <vector>

using namespace GInt;
typedef Accurate2D::GridType GridType;
typedef Accurate2D::GridFuncType GridFuncType;
typedef Accurate2D::MeshType MeshType;
typedef Accurate2D::MaxVLType MaxVLType;
typedef Accurate2D::GradType GradType;
typedef Accurate2D::MeshFuncType MeshFuncType;
typedef Accurate2D::StreamlineIntegratorTypeASC IntASC;
typedef Accurate2D::StreamlineIntegratorTypeDSC IntDSC;

static double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

struct Scene {
    int X, Y;
    std::vector<float> values;
    std::vector<uint8_t> grad;
    GridType* grid = NULL;
    GridFuncType* func = NULL;
    MeshType* mesh = NULL;
    MaxVLType* maxv = NULL;
    GradType* gradObj = NULL;
    MeshFuncType* meshfunc = NULL;
    ~Scene() { delete meshfunc; delete gradObj; delete maxv; delete mesh; delete func; delete grid; }
};

static bool BuildScene(Scene& s, int X, int Y, const std::function<float(int, int)>& f,
                       const gpu::Dgrad2DTables& tables) {
    s.X = X; s.Y = Y;
    s.values.resize((size_t)X * Y);
    for (int gy = 0; gy < Y; gy++)
        for (int gx = 0; gx < X; gx++) s.values[(size_t)gy * X + gx] = f(gx, gy);
    const size_t n_cells = (size_t)(2 * X - 1) * (size_t)(2 * Y - 1);
    s.grad.resize(n_cells);
    if (!gpu::ComputeDiscreteGradient2DFused(s.values.data(), X, Y, tables, s.grad.data(), NULL))
        return false;
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
    return true;
}

// Serial replica of ReIntegrateUpFromSaddles (asc) / ReIntegrateDownFromSaddles.
static void CpuLinesSerial(Scene& s, bool asc, std::vector<uint8_t>& out) {
    const INDEX_TYPE n_cells = s.mesh->numCells();
    DenseLabeling<char>* labels = new DenseLabeling<char>(n_cells);
    labels->SetAll(0);
    std::vector<std::pair<float, INDEX_TYPE> > criticals;
    MeshType::DCellsIterator it(s.mesh, 1, 0, n_cells);
    for (it.begin(); it.valid(); it.advance()) {
        const INDEX_TYPE c = it.value();
        if (s.gradObj->getCritical(c))
            criticals.push_back(std::make_pair((asc ? 1.f : -1.f) * (float)s.meshfunc->cellValue(c), c));
    }
    std::sort(criticals.begin(), criticals.end());
    IntASC* ia = asc ? new IntASC(s.grid, s.func, s.mesh, 0.0001f, 0.0f, 500) : NULL;
    IntDSC* id = asc ? NULL : new IntDSC(s.grid, s.func, s.mesh, 0.0001f, 0.0f, 500);
    if (ia) ia->SetDigitizingTarget(labels, 1);
    if (id) id->SetDigitizingTarget(labels, 2);
    for (long long k = (long long)criticals.size() - 1; k >= 0; k--) {
        const INDEX_TYPE sad = criticals[(size_t)k].second;
        std::vector<INDEX_TYPE> result;
        std::queue<INDEX_TYPE> q;
        std::set<INDEX_TYPE> visited;
        q.push(sad);
        int counter = 10;
        while (!q.empty() && counter >= 0) {
            const INDEX_TYPE cur = q.front(); q.pop();
            visited.insert(cur);
            if (asc) {
                MeshType::CofacetsIterator cf(s.mesh);
                for (cf.begin(cur); cf.valid(); cf.advance()) {
                    const INDEX_TYPE t = cf.value();
                    if (s.gradObj->getCritical(t) && visited.count(t) == 0) { result.push_back(t); visited.insert(t); }
                    else if (visited.count(t) == 0) {
                        const INDEX_TYPE p = s.gradObj->getPair(t);
                        result.push_back(t); visited.insert(t); visited.insert(p); q.push(p);
                    }
                }
            } else {
                MeshType::FacetsIterator fc(s.mesh);
                for (fc.begin(cur); fc.valid(); fc.advance()) {
                    const INDEX_TYPE t = fc.value();
                    if (s.gradObj->getCritical(t) && visited.count(t) == 0) { result.push_back(t); visited.insert(t); }
                    else if (visited.count(t) == 0) {
                        const INDEX_TYPE p = s.gradObj->getPair(t);
                        result.push_back(t); visited.insert(t); visited.insert(p); q.push(p);
                    }
                }
            }
            counter--;
        }
        for (size_t r = 0; r < result.size(); r++) {
            const INDEX_TYPE c = result[r];
            if (s.mesh->dimension(c) != (asc ? 2 : 0)) continue;
            std::vector<Vec2d> pts;
            std::vector<INDEX_TYPE> dline;
            Vec2d seed;
            s.mesh->centroid(c, seed);
            seed = seed * 0.5;
            if (asc) ia->IntegrateStreamline(seed, pts, dline);
            else id->IntegrateStreamline(seed, pts, dline);
        }
        if (asc) ia->set_label(sad); else id->set_label(sad);
    }
    out.assign((size_t)n_cells, 0);
    for (INDEX_TYPE c = 0; c < n_cells; c++) out[(size_t)c] = (uint8_t)labels->GetLabel(c);
    delete ia; delete id; delete labels;
}

struct ClassStats { long long cells = 0, thick = 0, inter = 0, uni = 0; };

// bits: class-bit field per cell (GPU style); cpuCls/cpuVal: the reference labeling
static ClassStats Compare(const Scene& s, const std::vector<uint8_t>& bits, uint8_t bit,
                          const std::vector<uint8_t>& cpuCls, uint8_t cpuVal) {
    ClassStats cs;
    const int64_t W = 2LL * s.X - 1, H = 2LL * s.Y - 1;
    for (int64_t cy = 0; cy < H; cy++)
        for (int64_t cx = 0; cx < W; cx++) {
            const size_t c = (size_t)(cx + cy * W);
            const bool g = (bits[c] & bit) != 0, p = cpuCls[c] == cpuVal;
            if (g) cs.cells++;
            if (g && p) cs.inter++;
            if (g || p) cs.uni++;
            // thickness: a line vertex (bit 2) / line quad (bit 1) with all 4 incident edges marked
            const int dim = (int)((cx & 1) + (cy & 1));
            if (g && ((bit == 2 && dim == 0) || (bit == 1 && dim == 2))) {
                int n = 0;
                const int64_t ex[4] = {cx + 1, cx - 1, cx, cx}, ey[4] = {cy, cy, cy + 1, cy - 1};
                for (int k = 0; k < 4; k++) {
                    if (ex[k] < 0 || ey[k] < 0 || ex[k] >= W || ey[k] >= H) continue;
                    if (bits[(size_t)(ex[k] + ey[k] * W)] & bit) n++;
                }
                if (n >= 4) cs.thick++;
            }
        }
    return cs;
}

static long long RunCase(const char* name, int X, int Y, const std::function<float(int, int)>& f,
                         const gpu::Dgrad2DTables& tables, bool cpuRef) {
    printf("\n=== case %s (%dx%d) ===\n", name, X, Y);
    Scene s;
    if (!BuildScene(s, X, Y, f, tables)) { printf("SCENE FAIL\n"); return 1; }
    const size_t n_cells = (size_t)(2 * X - 1) * (size_t)(2 * Y - 1);
    const int64_t W = 2LL * X - 1;

    auto run = [&](int T, int slots, int tpl, std::vector<uint8_t>& out, gpu::LineStats2D& st) {
        const gpu::TileGrid2D tg = gpu::MakeTileGrid2D(X, Y, T);
        gpu::LineOptions2D o;
        o.slots = slots; o.tilesPerLaunch = tpl;
        if (getenv("LINES_BFS")) o.bfsDepth = atoi(getenv("LINES_BFS"));
        out.assign(n_cells, 0);
        gpu::ByteSink2D sink = [&](int64_t cx0, int64_t cy0, int w, int h, const uint8_t* src, int pitch) {
            for (int r = 0; r < h; r++)
                for (int c = 0; c < w; c++) out[(size_t)((cy0 + r) * W + cx0 + c)] |= src[(size_t)r * pitch + c];
            return true;
        };
        return gpu::ComputeLineClasses2DTiled(gpu::HostImageSource2D(s.values.data(), X),
                                              gpu::HostGradSource2D(s.grad.data(), X), tg, tables, sink, o, &st);
    };

    long long bad = 0;
    const int Tin = X > Y ? X : Y;
    std::vector<uint8_t> ref, ref2, t64a, t64a2, t64b, t64b2, t256;
    gpu::LineStats2D st, st2, s64a, s64a2, s64b, s64b2, s256;
    if (!run(Tin, 1, 1, ref, st) || !run(Tin, 1, 1, ref2, st2)) { printf("GPU FAIL in-core\n"); return 1; }
    long long det = 0;
    for (size_t i = 0; i < n_cells; i++) if (ref[i] != ref2[i]) det++;
    printf("  in-core: rounds %lld, resumes %lld, dropped seeds %lld, walk ends terminal/boundary/exit/cap/merged = %lld/%lld/%lld/%lld/%lld, owner keys %.1f MB, %.1f ms; run-twice diff %lld %s\n",
           (long long)st.rounds, (long long)st.resumes, (long long)st.droppedSeeds,
           (long long)st.walkEnds[0], (long long)st.walkEnds[1], (long long)st.walkEnds[2], (long long)st.walkEnds[3],
           (long long)st.walkEnds[4], st.hostOwnerBytes / 1e6, st.wall_ms, det, det == 0 ? "PASS" : "FAIL");
    bad += det + st.droppedSeeds;

    if (!run(64, 1, 1, t64a, s64a) || !run(64, 1, 1, t64a2, s64a2) || !run(64, 3, 4, t64b, s64b) ||
        !run(64, 3, 4, t64b2, s64b2) || !run(256, 3, 4, t256, s256)) { printf("GPU FAIL tiled\n"); return 1; }
    {   // factor test for residency differences
        std::vector<uint8_t> t14, t31; gpu::LineStats2D s14, s31;
        run(64, 1, 4, t14, s14); run(64, 3, 1, t31, s31);
        long long d14 = 0, d31 = 0;
        for (size_t i = 0; i < n_cells; i++) { if (t14[i] != t64a[i]) d14++; if (t31[i] != t64a[i]) d31++; }
        printf("  T=64 factor test vs (1,1): (1,4) %lld, (3,1) %lld\n", d14, d31);
    }
    long long res = 0, rep11 = 0, rep34 = 0, d64 = 0, d256 = 0, shown = 0;
    for (size_t i = 0; i < n_cells; i++) {
        if (t64a[i] != t64a2[i]) rep11++;
        if (t64b[i] != t64b2[i]) rep34++;
        if (t64a[i] != t64b[i]) {
            res++;
            if (shown++ < 4)
                printf("    residency diff at cell (%lld,%lld) dim %d: (1,1)=%d (3,4)=%d\n", (long long)(i % W),
                       (long long)(i / W), (int)((i % W) % 2 + (i / W) % 2), t64a[i], t64b[i]);
        }
        if (t64a[i] != ref[i]) d64++;
        if (t256[i] != ref[i]) d256++;
    }
    printf("  T=64: residency diff %lld, repeat diff (1,1) %lld, (3,4) %lld %s; rounds %lld, resumes %lld, %.1f ms | cross-T diff vs in-core (info): T=64 %lld, T=256 %lld cells\n",
           res, rep11, rep34, (res == 0 && rep11 == 0 && rep34 == 0) ? "PASS" : "FAIL", (long long)s64b.rounds,
           (long long)s64b.resumes, s64b.wall_ms, d64, d256);
    bad += res + rep11 + rep34;

    // coverage: every critical edge carries both classes
    long long saddles = 0, missing = 0;
    for (int64_t cy = 0; cy < 2 * Y - 1; cy++)
        for (int64_t cx = 0; cx < W; cx++) {
            if (((cx ^ cy) & 1) == 0) continue;
            const size_t c = (size_t)(cx + cy * W);
            if (((s.grad[c] >> 2) & 7) != 7) continue;
            saddles++;
            if ((ref[c] & 3) != 3) missing++;
        }
    printf("  coverage: %lld saddles, %lld without both classes %s\n", saddles, missing, missing == 0 ? "PASS" : "FAIL");
    bad += missing;

    if (cpuRef) {
        std::vector<uint8_t> cpuA, cpuD;
        auto t0 = std::chrono::steady_clock::now();
        CpuLinesSerial(s, true, cpuA);
        CpuLinesSerial(s, false, cpuD);
        const double cpuMs = ms_since(t0);
        const ClassStats a = Compare(s, ref, 1, cpuA, 1), d = Compare(s, ref, 2, cpuD, 2);
        // the CPU's own thickness under the same metric (its lines are the product)
        std::vector<uint8_t> cpuAbits(n_cells, 0), cpuDbits(n_cells, 0);
        for (size_t i = 0; i < n_cells; i++) { if (cpuA[i] == 1) cpuAbits[i] = 1; if (cpuD[i] == 2) cpuDbits[i] = 2; }
        const ClassStats ca = Compare(s, cpuAbits, 1, cpuA, 1), cd = Compare(s, cpuDbits, 2, cpuD, 2);
        const double thickA = a.cells ? 100.0 * a.thick / (double)a.cells : 0.0;
        const double thickD = d.cells ? 100.0 * d.thick / (double)d.cells : 0.0;
        const double cthickA = ca.cells ? 100.0 * ca.thick / (double)ca.cells : 0.0;
        const double cthickD = cd.cells ? 100.0 * cd.thick / (double)cd.cells : 0.0;
        const bool thinOk = thickA <= 1.5 * cthickA + 1.0 && thickD <= 1.5 * cthickD + 1.0;
        printf("  asc lines: gpu %lld cells, cpu %lld, Jaccard %.3f, thick gpu %.2f%% cpu %.2f%% | dsc lines: gpu %lld, cpu %lld, Jaccard %.3f, thick gpu %.2f%% cpu %.2f%% (cpu serial %.0f ms) %s\n",
               a.cells, ca.cells, a.uni ? (double)a.inter / a.uni : 0.0, thickA, cthickA,
               d.cells, cd.cells, d.uni ? (double)d.inter / d.uni : 0.0, thickD, cthickD, cpuMs,
               thinOk ? "PASS" : "FAIL");
        if (!thinOk) bad++;
    } else {
        const ClassStats a = Compare(s, ref, 1, ref, 0xFF), d = Compare(s, ref, 2, ref, 0xFF);
        printf("  asc lines: %lld cells, thick %.2f%% | dsc lines: %lld cells, thick %.2f%%\n", a.cells,
               a.cells ? 100.0 * a.thick / (double)a.cells : 0.0, d.cells,
               d.cells ? 100.0 * d.thick / (double)d.cells : 0.0);
    }
    return bad;
}

int main(int argc, char** argv) {
    setvbuf(stdout, NULL, _IONBF, 0);
    const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
    gpu::Dgrad2DTables tables;
    {
        GridType grid(Vec2i(8, 8), Vec2b(false, false));
        MeshType mesh(&grid);
        if (!gpu::BuildDgrad2DTablesFromMesh(mesh, tables)) { printf("TABLE FAIL\n"); return 2; }
    }
    std::mt19937 rng(12345);
    std::uniform_real_distribution<float> uni(0.f, 1.f);
    auto noise = [&](int, int) { return uni(rng); };
    auto bumps = [](int gx, int gy) {
        auto g = [&](float cx, float cy, float s) {
            const float dx = gx - cx, dy = gy - cy;
            return std::exp(-(dx * dx + dy * dy) / (2 * s * s));
        };
        return g(50, 60, 25) + 0.8f * g(300, 150, 60) + 0.5f * g(120, 300, 40);
    };
    auto multi = [&](int gx, int gy) {
        const float x = gx * 0.02f, y = gy * 0.02f;
        return std::sin(x) * std::cos(y * 1.3f) + 0.5f * std::sin(3.1f * x + 0.7f) * std::sin(2.7f * y) +
               0.05f * uni(rng);
    };
    long long total = 0;
    total += RunCase("bumps-512x384", 512, 384, bumps, tables, true);
    total += RunCase("multi-512", 512, 512, multi, tables, true);
    total += RunCase("noise-256", 256, 256, noise, tables, true);
    total += RunCase("tiny-7x5", 7, 5, noise, tables, false);
    if (!quick) total += RunCase("multi-2048", 2048, 2048, multi, tables, false);
    printf("\n=== TOTAL failures: %lld -> %s ===\n", total, total == 0 ? "ALL PASS" : "FAIL");
    return total == 0 ? 0 : 1;
}
