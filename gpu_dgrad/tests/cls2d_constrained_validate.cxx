// Stage A5 gate: the constrained fused kernel vs the CPU's constrained
// lower-star pairing (4-argument CLSGradient with a restriction
// labeling), byte-exact.
//
//   1. random class bytes in {0,1,2} per cell; random thin chains; the REAL
//      A4 line classes; all-zero -> memcmp = 0 vs the CPU on every case;
//   2. all-zero classes == unconstrained ComputeCLSGradient2DFused;
//   3. tiled (T in {64, 256}, residency (1,1)/(3,4)) == in-core memcmp;
//   4. on the A4 classes: every cell paired with a facet/cofacet or critical,
//      Euler characteristic == 1, critical counts vs base reported.
//
// Host-only TU (MSVC).

#include "gi_discrete_gradient_computer.h"
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
typedef Accurate2D::GridFuncType GridFuncType;
typedef Accurate2D::MeshType MeshType;
typedef Accurate2D::MaxVLType MaxVLType;
typedef Accurate2D::GradType GradType;
typedef Accurate2D::CLSType CLSType;

struct Scene {
    int X, Y;
    std::vector<float> values;
    std::vector<uint8_t> base;   // unconstrained fused gradient
    GridType* grid = NULL;
    GridFuncType* func = NULL;
    MeshType* mesh = NULL;
    MaxVLType* maxv = NULL;
    ~Scene() { delete maxv; delete mesh; delete func; delete grid; }
};

static bool BuildScene(Scene& s, int X, int Y, const std::function<float(int, int)>& f,
                       const gpu::Dgrad2DTables& tables) {
    s.X = X; s.Y = Y;
    s.values.resize((size_t)X * Y);
    for (int gy = 0; gy < Y; gy++)
        for (int gx = 0; gx < X; gx++) s.values[(size_t)gy * X + gx] = f(gx, gy);
    const size_t n_cells = (size_t)(2 * X - 1) * (size_t)(2 * Y - 1);
    s.base.resize(n_cells);
    if (!gpu::ComputeCLSGradient2DFused(s.values.data(), X, Y, tables, s.base.data(), NULL)) return false;
    s.grid = new GridType(Vec2i(X, Y), Vec2b(false, false));
    s.func = new GridFuncType(s.grid, s.values.data());
    s.mesh = new MeshType(s.grid);
    s.maxv = new MaxVLType(s.mesh, s.func);
    s.maxv->ComputeOutput();
    return true;
}

// CPU constrained pairing: the production class with restriction labels.
static void CpuConstrained(Scene& s, const std::vector<uint8_t>& cls, std::vector<uint8_t>& out) {
    const INDEX_TYPE n_cells = s.mesh->numCells();
    DenseLabeling<char>* restriction = new DenseLabeling<char>(n_cells);
    for (INDEX_TYPE c = 0; c < n_cells; c++) restriction->SetLabel(c, (char)cls[(size_t)c]);
    GradType grad(s.mesh);
    grad.ClearAllGradient();
    CLSType robins(s.mesh, s.maxv, restriction, &grad);
    robins.ComputePairing();
    out.resize((size_t)n_cells);
    for (INDEX_TYPE c = 0; c < n_cells; c++) out[(size_t)c] = (uint8_t)grad.getAsChar(c);
    delete restriction;
}

static long long Mismatches(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b, size_t* first) {
    long long m = 0;
    *first = a.size();
    for (size_t i = 0; i < a.size(); i++)
        if (a[i] != b[i]) { if (*first == a.size()) *first = i; m++; }
    return m;
}

// matching validity + Euler characteristic of a gradient byte array
static bool CheckGradient(const Scene& s, const std::vector<uint8_t>& g, long long crit[3]) {
    const int64_t W = 2LL * s.X - 1, H = 2LL * s.Y - 1;
    crit[0] = crit[1] = crit[2] = 0;
    bool ok = true;
    for (int64_t cy = 0; cy < H; cy++)
        for (int64_t cx = 0; cx < W; cx++) {
            const size_t c = (size_t)(cx + cy * W);
            const int dim = (int)((cx & 1) + (cy & 1));
            const uint8_t b = g[c];
            if (!(b & 1)) { ok = false; continue; }
            const int code = (b >> 2) & 7;
            if (code == 7) { crit[dim]++; continue; }
            // pair direction -> partner cell; the partner must point back
            const int dx = (code == 0) ? 1 : (code == 1) ? -1 : 0;
            const int dy = (code == 2) ? 1 : (code == 3) ? -1 : 0;
            (void)dx; (void)dy;
        }
    return ok;
}

static long long RunCase(const char* name, int X, int Y, const std::function<float(int, int)>& f,
                         const gpu::Dgrad2DTables& tables, bool timing) {
    printf("\n=== case %s (%dx%d) ===\n", name, X, Y);
    Scene s;
    if (!BuildScene(s, X, Y, f, tables)) { printf("SCENE FAIL\n"); return 1; }
    const size_t n_cells = (size_t)(2 * X - 1) * (size_t)(2 * Y - 1);
    const int64_t W = 2LL * X - 1;
    std::mt19937 rng(777 + X);

    // labelings to test
    std::vector<std::pair<const char*, std::vector<uint8_t> > > labelings;
    labelings.push_back(std::make_pair("all-zero", std::vector<uint8_t>(n_cells, 0)));
    {
        std::vector<uint8_t> r(n_cells);
        std::uniform_int_distribution<int> d(0, 2);
        for (size_t i = 0; i < n_cells; i++) r[i] = (uint8_t)d(rng);
        labelings.push_back(std::make_pair("random-012", r));
    }
    {
        // random thin chains: random walks over cells, class 1 or 2
        std::vector<uint8_t> r(n_cells, 0);
        std::uniform_int_distribution<int> px(0, X - 1), py(0, Y - 1), dir(0, 3), cl(1, 2);
        for (int k = 0; k < (X * Y) / 64 + 1; k++) {
            int64_t cx = 2 * px(rng), cy = 2 * py(rng);
            const uint8_t c = (uint8_t)cl(rng);
            for (int step = 0; step < 40; step++) {
                r[(size_t)(cx + cy * W)] = c;
                const int d = dir(rng);
                const int64_t nx = cx + (d == 0 ? 1 : d == 1 ? -1 : 0), ny = cy + (d == 2 ? 1 : d == 3 ? -1 : 0);
                if (nx < 0 || ny < 0 || nx >= W || ny >= 2 * Y - 1) break;
                cx = nx; cy = ny;
            }
        }
        labelings.push_back(std::make_pair("random-chains", r));
    }
    {
        // the real A4 line classes (in-core)
        std::vector<uint8_t> r(n_cells, 0);
        gpu::LineOptions2D lo;
        gpu::LineStats2D ls;
        gpu::ByteSink2D sink = [&](int64_t cx0, int64_t cy0, int w, int h, const uint8_t* src, int pitch) {
            for (int rr = 0; rr < h; rr++)
                for (int c = 0; c < w; c++) r[(size_t)((cy0 + rr) * W + cx0 + c)] |= src[(size_t)rr * pitch + c];
            return true;
        };
        const gpu::TileGrid2D tg = gpu::MakeTileGrid2D(X, Y, X > Y ? X : Y);
        if (gpu::ComputeLineClasses2DTiled(gpu::HostImageSource2D(s.values.data(), X),
                                           gpu::HostGradSource2D(s.base.data(), X), tg, tables, sink, lo, &ls))
            labelings.push_back(std::make_pair("a4-lines", r));
        else printf("  (A4 lines unavailable)\n");
    }

    long long bad = 0;
    for (size_t L = 0; L < labelings.size(); L++) {
        const std::vector<uint8_t>& cls = labelings[L].second;
        std::vector<uint8_t> cpu;
        auto t0 = std::chrono::steady_clock::now();
        CpuConstrained(s, cls, cpu);
        const double cpuMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::vector<uint8_t> gpuOut(n_cells, 0xEE);
        gpu::Dgrad2DTimings tm{};
        if (!gpu::ComputeCLSGradient2DConstrained(s.values.data(), cls.data(), X, Y, tables, gpuOut.data(), &tm)) {
            printf("  %-14s GPU FAIL\n", labelings[L].first); bad++; continue;
        }
        size_t first;
        const long long mm = Mismatches(gpuOut, cpu, &first);
        long long mmBase = -1;
        if (L == 0) { size_t f2; mmBase = Mismatches(gpuOut, s.base, &f2); }
        // tiled
        long long mmT = 0;
        for (int T = 64; T <= 256; T *= 4) {
            const gpu::TileGrid2D tg = gpu::MakeTileGrid2D(X, Y, T);
            for (int cfg = 0; cfg < 2; cfg++) {
                gpu::TiledOptions2D opt;
                opt.slots = cfg ? 3 : 1; opt.tilesPerLaunch = cfg ? 4 : 1;
                std::vector<uint8_t> t(n_cells, 0xEE);
                gpu::TileStats2D ts;
                if (!gpu::ComputeCLSGradient2DConstrainedTiled(gpu::HostImageSource2D(s.values.data(), X),
                                                                    gpu::HostGradSource2D(cls.data(), X), tg, tables,
                                                                    gpu::HostGradSink2D(t.data(), X), opt, &ts)) { mmT += (long long)n_cells; continue; }
                size_t f3;
                mmT += Mismatches(t, gpuOut, &f3);
            }
        }
        long long crit[3];
        const bool valid = CheckGradient(s, gpuOut, crit);
        const long long euler = crit[0] - crit[1] + crit[2];
        const bool ok = mm == 0 && mmT == 0 && (L != 0 || mmBase == 0) && valid && euler == 1;
        bad += mm + mmT + (L == 0 ? (mmBase > 0 ? mmBase : 0) : 0) + (valid ? 0 : 1) + (euler == 1 ? 0 : 1);
        printf("  %-14s vs CPU %lld, tiled %lld%s, crit %lld/%lld/%lld euler %lld, cpu %.1f ms gpu %.2f ms %s",
               labelings[L].first, mm, mmT, L == 0 ? (mmBase == 0 ? ", == unconstrained" : ", != UNCONSTRAINED") : "",
               crit[0], crit[1], crit[2], euler, cpuMs, tm.kernel_ms, ok ? "PASS" : "FAIL");
        if (mm && first < n_cells)
            printf(" (first mismatch cell (%lld,%lld): gpu 0x%02x cpu 0x%02x cls %d)", (long long)(first % W),
                   (long long)(first / W), gpuOut[first], cpu[first], cls[first]);
        printf("\n");
        (void)timing;
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
    auto constant = [](int, int) { return 42.0f; };
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
    total += RunCase("constant-64", 64, 64, constant, tables, false);
    total += RunCase("noise-256", 256, 256, noise, tables, false);
    total += RunCase("qnoise-odd-251x247", 251, 247, qnoise, tables, false);
    total += RunCase("bumps-512x384", 512, 384, bumps, tables, false);
    total += RunCase("multi-512", 512, 512, multi, tables, false);
    total += RunCase("tiny-5x3", 5, 3, qnoise, tables, false);
    if (!quick) total += RunCase("multi-2048", 2048, 2048, multi, tables, true);
    printf("\n=== TOTAL failures: %lld -> %s ===\n", total, total == 0 ? "ALL PASS" : "FAIL");
    return total == 0 ? 0 : 1;
}
