// Stage 1 validation harness: 2D GPU discrete gradient vs the CPU reference.
//
// CPU reference = the pure combinatorial path of Accurate2D:
//   RegularGridMaxMinVertexLabeling2D::ComputeOutput()
//   + DiscreteGradientLabeling::ClearAllGradient()
//   + MyRobinsNoalloc<TopologicalRegularGrid2D,...,4,6>::ComputePairing()
// (no numeric re-integration, no setAscendingManifoldDimensions - the compare
// happens on the gradient bytes exactly as ComputePairing leaves them).
//
// The GPU kernel's 9-slot topology tables are generated HERE from the real
// mesh's iterators, so both implementations share one source of truth for
// iteration order; the strict byte compare then binds them.
//
// Note on the CPU baseline: MyRobinsNoalloc::HomotopyExpand contains a known
// dead per-lower-star heap allocation (vector<Vec2l> coords,
// gi_modified_robins.h:1005-1013) that only compiles on 2D meshes. It is part
// of today's shipping 2D path, so the baseline measured here includes it;
// removing it is a separate CPU-side fix.
//
// Exit code: 0 iff every case matches bit-exactly.

#include "gi_discrete_gradient_computer.h"
#include "dgrad_gpu_api.h"

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
typedef Accurate2D::RobinsType RobinsType;

static double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - t0)
        .count();
}

// Generate the 9-slot tables from a real TopologicalRegularGrid2D, using the
// same iterators the CPU algorithm runs. Aborts on any mismatch with the
// slot conventions the kernel assumes.
static bool BuildTablesFromMesh(gpu::Dgrad2DTables& t) {
    GridType grid(Vec2i(8, 8), Vec2b(false, false));
    MeshType mesh(&grid);
    const INDEX_TYPE W = mesh.numCellsAxis(0);
    const INDEX_TYPE base = mesh.coords2Cellid(Vec2l(4, 4)); // interior vertex

    INDEX_TYPE off9[9];
    for (int s = 0; s < 9; s++) {
        off9[s] = mesh.get8NeighborOffset(s);
        const INDEX_TYPE expect = (s % 3 - 1) + (INDEX_TYPE)(s / 3 - 1) * W;
        if (off9[s] != expect) {
            printf("TABLE FAIL: adjacent offset %d is %lld, expected %lld\n", s,
                   (long long)off9[s], (long long)expect);
            return false;
        }
    }
    for (int s = 0; s < 9; s++) {
        if (off9[s] != -off9[8 - s]) { // kernel's maxbyte==8-s membership test
            printf("TABLE FAIL: offset table not symmetric at %d\n", s);
            return false;
        }
    }

    std::memset(&t, 0, sizeof(t));
    for (int s = 0; s < 9; s++) {
        const INDEX_TYPE c = base + off9[s];
        t.dim[s] = (uint8_t)mesh.dimension(c);

        MeshType::FacetsIterator fit(&mesh);
        int nf = 0;
        for (fit.begin(c); fit.valid(); fit.advance()) {
            const INDEX_TYPE d = fit.value() - base;
            for (int q = 0; q < 9; q++)
                if (off9[q] == d) t.fac[s][nf++] = (uint8_t)q;
            // facets outside the 3x3 can never be in a lower star; skipped
        }
        t.fac_count[s] = (uint8_t)nf;

        MeshType::CofacetsIterator cfit(&mesh);
        int nc = 0;
        for (cfit.begin(c); cfit.valid(); cfit.advance()) {
            const INDEX_TYPE d = cfit.value() - base;
            for (int q = 0; q < 9; q++)
                if (off9[q] == d) t.cof[s][nc++] = (uint8_t)q;
        }
        t.cof_count[s] = (uint8_t)nc;
    }
    t.dir_px = mesh.Compress6NeighborOffsetToByte(base, base + 1);
    t.dir_mx = mesh.Compress6NeighborOffsetToByte(base, base - 1);
    t.dir_py = mesh.Compress6NeighborOffsetToByte(base, base + W);
    t.dir_my = mesh.Compress6NeighborOffsetToByte(base, base - W);

    printf("tables: dirs px=%d mx=%d py=%d my=%d\n", t.dir_px, t.dir_mx, t.dir_py,
           t.dir_my);
    for (int s = 0; s < 9; s++) {
        printf("  slot %d dim %d fac[", s, t.dim[s]);
        for (int f = 0; f < t.fac_count[s]; f++) printf(" %d", t.fac[s][f]);
        printf(" ] cof[");
        for (int f = 0; f < t.cof_count[s]; f++) printf(" %d", t.cof[s][f]);
        printf(" ]\n");
    }
    return true;
}

// Runs one image through CPU and GPU; returns the number of mismatched cells.
static long long RunCase(const char* name, int X, int Y,
                         const std::function<float(int, int)>& f,
                         const gpu::Dgrad2DTables& tables, bool report_timing) {
    printf("\n=== case %s (%dx%d) ===\n", name, X, Y);
    std::vector<float> values((size_t)X * Y);
    for (int gy = 0; gy < Y; gy++)
        for (int gx = 0; gx < X; gx++)
            values[(size_t)gy * X + gx] = f(gx, gy);

    GridType* grid = new GridType(Vec2i(X, Y), Vec2b(false, false));
    GridFuncType* func = new GridFuncType(grid, values.data());
    MeshType* mesh = new MeshType(grid);
    const INDEX_TYPE n_cells = mesh->numCells();

    MaxVLType* maxv = new MaxVLType(mesh, func);
    auto t0 = std::chrono::steady_clock::now();
    maxv->ComputeOutput();
    const double label_ms = ms_since(t0);

    GradType* grad = new GradType(mesh);
    grad->ClearAllGradient();
    RobinsType* robins = new RobinsType(mesh, maxv, grad);
    t0 = std::chrono::steady_clock::now();
    robins->ComputePairing();
    const double robins_ms = ms_since(t0);

    // hand the CPU-computed labels to the GPU (fusion removes these in Stage 2)
    std::vector<uint8_t> maxb((size_t)n_cells), minb((size_t)n_cells);
    for (INDEX_TYPE c = 0; c < n_cells; c++) {
        maxb[(size_t)c] = maxv->GetUncompressedMaxVal(c);
        minb[(size_t)c] = maxv->GetUncompressedMinVal(c);
    }

    // both GPU variants: Stage 1 (CPU labels, L2-only) and Stage 2 (fused)
    struct Variant {
        const char* name;
        bool fused;
        gpu::Dgrad2DTimings tm;
        std::vector<uint8_t> out;
    } variants[2] = {{"labels", false, {}, {}}, {"fused ", true, {}, {}}};

    for (auto& v : variants) {
        v.out.assign((size_t)n_cells, 0);
        bool ok;
        for (int run = 0; run < 2; run++) { // 2nd run: warm timings
            std::fill(v.out.begin(), v.out.end(), 0);
            ok = v.fused
                     ? gpu::ComputeDiscreteGradient2DFused(values.data(), X, Y,
                                                           tables, v.out.data(),
                                                           true, &v.tm)
                     : gpu::ComputeDiscreteGradient2D(values.data(), X, Y,
                                                      maxb.data(), minb.data(),
                                                      tables, v.out.data(), true,
                                                      &v.tm);
            if (!ok) break;
        }
        if (!ok) {
            printf("GPU FAIL: %s variant returned false\n", v.name);
            delete robins; delete grad; delete maxv; delete mesh; delete func;
            delete grid;
            return -1;
        }
    }

    long long mismatches = 0, checked = 0, shown = 0;
    for (INDEX_TYPE c = 0; c < n_cells; c++) {
        const INDEX_TYPE hv = maxv->Cell2HighestVertex(c);
        Vec2l hc;
        mesh->cellid2Coords(hv, hc);
        const INDEX_TYPE gx = hc[0] / 2, gy = hc[1] / 2;
        const bool interior = gx >= 1 && gx <= X - 2 && gy >= 1 && gy <= Y - 2;
        const uint8_t expect =
            interior ? (uint8_t)grad->getAsChar(c) : (uint8_t)0;
        if (interior) checked++;
        for (int v = 0; v < 2; v++) {
            if (variants[v].out[(size_t)c] != expect) {
                mismatches++;
                if (shown < 5) {
                    Vec2l cc;
                    mesh->cellid2Coords(c, cc);
                    printf("  MISMATCH [%s] cell %lld (%lld,%lld) dim %d "
                           "interior=%d: cpu=0x%02x gpu=0x%02x\n",
                           variants[v].name, (long long)c, (long long)cc[0],
                           (long long)cc[1], (int)mesh->dimension(c),
                           (int)interior, expect, variants[v].out[(size_t)c]);
                    shown++;
                }
            }
        }
    }

    printf("  interior cells checked : %lld (of %lld total)\n", checked,
           (long long)n_cells);
    printf("  mismatches             : %lld  -> %s\n", mismatches,
           mismatches == 0 ? "PASS" : "FAIL");
    if (report_timing) {
        printf("  CPU maxmin labeling    : %8.2f ms\n", label_ms);
        printf("  CPU robins pairing     : %8.2f ms  (%d omp threads)\n",
               robins_ms, omp_get_max_threads());
        for (const auto& v : variants)
            printf("  GPU %s h2d/kern/d2h : %8.2f / %.2f / %.2f ms\n", v.name,
                   v.tm.h2d_ms, v.tm.kernel_ms, v.tm.d2h_ms);
        // fused-kernel effective traffic: values once per 18x18 tile per 16x16
        // block (1.27x halo amplification) + one gradient byte per cell
        const double fused_bytes =
            (double)X * Y * 4.0 * (18.0 * 18.0) / (16.0 * 16.0) +
            (double)n_cells;
        printf("  fused kernel traffic   : %8.1f MB model -> %.1f GB/s effective\n",
               fused_bytes / 1e6, fused_bytes / (variants[1].tm.kernel_ms * 1e6));
    }

    delete robins; delete grad; delete maxv; delete mesh; delete func; delete grid;
    return mismatches;
}

int main() {
    gpu::Dgrad2DTables tables;
    if (!BuildTablesFromMesh(tables)) return 2;

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
    // correctness set: heavy tie-breaking (const, ramp, quantized) + odd/even dims
    total += RunCase("constant", 64, 64, constant, tables, false);
    total += RunCase("xramp", 128, 96, xramp, tables, false);
    total += RunCase("noise", 256, 256, noise, tables, false);
    total += RunCase("qnoise-odd", 251, 247, qnoise, tables, false);
    total += RunCase("bumps", 512, 384, bumps, tables, false);
    // baseline/timing set
    total += RunCase("noise-1024", 1024, 1024, noise, tables, true);
    total += RunCase("qnoise-4096", 4096, 4096, qnoise, tables, true);

    printf("\n=== TOTAL mismatches: %lld -> %s ===\n", total,
           total == 0 ? "ALL PASS" : "FAIL");
    return total == 0 ? 0 : 1;
}
