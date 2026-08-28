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
#include "gi_morse_smale_complex_basic.h"
#include "dgrad_gpu_api.h"
#include "dgrad2d_tables.h"

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

// Generate the 9-slot tables from a real TopologicalRegularGrid2D via the
// shared builder (the same one msc_2d_lib uses), and print them.
static bool BuildTablesFromMesh(gpu::Dgrad2DTables& t) {
    GridType grid(Vec2i(8, 8), Vec2b(false, false));
    MeshType mesh(&grid);
    if (!gpu::BuildDgrad2DTablesFromMesh(mesh, t)) {
        printf("TABLE FAIL: mesh offsets do not match kernel slot conventions\n");
        return false;
    }

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

    // GPU variants: Stage 1 (labels, interior-only, A/B) and Stage 3 (fused,
    // full domain incl. boundary). The fused variant is compared on EVERY cell.
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
                                                           &v.tm)
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

    long long mismatches = 0, shown = 0;
    long long cpu_crit[3] = {0, 0, 0}, gpu_crit[3] = {0, 0, 0};
    for (INDEX_TYPE c = 0; c < n_cells; c++) {
        const uint8_t cpu_byte = (uint8_t)grad->getAsChar(c);
        const int d = (int)mesh->dimension(c);
        if ((cpu_byte >> 2 & 7) == 7) cpu_crit[d]++;
        if ((variants[1].out[(size_t)c] >> 2 & 7) == 7) gpu_crit[d]++;

        // fused: full-array compare
        if (variants[1].out[(size_t)c] != cpu_byte) {
            mismatches++;
            if (shown < 5) {
                Vec2l cc;
                mesh->cellid2Coords(c, cc);
                printf("  MISMATCH [fused] cell %lld (%lld,%lld) dim %d bv %d: "
                       "cpu=0x%02x gpu=0x%02x\n",
                       (long long)c, (long long)cc[0], (long long)cc[1], d,
                       (int)mesh->boundaryValue(c), cpu_byte,
                       variants[1].out[(size_t)c]);
                shown++;
            }
        }
        // labels variant: interior-owned cells only
        const INDEX_TYPE hv = maxv->Cell2HighestVertex(c);
        Vec2l hc;
        mesh->cellid2Coords(hv, hc);
        const INDEX_TYPE gx = hc[0] / 2, gy = hc[1] / 2;
        const bool interior = gx >= 1 && gx <= X - 2 && gy >= 1 && gy <= Y - 2;
        const uint8_t expect_lab = interior ? cpu_byte : (uint8_t)0;
        if (variants[0].out[(size_t)c] != expect_lab) {
            mismatches++;
            if (shown < 5) {
                printf("  MISMATCH [labels] cell %lld: cpu=0x%02x gpu=0x%02x\n",
                       (long long)c, expect_lab, variants[0].out[(size_t)c]);
                shown++;
            }
        }
    }

    const long long euler = cpu_crit[0] - cpu_crit[1] + cpu_crit[2];
    const bool crit_match = cpu_crit[0] == gpu_crit[0] &&
                            cpu_crit[1] == gpu_crit[1] &&
                            cpu_crit[2] == gpu_crit[2];
    printf("  cells compared (full)  : %lld\n", (long long)n_cells);
    printf("  criticals cpu/gpu      : %lld/%lld min, %lld/%lld sad, %lld/%lld "
           "max; euler=%lld %s\n",
           cpu_crit[0], gpu_crit[0], cpu_crit[1], gpu_crit[1], cpu_crit[2],
           gpu_crit[2], euler,
           (crit_match && euler == 1) ? "(match, chi==1)" : "(CHECK)");
    printf("  mismatches             : %lld  -> %s\n", mismatches,
           mismatches == 0 ? "PASS" : "FAIL");
    if (report_timing) {
        printf("  CPU maxmin labeling    : %8.2f ms\n", label_ms);
        printf("  CPU robins pairing     : %8.2f ms  (%d omp threads)\n",
               robins_ms, omp_get_max_threads());
        for (const auto& v : variants)
            printf("  GPU %s h2d/kern/d2h : %8.2f / %.2f / %.2f ms\n", v.name,
                   v.tm.h2d_ms, v.tm.kernel_ms, v.tm.d2h_ms);
        const double fused_bytes =
            (double)X * Y * 4.0 * (18.0 * 18.0) / (16.0 * 16.0) +
            (double)n_cells;
        printf("  fused kernel traffic   : %8.1f MB model -> %.1f GB/s effective\n",
               fused_bytes / 1e6, fused_bytes / (variants[1].tm.kernel_ms * 1e6));
        // block-shape tuning matrix on this case
        for (int shape = 0; shape < 3; shape++) {
            static const char* shape_names[3] = {"16x16", "32x8 ", "8x8  "};
            gpu::Dgrad2DTimings ts{};
            std::vector<uint8_t> tout((size_t)n_cells);
            gpu::ComputeDiscreteGradient2DBatched(values.data(), X, Y, 1, tables,
                                                  tout.data(), shape, &ts);
            gpu::ComputeDiscreteGradient2DBatched(values.data(), X, Y, 1, tables,
                                                  tout.data(), shape, &ts);
            const bool same = std::memcmp(tout.data(), variants[1].out.data(),
                                          (size_t)n_cells) == 0;
            printf("  block %s kernel     : %8.2f ms %s\n", shape_names[shape],
                   ts.kernel_ms, same ? "" : "(OUTPUT DIFFERS!)");
        }
    }

    delete robins; delete grad; delete maxv; delete mesh; delete func; delete grid;
    return mismatches;
}

// Stage 4 labeling gate: flow-terminal labeling by pointer doubling over the
// offset-doubled successor (vertex flow: succ(v)=v+2*pairdir; quad flow:
// succ(q)=q+2*pairdir) must reproduce the per-cell base manifold identity the
// serial fillGeometry recursion paints - rec_man_trace_up/down is exactly the
// reverse of the successor, so the partitions coincide cell for cell.
typedef GInt::MorseSmaleComplexBasic<float, MeshType, Accurate2D::MeshFuncType,
                                     GradType>
    MscT;

static long long RunLabelCase(const char* name, int X, int Y,
                              const std::function<float(int, int)>& f,
                              const gpu::Dgrad2DTables& tables) {
    printf("\n=== labeling case %s (%dx%d) ===\n", name, X, Y);
    std::vector<float> values((size_t)X * Y);
    for (int gy = 0; gy < Y; gy++)
        for (int gx = 0; gx < X; gx++)
            values[(size_t)gy * X + gx] = f(gx, gy);

    GridType* grid = new GridType(Vec2i(X, Y), Vec2b(false, false));
    GridFuncType* func = new GridFuncType(grid, values.data());
    MeshType* mesh = new MeshType(grid);
    const INDEX_TYPE n_cells = mesh->numCells();

    MaxVLType* maxv = new MaxVLType(mesh, func);
    maxv->ComputeOutput();
    GradType* grad = new GradType(mesh);
    grad->ClearAllGradient();
    RobinsType* robins = new RobinsType(mesh, maxv, grad);
    robins->ComputePairing();

    // asc-manifold dims are a prerequisite for arc tracing in ComputeFromGrad;
    // they write only the ldir bits, which the flow labeling ignores.
    Accurate2D::MeshFuncType* mf = new Accurate2D::MeshFuncType();
    mf->setMeshAndFuncAndMaxVLabeling(mesh, func, maxv);
    TopologicalGradientUsingAlgorithms<MeshType, Accurate2D::MeshFuncType,
                                       GradType>
        talg(mf, mesh, grad);
    talg.setAscendingManifoldDimensions();

    MscT* msc = new MscT(grad, mesh, mf);
    msc->ComputeFromGrad();
    msc->ComputeHierarchy(0.0f);
    msc->SetSelectPersAbs(-1);

    // CPU reference: per-vertex terminal minimum (as grid vertex number) and
    // per-quad terminal maximum (as quad lattice index) via fillGeometry.
    const int QX = X - 1, QY = Y - 1;
    std::vector<int32_t> ref_v((size_t)X * Y, -1), ref_q((size_t)QX * QY, -1);
    auto t0 = std::chrono::steady_clock::now();
    MscT::LivingNodesIterator nit(msc);
    for (nit.begin(); nit.valid(); nit.advance()) {
        const INT_TYPE nid = nit.value();
        const node<float>& nd = msc->getNode(nid);
        if (nd.dim == 0) {
            std::set<INDEX_TYPE> man;
            msc->fillGeometry(nid, man, true);
            const int32_t lab = (int32_t)mesh->VertexNumberFromCellID(nd.cellindex);
            for (std::set<INDEX_TYPE>::const_iterator it = man.begin();
                 it != man.end(); ++it) {
                if (mesh->dimension(*it) != 0) continue;
                ref_v[(size_t)mesh->VertexNumberFromCellID(*it)] = lab;
            }
        } else if (nd.dim == 2) {
            std::set<INDEX_TYPE> man;
            msc->fillGeometry(nid, man, false);
            Vec2l nc;
            mesh->cellid2Coords(nd.cellindex, nc);
            const int32_t lab = (int32_t)((nc[0] - 1) / 2 + ((nc[1] - 1) / 2) * QX);
            for (std::set<INDEX_TYPE>::const_iterator it = man.begin();
                 it != man.end(); ++it) {
                if (mesh->dimension(*it) != 2) continue;
                Vec2l cc;
                mesh->cellid2Coords(*it, cc);
                ref_q[(size_t)((cc[0] - 1) / 2 + ((cc[1] - 1) / 2) * QX)] = lab;
            }
        }
    }
    const double cpu_fill_ms = ms_since(t0);

    // GPU flow labeling on the same gradient bytes
    std::vector<uint8_t> gbytes((size_t)n_cells);
    for (INDEX_TYPE c = 0; c < n_cells; c++)
        gbytes[(size_t)c] = (uint8_t)grad->getAsChar(c);
    std::vector<int32_t> gv((size_t)X * Y), gq((size_t)QX * QY);
    gpu::Dgrad2DTimings tm{};
    bool ok = gpu::Label2DFlowTerminals(gbytes.data(), X, Y, tables, gv.data(),
                                        gq.data(), &tm);
    ok = ok && gpu::Label2DFlowTerminals(gbytes.data(), X, Y, tables, gv.data(),
                                         gq.data(), &tm); // warm timings
    if (!ok) {
        printf("GPU FAIL: Label2DFlowTerminals returned false\n");
        delete msc; delete mf; delete robins; delete grad; delete maxv;
        delete mesh; delete func; delete grid;
        return -1;
    }

    long long mism = 0, unref = 0, shown = 0;
    for (size_t i = 0; i < ref_v.size(); i++) {
        if (ref_v[i] < 0) unref++;
        if (gv[i] != ref_v[i]) {
            if (shown++ < 5)
                printf("  VERT MISMATCH pix %zu: cpu=%d gpu=%d\n", i, ref_v[i],
                       gv[i]);
            mism++;
        }
    }
    for (size_t i = 0; i < ref_q.size(); i++) {
        if (ref_q[i] < 0) unref++;
        if (gq[i] != ref_q[i]) {
            if (shown++ < 5)
                printf("  QUAD MISMATCH q %zu: cpu=%d gpu=%d\n", i, ref_q[i],
                       gq[i]);
            mism++;
        }
    }
    printf("  vertex+quad labels     : %zu (+%zu quads), cpu-unlabeled=%lld\n",
           ref_v.size(), ref_q.size(), unref);
    printf("  mismatches             : %lld  -> %s\n", mism,
           mism == 0 ? "PASS" : "FAIL");
    printf("  CPU fillGeometry paint : %8.2f ms (serial recursion + std::set)\n",
           cpu_fill_ms);
    printf("  GPU flow h2d/kern/d2h  : %8.2f / %.2f / %.2f ms\n", tm.h2d_ms,
           tm.kernel_ms, tm.d2h_ms);

    delete msc; delete mf; delete robins; delete grad; delete maxv; delete mesh;
    delete func; delete grid;
    return mism + (unref ? 1 : 0);
}

// Batched-slices gate: S independent images in one launch, each slice's full
// gradient memcmp-equal to its per-slice CPU reference; kernel time compared
// against S x the single-slice kernel time.
static long long RunBatchedCase(int X, int Y, int S,
                                const gpu::Dgrad2DTables& tables) {
    printf("\n=== batched case (%dx%d x %d slices) ===\n", X, Y, S);
    std::mt19937 rng(777);
    std::uniform_real_distribution<float> uni(0.f, 1.f);
    std::vector<float> values((size_t)X * Y * S);
    for (auto& v : values) v = uni(rng);

    GridType* grid = new GridType(Vec2i(X, Y), Vec2b(false, false));
    MeshType* mesh = new MeshType(grid);
    const INDEX_TYPE n_cells = mesh->numCells();

    // per-slice CPU references
    std::vector<uint8_t> cpu_ref((size_t)n_cells * S);
    auto t0 = std::chrono::steady_clock::now();
    for (int s = 0; s < S; s++) {
        GridFuncType func(grid, values.data() + (size_t)s * X * Y);
        MaxVLType maxv(mesh, &func);
        maxv.ComputeOutput();
        GradType grad(mesh);
        grad.ClearAllGradient();
        RobinsType robins(mesh, &maxv, &grad);
        robins.ComputePairing();
        for (INDEX_TYPE c = 0; c < n_cells; c++)
            cpu_ref[(size_t)s * n_cells + c] = (uint8_t)grad.getAsChar(c);
    }
    const double cpu_ms = ms_since(t0);

    // one batched GPU launch (2nd run for warm timings), plus a single-slice
    // reference timing for the linear-scaling check
    std::vector<uint8_t> gpu_out((size_t)n_cells * S);
    gpu::Dgrad2DTimings tm{}, tm1{};
    bool ok = gpu::ComputeDiscreteGradient2DBatched(values.data(), X, Y, S,
                                                    tables, gpu_out.data(), 0,
                                                    &tm);
    ok = ok && gpu::ComputeDiscreteGradient2DBatched(values.data(), X, Y, S,
                                                     tables, gpu_out.data(), 0,
                                                     &tm);
    std::vector<uint8_t> one((size_t)n_cells);
    ok = ok && gpu::ComputeDiscreteGradient2DBatched(values.data(), X, Y, 1,
                                                     tables, one.data(), 0, &tm1);
    if (!ok) {
        printf("GPU FAIL: batched call returned false\n");
        delete mesh; delete grid;
        return -1;
    }

    long long mismatches = 0;
    for (size_t i = 0; i < cpu_ref.size(); i++)
        if (cpu_ref[i] != gpu_out[i]) {
            if (mismatches < 5)
                printf("  MISMATCH slice %lld cell %lld: cpu=0x%02x gpu=0x%02x\n",
                       (long long)(i / n_cells), (long long)(i % n_cells),
                       cpu_ref[i], gpu_out[i]);
            mismatches++;
        }

    printf("  mismatches             : %lld  -> %s\n", mismatches,
           mismatches == 0 ? "PASS" : "FAIL");
    printf("  CPU %d slices          : %8.2f ms (label+pairing, serial loop)\n",
           S, cpu_ms);
    printf("  GPU batch h2d/kern/d2h : %8.2f / %.2f / %.2f ms\n", tm.h2d_ms,
           tm.kernel_ms, tm.d2h_ms);
    printf("  scaling: batch kernel %.2f ms vs %d x single %.2f = %.2f ms "
           "(%.2fx of linear)\n",
           tm.kernel_ms, S, tm1.kernel_ms, S * tm1.kernel_ms,
           tm.kernel_ms / (S * tm1.kernel_ms));

    delete mesh; delete grid;
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
    // tiny images: boundary classes dominate (corners, 1-wide strips)
    total += RunCase("tiny-3x3", 3, 3, noise, tables, false);
    total += RunCase("tiny-5x3", 5, 3, qnoise, tables, false);
    // batched slices
    total += RunBatchedCase(512, 512, 64, tables);
    // Stage 4: flow-terminal labeling vs fillGeometry base manifolds
    total += RunLabelCase("label-qnoise-512", 512, 512, qnoise, tables);
    total += RunLabelCase("label-bumps-256", 256, 192, bumps, tables);
    total += RunLabelCase("label-noise-1024", 1024, 1024, noise, tables);
    total += RunLabelCase("label-qnoise-4096", 4096, 4096, qnoise, tables);
    // baseline/timing set
    total += RunCase("noise-1024", 1024, 1024, noise, tables, true);
    total += RunCase("qnoise-4096", 4096, 4096, qnoise, tables, true);

    printf("\n=== TOTAL mismatches: %lld -> %s ===\n", total,
           total == 0 ? "ALL PASS" : "FAIL");
    return total == 0 ? 0 : 1;
}
