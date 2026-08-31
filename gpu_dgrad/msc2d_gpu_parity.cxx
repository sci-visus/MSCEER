// Stage 4 end-to-end parity gate: the full msc_2d_lib pipeline run with the
// CPU gradient and with the GPU gradient (ComputeOptions::useGpuGradient)
// must produce IDENTICAL results - critical points, arc geometry, and
// ascending/descending manifold label images, in both builder modes, at
// multiple persistence values. The gradient bytes are bit-equal (proven by
// dgrad2d_validate), and the downstream pipeline is deterministic, so any
// difference here is an integration bug.
//
// accurateAsc/Dsc are off: the numeric re-integration path is not part of the
// combinatorial parity contract.
//
// Exit code 0 iff everything matches.

#include "msc_2d_lib.h"
#include "dgrad_gpu_api.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

using GInt::Msc2D::Msc2D;
using GInt::Msc2D::LabelImage;
using GInt::Msc2D::CriticalPoint;
using GInt::Msc2D::ArcGeometry;

static double ms_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(
               std::chrono::steady_clock::now() - t0)
        .count();
}

static long long compareCrits(const std::vector<CriticalPoint>& a,
                              const std::vector<CriticalPoint>& b) {
    if (a.size() != b.size()) {
        printf("    criticalPoints size %zu vs %zu\n", a.size(), b.size());
        return (long long)(a.size() > b.size() ? a.size() - b.size()
                                               : b.size() - a.size());
    }
    long long bad = 0;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i].id != b[i].id || a[i].x != b[i].x || a[i].y != b[i].y ||
            a[i].dim != b[i].dim || a[i].value != b[i].value)
            bad++;
    }
    return bad;
}

static long long compareLabels(const char* what, const LabelImage& a,
                               const LabelImage& b) {
    if (a.width != b.width || a.height != b.height ||
        a.labels.size() != b.labels.size()) {
        printf("    %s shape mismatch\n", what);
        return 1;
    }
    long long bad = 0;
    for (size_t i = 0; i < a.labels.size(); i++)
        if (a.labels[i] != b.labels[i]) bad++;
    if (bad) printf("    %s: %lld differing pixels\n", what, bad);
    return bad;
}

// Arcs are compared as an order-independent multiset keyed by endpoints, dim
// and full polyline, with arc ids excluded: in serial mode with geometry on,
// createArc runs under "#pragma omp critical" inside a parallel node loop
// (gi_morse_smale_complex_basic.h:512), so arc creation ORDER - and thus arc
// ids - is nondeterministic run to run on the CPU itself. Node ids ARE
// deterministic (serial creation scan), so keys are stable.
static long long compareArcs(std::vector<ArcGeometry> a,
                             std::vector<ArcGeometry> b) {
    if (a.size() != b.size()) {
        printf("    arcGeometry size %zu vs %zu\n", a.size(), b.size());
        return 1;
    }
    const auto keyLess = [](const ArcGeometry& p, const ArcGeometry& q) {
        if (p.lowerNodeId != q.lowerNodeId) return p.lowerNodeId < q.lowerNodeId;
        if (p.upperNodeId != q.upperNodeId) return p.upperNodeId < q.upperNodeId;
        if (p.dim != q.dim) return p.dim < q.dim;
        if (p.line.size() != q.line.size()) return p.line.size() < q.line.size();
        for (size_t j = 0; j < p.line.size(); j++) {
            if (p.line[j].x != q.line[j].x) return p.line[j].x < q.line[j].x;
            if (p.line[j].y != q.line[j].y) return p.line[j].y < q.line[j].y;
        }
        return false;
    };
    std::sort(a.begin(), a.end(), keyLess);
    std::sort(b.begin(), b.end(), keyLess);
    long long bad = 0;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i].lowerNodeId != b[i].lowerNodeId ||
            a[i].upperNodeId != b[i].upperNodeId || a[i].dim != b[i].dim ||
            a[i].line.size() != b[i].line.size()) {
            bad++;
            continue;
        }
        for (size_t j = 0; j < a[i].line.size(); j++)
            if (a[i].line[j].x != b[i].line[j].x ||
                a[i].line[j].y != b[i].line[j].y) {
                bad++;
                break;
            }
    }
    if (bad) printf("    arcGeometry: %lld differing arcs (order-independent)\n", bad);
    return bad;
}

static long long runMode(const char* modename, Msc2D::BuilderMode mode,
                         const std::vector<float>& data, int rows, int cols,
                         float range) {
    printf("  mode %s:\n", modename);
    Msc2D::ComputeOptions opt;
    opt.builderMode = mode;
    opt.requestedParallelism = 4;
    opt.accurateAsc = false;
    opt.accurateDsc = false;
    opt.buildArcGeometry[1] = true; // exercise arc geometry through the pipeline

    long long bad = 0;
    Msc2D cpu, gpu;

    opt.useGpuGradient = false;
    auto t0 = std::chrono::steady_clock::now();
    cpu.compute(data.data(), rows, cols, opt);
    const double cpu_ms = ms_since(t0);

    opt.useGpuGradient = true;
    t0 = std::chrono::steady_clock::now();
    gpu.compute(data.data(), rows, cols, opt);
    const double gpu_ms = ms_since(t0);

    bad += compareCrits(cpu.criticalPoints(), gpu.criticalPoints());
    bad += compareArcs(cpu.arcGeometry(), gpu.arcGeometry());
    bad += compareLabels("asc2 @cancel", cpu.ascending2Manifolds(),
                         gpu.ascending2Manifolds());
    bad += compareLabels("dsc2 @cancel", cpu.descending2Manifolds(),
                         gpu.descending2Manifolds());

    // re-threshold and compare again
    const float p2 = 0.02f * range;
    cpu.setPersistence(p2);
    gpu.setPersistence(p2);
    bad += compareCrits(cpu.criticalPoints(), gpu.criticalPoints());
    bad += compareLabels("asc2 @p2", cpu.ascending2Manifolds(),
                         gpu.ascending2Manifolds());
    bad += compareLabels("dsc2 @p2", cpu.descending2Manifolds(),
                         gpu.descending2Manifolds());

    printf("    compute() wall: cpu-gradient %.1f ms, gpu-gradient %.1f ms\n",
           cpu_ms, gpu_ms);
    printf("    -> %s\n", bad == 0 ? "PASS" : "FAIL");
    return bad;
}

// ---------------------------------------------------------------------------
// Stage 6: persistent device label context (GInt::gpu::LabelCtx2D).
//
// The gather is integer and exact, so the gate is memcmp against a CPU
// reference over random compact images and random remaps - including the
// degenerate shapes a caller can legitimately hand it (m == 0, an
// all-unlabeled image, an all -1 remap, a re-upload onto a live context) and
// the out-of-range compact ids the contract pins to -1.
// ---------------------------------------------------------------------------

static void cpuGather(const std::vector<int32_t>& base,
                      const std::vector<int32_t>& remap, int32_t m,
                      std::vector<int32_t>& out) {
    out.assign(base.size(), -1);
    for (size_t i = 0; i < base.size(); i++) {
        const int32_t b = base[i];
        out[i] = (b < 0 || b >= m) ? -1 : remap[(size_t)b];
    }
}

static long long runLabelCtxSection() {
    printf("\n=== label ctx (LabelCtx2D) ===\n");
    long long bad = 0;

    GInt::gpu::LabelCtx2D* probe = GInt::gpu::CreateLabelCtx2D(8, 8);
    if (!probe) {
        printf("    -> SKIP (no device)\n");
        return 0;
    }
    GInt::gpu::DestroyLabelCtx2D(probe);

    struct Shape {
        const char* name;
        int X, Y;
        int32_t m;      // remap size
        int32_t idmax;  // compact ids drawn from [-1, idmax)
    } shapes[] = {
        {"257x129 dense", 257, 129, 400, 400},
        {"512x384 sparse", 512, 384, 5000, 5000},
        {"64x64 m=0", 64, 64, 0, 8},
        {"64x64 all-unlabeled", 64, 64, 32, 0},
        {"128x97 out-of-range", 128, 97, 16, 64},
        {"1x1 degenerate", 1, 1, 4, 4},
    };

    std::mt19937 rng(60606);
    for (const Shape& sh : shapes) {
        const size_t n = (size_t)sh.X * sh.Y;
        std::vector<int32_t> base(n), remap((size_t)(sh.m > 0 ? sh.m : 1)), ref, got(n);
        for (size_t i = 0; i < n; i++)
            base[i] = sh.idmax > 0
                          ? (int32_t)((int)(rng() % (unsigned)(sh.idmax + 1)) - 1)
                          : -1;
        for (int32_t j = 0; j < sh.m; j++)
            remap[(size_t)j] = (rng() % 5u == 0u) ? -1 : (int32_t)(1000000 - 7 * j);

        GInt::gpu::LabelCtx2D* ctx = GInt::gpu::CreateLabelCtx2D(sh.X, sh.Y);
        if (!ctx) {
            printf("    %s: CreateLabelCtx2D failed\n", sh.name);
            bad++;
            continue;
        }

        long long shape_bad = 0;
        // A relabel before any base upload must fail, not scribble.
        if (GInt::gpu::CtxRelabel(ctx, true, remap.data(), sh.m, got.data(), NULL))
            shape_bad++;
        if (!GInt::gpu::CtxSetBaseLabels(ctx, true, base.data())) shape_bad++;
        if (!GInt::gpu::CtxBaseLabelsDevice(ctx, true)) shape_bad++;
        if (GInt::gpu::CtxBaseLabelsDevice(ctx, false)) shape_bad++;  // dsc not uploaded

        GInt::gpu::Dgrad2DTimings tm;
        if (!GInt::gpu::CtxRelabel(ctx, true, remap.data(), sh.m, got.data(), &tm)) {
            shape_bad++;
        } else {
            cpuGather(base, remap, sh.m, ref);
            if (std::memcmp(got.data(), ref.data(), n * sizeof(int32_t)) != 0)
                shape_bad++;
        }

        // An all -1 remap must blank the image; the output buffer is reused, so
        // this also proves the previous result is fully overwritten.
        if (sh.m > 0) {
            std::vector<int32_t> blank((size_t)sh.m, -1);
            if (!GInt::gpu::CtxRelabel(ctx, true, blank.data(), sh.m, got.data(), NULL)) {
                shape_bad++;
            } else {
                cpuGather(base, blank, sh.m, ref);
                if (std::memcmp(got.data(), ref.data(), n * sizeof(int32_t)) != 0)
                    shape_bad++;
            }
        }

        // Re-upload a different base labeling onto the live context, and use the
        // descending direction to prove the two directions are independent.
        std::vector<int32_t> base2(n);
        for (size_t i = 0; i < n; i++)
            base2[i] = sh.idmax > 0
                           ? (int32_t)((int)(rng() % (unsigned)(sh.idmax + 1)) - 1)
                           : -1;
        if (!GInt::gpu::CtxSetBaseLabels(ctx, true, base2.data())) shape_bad++;
        if (!GInt::gpu::CtxSetBaseLabels(ctx, false, base.data())) shape_bad++;
        if (!GInt::gpu::CtxRelabel(ctx, true, remap.data(), sh.m, got.data(), NULL)) {
            shape_bad++;
        } else {
            cpuGather(base2, remap, sh.m, ref);
            if (std::memcmp(got.data(), ref.data(), n * sizeof(int32_t)) != 0)
                shape_bad++;
        }

        // The no-download variant must hand back a device pointer and leave the
        // same image behind: the following downloading call reads that buffer.
        const void* dev = NULL;
        if (!GInt::gpu::CtxRelabelDevice(ctx, false, remap.data(), sh.m, &dev, NULL) ||
            dev == NULL) {
            shape_bad++;
        } else if (!GInt::gpu::CtxRelabel(ctx, false, remap.data(), sh.m, got.data(),
                                          NULL)) {
            shape_bad++;
        } else {
            cpuGather(base, remap, sh.m, ref);
            if (std::memcmp(got.data(), ref.data(), n * sizeof(int32_t)) != 0)
                shape_bad++;
        }

        GInt::gpu::DestroyLabelCtx2D(ctx);
        printf("    %-22s m=%-6d -> %s\n", sh.name, (int)sh.m,
               shape_bad == 0 ? "ok" : "MISMATCH");
        bad += shape_bad;
    }

    // Leak check: a create/upload/relabel/destroy loop must return the device to
    // its starting free-byte count. The baseline is taken after the first
    // iteration, which also pays one-time runtime/context allocations.
    const int X = 512, Y = 512, m = 4096;
    std::vector<int32_t> base((size_t)X * Y), remap(m), out((size_t)X * Y);
    for (size_t i = 0; i < base.size(); i++) base[i] = (int32_t)(i % (m + 1)) - 1;
    for (int j = 0; j < m; j++) remap[(size_t)j] = j * 3 - 1;

    int64_t before = 0, after = 0;
    for (int it = 0; it < 200; it++) {
        GInt::gpu::LabelCtx2D* ctx = GInt::gpu::CreateLabelCtx2D(X, Y);
        if (!ctx) {
            bad++;
            break;
        }
        if (!GInt::gpu::CtxSetBaseLabels(ctx, true, base.data())) bad++;
        if (!GInt::gpu::CtxSetBaseLabels(ctx, false, base.data())) bad++;
        if (!GInt::gpu::CtxRelabel(ctx, true, remap.data(), m, out.data(), NULL)) bad++;
        if (!GInt::gpu::CtxRelabel(ctx, false, remap.data(), m, out.data(), NULL)) bad++;
        GInt::gpu::DestroyLabelCtx2D(ctx);
        if (it == 0 && !GInt::gpu::QueryDeviceMemory(&before, NULL)) bad++;
    }
    if (!GInt::gpu::QueryDeviceMemory(&after, NULL)) bad++;
    const long long leaked = (long long)(before - after);
    if (leaked != 0) bad++;
    printf("    leak loop 200x create/upload/relabel/destroy: delta=%lld bytes -> %s\n",
           leaked, leaked == 0 ? "ok" : "LEAK");

    printf("    -> %s\n", bad == 0 ? "PASS" : "FAIL");
    return bad;
}

int main() {
    long long total = 0;

    struct Case {
        const char* name;
        int rows, cols;
        int kind; // 0 noise, 1 bumps
    } cases[] = {{"noise-512", 512, 512, 0}, {"bumps-512x384", 512, 384, 1}};

    for (const Case& c : cases) {
        printf("\n=== parity case %s (%dx%d) ===\n", c.name, c.rows, c.cols);
        std::mt19937 rng(2024);
        std::uniform_real_distribution<float> uni(0.f, 1.f);
        std::vector<float> data((size_t)c.rows * c.cols);
        float mn = 1e30f, mx = -1e30f;
        for (int i = 0; i < c.rows; i++)
            for (int j = 0; j < c.cols; j++) {
                float v;
                if (c.kind == 0) {
                    v = uni(rng);
                } else {
                    auto g = [&](float cx, float cy, float s) {
                        const float dx = j - cx, dy = i - cy;
                        return std::exp(-(dx * dx + dy * dy) / (2 * s * s));
                    };
                    v = g(60, 70, 30) + 0.8f * g(300, 200, 70) +
                        0.5f * g(120, 400, 50) + 0.01f * uni(rng);
                }
                data[(size_t)i * c.cols + j] = v;
                mn = v < mn ? v : mn;
                mx = v > mx ? v : mx;
            }

        total += runMode("serial     ", Msc2D::BuilderMode::Serial, data, c.rows,
                         c.cols, mx - mn);
        total += runMode("partitioned", Msc2D::BuilderMode::Partitioned, data,
                         c.rows, c.cols, mx - mn);
    }

    total += runLabelCtxSection();

    printf("\n=== PARITY TOTAL: %lld -> %s ===\n", total,
           total == 0 ? "ALL PASS" : "FAIL");
    return total == 0 ? 0 : 1;
}
