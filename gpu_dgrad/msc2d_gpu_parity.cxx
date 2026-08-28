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

    printf("\n=== PARITY TOTAL: %lld -> %s ===\n", total,
           total == 0 ? "ALL PASS" : "FAIL");
    return total == 0 ? 0 : 1;
}
