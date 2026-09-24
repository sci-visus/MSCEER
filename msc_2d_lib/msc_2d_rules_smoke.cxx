// Custom simplification rules smoke test (ComputeOptions::rules, max-area).
//
// The field is a noisy egg-crate: a regular grid of wells, so every basin has
// roughly the same area and, with the cancel persistence above the whole value
// range, an UNRESTRICTED simplification merges nearly everything. With a
// max-area rule on the minima (or maxima) the result is known in shape if not
// in detail: no region formed by merging two or more base regions may exceed
// the cap, and many more regions survive than without the rule.
//
// Checked for both Simplification modes and both directions:
//   1. neutrality - a cap larger than the image changes nothing: labels are
//      identical to the rule-free run (the MSC subclass and the forest's gated
//      serial build reproduce the ungated order exactly);
//   2. the cap - every merged living region has <= cap lattice points;
//   3. the rule bites - strictly more living regions than the rule-free run.

#include "msc_2d_lib.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <random>
#include <vector>

namespace {

typedef GInt::Msc2D::Msc2D M2;

const int kDim = 160;
const float kPeriod = 16.0f;       // wells every 16 px -> ~100 basins of ~256 px
const long long kCap = 1200;       // a merged region may hold <= this many points

std::vector<float> buildField(bool negate) {
    std::vector<float> field(static_cast<size_t>(kDim) * kDim);
    std::mt19937 rng(7u);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    const float w = 2.0f * 3.14159265f / kPeriod;
    for (int r = 0; r < kDim; ++r) {
        for (int c = 0; c < kDim; ++c) {
            // slow tilt so well depths are distinct and the merge order has
            // somewhere to go, small noise as a tie-break
            float v = std::sin(w * c) * std::sin(w * r) + 0.002f * (c + 0.7f * r) +
                      0.01f * dist(rng);
            field[static_cast<size_t>(r) * kDim + c] = negate ? -v : v;
        }
    }
    return field;
}

struct RunResult {
    std::vector<int> labels;        // painted living labels at the selected persistence
    std::vector<int> base;          // base labeling
    std::vector<int> baseToLiving;
    int living = 0;
};

RunResult run(M2::Simplification mode, bool ascending, bool rulesOn, long long cap) {
    const std::vector<float> field = buildField(!ascending);
    M2::ComputeOptions opt;
    opt.simplification = mode;
    opt.basePersistenceAbs = 0.0f;
    opt.cancelPersistenceAbs = 100.0f;  // far above the value range
    opt.requestedParallelism = 4;
    if (rulesOn) {
        opt.rules.restrictMinima = ascending;
        opt.rules.restrictMaxima = !ascending;
        opt.rules.maxRegionCells = cap;
    }
    M2 msc;
    msc.compute(field.data(), kDim, kDim, opt);
    msc.setPersistence(100.0f);
    RunResult out;
    GInt::Msc2D::LabelImage img = ascending ? msc.ascending2Manifolds() : msc.descending2Manifolds();
    out.labels = img.labels;
    out.base = msc.baseLabeling(ascending);
    out.baseToLiving = msc.baseToLiving(ascending);
    std::map<int, int> ids;
    for (size_t i = 0; i < out.labels.size(); ++i)
        if (out.labels[i] >= 0) ids[out.labels[i]] = 1;
    out.living = static_cast<int>(ids.size());
    return out;
}

// Largest living region that absorbed >= 2 base regions, in base-labeling
// points (base region sizes summed through baseToLiving).
long long largestMergedRegion(const RunResult& r) {
    std::vector<long long> baseSize;
    for (size_t i = 0; i < r.base.size(); ++i) {
        const int b = r.base[i];
        if (b < 0) continue;
        if (static_cast<size_t>(b) >= baseSize.size()) baseSize.resize(b + 1, 0);
        baseSize[b]++;
    }
    std::map<int, long long> size, parts;
    for (size_t b = 0; b < baseSize.size() && b < r.baseToLiving.size(); ++b) {
        const int l = r.baseToLiving[b];
        if (l < 0) continue;
        size[l] += baseSize[b];
        parts[l] += 1;
    }
    long long worst = 0;
    for (std::map<int, long long>::const_iterator it = size.begin(); it != size.end(); ++it)
        if (parts[it->first] >= 2 && it->second > worst) worst = it->second;
    return worst;
}

bool check(M2::Simplification mode, bool ascending) {
    const char* tag = mode == M2::Simplification::MergeForest ? "forest" : "msc";
    const char* dir = ascending ? "min/asc" : "max/dsc";
    bool ok = true;

    const RunResult free_ = run(mode, ascending, false, -1);
    const RunResult huge = run(mode, ascending, true, 10LL * kDim * kDim);
    const RunResult capped = run(mode, ascending, true, kCap);
    const long long worst = largestMergedRegion(capped);

    std::cout << "rules " << tag << " " << dir << ": living free=" << free_.living
              << " hugeCap=" << huge.living << " cap" << kCap << "=" << capped.living
              << " largestMerged=" << worst << " (free largestMerged="
              << largestMergedRegion(free_) << ")" << std::endl;

    if (huge.labels != free_.labels) {
        std::cerr << "FAIL[" << tag << " " << dir << "]: an unreachable cap changed the labels"
                  << std::endl;
        ok = false;
    }
    if (worst > kCap) {
        std::cerr << "FAIL[" << tag << " " << dir << "]: merged region of " << worst
                  << " points exceeds the cap " << kCap << std::endl;
        ok = false;
    }
    if (capped.living <= free_.living) {
        std::cerr << "FAIL[" << tag << " " << dir << "]: the cap did not keep more regions alive"
                  << std::endl;
        ok = false;
    }
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    const M2::Simplification modes[2] = {
        M2::Simplification::MscHierarchy, M2::Simplification::MergeForest };
    for (int m = 0; m < 2; ++m) {
        ok = check(modes[m], true) && ok;
        ok = check(modes[m], false) && ok;
    }
    std::cout << (ok ? "PASS" : "FAIL") << std::endl;
    return ok ? 0 : 1;
}
