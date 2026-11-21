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
//
// And for the parallel builders (builderMode = Partitioned, parallelism 4):
//   4. forest, banded: labels identical to the serial gated build (the banded
//      build is order-consistent with the serial one under a gate);
//   5. MSC, partitioned, with a base persistence high enough that the local
//      phase cancels: the cap holds and the rule bites. The cap is checked
//      against RAW flow basins (a serial forest at persistence 0), since the
//      partitioned base regions are themselves merged;
//   6. MSC, partitioned, allowParallel = false: identical to the serial run.

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

struct RunConfig {
    bool partitioned = false;
    bool allowParallel = true;
    float basePersistence = 0.0f;
    float selectPersistence = 100.0f;  // far above the value range
};

RunResult run(M2::Simplification mode, bool ascending, bool rulesOn, long long cap,
              const RunConfig& cfg = RunConfig()) {
    const std::vector<float> field = buildField(!ascending);
    M2::ComputeOptions opt;
    opt.simplification = mode;
    opt.builderMode = cfg.partitioned ? M2::BuilderMode::Partitioned : M2::BuilderMode::Serial;
    opt.basePersistenceAbs = cfg.basePersistence;
    opt.cancelPersistenceAbs = cfg.selectPersistence;
    opt.requestedParallelism = 4;
    // The CPU accurate gradient is not reproducible run to run, which would
    // make every exact cross-run comparison below meaningless.
    opt.accurateAsc = false;
    opt.accurateDsc = false;
    if (rulesOn) {
        opt.rules.restrictMinima = ascending;
        opt.rules.restrictMaxima = !ascending;
        opt.rules.maxRegionCells = cap;
        opt.rules.allowParallel = cfg.allowParallel;
    }
    M2 msc;
    msc.compute(field.data(), kDim, kDim, opt);
    msc.setPersistence(cfg.selectPersistence);
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

// Largest living region that covers >= 2 raw flow basins (a basin counts
// when the majority of its points carry that living label), in points.
long long largestMergedByRaw(const std::vector<int>& labels, const std::vector<int>& raw) {
    std::map<int, long long> rawSize, livingSize;
    std::map<std::pair<int, int>, long long> overlap;
    for (size_t i = 0; i < labels.size() && i < raw.size(); ++i) {
        if (raw[i] >= 0) rawSize[raw[i]]++;
        if (labels[i] >= 0) livingSize[labels[i]]++;
        if (raw[i] >= 0 && labels[i] >= 0) overlap[std::make_pair(labels[i], raw[i])]++;
    }
    std::map<int, int> basins;
    for (std::map<std::pair<int, int>, long long>::const_iterator it = overlap.begin();
         it != overlap.end(); ++it)
        if (2 * it->second > rawSize[it->first.second]) basins[it->first.first]++;
    long long worst = 0;
    for (std::map<int, long long>::const_iterator it = livingSize.begin(); it != livingSize.end(); ++it)
        if (basins[it->first] >= 2 && it->second > worst) worst = it->second;
    return worst;
}

// Same partition of the image up to a renaming of the ids (node ids are not
// stable across builds whose threaded stages differ).
bool samePartition(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size()) return false;
    std::map<int, int> ab, ba;
    for (size_t i = 0; i < a.size(); ++i) {
        std::map<int, int>::iterator x = ab.find(a[i]);
        std::map<int, int>::iterator y = ba.find(b[i]);
        if (x == ab.end()) ab[a[i]] = b[i]; else if (x->second != b[i]) return false;
        if (y == ba.end()) ba[b[i]] = a[i]; else if (y->second != a[i]) return false;
    }
    return true;
}

// Base labeling of a serial forest at persistence 0: the raw flow basins.
std::vector<int> rawBasins(bool ascending) {
    RunConfig cfg;
    cfg.selectPersistence = 0.0f;
    return run(M2::Simplification::MergeForest, ascending, false, -1, cfg).base;
}

bool checkParallel(bool ascending) {
    const char* dir = ascending ? "min/asc" : "max/dsc";
    bool ok = true;

    // 4. banded forest == serial forest, both gated
    {
        RunConfig par;
        par.partitioned = true;
        const RunResult serial = run(M2::Simplification::MergeForest, ascending, true, kCap);
        const RunResult banded = run(M2::Simplification::MergeForest, ascending, true, kCap, par);
        std::cout << "rules forest-banded " << dir << ": living serial=" << serial.living
                  << " banded=" << banded.living << std::endl;
        if (banded.labels != serial.labels) {
            std::cerr << "FAIL[forest-banded " << dir
                      << "]: gated banded build differs from the gated serial build" << std::endl;
            ok = false;
        }
    }

    // 5./6. partitioned MSC, local phase active
    {
        RunConfig par;
        par.partitioned = true;
        par.basePersistence = 1.0f;
        RunConfig strict = par;
        strict.allowParallel = false;
        RunConfig serialCfg;
        serialCfg.basePersistence = par.basePersistence;

        const std::vector<int> raw = rawBasins(ascending);
        const RunResult free_ = run(M2::Simplification::MscHierarchy, ascending, false, -1, par);
        const RunResult capped = run(M2::Simplification::MscHierarchy, ascending, true, kCap, par);
        const RunResult strictRun = run(M2::Simplification::MscHierarchy, ascending, true, kCap, strict);
        const RunResult serial = run(M2::Simplification::MscHierarchy, ascending, true, kCap, serialCfg);
        const long long worst = largestMergedByRaw(capped.labels, raw);
        std::cout << "rules msc-partitioned " << dir << ": living free=" << free_.living
                  << " cap" << kCap << "=" << capped.living << " largestMerged=" << worst
                  << " (free largestMerged=" << largestMergedByRaw(free_.labels, raw)
                  << ") serial=" << serial.living << " strict=" << strictRun.living << std::endl;
        if (worst > kCap) {
            std::cerr << "FAIL[msc-partitioned " << dir << "]: merged region of " << worst
                      << " points exceeds the cap " << kCap << std::endl;
            ok = false;
        }
        if (capped.living <= free_.living) {
            std::cerr << "FAIL[msc-partitioned " << dir
                      << "]: the cap did not keep more regions alive" << std::endl;
            ok = false;
        }
        if (!samePartition(strictRun.labels, serial.labels)) {
            std::cerr << "FAIL[msc-partitioned " << dir
                      << "]: allowParallel=false differs from the serial build" << std::endl;
            ok = false;
        }
    }
    return ok;
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
    ok = checkParallel(true) && ok;
    ok = checkParallel(false) && ok;
    std::cout << (ok ? "PASS" : "FAIL") << std::endl;
    return ok ? 0 : 1;
}
