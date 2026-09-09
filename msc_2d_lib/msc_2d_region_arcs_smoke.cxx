// livingRegionArcs() smoke test.
//
// The field is four Gaussian wells at the corners of a square. That geometry
// makes the ANSWER known in advance, which is the whole point: the four wells
// form a 4-cycle -- each is separated from its two SIDE neighbours by a saddle
// on the segment between their centres, while the DIAGONAL pairs meet only at
// the centre and are therefore not adjacent at all. A query that invented a
// diagonal, dropped a side, or read its saddle value off the wrong cell would
// be wrong in a way no structural invariant could catch.
//
// The assertions are about the four wells specifically, not about the whole
// image, because the whole image has more regions than the wells: GInt's 2D
// complex makes the domain boundary critical, so every local minimum of the
// field ALONG the boundary -- one per edge here, plus the corners -- is a living
// region of its own at any threshold. Those are a property of the library
// boundary handling, not of the adjacency query, so the test states what it
// actually knows: the four well regions, all four sides, neither diagonal, and
// the saddle values.
//
// Those assertions are geometric -- which wells touch, and at what value -- so
// they survive the banded configuration below running at parallelism 4, where
// the threaded gradient build makes the region IDS differ run to run.

#include "msc_2d_lib.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

const int kDim = 192;          // square field, kDim x kDim
// sigma is load-bearing. At 35 the four basins tile the whole domain -- the
// corners drain into the nearest well, so there is no background region and no
// boundary minimum -- while the diagonals still meet only at the centre. Widen
// it to 40 and wells 1 and 2 acquire a shared boundary; narrow it enough and
// the flat corners grow minima of their own. Both were observed.
const float kSigma = 35.0f;
const float kNoise = 1e-4f;    // tie-break only

struct Well {
    float cx;
    float cy;
    float depth;
};

// Distinct depths so the merge order is unambiguous; corners of a square with
// enough clearance that the four barriers are well separated from the wells.
const Well kWells[4] = {
    { 48.0f,  48.0f, 1.00f },
    { 144.0f, 48.0f, 0.90f },
    { 48.0f,  144.0f, 0.80f },
    { 144.0f, 144.0f, 0.70f }
};

// Side pairs of the square, by index into kWells. The diagonals (0,3) and (1,2)
// are deliberately absent.
const int kExpectedPairs[4][2] = { { 0, 1 }, { 0, 2 }, { 1, 3 }, { 2, 3 } };

float wellField(float x, float y) {
    float v = 0.0f;
    for (int w = 0; w < 4; ++w) {
        const float dx = x - kWells[w].cx;
        const float dy = y - kWells[w].cy;
        v -= kWells[w].depth *
             std::exp(-(dx * dx + dy * dy) / (2.0f * kSigma * kSigma));
    }
    return v;
}

// Row-major (row = y, col = x), matching compute()'s rowMajorValues contract.
std::vector<float> buildField(bool negate) {
    std::vector<float> field(static_cast<size_t>(kDim) * static_cast<size_t>(kDim), 0.0f);
    std::mt19937 rng(2024u);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    for (int r = 0; r < kDim; ++r) {
        for (int c = 0; c < kDim; ++c) {
            float v = wellField(static_cast<float>(c), static_cast<float>(r)) +
                      kNoise * dist(rng);
            if (negate) v = -v;
            field[static_cast<size_t>(r) * static_cast<size_t>(kDim) +
                  static_cast<size_t>(c)] = v;
        }
    }
    return field;
}

// The analytic barrier between two wells: the highest point on the segment
// joining their centres. NOT the midpoint -- unequal depths shift the saddle
// toward the shallower well by enough (~0.01) to matter at this tolerance.
float analyticBarrier(int wa, int wb) {
    float best = -1e30f;
    const int kSamples = 2001;
    for (int i = 0; i < kSamples; ++i) {
        const float t = static_cast<float>(i) / static_cast<float>(kSamples - 1);
        const float v = wellField(kWells[wa].cx + t * (kWells[wb].cx - kWells[wa].cx),
                                  kWells[wa].cy + t * (kWells[wb].cy - kWells[wa].cy));
        if (v > best) best = v;
    }
    return best;
}

typedef std::pair<int, int> Pair;

std::set<Pair> pairSet(const std::vector<GInt::Msc2D::RegionArc>& arcs) {
    std::set<Pair> out;
    for (size_t i = 0; i < arcs.size(); ++i) {
        out.insert(std::make_pair(arcs[i].a, arcs[i].b));
    }
    return out;
}

void printArcs(const char* tag, const std::vector<GInt::Msc2D::RegionArc>& arcs) {
    std::cout << "region_arcs " << tag << " count=" << arcs.size();
    for (size_t i = 0; i < arcs.size(); ++i) {
        std::cout << " (" << arcs[i].a << "," << arcs[i].b << ",sv="
                  << arcs[i].saddleValue << ",n=" << arcs[i].count << ")";
    }
    std::cout << std::endl;
}

// --- invariants that must hold in every mode, direction and threshold ------

// a < b, no self-loops, no duplicate unordered pair, sorted by (a, b), and a
// positive collapse count.
bool checkStructure(const char* tag, const std::vector<GInt::Msc2D::RegionArc>& arcs) {
    std::set<Pair> seen;
    for (size_t i = 0; i < arcs.size(); ++i) {
        const GInt::Msc2D::RegionArc& e = arcs[i];
        if (e.a >= e.b) {
            std::cerr << "FAIL[" << tag << "]: entry " << i << " is not ordered a<b ("
                      << e.a << "," << e.b << ")" << std::endl;
            return false;
        }
        if (e.count < 1) {
            std::cerr << "FAIL[" << tag << "]: entry " << i << " has count " << e.count
                      << std::endl;
            return false;
        }
        if (!seen.insert(std::make_pair(e.a, e.b)).second) {
            std::cerr << "FAIL[" << tag << "]: duplicate pair (" << e.a << "," << e.b
                      << ")" << std::endl;
            return false;
        }
        if (i > 0) {
            const GInt::Msc2D::RegionArc& p = arcs[i - 1];
            if (p.a > e.a || (p.a == e.a && p.b >= e.b)) {
                std::cerr << "FAIL[" << tag << "]: output not sorted at entry " << i
                          << std::endl;
                return false;
            }
        }
    }
    return true;
}

// Every id handed out must be reachable through baseToLiving at this same
// persistence -- otherwise a caller holding the remap cannot address it.
bool checkIdSpace(const char* tag, const std::vector<GInt::Msc2D::RegionArc>& arcs,
                  const std::vector<int>& baseToLiving) {
    std::set<int> living;
    for (size_t i = 0; i < baseToLiving.size(); ++i) {
        if (baseToLiving[i] >= 0) living.insert(baseToLiving[i]);
    }
    for (size_t i = 0; i < arcs.size(); ++i) {
        for (int k = 0; k < 2; ++k) {
            const int id = (k == 0) ? arcs[i].a : arcs[i].b;
            if (living.find(id) == living.end()) {
                std::cerr << "FAIL[" << tag << "]: region id " << id
                          << " is not in baseToLiving (living regions="
                          << living.size() << ")" << std::endl;
                return false;
            }
        }
    }
    return true;
}

bool sameArcs(const std::vector<GInt::Msc2D::RegionArc>& a,
              const std::vector<GInt::Msc2D::RegionArc>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].a != b[i].a || a[i].b != b[i].b || a[i].count != b[i].count) return false;
        if (a[i].saddleValue != b[i].saddleValue) return false;
    }
    return true;
}

// The configurations under test. The banded one is not redundant: at
// parallelism 1 ExtractArcsBanded2D emits a single band and the cross-arc list
// is empty, so nothing else here would ever exercise the retention of arcs that
// span bands -- and those are exactly the arcs whose endpoints the forest build
// freezes and replays out of the residue heap.
struct Config {
    const char* name;
    GInt::Msc2D::Msc2D::Simplification simplification;
    GInt::Msc2D::Msc2D::BuilderMode builderMode;
    int parallelism;
    bool gpuGradient;
};

// The GPU one is forest mode production configuration downstream, and it is the
// path that produces the flow-terminal maps the retained arcs are read against.
// In a build without CUDA it prints a notice and falls back to the CPU, which
// makes it a harmless duplicate rather than a skipped case.
const int kNumConfigs = 4;
const Config kConfigs[kNumConfigs] = {
    { "msc", GInt::Msc2D::Msc2D::Simplification::MscHierarchy,
      GInt::Msc2D::Msc2D::BuilderMode::Serial, 1, false },
    { "forest", GInt::Msc2D::Msc2D::Simplification::MergeForest,
      GInt::Msc2D::Msc2D::BuilderMode::Serial, 1, false },
    { "forest-banded", GInt::Msc2D::Msc2D::Simplification::MergeForest,
      GInt::Msc2D::Msc2D::BuilderMode::Partitioned, 4, false },
    { "forest-gpu", GInt::Msc2D::Msc2D::Simplification::MergeForest,
      GInt::Msc2D::Msc2D::BuilderMode::Serial, 1, true }
};

// Which living region owns each well centre, in the id space livingRegionArcs
// uses. This is the bridge from "well 0" to whatever id the pipeline assigned.
std::vector<int> ownersAtWells(GInt::Msc2D::Msc2D& msc, bool ascending) {
    const GInt::Msc2D::LabelImage img =
        ascending ? msc.ascending2Manifolds() : msc.descending2Manifolds();
    std::vector<int> owners(4, -1);
    for (int w = 0; w < 4; ++w) {
        const int x = static_cast<int>(kWells[w].cx);
        const int y = static_cast<int>(kWells[w].cy);
        owners[w] = img.labels[static_cast<size_t>(y) * static_cast<size_t>(img.width) +
                               static_cast<size_t>(x)];
    }
    return owners;
}

// The headline check: at a threshold that has killed the tie-break noise but is
// far below the inter-well barriers, the adjacency must be exactly the 4-cycle,
// with each saddle value close to the analytic barrier.
// Called only on failure. An unexpected pair count almost always means the
// FIELD grew an extremum somewhere rather than that the adjacency logic is
// wrong, and the only way to tell is to look at them. Kept off the passing path
// because criticalPoints() forces the lazy MSC build in MergeForest mode, and
// the point of that mode is that a segmentation-only caller never pays for one.
void dumpExtrema(const char* tag, GInt::Msc2D::Msc2D& msc, bool ascending) {
    const std::vector<GInt::Msc2D::CriticalPoint> cps = msc.criticalPoints();
    const int extremalDim = ascending ? 0 : 2;
    for (size_t i = 0; i < cps.size(); ++i) {
        if (cps[i].dim != extremalDim) continue;
        std::cerr << "  living extremum " << tag << " id=" << cps[i].id
                  << " x=" << cps[i].x << " y=" << cps[i].y
                  << " value=" << cps[i].value << std::endl;
    }
}

bool checkWellAdjacencyInner(const char* tag, GInt::Msc2D::Msc2D& msc, bool ascending,
                             float threshold, bool negated) {
    msc.setPersistence(threshold);
    const std::vector<int> owners = ownersAtWells(msc, ascending);
    const std::vector<GInt::Msc2D::RegionArc>& arcs = msc.livingRegionArcs(ascending);
    printArcs(tag, arcs);

    // The wells must still be four separate regions at this threshold.
    std::set<int> distinct;
    for (int w = 0; w < 4; ++w) {
        if (owners[w] < 0) {
            std::cerr << "FAIL[" << tag << "]: well " << w << " centre is unlabeled"
                      << std::endl;
            return false;
        }
        distinct.insert(owners[w]);
    }
    if (distinct.size() != 4) {
        std::cerr << "FAIL[" << tag << "]: the four wells collapsed into "
                  << distinct.size() << " regions at t=" << threshold << std::endl;
        return false;
    }

    const std::set<Pair> present = pairSet(arcs);

    // The two diagonals must NOT be adjacent: they meet only at the centre.
    const int kDiagonals[2][2] = { { 0, 3 }, { 1, 2 } };
    for (int i = 0; i < 2; ++i) {
        const int a = owners[kDiagonals[i][0]];
        const int b = owners[kDiagonals[i][1]];
        const Pair key(a < b ? a : b, a < b ? b : a);
        if (present.find(key) != present.end()) {
            std::cerr << "FAIL[" << tag << "]: wells " << kDiagonals[i][0] << " and "
                      << kDiagonals[i][1] << " are diagonal but reported adjacent"
                      << std::endl;
            return false;
        }
    }

    // All four sides must be adjacent, with the saddle value on the barrier.
    // The mesh function samples vertices and the discrete saddle lands on a
    // nearby grid edge, so a couple of thousandths of slack is expected; the
    // tolerance is here to catch a value read from the wrong cell (which is off
    // by tenths), not to pin interpolation.
    const float tol = 0.02f;
    for (int i = 0; i < 4; ++i) {
        const int wa = kExpectedPairs[i][0];
        const int wb = kExpectedPairs[i][1];
        const int a = owners[wa];
        const int b = owners[wb];
        const Pair key(a < b ? a : b, a < b ? b : a);
        if (present.find(key) == present.end()) {
            std::cerr << "FAIL[" << tag << "]: wells " << wa << " and " << wb
                      << " share a saddle but are not reported adjacent" << std::endl;
            return false;
        }
        float barrier = analyticBarrier(wa, wb);
        if (negated) barrier = -barrier;
        for (size_t k = 0; k < arcs.size(); ++k) {
            if (arcs[k].a != key.first || arcs[k].b != key.second) continue;
            const float err = std::fabs(arcs[k].saddleValue - barrier);
            std::cout << "region_arcs_barrier " << tag << " wells=" << wa << "-" << wb
                      << " reported=" << arcs[k].saddleValue << " analytic=" << barrier
                      << " err=" << err << std::endl;
            if (err > tol) {
                std::cerr << "FAIL[" << tag << "]: saddle value off by " << err
                          << " (tolerance " << tol << ")" << std::endl;
                return false;
            }
        }
    }
    std::cout << "region_arcs_wells " << tag << " total_pairs=" << arcs.size()
              << " (4 well-to-well, the rest are boundary regions)" << std::endl;
    return true;
}

bool checkWellAdjacency(const char* tag, GInt::Msc2D::Msc2D& msc, bool ascending,
                        float threshold, bool negated) {
    if (checkWellAdjacencyInner(tag, msc, ascending, threshold, negated)) return true;
    dumpExtrema(tag, msc, ascending);
    return false;
}

}  // namespace

int main() {
    const std::vector<float> field = buildField(false);
    const std::vector<float> negated = buildField(true);

    // The shallowest inter-well merge is |-0.630 - (-0.740)| = 0.11 (wells 1
    // and 3), so 0.05 keeps all four alive; the tie-break noise is 1e-4, so it
    // kills every spurious minimum. Anything between those two works.
    const float kFourCycleT = 0.05f;

    // Per-config arc counts at two thresholds, compared (not gated) at the end.
    size_t ascCountByMode[kNumConfigs][2] = { { 0, 0 } };
    const float kCompareT[2] = { 0.01f, 0.05f };

    for (int m = 0; m < kNumConfigs; ++m) {
        const char* mode = kConfigs[m].name;

        GInt::Msc2D::Msc2D::ComputeOptions options;
        options.requestedParallelism = kConfigs[m].parallelism;
        options.builderMode = kConfigs[m].builderMode;
        options.basePersistenceAbs = 0.0f;
        options.cancelPersistenceAbs = 2.0f;   // hierarchy cap above the full range
        options.simplification = kConfigs[m].simplification;
        options.useGpuGradient = kConfigs[m].gpuGradient;

        GInt::Msc2D::Msc2D msc;
        msc.compute(field.data(), kDim, kDim, options);

        // (1) the known well-to-well adjacency, ascending.
        std::string tag = std::string(mode) + "/asc";
        if (!checkWellAdjacency(tag.c_str(), msc, true, kFourCycleT, false)) return 1;
        for (int k = 0; k < 2; ++k) {
            msc.setPersistence(kCompareT[k]);
            ascCountByMode[m][k] = msc.livingRegionArcs(true).size();
        }
        msc.setPersistence(kFourCycleT);
        const size_t pairsAtFourCycle = msc.livingRegionArcs(true).size();

        // (2) Past the last merge the four wells are one region, so no
        // well-to-well adjacency is left. The output does not go empty: the four
        // domain CORNERS are boundary minima the hierarchy never cancels, and
        // they stay adjacent to each other forever. Asserting emptiness here
        // would be asserting something about GInt boundary handling, not about
        // this query, so the check is the one that is actually about the wells.
        msc.setPersistence(2.0f);
        const std::vector<int> collapsedOwners = ownersAtWells(msc, true);
        const std::vector<GInt::Msc2D::RegionArc>& collapsed = msc.livingRegionArcs(true);
        printArcs((std::string(mode) + "/collapsed").c_str(), collapsed);
        std::set<int> collapsedDistinct(collapsedOwners.begin(), collapsedOwners.end());
        std::cout << "region_arcs_collapsed mode=" << mode << " count=" << collapsed.size()
                  << " well_regions=" << collapsedDistinct.size() << std::endl;
        if (collapsedDistinct.size() != 1) {
            std::cerr << "FAIL[" << mode << "]: the four wells did not collapse into one "
                         "region past the last merge (" << collapsedDistinct.size()
                      << " remain)" << std::endl;
            return 2;
        }
        if (!(collapsed.size() < pairsAtFourCycle)) {
            std::cerr << "FAIL[" << mode << "]: adjacency count did not shrink from "
                      << pairsAtFourCycle << " at t=" << kFourCycleT << " to "
                      << collapsed.size() << " past the last merge" << std::endl;
            return 2;
        }

        // (3)+(4) id space and structure, both directions, several thresholds.
        const float sweep[5] = { 0.0f, 0.02f, 0.05f, 0.15f, 0.01f };
        for (int k = 0; k < 5; ++k) {
            msc.setPersistence(sweep[k]);
            for (int dir = 0; dir < 2; ++dir) {
                const bool ascending = (dir == 0);
                char buf[128];
                std::snprintf(buf, sizeof(buf), "%s/%s@%.3f", mode,
                              ascending ? "asc" : "dsc", sweep[k]);
                const std::vector<GInt::Msc2D::RegionArc>& arcs =
                    msc.livingRegionArcs(ascending);
                const std::vector<int>& b2l = msc.baseToLiving(ascending);
                std::cout << "region_arcs_sweep tag=" << buf << " pairs=" << arcs.size()
                          << " base_regions=" << b2l.size() << std::endl;
                if (!checkStructure(buf, arcs)) return 3;
                if (!checkIdSpace(buf, arcs, b2l)) return 4;
            }
        }

        // (5) determinism. The repeat call is a cache HIT; the round trip
        // through another persistence forces a genuine rebuild, which is the
        // stronger of the two.
        msc.setPersistence(kFourCycleT);
        const std::vector<GInt::Msc2D::RegionArc> first = msc.livingRegionArcs(true);
        const std::vector<GInt::Msc2D::RegionArc> repeat = msc.livingRegionArcs(true);
        if (!sameArcs(first, repeat)) {
            std::cerr << "FAIL[" << mode << "]: repeated call at one persistence differs"
                      << std::endl;
            return 5;
        }
        msc.setPersistence(0.15f);
        msc.livingRegionArcs(true);
        msc.setPersistence(kFourCycleT);
        const std::vector<GInt::Msc2D::RegionArc> rebuilt = msc.livingRegionArcs(true);
        if (!sameArcs(first, rebuilt)) {
            std::cerr << "FAIL[" << mode << "]: rebuild after a persistence round trip "
                         "differs from the first build" << std::endl;
            return 6;
        }

        // (7) descending, on the negated field: the maxima mirror the minima, so
        // the same 4-cycle must come back with the barriers negated.
        GInt::Msc2D::Msc2D mscNeg;
        mscNeg.compute(negated.data(), kDim, kDim, options);
        tag = std::string(mode) + "/dsc";
        if (!checkWellAdjacency(tag.c_str(), mscNeg, false, kFourCycleT, true)) return 7;
    }

    // (6) The two modes are compared, not gated. Their ids live in different
    // spaces (MSC node ids vs base region ids), so only the counts are
    // comparable -- and even those may legitimately differ, because the forest
    // arc set is fixed while the MSC rewires arcs as it cancels. The exact
    // 4-cycle assertion above is what pins correctness in each mode; this line
    // exists so a divergence in the noisy regime is visible rather than silent.
    for (int k = 0; k < 2; ++k) {
        std::cout << "region_arcs_mode_compare t=" << kCompareT[k];
        for (int m = 0; m < kNumConfigs; ++m) {
            std::cout << " " << kConfigs[m].name << "=" << ascCountByMode[m][k];
        }
        std::cout << std::endl;
    }

    std::cout << "region_arcs_smoke ok" << std::endl;
    return 0;
}
