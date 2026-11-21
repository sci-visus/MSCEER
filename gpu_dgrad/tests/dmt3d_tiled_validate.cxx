// Gate for the tiled GPU DMT drivers (dmt3d_tiled.cu): ComputeDMFGradient3DTiled
// (multi-pass DMF, with and without the second pass) and
// ComputeLowerStarGradient3DTiled (PLS, PLS-BR).
//
// Gates (exit code 0 iff all pass):
//   1. tiled == in-core GPU driver (ComputeDMFGradient3DGPU /
//      ComputeLowerStarGradient3DGPU) == CPU raw driver (on volumes up to 2^22
//      cells), memcmp on the gradient (and the DMF stage), equal stats, for
//      degenerate, odd and cubic volumes and many tile shapes: tiles smaller
//      than the halo, anisotropic tiles, slabs, tiles touching the boundary.
//   2. halo sweep: tiles staged with halo (lo, hi) for lo, hi in 0..2 on
//      random 22^3 fields; the default halos (DMF 2/1, the others 1/1) must
//      match everywhere, and the report shows which smaller halos fail
//      (each tile is a truncated complex of its own, so this measures the
//      halo of the actual blocking, not only of the value dependence).
//   3. fixed footprint: DMF on procedural 512^3 and 1024^3 ranks streamed
//      tile by tile (no host volume or gradient): every cell written, Euler
//      characteristic 1, 512^3 hash == in-core hash, device free memory flat,
//      0-byte device leak.
// Usage: dmt3d_tiled_validate [--quick]   (--quick skips gate 3 and the big volumes)

#include "gi_dmt_dmf_gradient.h"
#include "gi_dmt_lowerstar_gradient.h"
#include "dmt_core/dmt_cubical3d.h"
#include "dmt_core/dmt_pls.h"
#include "dmt_core/dmt_pls_br.h"
#include "dmt_core/dmt_star27_tables.h"
#include "dgrad_gpu_api.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

namespace hd = GInt::dmt;
namespace gpu = GInt::gpu;

namespace {

double MsSince(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

std::vector<uint32_t> MakeRanks(int64_t X, int64_t Y, int64_t Z, unsigned seed, bool sparse) {
    const int64_t nv = X * Y * Z;
    std::vector<uint32_t> r((size_t)nv);
    for (int64_t i = 0; i < nv; i++) r[(size_t)i] = (uint32_t)i;
    std::mt19937 rng(seed);
    std::shuffle(r.begin(), r.end(), rng);
    if (sparse)
        for (auto& v : r) v = 37u * v + 11u;
    return r;
}

enum class R { PLS, PLSBR, DMF2, DMF1 };
const char* Name(R r) {
    switch (r) {
        case R::PLS: return "PLS";
        case R::PLSBR: return "PLS-BR";
        case R::DMF2: return "DMF";
        default: return "DMF(alg1)";
    }
}

struct Out {
    std::vector<uint8_t> grad, stage;
    int64_t s[8] = {0, 0, 0, 0, 0, 0, 0, 0};   // crit0..3, then rule-specific
    bool ok = false;
};

bool SameStats(const Out& a, const Out& b) {
    for (int i = 0; i < 8; i++)
        if (a.s[i] != b.s[i]) return false;
    return true;
}

Out RunInCore(R r, const std::vector<uint32_t>& ranks, int64_t X, int64_t Y, int64_t Z) {
    Out o;
    const int64_t nc = (2 * X - 1) * (2 * Y - 1) * (2 * Z - 1);
    o.grad.assign((size_t)nc, 0xCD);
    if (r == R::DMF2 || r == R::DMF1) {
        o.stage.assign((size_t)nc, 0xCD);
        gpu::DmtDMFStats st;
        o.ok = gpu::ComputeDMFGradient3DGPU(ranks.data(), X, Y, Z, r == R::DMF2, o.grad.data(), o.stage.data(),
                                            nullptr, &st, nullptr);
        for (int d = 0; d < 4; d++) o.s[d] = st.critical[d];
        o.s[4] = st.guard1; o.s[5] = st.guard2; o.s[6] = st.alg1; o.s[7] = st.alg2;
    } else {
        gpu::DmtLowerStarStats st;
        o.ok = gpu::ComputeLowerStarGradient3DGPU(ranks.data(), X, Y, Z,
                                                  r == R::PLS ? gpu::DmtLowerStarRule::PLS : gpu::DmtLowerStarRule::PLS_BR,
                                                  o.grad.data(), nullptr, nullptr, &st, nullptr);
        for (int d = 0; d < 4; d++) o.s[d] = st.critical[d];
        o.s[4] = st.orphans; o.s[5] = st.unassigned;
    }
    return o;
}

Out RunCpu(R r, const std::vector<uint32_t>& ranks, int64_t X, int64_t Y, int64_t Z) {
    Out o;
    const int64_t nc = (2 * X - 1) * (2 * Y - 1) * (2 * Z - 1);
    o.grad.assign((size_t)nc, 0xCD);
    if (r == R::DMF2 || r == R::DMF1) {
        o.stage.assign((size_t)nc, 0xCD);
        hd::DMFOptions opt;
        opt.second_pass = r == R::DMF2;
        hd::DMFStats st;
        hd::ComputeDMFGradient3D(ranks.data(), X, Y, Z, opt, o.grad.data(), o.stage.data(), nullptr, &st);
        for (int d = 0; d < 4; d++) o.s[d] = st.critical[d];
        o.s[4] = st.guard1; o.s[5] = st.guard2; o.s[6] = st.alg1; o.s[7] = st.alg2;
    } else {
        hd::LowerStarStats st;
        hd::ComputeLowerStarGradient3D(ranks.data(), X, Y, Z, r == R::PLS ? hd::LowerStarRule::PLS : hd::LowerStarRule::PLS_BR,
                                       o.grad.data(), nullptr, nullptr, &st);
        for (int d = 0; d < 4; d++) o.s[d] = st.critical[d];
        o.s[4] = st.orphans; o.s[5] = st.unassigned;
    }
    o.ok = true;
    return o;
}

Out RunTiled(R r, const std::vector<uint32_t>& ranks, int64_t X, int64_t Y, int64_t Z, int Tx, int Ty, int Tz,
             int slots, int lo, int hi, gpu::TileStats3D* ts = nullptr) {
    Out o;
    const int64_t nc = (2 * X - 1) * (2 * Y - 1) * (2 * Z - 1);
    o.grad.assign((size_t)nc, 0xEE);
    const gpu::TileGrid3D tg = gpu::MakeTileGrid3D(X, Y, Z, Tx, Ty, Tz);
    gpu::DmtTiledOptions3D opt;
    opt.slots = slots;
    opt.haloLo = lo;
    opt.haloHi = hi;
    const auto src = gpu::HostRankSource3D(ranks.data(), X, Y);
    const auto gs = gpu::HostGradSink3D(o.grad.data(), X, Y);
    if (r == R::DMF2 || r == R::DMF1) {
        o.stage.assign((size_t)nc, 0xEE);
        const auto ss = gpu::HostGradSink3D(o.stage.data(), X, Y);
        gpu::DmtDMFStats st;
        o.ok = gpu::ComputeDMFGradient3DTiled(src, tg, r == R::DMF2, gs, &ss, opt, &st, ts);
        for (int d = 0; d < 4; d++) o.s[d] = st.critical[d];
        o.s[4] = st.guard1; o.s[5] = st.guard2; o.s[6] = st.alg1; o.s[7] = st.alg2;
    } else {
        gpu::DmtLowerStarStats st;
        o.ok = gpu::ComputeLowerStarGradient3DTiled(src, tg, r == R::PLS ? gpu::DmtLowerStarRule::PLS : gpu::DmtLowerStarRule::PLS_BR,
                                                    gs, opt, &st, ts);
        for (int d = 0; d < 4; d++) o.s[d] = st.critical[d];
        o.s[4] = st.orphans; o.s[5] = st.unassigned;
    }
    return o;
}

// Streamed entry points (dgrad_gpu_api.h): mode 2 = direct slab DMA with the
// arrays pinned for the call, mode 1 = copy through the tiled ring.
Out RunStreamed(R r, std::vector<uint32_t>& ranks, int64_t X, int64_t Y, int64_t Z, int mode, int tx, int ty,
                int tz, int slots) {
    Out o;
    const int64_t nc = (2 * X - 1) * (2 * Y - 1) * (2 * Z - 1);
    o.grad.assign((size_t)nc, 0xEE);
    gpu::DmtStreamOptions3D opt;
    opt.hostMode = mode;
    opt.tileX = tx;
    opt.tileY = ty;
    opt.slabDepth = tz;
    opt.slots = slots;
    if (r == R::DMF2 || r == R::DMF1) {
        o.stage.assign((size_t)nc, 0xEE);
        gpu::DmtDMFStats st;
        o.ok = gpu::ComputeDMFGradient3DStreamed(ranks.data(), X, Y, Z, r == R::DMF2, o.grad.data(), o.stage.data(),
                                                 opt, &st, nullptr);
        for (int d = 0; d < 4; d++) o.s[d] = st.critical[d];
        o.s[4] = st.guard1; o.s[5] = st.guard2; o.s[6] = st.alg1; o.s[7] = st.alg2;
    } else {
        gpu::DmtLowerStarStats st;
        o.ok = gpu::ComputeLowerStarGradient3DStreamed(ranks.data(), X, Y, Z,
                                                       r == R::PLS ? gpu::DmtLowerStarRule::PLS : gpu::DmtLowerStarRule::PLS_BR,
                                                       o.grad.data(), opt, &st, nullptr);
        for (int d = 0; d < 4; d++) o.s[d] = st.critical[d];
        o.s[4] = st.orphans; o.s[5] = st.unassigned;
    }
    return o;
}

int64_t CountDiff(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
    if (a.size() != b.size()) return -1;
    int64_t n = 0;
    for (size_t i = 0; i < a.size(); i++) n += a[i] != b[i];
    return n;
}

// grad + stage + stats; returns mismatching cells (or 1 for a stats/run failure)
int64_t Compare(const Out& ref, const Out& t) {
    if (!t.ok) return 1;
    int64_t d = CountDiff(ref.grad, t.grad);
    if (!ref.stage.empty()) d += CountDiff(ref.stage, t.stage);
    if (d == 0 && !SameStats(ref, t)) d = 1;
    return d;
}

struct Cfg { int Tx, Ty, Tz, slots; };

std::vector<Cfg> Configs(int64_t X, int64_t Y, int64_t Z) {
    const int64_t nv = X * Y * Z;
    std::vector<Cfg> c;
    if (nv <= 5000) { c.push_back({1, 1, 1, 3}); c.push_back({2, 2, 2, 2}); }
    if (nv <= 60000) { c.push_back({3, 3, 3, 3}); c.push_back({5, 3, 2, 1}); }
    if (nv <= 300000) { c.push_back({7, 7, 7, 3}); c.push_back({16, 16, 16, 2}); c.push_back({(int)X, 4, 3, 3}); }
    c.push_back({(int)X, (int)Y, 4, 3});
    c.push_back({32, 32, 32, 3});
    if (nv > 300000) { c.push_back({64, 64, 64, 3}); c.push_back({100, 37, 61, 2}); c.push_back({128, 128, 16, 3}); }
    return c;
}

// ---------------------------------------------------------------------------
// fixed footprint
// ---------------------------------------------------------------------------

inline uint32_t Fmix32(uint32_t h) {   // a bijection on uint32: distinct ranks
    h ^= h >> 16; h *= 0x85ebca6bu; h ^= h >> 13; h *= 0xc2b2ae35u; h ^= h >> 16;
    return h;
}

inline uint64_t Mix(uint64_t i, uint8_t b) {
    uint64_t h = (i * 0x9E3779B97F4A7C15ull) ^ ((uint64_t)b * 0xC2B2AE3D27D4EB4Full);
    h ^= h >> 31; h *= 0xBF58476D1CE4E5B9ull; h ^= h >> 29;
    return h;
}

struct FootResult {
    bool ok = false;
    uint64_t hash = 0;
    int64_t cells = 0;
    gpu::TileStats3D ts;
    gpu::DmtDMFStats st;
    int64_t freeBefore = 0, freeAfter = 0;
};

FootResult RunFootprint(int64_t N, int Tx, int Ty, int Tz) {
    FootResult r;
    printf("\n=== DMF procedural %lld^3 ranks, T=(%d,%d,%d), 3 slots (no host volume, no host gradient) ===\n",
           (long long)N, Tx, Ty, Tz);
    gpu::QueryDeviceMemory(&r.freeBefore, nullptr);
    gpu::RankSource3D src = [N](const gpu::Box3D& b, uint32_t* dst, int64_t pitch, int64_t spitch) {
#pragma omp parallel for schedule(static)
        for (int z = 0; z < b.d; z++)
            for (int y = 0; y < b.h; y++) {
                uint32_t* d = dst + z * spitch + y * pitch;
                const int64_t base = ((b.z0 + z) * N + (b.y0 + y)) * N + b.x0;
                for (int x = 0; x < b.w; x++) d[x] = Fmix32((uint32_t)(base + x));
            }
        return true;
    };
    const int64_t W = 2 * N - 1;
    gpu::ByteSink3D sink = [&](const gpu::Box3D& b, const uint8_t* s, int64_t pitch, int64_t spitch) {
        uint64_t hs = 0;
#pragma omp parallel for schedule(static) reduction(+ : hs)
        for (int z = 0; z < b.d; z++)
            for (int y = 0; y < b.h; y++) {
                const uint8_t* row = s + z * spitch + y * pitch;
                const int64_t base = ((b.z0 + z) * W + (b.y0 + y)) * W + b.x0;
                for (int x = 0; x < b.w; x++) hs += Mix((uint64_t)(base + x), row[x]);
            }
        r.hash += hs;
        r.cells += b.Count();
        return true;
    };
    const gpu::TileGrid3D tg = gpu::MakeTileGrid3D(N, N, N, Tx, Ty, Tz);
    gpu::DmtTiledOptions3D opt;
    opt.slots = 3;
    r.ok = gpu::ComputeDMFGradient3DTiled(src, tg, true, sink, nullptr, opt, &r.st, &r.ts);
    gpu::QueryDeviceMemory(&r.freeAfter, nullptr);
    const int64_t euler = r.st.critical[0] - r.st.critical[1] + r.st.critical[2] - r.st.critical[3];
    printf("  ok=%d cells=%lld (expected %lld) crit=%lld/%lld/%lld/%lld euler=%lld guards %lld/%lld claims %lld/%lld\n",
           (int)r.ok, (long long)r.cells, (long long)(W * W * W), (long long)r.st.critical[0],
           (long long)r.st.critical[1], (long long)r.st.critical[2], (long long)r.st.critical[3], (long long)euler,
           (long long)r.st.guard1, (long long)r.st.guard2, (long long)r.st.alg1, (long long)r.st.alg2);
    printf("  tiles=%lld wall %.1f ms (kernel %.1f, h2d %.1f, d2h %.1f summed; source %.1f, sink %.1f host ms)\n",
           (long long)r.ts.tiles, r.ts.wall_ms, r.ts.kernel_ms, r.ts.h2d_ms, r.ts.d2h_ms, r.ts.source_ms, r.ts.sink_ms);
    printf("  device: allocated %.1f MB, free min/max %.1f / %.1f MB, before/after %.1f / %.1f MB (delta %lld B)\n",
           r.ts.deviceBytesAllocated / 1e6, r.ts.deviceFreeMin / 1e6, r.ts.deviceFreeMax / 1e6, r.freeBefore / 1e6,
           r.freeAfter / 1e6, (long long)(r.freeAfter - r.freeBefore));
    r.ok = r.ok && r.cells == W * W * W && euler == 1;
    printf("  -> %s\n", r.ok ? "PASS" : "FAIL");
    return r;
}

// Gate 0: the lower-star core with rank-mask keys (PLSStar / PLSBRStar,
// StarResult recorder; PLSStarCodes / PLSBRStarCodes, pair-code recorder) against the
// original packed 64-bit keys (PLSStarPacked / PLSBRStarPacked), star by star:
// partner, step, part and orphan count, and the pair codes. Random stars with
// distinct keys, every clipping pattern of the 27 slots (domain faces, edges,
// corners and degenerate axes) and random boundary masks, plus every vertex
// of a few volumes. Returns the number of differing stars.
int64_t CheckStarCore(int64_t nRandom) {
    dmt::Star27Tables T;
    dmt::BuildStar27TablesCanonical(T);
    std::mt19937 rng(424242);
    int64_t bad = 0, n = 0;
    auto one = [&](const uint32_t key[27], uint32_t vvalid, uint32_t bnd) {
        for (int br = 0; br < 2; br++) {
            dmt::StarResult a, b;
            if (br) { dmt::PLSBRStar(T, key, vvalid, bnd, a); dmt::PLSBRStarPacked(T, key, vvalid, bnd, b); }
            else { dmt::PLSStar(T, key, vvalid, a); dmt::PLSStarPacked(T, key, vvalid, b); }
            bool same = a.n_orphans == b.n_orphans;
            for (int q = 0; q < 27 && same; q++)
                same = a.pair[q] == b.pair[q] && a.step[q] == b.step[q] && a.part[q] == b.part[q];
            uint32_t L;
            dmt::StarCodeRecorder rec;
            if (br) dmt::PLSBRStarCodes(T, key, vvalid, bnd, L, rec);
            else dmt::PLSStarCodes(T, key, vvalid, L, rec);
            for (int q = 0; q < 27 && same; q++) {
                const bool inL = ((rec.done >> q) & 1u) != 0 && ((L >> q) & 1u) != 0;
                if (inL != (b.pair[q] != -2)) { same = false; break; }
                if (!inL) continue;
                const int code = b.pair[q] < 0 ? dmt::kPairCritical : dmt::StarDir(T, q, b.pair[q]);
                same = (int)rec.Code(q) == code;
            }
            same = same && rec.n_orphans == b.n_orphans;
            if (!same) bad++;
            n++;
        }
    };
    // clip patterns: per axis the vertex offsets present, {-1,0,1} minus a side
    const int clipLo[3] = {-1, 0, -1}, clipHi[3] = {1, 1, 0};
    std::uniform_int_distribution<uint32_t> any;
    for (int64_t it = 0; it < nRandom; it++) {
        const int cx = (int)(it % 4), cy = (int)((it / 4) % 4), cz = (int)((it / 16) % 4);
        // pattern 3 = degenerate axis (only offset 0)
        auto lo = [&](int c) { return c == 3 ? 0 : clipLo[c]; };
        auto hi = [&](int c) { return c == 3 ? 0 : clipHi[c]; };
        uint32_t key[27], vvalid = 0;
        std::vector<uint32_t> vals(27);
        for (int q = 0; q < 27; q++) vals[q] = (uint32_t)q;
        std::shuffle(vals.begin(), vals.end(), rng);
        const bool spread = (it & 1) != 0;
        for (int q = 0; q < 27; q++) {
            const int dx = q % 3 - 1, dy = (q / 3) % 3 - 1, dz = q / 9 - 1;
            const bool in = dx >= lo(cx) && dx <= hi(cx) && dy >= lo(cy) && dy <= hi(cy) && dz >= lo(cz) && dz <= hi(cz);
            key[q] = in ? (spread ? vals[q] * 1000003u + (any(rng) % 1000u) : vals[q]) : 0u;
            if (in) vvalid |= 1u << q;
        }
        // keep the centre high now and then so the lower stars are large
        if (it % 3 == 0) key[13] = 0xFFFFFFF0u;
        const uint32_t bnd = any(rng) & ((1u << 27) - 1);
        one(key, vvalid, bnd);
    }
    // every vertex of a few volumes (true boundary masks)
    const int64_t dims[4][3] = {{9, 8, 7}, {1, 6, 5}, {12, 1, 1}, {16, 16, 16}};
    for (int v = 0; v < 4; v++) {
        const int64_t X = dims[v][0], Y = dims[v][1], Z = dims[v][2];
        const auto ranks = MakeRanks(X, Y, Z, 77u + (unsigned)v, v == 1);
        const dmt::Cubical3D G = dmt::Cubical3D::Make(X, Y, Z);
        for (int64_t z = 0; z < Z; z++)
            for (int64_t y = 0; y < Y; y++)
                for (int64_t x = 0; x < X; x++) {
                    uint32_t key[27], vvalid;
                    dmt::GatherStarKeys(G, ranks.data(), x, y, z, key, vvalid);
                    one(key, vvalid, dmt::StarBoundaryMask(T, G, x, y, z));
                }
    }
    printf("\n=== lower-star core: rank-mask keys vs packed keys, %lld stars x {PLS, PLS-BR}: %lld differ -> %s ===\n",
           (long long)(n / 2), (long long)bad, bad ? "FAIL" : "PASS");
    return bad;
}

} // namespace

int main(int argc, char** argv) {
    const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
    gpu::DeviceInfo di;
    if (!gpu::QueryDevice(di)) { printf("no CUDA device\n"); return 2; }
    printf("device: %s\n", di.name);

    struct Shape { int64_t X, Y, Z; bool sparse; };
    std::vector<Shape> shapes = {
        {1, 1, 1, false}, {2, 2, 2, false}, {3, 3, 3, true}, {1, 5, 7, false}, {6, 1, 4, false},
        {5, 6, 1, true}, {2, 9, 3, false}, {9, 2, 2, false}, {17, 13, 11, false}, {24, 20, 18, true},
        {33, 7, 29, false}, {40, 41, 39, false},
    };
    if (!quick) { shapes.push_back({96, 80, 72, false}); shapes.push_back({200, 130, 90, true}); }

    int64_t bad = 0, runs = 0;
    bad += CheckStarCore(quick ? 200000 : 2000000);
    const R rules[] = {R::DMF2, R::DMF1, R::PLS, R::PLSBR};
    unsigned seed = 9001;
    for (const Shape& sh : shapes) {
        const auto ranks = MakeRanks(sh.X, sh.Y, sh.Z, seed++, sh.sparse);
        const int64_t nc = (2 * sh.X - 1) * (2 * sh.Y - 1) * (2 * sh.Z - 1);
        printf("\n=== %lldx%lldx%lld%s ===\n", (long long)sh.X, (long long)sh.Y, (long long)sh.Z, sh.sparse ? " sparse" : "");
        for (R r : rules) {
            const Out inc = RunInCore(r, ranks, sh.X, sh.Y, sh.Z);
            if (!inc.ok) { printf("  %s in-core GPU failed\n", Name(r)); bad++; continue; }
            int64_t dc = 0;
            if (nc <= ((int64_t)1 << 22)) {
                const Out cpu = RunCpu(r, ranks, sh.X, sh.Y, sh.Z);
                dc = Compare(cpu, inc);
            }
            int64_t worst = 0;
            int ncfg = 0;
            for (const Cfg& c : Configs(sh.X, sh.Y, sh.Z)) {
                const Out t = RunTiled(r, ranks, sh.X, sh.Y, sh.Z, c.Tx, c.Ty, c.Tz, c.slots, -1, -1);
                const int64_t d = Compare(inc, t);
                if (d) {
                    printf("    %s tiled T=(%d,%d,%d) slots=%d: %lld mismatching cells/stats\n", Name(r), c.Tx, c.Ty,
                           c.Tz, c.slots, (long long)d);
                    bad++;
                }
                worst = std::max(worst, d);
                ncfg++;
                runs++;
            }
            if (dc) bad++;
            // streamed: direct (pinned for the call) and copy, tiny and default tiles / slabs
            int64_t sworst = 0;
            {
                auto& rk = const_cast<std::vector<uint32_t>&>(ranks);
                const int cfgs[4][5] = {{2, 3, 2, 2, 2}, {2, 128, 128, 64, 2}, {1, 5, 4, 3, 3}, {2, 7, 9, 1, 3}};
                for (const auto& c : cfgs) {
                    const Out t = RunStreamed(r, rk, sh.X, sh.Y, sh.Z, c[0], c[1], c[2], c[3], c[4]);
                    const int64_t d = Compare(inc, t);
                    if (d) {
                        printf("    %s streamed mode=%d tile=(%d,%d) slab=%d slots=%d: %lld mismatching cells/stats\n",
                               Name(r), c[0], c[1], c[2], c[3], c[4], (long long)d);
                        bad++;
                    }
                    sworst = std::max(sworst, d);
                    runs++;
                }
            }
            printf("  %-9s in-core vs CPU: %s   tiled vs in-core (%d tile shapes): %s   streamed (4): %s\n", Name(r),
                   nc <= ((int64_t)1 << 22) ? (dc ? "FAIL" : "PASS") : "skipped", ncfg, worst ? "FAIL" : "PASS",
                   sworst ? "FAIL" : "PASS");
        }
    }

    // gate 2: halo sweep
    printf("\n=== halo sweep: 22^3 random fields, T=6 (27 tiles, 8 interior-ish), staged [o-lo, o+t+1+hi) ===\n");
    printf("  runs with any mismatch out of 24 fields (mismatching cells summed):\n");
    bool sweepOk = true;
    for (R r : rules) {
        printf("  %-9s", Name(r));
        for (int lo = 0; lo <= 2; lo++)
            for (int hi = 0; hi <= 2; hi++) {
                int failRuns = 0;
                int64_t cells = 0;
                for (int f = 0; f < 24; f++) {
                    const auto ranks = MakeRanks(22, 22, 22, 5000u + (unsigned)f, false);
                    const Out inc = RunInCore(r, ranks, 22, 22, 22);
                    const Out t = RunTiled(r, ranks, 22, 22, 22, 6, 6, 6, 3, lo, hi);
                    const int64_t d = Compare(inc, t);
                    if (d) { failRuns++; cells += d; }
                }
                const int defLo = r == R::DMF2 ? 2 : 1;
                const bool atDefault = lo >= defLo && hi >= 1;
                if (atDefault && failRuns) sweepOk = false;
                printf("  (%d,%d): %2d/24 %6lld%s", lo, hi, failRuns, (long long)cells, atDefault ? "*" : " ");
            }
        printf("\n");
    }
    printf("  (* = at or above the rule's default halo; must be 0)  -> %s\n", sweepOk ? "PASS" : "FAIL");
    if (!sweepOk) bad++;

    bool footOk = true;
    if (!quick) {
        const int64_t N0 = 512;
        uint64_t refHash = 0;
        {
            std::vector<uint32_t> ranks((size_t)(N0 * N0 * N0));
#pragma omp parallel for schedule(static)
            for (long long i = 0; i < (long long)ranks.size(); i++) ranks[(size_t)i] = Fmix32((uint32_t)i);
            const int64_t nc = (2 * N0 - 1) * (2 * N0 - 1) * (2 * N0 - 1);
            std::vector<uint8_t> g((size_t)nc);
            gpu::DmtDMFStats st;
            gpu::Dmt3DTimings tm{};
            if (!gpu::ComputeDMFGradient3DGPU(ranks.data(), N0, N0, N0, true, g.data(), nullptr, nullptr, &st, &tm)) {
                printf("in-core 512^3 DMF failed\n");
                footOk = false;
            }
            uint64_t hs = 0;
#pragma omp parallel for schedule(static) reduction(+ : hs)
            for (long long i = 0; i < (long long)nc; i++) hs += Mix((uint64_t)i, g[(size_t)i]);
            refHash = hs;
            printf("\nin-core GPU DMF 512^3: h2d %.1f, kernel %.1f, d2h %.1f ms\n", tm.h2d_ms, tm.kernel_ms, tm.d2h_ms);
        }
        FootResult a = RunFootprint(512, 128, 128, 64);
        FootResult b = RunFootprint(1024, 128, 128, 64);
        const bool hashOk = a.hash == refHash;
        const int64_t fd = a.ts.deviceFreeMin - b.ts.deviceFreeMin;
        const int64_t la = a.freeAfter - a.freeBefore, lb = b.freeAfter - b.freeBefore;
        printf("\n  tiled 512^3 hash == in-core hash: %s; device free-min 512^3 vs 1024^3 differ by %lld B; leak %lld / %lld B\n",
               hashOk ? "PASS" : "FAIL", (long long)fd, (long long)la, (long long)lb);
        footOk = footOk && a.ok && b.ok && hashOk && std::llabs(fd) <= (64ll << 20) && la == 0 && lb == 0;
        printf("  footprint gate -> %s\n", footOk ? "PASS" : "FAIL");
    }

    const bool pass = bad == 0 && footOk;
    printf("\n=== %lld tiled runs, failures %lld, footprint %s -> %s ===\n", (long long)runs, (long long)bad,
           quick ? "skipped" : (footOk ? "ok" : "FAIL"), pass ? "ALL PASS" : "FAIL");
    return pass ? 0 : 1;
}
