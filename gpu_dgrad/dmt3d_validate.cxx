// Gate for the 3D discrete Morse theory GPU drivers (dmt3d_gradient.cu):
// ComputeLowerStarGradient3DGPU (PLS, PLS-BR) and ComputeDMFGradient3DGPU
// against the CPU raw drivers ComputeLowerStarGradient3D and
// ComputeDMFGradient3D, on random-permutation rank fields of many shapes
// (degenerate axes, slabs, cubes up to a few hundred per axis).
//
// Compared per case, memcmp-exact: gradient bytes; step and part journals
// (lower star); stage and (dimension, key) order (DMF, on cases small enough
// for the CPU order sort); and every stats field. Then a timing table:
// CPU (all OpenMP threads) vs GPU end-to-end and its h2d / kernel / d2h split.
//
// Also checked on the host for the small and medium shapes: the key-cache
// functions of dmt_core/dmt_dmf_keycache.h (used by the GPU DMF kernels)
// against the recomputing ones of dmt_core/dmt_dmf_cell.h, cell by cell.
//
// Host-only translation unit (MSVC); it may include the gi_* stack.
// Exit code: 0 iff every case matches, 1 on a mismatch or GPU failure, 2 if
// there is no CUDA device.

#include "gi_dmt_dmf_gradient.h"
#include "gi_dmt_lowerstar_gradient.h"
#include "dmt_core/dmt_cubical3d.h"
#include "dmt_core/dmt_dmf_keycache.h"
#include "dgrad_gpu_api.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

namespace hd = GInt::dmt;   // host drivers
namespace gpu = GInt::gpu;

namespace {

double MsSince(std::chrono::steady_clock::time_point t0) {
	return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

struct Case {
	int64_t X, Y, Z;
	bool sparse;   // ranks = 37 * perm + 11 instead of a permutation of 0..nv-1
};

std::vector<uint32_t> MakeRanks(const Case& c, unsigned seed) {
	const int64_t nv = c.X * c.Y * c.Z;
	std::vector<uint32_t> r((size_t)nv);
	for (int64_t i = 0; i < nv; i++) r[(size_t)i] = (uint32_t)i;
	std::mt19937 rng(seed);
	std::shuffle(r.begin(), r.end(), rng);
	if (c.sparse)
		for (auto& v : r) v = 37u * v + 11u;
	return r;
}

// Prints the first differing cell of two byte arrays; returns true if equal.
template <class T>
bool SameArrays(const char* what, const ::dmt::Cubical3D& G, const T* a, const T* b, int64_t n) {
	if (std::memcmp(a, b, sizeof(T) * (size_t)n) == 0) return true;
	int64_t diff = 0, first = -1;
	for (int64_t i = 0; i < n; i++)
		if (a[i] != b[i]) { if (first < 0) first = i; diff++; }
	int64_t q[3];
	G.Coords(first, q);
	printf("    %s: %lld cells differ, first at cell (%lld,%lld,%lld) cpu=%lld gpu=%lld\n", what, (long long)diff,
		(long long)q[0], (long long)q[1], (long long)q[2], (long long)a[first], (long long)b[first]);
	return false;
}

struct Timing {
	std::string label;
	double cpu_ms, gpu_ms;
	gpu::Dmt3DTimings split;
};

bool CheckLowerStar(const Case& c, const std::vector<uint32_t>& ranks, hd::LowerStarRule rule, bool journal,
	std::vector<Timing>* timings) {
	const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(c.X, c.Y, c.Z);
	const int64_t nc = G.NumCells();
	std::vector<uint8_t> cg((size_t)nc), cs, cp, gg((size_t)nc, 0xCD), gs, gp;
	if (journal) {
		cs.assign((size_t)nc, 0); cp.assign((size_t)nc, 0);
		gs.assign((size_t)nc, 0xCD); gp.assign((size_t)nc, 0xCD);
	}
	hd::LowerStarStats cst;
	auto t0 = std::chrono::steady_clock::now();
	const bool cok = hd::ComputeLowerStarGradient3D(ranks.data(), c.X, c.Y, c.Z, rule, cg.data(),
		journal ? cs.data() : nullptr, journal ? cp.data() : nullptr, &cst);
	const double cpu_ms = MsSince(t0);

	gpu::DmtLowerStarStats gst;
	std::memset(&gst, 0xCD, sizeof(gst));
	gpu::Dmt3DTimings tm = {};
	t0 = std::chrono::steady_clock::now();
	const bool gok = gpu::ComputeLowerStarGradient3DGPU(ranks.data(), c.X, c.Y, c.Z,
		rule == hd::LowerStarRule::PLS ? gpu::DmtLowerStarRule::PLS : gpu::DmtLowerStarRule::PLS_BR,
		gg.data(), journal ? gs.data() : nullptr, journal ? gp.data() : nullptr, &gst, &tm);
	const double gpu_ms = MsSince(t0);
	if (!gok) { printf("    GPU call failed\n"); return false; }

	bool ok = cok && cst.unassigned == 0;
	if (!ok) printf("    CPU driver reported unassigned cells\n");
	ok = SameArrays("grad", G, cg.data(), gg.data(), nc) && ok;
	if (journal) {
		ok = SameArrays("step", G, cs.data(), gs.data(), nc) && ok;
		ok = SameArrays("part", G, cp.data(), gp.data(), nc) && ok;
	}
	for (int d = 0; d < 4; d++) ok = ok && cst.critical[d] == gst.critical[d];
	ok = ok && cst.orphans == gst.orphans && cst.unassigned == gst.unassigned;
	if (!ok) {
		printf("    stats cpu crit=(%lld,%lld,%lld,%lld) orph=%lld unas=%lld | gpu crit=(%lld,%lld,%lld,%lld) orph=%lld unas=%lld\n",
			(long long)cst.critical[0], (long long)cst.critical[1], (long long)cst.critical[2], (long long)cst.critical[3],
			(long long)cst.orphans, (long long)cst.unassigned,
			(long long)gst.critical[0], (long long)gst.critical[1], (long long)gst.critical[2], (long long)gst.critical[3],
			(long long)gst.orphans, (long long)gst.unassigned);
	}
	if (timings) {
		char buf[128];
		std::snprintf(buf, sizeof(buf), "%s %lldx%lldx%lld%s", rule == hd::LowerStarRule::PLS ? "PLS" : "PLS-BR",
			(long long)c.X, (long long)c.Y, (long long)c.Z, journal ? " +journal" : "");
		timings->push_back({ buf, cpu_ms, gpu_ms, tm });
	}
	return ok;
}

bool CheckDMF(const Case& c, const std::vector<uint32_t>& ranks, bool second_pass, bool with_order,
	std::vector<Timing>* timings) {
	const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(c.X, c.Y, c.Z);
	const int64_t nc = G.NumCells();
	std::vector<uint8_t> cg((size_t)nc), cs((size_t)nc), gg((size_t)nc, 0xCD), gs((size_t)nc, 0xCD);
	std::vector<int64_t> co, go;
	if (with_order) { co.assign((size_t)nc, -1); go.assign((size_t)nc, -7); }
	hd::DMFOptions opt;
	opt.second_pass = second_pass;
	hd::DMFStats cst;
	auto t0 = std::chrono::steady_clock::now();
	hd::ComputeDMFGradient3D(ranks.data(), c.X, c.Y, c.Z, opt, cg.data(), cs.data(),
		with_order ? co.data() : nullptr, &cst);
	const double cpu_ms = MsSince(t0);

	gpu::DmtDMFStats gst;
	std::memset(&gst, 0xCD, sizeof(gst));
	gpu::Dmt3DTimings tm = {};
	t0 = std::chrono::steady_clock::now();
	const bool gok = gpu::ComputeDMFGradient3DGPU(ranks.data(), c.X, c.Y, c.Z, second_pass, gg.data(), gs.data(),
		with_order ? go.data() : nullptr, &gst, &tm);
	const double gpu_ms = MsSince(t0);
	if (!gok) { printf("    GPU call failed\n"); return false; }

	bool ok = SameArrays("grad", G, cg.data(), gg.data(), nc);
	ok = SameArrays("stage", G, cs.data(), gs.data(), nc) && ok;
	if (with_order) ok = SameArrays("order", G, co.data(), go.data(), nc) && ok;
	bool sok = cst.guard1 == gst.guard1 && cst.guard2 == gst.guard2 && cst.alg1 == gst.alg1 && cst.alg2 == gst.alg2;
	for (int d = 0; d < 4; d++) sok = sok && cst.critical[d] == gst.critical[d];
	if (!sok) {
		printf("    stats cpu g1=%lld g2=%lld a1=%lld a2=%lld crit=(%lld,%lld,%lld,%lld) | gpu g1=%lld g2=%lld a1=%lld a2=%lld crit=(%lld,%lld,%lld,%lld)\n",
			(long long)cst.guard1, (long long)cst.guard2, (long long)cst.alg1, (long long)cst.alg2,
			(long long)cst.critical[0], (long long)cst.critical[1], (long long)cst.critical[2], (long long)cst.critical[3],
			(long long)gst.guard1, (long long)gst.guard2, (long long)gst.alg1, (long long)gst.alg2,
			(long long)gst.critical[0], (long long)gst.critical[1], (long long)gst.critical[2], (long long)gst.critical[3]);
	}
	ok = ok && sok;
	if (timings) {
		char buf[128];
		std::snprintf(buf, sizeof(buf), "DMF%s %lldx%lldx%lld%s", second_pass ? "" : " (alg1 only)",
			(long long)c.X, (long long)c.Y, (long long)c.Z, with_order ? " +order" : "");
		timings->push_back({ buf, cpu_ms, gpu_ms, tm });
	}
	return ok;
}

// dmt_dmf_keycache.h vs dmt_dmf_cell.h on every cell: the top table filled
// one dimension at a time, every key and key component, and the lowest
// selecting cofacet for both selectors.
bool CheckKeyCache(const Case& c, const std::vector<uint32_t>& ranks) {
	const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(c.X, c.Y, c.Z);
	const int64_t nc = G.NumCells();
	std::vector<uint8_t> top_ref((size_t)nc), top((size_t)nc, 0x77);
	for (int64_t id = 0; id < nc; id++) {
		uint8_t g0 = 7, h2 = 7;
		::dmt::DMFTopFacets(G, ranks.data(), id, g0, h2);
		top_ref[(size_t)id] = (uint8_t)(g0 | (h2 << 4));
	}
	for (int d = 1; d <= 3; d++)
		for (int64_t id = 0; id < nc; id++) {
			if (G.Dim(id) != d) continue;
			uint8_t g0 = 7, h2 = 7;
			::dmt::DMFTopFacetsFromTop(G, ranks.data(), top.data(), id, g0, h2);
			top[(size_t)id] = (uint8_t)(g0 | (h2 << 4));
		}
	if (!SameArrays("keycache top", G, top_ref.data(), top.data(), nc)) return false;
	for (int64_t id = 0; id < nc; id++) {
		uint32_t k[8], kc[8];
		const int len = ::dmt::DMFKey(G, ranks.data(), id, k);
		const int lenc = ::dmt::DMFKeyFromTop(G, ranks.data(), top.data(), id, kc);
		bool same = len == lenc && ::dmt::DMFKeyCmp(k, kc, len) == 0;
		for (int i = 0; i < len && same; i++) same = ::dmt::DMFKeyComponentFromTop(G, ranks.data(), top.data(), id, i) == k[i];
		for (int w = 0; w < 2 && same; w++)
			same = ::dmt::DMFLowestSelectingCofacet(G, ranks.data(), top.data(), id, w) ==
				::dmt::DMFLowestSelectingCofacetFromTop(G, ranks.data(), top.data(), id, w);
		if (!same) {
			int64_t q[3];
			G.Coords(id, q);
			printf("    keycache: cell (%lld,%lld,%lld) differs\n", (long long)q[0], (long long)q[1], (long long)q[2]);
			return false;
		}
	}
	return true;
}

} // namespace

int main() {
	gpu::DeviceInfo info;
	if (!gpu::QueryDevice(info)) {
		printf("dmt3d_validate: no CUDA device\n");
		return 2;
	}
	printf("dmt3d_validate on %s (sm_%d%d, %d SMs), CPU threads %d\n", info.name, info.cc_major, info.cc_minor,
		info.sm_count, hd::ResolveThreads(0));

	std::vector<Case> cases = {
		{ 1, 1, 1, false }, { 2, 1, 1, false }, { 1, 1, 7, false }, { 2, 2, 2, false }, { 3, 3, 4, false },
		{ 1, 9, 13, false }, { 11, 1, 6, false }, { 33, 40, 1, false }, { 5, 6, 7, false }, { 9, 8, 7, true },
		{ 16, 17, 18, false }, { 31, 2, 29, false }, { 64, 48, 40, false }, { 64, 64, 64, true },
		{ 129, 97, 65, false },
	};
	std::mt19937 pick(2026u);
	for (int k = 0; k < 24; k++)
		cases.push_back({ 1 + (int64_t)(pick() % 14), 1 + (int64_t)(pick() % 14), 1 + (int64_t)(pick() % 14), (k % 5) == 4 });
	// Large cases: correctness and the timing table.
	const std::vector<Case> big = { { 200, 200, 200, false }, { 256, 256, 256, false }, { 384, 256, 192, false } };
	const int64_t kOrderCellLimit = (int64_t)1 << 24;   // CPU order sort holds 32 B of keys per cell

	// Warm-up: context creation and module load outside every timing.
	{
		const Case w = { 8, 8, 8, false };
		std::vector<uint32_t> r = MakeRanks(w, 1u);
		std::vector<uint8_t> g(15 * 15 * 15);
		gpu::ComputeLowerStarGradient3DGPU(r.data(), 8, 8, 8, gpu::DmtLowerStarRule::PLS, g.data(), nullptr, nullptr,
			nullptr, nullptr);
		gpu::ComputeDMFGradient3DGPU(r.data(), 8, 8, 8, true, g.data(), nullptr, nullptr, nullptr, nullptr);
	}

	int runs = 0, bad = 0;
	auto report = [&](bool ok, const char* what, const Case& c) {
		runs++;
		if (!ok) {
			bad++;
			printf("  FAIL %s %lldx%lldx%lld%s\n", what, (long long)c.X, (long long)c.Y, (long long)c.Z,
				c.sparse ? " (sparse ranks)" : "");
		}
	};

	for (size_t i = 0; i < cases.size(); i++) {
		const Case& c = cases[i];
		const std::vector<uint32_t> ranks = MakeRanks(c, 1000u + (unsigned)i);
		const int64_t nc = (2 * c.X - 1) * (2 * c.Y - 1) * (2 * c.Z - 1);
		const bool order = nc <= kOrderCellLimit;
		report(CheckLowerStar(c, ranks, hd::LowerStarRule::PLS, true, nullptr), "PLS", c);
		report(CheckLowerStar(c, ranks, hd::LowerStarRule::PLS_BR, true, nullptr), "PLS-BR", c);
		report(CheckLowerStar(c, ranks, hd::LowerStarRule::PLS, false, nullptr), "PLS (no journal)", c);
		report(CheckDMF(c, ranks, true, order, nullptr), "DMF", c);
		report(CheckDMF(c, ranks, false, order, nullptr), "DMF alg1", c);
		if (nc <= ((int64_t)1 << 21)) report(CheckKeyCache(c, ranks), "DMF key cache (host)", c);
	}
	printf("small/medium: %zu shapes, %d runs, %d mismatches\n", cases.size(), runs, bad);

	std::vector<Timing> timings;
	for (size_t i = 0; i < big.size(); i++) {
		const Case& c = big[i];
		const std::vector<uint32_t> ranks = MakeRanks(c, 5000u + (unsigned)i);
		report(CheckLowerStar(c, ranks, hd::LowerStarRule::PLS, false, &timings), "PLS", c);
		report(CheckLowerStar(c, ranks, hd::LowerStarRule::PLS_BR, false, &timings), "PLS-BR", c);
		report(CheckLowerStar(c, ranks, hd::LowerStarRule::PLS_BR, true, &timings), "PLS-BR journal", c);
		report(CheckDMF(c, ranks, true, false, &timings), "DMF", c);
		if (i == 0) report(CheckDMF(c, ranks, true, true, &timings), "DMF order", c);
	}

	printf("\n%-34s %10s %10s %8s %8s %8s %8s\n", "case", "CPU ms", "GPU ms", "speedup", "h2d", "kernel", "d2h");
	for (const Timing& t : timings)
		printf("%-34s %10.1f %10.1f %7.1fx %8.2f %8.2f %8.2f\n", t.label.c_str(), t.cpu_ms, t.gpu_ms,
			t.cpu_ms / (t.gpu_ms > 0 ? t.gpu_ms : 1e-9), t.split.h2d_ms, t.split.kernel_ms, t.split.d2h_ms);

	printf("\ndmt3d_validate: %d runs, %d mismatches -> %s\n", runs, bad, bad ? "FAIL" : "ok");
	return bad ? 1 : 0;
}
