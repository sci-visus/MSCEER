// Alternating-cycle flips and the hunts built on them, on random fields.
//
//   - every 4-cycle move and every random alternating cycle (k = 2..6) that
//     CycleOk accepts keeps the critical set, and undoing it restores the
//     gradient bytes exactly;
//   - the local acyclicity gate is sound: whenever CreatesCycleAt says no
//     closed V-path was created, the full IsAcyclic agrees (and vice versa);
//   - every k = 2 cycle is one of the listed 4-cycle moves;
//   - Hunt, HuntLocal, HuntCycles and HuntAll leave an acyclic gradient with
//     the same critical cells, d~ o d~ = 0 and Betti [1, 0, 0, 0].
//
// Exit code: 0 iff every check passes.

#include "gi_dmt_acyclic.h"
#include "gi_dmt_alternating_cycles.h"
#include "gi_dmt_cycle_hunts.h"
#include "gi_dmt_lowerstar_gradient.h"
#include "gi_dmt_morse_complex.h"

#include <cstdio>
#include <cstring>
#include <random>
#include <set>
#include <tuple>
#include <vector>

static int g_fail = 0;
static int g_checks = 0;
#define CHECK(cond, ...) do { g_checks++; if (!(cond)) { g_fail++; if (g_fail <= 20) { printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while (0)

namespace gd = GInt::dmt;

static std::vector<uint8_t> CritMask(const gd::GradientView3D& V) {
	std::vector<uint8_t> m((size_t)V.G.NumCells());
	for (int64_t c = 0; c < V.G.NumCells(); c++) m[(size_t)c] = V.IsCritical(c);
	return m;
}

static void ValidMorse(const char* what, gd::GradientView3D& V, const std::vector<uint8_t>& crit0) {
	CHECK(gd::IsAcyclic(V), "%s: closed V-path", what);
	CHECK(gd::CheckConsistency(V) == 0, "%s: inconsistent matching", what);
	CHECK(CritMask(V) == crit0, "%s: critical set changed", what);
	const gd::MorseComplexCheck R = gd::CheckMorseComplex(V, 3, 1);
	CHECK(R.bad.empty() && R.exact, "%s: d~ o d~ != 0", what);
	const auto b = gd::BettiModP(V);
	CHECK(b[0] == 1 && b[1] == 0 && b[2] == 0 && b[3] == 0, "%s: Betti %lld %lld %lld %lld", what,
		(long long)b[0], (long long)b[1], (long long)b[2], (long long)b[3]);
}

int main(int argc, char** argv) {
	const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
	printf("dmt_cycles_smoke%s\n", quick ? " --quick" : "");
	std::mt19937 pick(77u);
	int64_t ncycles = 0, nmoves = 0;
	const int nfields = quick ? 12 : 60;
	for (int i = 0; i < nfields; i++) {
		const int X = 3 + (int)(pick() % 6), Y = 3 + (int)(pick() % 6), Z = 3 + (int)(pick() % 6);
		const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(X, Y, Z);
		std::vector<uint32_t> ranks((size_t)G.NumVertices());
		for (size_t k = 0; k < ranks.size(); k++) ranks[k] = (uint32_t)k;
		std::mt19937 rng(500u + (unsigned)i);
		std::shuffle(ranks.begin(), ranks.end(), rng);
		const gd::LowerStarRule rule = (i % 2) ? gd::LowerStarRule::PLS_BR : gd::LowerStarRule::PLS;
		std::vector<uint8_t> bytes((size_t)G.NumCells());
		gd::ComputeLowerStarGradient3D(ranks.data(), X, Y, Z, rule, bytes.data(), nullptr, nullptr, nullptr, 1);
		gd::GradientView3D V(G, bytes.data());
		const std::vector<uint8_t> crit0 = CritMask(V);

		// 4-cycle moves
		const std::vector<gd::Move> mv = gd::Moves(V);
		for (const gd::Move& m : mv) {
			const std::vector<uint8_t> before = bytes;
			gd::ApplyMove(V, m);
			nmoves++;
			CHECK(CritMask(V) == crit0, "%dx%dx%d: a move changed the critical set", X, Y, Z);
			const bool gate = gd::CreatesCycle(V, m);
			CHECK(gate == !gd::IsAcyclic(V), "%dx%dx%d: move gate %d disagrees with IsAcyclic", X, Y, Z, (int)gate);
			gd::UndoMove(V, m);
			CHECK(bytes == before, "%dx%dx%d: UndoMove does not restore", X, Y, Z);
		}

		// random alternating cycles over the whole grid, applied and undone
		gd::Zone all;
		{
			std::unordered_set<int64_t> s;
			for (int64_t c = 0; c < G.NumCells(); c++) s.insert(c);
			all.Assign(s);
		}
		std::mt19937_64 r64(900u + (unsigned)i);
		for (int t = 0; t < (quick ? 150 : 600); t++) {
			const int kmax = 2 + (int)(r64() % 5);
			const std::vector<int64_t> xs = gd::FindCycle(V, all, r64, 2, kmax);
			if (xs.empty()) continue;
			if (!gd::CycleOk(V, xs)) continue;
			if (xs.size() == 2) {
				// the same four cells form one of the CURRENT gradient's 4-cycle moves
				const int64_t c0 = xs[0], c1 = xs[1], p0 = V.Partner(c0), p1 = V.Partner(c1);
				const std::vector<int64_t> four = { c0, c1, p0, p1 };
				const std::set<int64_t> want(four.begin(), four.end());
				bool listed = false;
				for (const gd::Move& m : gd::LocalMoves(V, four))
					if (std::set<int64_t>{ m.b, m.v, m.a, m.q } == want) { listed = true; break; }
				CHECK(listed, "%dx%dx%d: a k=2 cycle is not a listed move", X, Y, Z);
			}
			const std::vector<uint8_t> before = bytes;
			const gd::CycleToken tok = gd::ApplyCycle(V, xs);
			ncycles++;
			CHECK(CritMask(V) == crit0, "%dx%dx%d: a cycle changed the critical set", X, Y, Z);
			CHECK(gd::CheckConsistency(V) == 0, "%dx%dx%d: a cycle broke the matching", X, Y, Z);
			const bool gate = gd::CycleCreatesClosedVPath(V, tok);
			CHECK(gate == !gd::IsAcyclic(V), "%dx%dx%d: cycle gate %d disagrees with IsAcyclic (k=%zu)", X, Y, Z, (int)gate, xs.size());
			if (gate || (r64() & 1)) {
				gd::UndoCycle(V, tok);
				CHECK(bytes == before, "%dx%dx%d: UndoCycle does not restore", X, Y, Z);
			}
		}
		ValidMorse("after random acyclic flips", V, crit0);

		// hunts on the most-strangled connection, if any
		const auto st = gd::Strangled(V, false);
		if (!st.empty()) {
			const int64_t tau = st[0].tau, sigma = st[0].sigma;
			gd::HuntOptions o;
			o.rounds = quick ? 300 : 1500;
			o.verify = true;
			std::vector<uint8_t> keep = bytes;
			gd::HuntResult R1 = gd::HuntCycles(V, tau, sigma, (o.kmin = 2, o.kmax = 6, o));
			CHECK(R1.critical_set_preserved, "HuntCycles changed the critical set");
			ValidMorse("HuntCycles", V, crit0);
			bytes = keep;
			gd::HuntResult R2 = gd::HuntLocal(V, tau, sigma, o);
			CHECK(R2.critical_set_preserved, "HuntLocal changed the critical set");
			ValidMorse("HuntLocal", V, crit0);
			bytes = keep;
			gd::HuntResult R3 = gd::Hunt(V, tau, sigma, o);
			CHECK(R3.critical_set_preserved, "Hunt changed the critical set");
			ValidMorse("Hunt", V, crit0);
			bytes = keep;
			o.rounds = quick ? 150 : 600;
			gd::HuntResult R4 = gd::HuntAll(V, o, false);
			CHECK(R4.critical_set_preserved, "HuntAll changed the critical set");
			ValidMorse("HuntAll", V, crit0);
			bytes = keep;
		}
	}
	printf("dmt_cycles_smoke: %lld moves, %lld cycles, %d checks, %d failures -> %s\n",
		(long long)nmoves, (long long)ncycles, g_checks, g_fail, g_fail ? "FAIL" : "ok");
	return g_fail ? 1 : 0;
}
