// Gradient-path counting and the Morse complex.
//
//   - U128 / S128: carries at 2^64, the 2^127 boundary, sign changes,
//     products, and the sticky overflow flag;
//   - d o d = 0 on the cubical complex;
//   - on random-permutation fields, for PLS, PLS-BR and DMF: d~ o d~ = 0
//     (exactly), Betti numbers [1, 0, 0, 0] over Z/p, the serial and the
//     level-parallel memo give identical boundary matrices, and the scalar
//     PairScorer equals the MorseRow entry on every connection;
//   - a hand-built closed V-path raises ClosedVPath.
//
// Exit code: 0 iff every check passes.

#include "dmt/gi_dmt_acyclic.h"
#include "gi_dmt_lowerstar_gradient.h"
#include "gradient/gi_dmf_gradient.h"
#include "dmt/gi_dmt_morse_complex.h"
#include "dmt/gi_dmt_vpath_counter.h"

#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

static int g_fail = 0;
static int g_checks = 0;
#define CHECK(cond, ...) do { g_checks++; if (!(cond)) { g_fail++; if (g_fail <= 20) { printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while (0)

namespace gd = GInt::dmt;
using ::dmt::U128;
using ::dmt::S128;

static void Count128Tests() {
	U128 a(~0ull), b(1);
	U128 c = ::dmt::Add(a, b);
	CHECK(c.lo == 0 && c.hi == 1 && !c.overflow, "U128 carry at 2^64");
	U128 top; top.lo = ~0ull; top.hi = ~0ull;
	CHECK(::dmt::Add(top, b).overflow, "U128 overflow at 2^128");
	CHECK(!::dmt::Add(top, U128()).overflow, "U128 + 0 at max is exact");

	S128 m1 = S128::FromInt(-1), p1 = S128::FromInt(1);
	CHECK(::dmt::Add(m1, p1).IsZero(), "S128 -1 + 1");
	S128 x = S128::FromInt(INT64_MIN);
	S128 y = ::dmt::Add(x, m1);
	CHECK(!y.overflow && !y.FitsI64() && y.Negative(), "S128 below INT64_MIN");
	S128 z = ::dmt::Neg(y);
	// y = -(2^63 + 1), so -y = 2^63 + 1
	CHECK(!z.overflow && !z.Negative() && z.lo == 0x8000000000000001ull && z.hi == 0, "S128 negate past 2^63");
	S128 maxp; maxp.lo = ~0ull; maxp.hi = ~0ull >> 1;   // 2^127 - 1
	CHECK(::dmt::Add(maxp, p1).overflow, "S128 overflow at 2^127");
	CHECK(!::dmt::Add(maxp, m1).overflow, "S128 2^127 - 2");
	S128 minv; minv.lo = 0; minv.hi = (uint64_t)1 << 63;  // -2^127
	CHECK(::dmt::Neg(minv).overflow, "S128 -(-2^127) overflows");
	U128 big; big.lo = 0; big.hi = (uint64_t)1 << 63;     // 2^127
	CHECK(::dmt::ToSigned(big).overflow, "U128 2^127 does not fit S128");

	S128 two64; two64.lo = 0; two64.hi = 1;
	S128 p = ::dmt::Mul(two64, S128::FromInt(-3));
	CHECK(!p.overflow && p.Negative() && ::dmt::Abs(p).hi == 3 && ::dmt::Abs(p).lo == 0, "S128 2^64 * -3");
	CHECK(::dmt::Mul(two64, two64).overflow, "S128 2^64 * 2^64 overflows");
	S128 q = ::dmt::Mul(S128::FromInt(-7), S128::FromInt(-6));
	CHECK(!q.overflow && q.FitsI64() && q.ToI64() == 42, "S128 -7 * -6");
	S128 sticky = y;
	sticky.overflow = true;
	CHECK(::dmt::Add(sticky, p1).overflow && ::dmt::Neg(sticky).overflow, "overflow is sticky");
}

static bool SameMatrix(const gd::MorseMatrix& A, const gd::MorseMatrix& B) {
	if (A.taus != B.taus || A.rowptr != B.rowptr || A.entries.size() != B.entries.size()) return false;
	for (size_t i = 0; i < A.entries.size(); i++) {
		const auto& a = A.entries[i];
		const auto& b = B.entries[i];
		if (a.sigma != b.sigma || a.u.lo != b.u.lo || a.u.hi != b.u.hi || a.coeff.lo != b.coeff.lo ||
			a.coeff.hi != b.coeff.hi || a.best.lo != b.best.lo || a.best.hi != b.best.hi)
			return false;
	}
	return true;
}

static void CheckGradient(const char* rule, int X, int Y, int Z, unsigned seed, std::vector<uint8_t>& bytes) {
	const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(X, Y, Z);
	gd::GradientView3D V(G, bytes.data());
	gd::MorseComplexCheck R = gd::CheckMorseComplex(V, 3, 1);
	CHECK(R.bad.empty() && R.exact, "%s %dx%dx%d seed %u: d~ o d~ has %zu nonzero entries", rule, X, Y, Z, seed, R.bad.size());
	const std::vector<int64_t> b = gd::BettiModP(V);
	CHECK(b.size() == 4 && b[0] == 1 && b[1] == 0 && b[2] == 0 && b[3] == 0,
		"%s %dx%dx%d seed %u: Betti %lld %lld %lld %lld", rule, X, Y, Z, seed,
		(long long)b[0], (long long)b[1], (long long)b[2], (long long)b[3]);
	gd::PairScorer scorer(G);
	for (int p = 0; p < 3; p++) {
		const gd::MorseMatrix serial = R.B[(size_t)p];
		const gd::MorseMatrix lev = gd::MorseBoundary(V, p, true, 4);
		CHECK(SameMatrix(serial, lev), "%s %dx%dx%d seed %u p=%d: level-parallel memo differs", rule, X, Y, Z, seed, p);
		for (size_t i = 0; i < serial.taus.size(); i++)
			for (int64_t k = serial.rowptr[i]; k < serial.rowptr[i + 1]; k++) {
				const auto& e = serial.entries[(size_t)k];
				const gd::PairScore s = scorer.Score(V, serial.taus[i], e.sigma);
				CHECK(s.paths.lo == e.u.lo && s.paths.hi == e.u.hi && s.coeff.lo == e.coeff.lo && s.coeff.hi == e.coeff.hi,
					"%s %dx%dx%d seed %u p=%d: PairScorer differs from MorseRow", rule, X, Y, Z, seed, p);
			}
	}
}

int main(int argc, char** argv) {
	const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
	printf("dmt_morse_smoke%s\n", quick ? " --quick" : "");
	Count128Tests();

	int dims[][3] = { {2, 2, 2}, {3, 4, 5}, {9, 8, 7} };
	for (auto& d : dims)
		CHECK(gd::CheckDD(::dmt::Cubical3D::Make(d[0], d[1], d[2])) == 0, "d o d != 0 on %dx%dx%d", d[0], d[1], d[2]);

	std::mt19937 pick(4242u);
	const int n = quick ? 40 : 250;
	for (int i = 0; i < n; i++) {
		const int X = 2 + (int)(pick() % 7), Y = 2 + (int)(pick() % 7), Z = 2 + (int)(pick() % 7);
		const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(X, Y, Z);
		std::vector<uint32_t> ranks((size_t)G.NumVertices());
		for (size_t k = 0; k < ranks.size(); k++) ranks[k] = (uint32_t)k;
		std::mt19937 rng(9000u + (unsigned)i);
		std::shuffle(ranks.begin(), ranks.end(), rng);
		std::vector<uint8_t> g((size_t)G.NumCells());
		gd::ComputeLowerStarGradient3D(ranks.data(), X, Y, Z, gd::LowerStarRule::PLS, g.data(), nullptr, nullptr, nullptr, 1);
		CheckGradient("PLS", X, Y, Z, 9000u + i, g);
		gd::ComputeLowerStarGradient3D(ranks.data(), X, Y, Z, gd::LowerStarRule::PLS_BR, g.data(), nullptr, nullptr, nullptr, 1);
		CheckGradient("PLS-BR", X, Y, Z, 9000u + i, g);
		gd::DMFOptions opt;
		gd::ComputeDMFGradient3D(ranks.data(), X, Y, Z, opt, g.data(), nullptr, nullptr, nullptr);
		CheckGradient("DMF", X, Y, Z, 9000u + i, g);
	}

	// A closed V-path (four quads around a vertex) must raise.
	{
		const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(3, 3, 3);
		std::vector<uint8_t> bytes((size_t)G.NumCells(), ::dmt::MakeGradByte(::dmt::kPairCritical));
		gd::GradientView3D V(G, bytes.data());
		const int64_t E[4] = { G.Id(3, 2, 2), G.Id(2, 3, 2), G.Id(1, 2, 2), G.Id(2, 1, 2) };
		const int64_t Q[4] = { G.Id(3, 3, 2), G.Id(1, 3, 2), G.Id(1, 1, 2), G.Id(3, 1, 2) };
		for (int i = 0; i < 4; i++) V.SetPair(E[i], Q[i]);
		gd::CounterConfig cfg;
		cfg.sinks = V.Critical(1);
		bool raised = false, raised_lev = false;
		try { gd::VPathCounter C(V, cfg); C.PathsFrom(E[0]); }
		catch (const gd::ClosedVPath&) { raised = true; }
		try { gd::VPathCounter C(V, cfg); C.PrecomputeLevels(std::vector<int64_t>(E, E + 4), 2); }
		catch (const gd::ClosedVPath&) { raised_lev = true; }
		CHECK(raised && raised_lev, "closed V-path not reported (serial %d, levels %d)", (int)raised, (int)raised_lev);
		CHECK(!gd::IsAcyclic(V), "IsAcyclic missed the closed V-path");
	}

	printf("dmt_morse_smoke: %d checks, %d failures -> %s\n", g_checks, g_fail, g_fail ? "FAIL" : "ok");
	return g_fail ? 1 : 0;
}
