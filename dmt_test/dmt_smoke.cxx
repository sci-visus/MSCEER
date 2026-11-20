// Invariants of the lower-star rules on random-permutation fields.
//
// For every field (sizes 2..9 per axis, many seeds) and both rules:
//   - every cell assigned, no orphans, a consistent matching, acyclic;
//   - Euler characteristic of the critical cells = 1 (a ball);
//   - PLS-BR: no pair joins a boundary cell to an interior cell, and the
//     boundary critical cells have Euler characteristic 2 (a sphere);
//   - the result is identical for 1, 2 and 8 threads;
//   - the raw driver and the MSCEER-shaped template driver (fed by
//     RankMaxMinVertexLabeling3D, and by RegularGridMaxMinVertexLabeling3D
//     over the same values as a float image) write identical bytes;
//   - MyRobinsNoalloc fed by RankMaxMinVertexLabeling3D writes the same bytes
//     as when fed by RegularGridMaxMinVertexLabeling3D.
//
// Exit code: 0 iff every check passes.

#include "gi_basic_types.h"
#include "gi_timing.h"   // gi_max_vertex_labeling.h uses ThreadedTimer without including it
#include "gi_regular_grid.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_topological_regular_grid.h"
#include "gi_max_vertex_labeling.h"
#include "gi_discrete_gradient_labeling.h"
#include "gi_modified_robins.h"
#include "gi_dmt_rank_field.h"
#include "gi_dmt_lowerstar_gradient.h"
#include "gi_dmt_acyclic.h"

#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

using namespace GInt;

typedef RegularGrid3D GridType;
typedef RegularGridTrilinearFunction GridFuncType;
typedef TopologicalRegularGrid3D MeshType;
typedef DiscreteGradientLabeling<MeshType> GradType;
typedef RegularGridMaxMinVertexLabeling3D<MeshType, GridFuncType> FloatMaxVL;
typedef GInt::dmt::RankMaxMinVertexLabeling3D<MeshType> RankMaxVL;

static int g_fail = 0;
static int g_checks = 0;
#define CHECK(cond, ...) do { g_checks++; if (!(cond)) { g_fail++; if (g_fail <= 20) { printf("  FAIL: "); printf(__VA_ARGS__); printf("\n"); } } } while (0)

static const char* RuleName(GInt::dmt::LowerStarRule r) { return r == GInt::dmt::LowerStarRule::PLS ? "PLS" : "PLS-BR"; }

static void CheckField(int X, int Y, int Z, unsigned seed, bool template_checks) {
	const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(X, Y, Z);
	const int64_t nv = G.NumVertices(), nc = G.NumCells();
	std::vector<uint32_t> ranks((size_t)nv);
	for (int64_t i = 0; i < nv; i++) ranks[(size_t)i] = (uint32_t)i;
	std::mt19937 rng(seed);
	std::shuffle(ranks.begin(), ranks.end(), rng);

	for (GInt::dmt::LowerStarRule rule : { GInt::dmt::LowerStarRule::PLS, GInt::dmt::LowerStarRule::PLS_BR }) {
		std::vector<uint8_t> g1((size_t)nc), g2((size_t)nc), g8((size_t)nc);
		GInt::dmt::LowerStarStats st;
		const bool ok = GInt::dmt::ComputeLowerStarGradient3D(ranks.data(), X, Y, Z, rule, g1.data(), nullptr, nullptr, &st, 1);
		CHECK(ok && st.unassigned == 0, "%s %dx%dx%d seed %u: %lld unassigned", RuleName(rule), X, Y, Z, seed, (long long)st.unassigned);
		CHECK(st.orphans == 0, "%s %dx%dx%d seed %u: %lld orphans", RuleName(rule), X, Y, Z, seed, (long long)st.orphans);
		GInt::dmt::ComputeLowerStarGradient3D(ranks.data(), X, Y, Z, rule, g2.data(), nullptr, nullptr, nullptr, 2);
		GInt::dmt::ComputeLowerStarGradient3D(ranks.data(), X, Y, Z, rule, g8.data(), nullptr, nullptr, nullptr, 8);
		CHECK(g1 == g2 && g1 == g8, "%s %dx%dx%d seed %u: thread count changes the result", RuleName(rule), X, Y, Z, seed);

		GInt::dmt::GradientView3D V(G, g1.data());
		CHECK(GInt::dmt::CheckConsistency(V) == 0, "%s %dx%dx%d seed %u: inconsistent matching", RuleName(rule), X, Y, Z, seed);
		CHECK(GInt::dmt::IsAcyclic(V), "%s %dx%dx%d seed %u: closed V-path", RuleName(rule), X, Y, Z, seed);
		const int64_t chi = st.critical[0] - st.critical[1] + st.critical[2] - st.critical[3];
		CHECK(chi == 1, "%s %dx%dx%d seed %u: Euler characteristic %lld", RuleName(rule), X, Y, Z, seed, (long long)chi);

		if (rule == GInt::dmt::LowerStarRule::PLS_BR) {
			int64_t cross = 0, bchi = 0;
			for (int64_t c = 0; c < nc; c++) {
				const bool b = G.OnBoundary(c);
				if (V.IsCritical(c)) {
					if (b) bchi += (G.Dim(c) % 2 == 0) ? 1 : -1;
				}
				else if (b != G.OnBoundary(V.Partner(c))) cross++;
			}
			CHECK(cross == 0, "PLS-BR %dx%dx%d seed %u: %lld pairs cross boundary/interior", X, Y, Z, seed, (long long)cross);
			CHECK(bchi == 2, "PLS-BR %dx%dx%d seed %u: boundary Euler characteristic %lld", X, Y, Z, seed, (long long)bchi);
		}

		if (!template_checks) continue;

		// Template driver on exact ranks.
		GridType grid(Vec3l(X, Y, Z), Vec3b(0, 0, 0));
		MeshType mesh(&grid);
		GInt::dmt::RankField3D field(X, Y, Z);
		field.ranks = ranks;
		RankMaxVL rmaxv(&mesh, &field);
		rmaxv.ComputeOutput();
		GradType gt(&mesh);
		gt.ClearAllGradient();
		if (rule == GInt::dmt::LowerStarRule::PLS) {
			GInt::dmt::PLSGradient<MeshType, RankMaxVL, GradType> r(&mesh, &rmaxv, &gt);
			r.ComputePairing();
		}
		else {
			GInt::dmt::PLSBRGradient<MeshType, RankMaxVL, GradType> r(&mesh, &rmaxv, &gt);
			r.ComputePairing();
		}
		CHECK(std::memcmp(gt.m_dgrad->LabelArray(), g1.data(), (size_t)nc) == 0,
			"%s %dx%dx%d seed %u: template (rank labeling) differs from raw driver", RuleName(rule), X, Y, Z, seed);

		// Template driver on MSCEER's float-image labeling of the same order.
		GridFuncType func(&grid);
		std::vector<float> img((size_t)nv);
		for (int64_t i = 0; i < nv; i++) img[(size_t)i] = (float)ranks[(size_t)i];
		func.DeepCopyImage(img.data());
		FloatMaxVL fmaxv(&mesh, &func);
		fmaxv.ComputeOutput();
		GradType gf(&mesh);
		gf.ClearAllGradient();
		if (rule == GInt::dmt::LowerStarRule::PLS) {
			GInt::dmt::PLSGradient<MeshType, FloatMaxVL, GradType> r(&mesh, &fmaxv, &gf);
			r.ComputePairing();
		}
		else {
			GInt::dmt::PLSBRGradient<MeshType, FloatMaxVL, GradType> r(&mesh, &fmaxv, &gf);
			r.ComputePairing();
		}
		CHECK(std::memcmp(gf.m_dgrad->LabelArray(), g1.data(), (size_t)nc) == 0,
			"%s %dx%dx%d seed %u: template (float labeling) differs from raw driver", RuleName(rule), X, Y, Z, seed);

		// The rank labeling is a faithful drop-in for MSCEER's own rule.
		if (rule == GInt::dmt::LowerStarRule::PLS) {
			GradType ga(&mesh), gb(&mesh);
			ga.ClearAllGradient();
			gb.ClearAllGradient();
			MyRobinsNoalloc<MeshType, RankMaxVL, GradType, 4, 6> ra(&mesh, &rmaxv, &ga);
			ra.ComputePairing();
			MyRobinsNoalloc<MeshType, FloatMaxVL, GradType, 4, 6> rb(&mesh, &fmaxv, &gb);
			rb.ComputePairing();
			CHECK(std::memcmp(ga.m_dgrad->LabelArray(), gb.m_dgrad->LabelArray(), (size_t)nc) == 0,
				"MyRobinsNoalloc %dx%dx%d seed %u: rank labeling differs from float labeling", X, Y, Z, seed);
		}
	}
}

int main(int argc, char** argv) {
	const bool quick = argc > 1 && std::strcmp(argv[1], "--quick") == 0;
	printf("dmt_smoke%s\n", quick ? " --quick" : "");
	int fields = 0;
	std::mt19937 pick(20261119u);
	const int nrandom = quick ? 120 : 1000;
	for (int i = 0; i < nrandom; i++) {
		const int X = 2 + (int)(pick() % 8), Y = 2 + (int)(pick() % 8), Z = 2 + (int)(pick() % 8);
		// The template drivers build MSCEER meshes and print; run them on a subset.
		CheckField(X, Y, Z, 1000u + (unsigned)i, i % 10 == 0);
		fields++;
	}
	// A couple of larger ones.
	CheckField(17, 13, 11, 7u, true);
	CheckField(24, 20, 16, 8u, false);
	fields += 2;
	printf("dmt_smoke: %d fields, %d checks, %d failures -> %s\n", fields, g_checks, g_fail, g_fail ? "FAIL" : "ok");
	return g_fail ? 1 : 0;
}
