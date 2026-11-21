// TopologicalGradientUsingAlgorithms::CheckGradientForLoops on three gradients:
//   1. a hand-built closed V-path (four quads around a vertex in the z=1 plane);
//   2. the same four pairs with one arrow reversed, which is an open path;
//   3. MSCEER's own 3D lower-star gradient (MyRobinsNoalloc) on random fields,
//      which is acyclic by construction.
//
// Exit code: 0 iff every case reports the expected number of cycles.

#include "gi_basic_types.h"
#include "gi_timing.h"  // gi_max_vertex_labeling.h uses ThreadedTimer without including it
#include "gi_regular_grid.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_topological_regular_grid.h"
#include "gi_max_vertex_labeling.h"
#include "gi_discrete_gradient_labeling.h"
#include "gi_modified_robins.h"
#include "gi_topological_gradient_using_algorithms.h"

#include <cstdio>
#include <random>
#include <vector>

using namespace GInt;

typedef RegularGrid3D GridType;
typedef RegularGridTrilinearFunction GridFuncType;
typedef TopologicalRegularGrid3D MeshType;
typedef DiscreteGradientLabeling<MeshType> GradType;
typedef RegularGridMaxMinVertexLabeling3D<MeshType, GridFuncType> MaxVLType;
typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType, 4, 6> RobinsType;
typedef TopologicalGradientUsingAlgorithms<MeshType, void, GradType> AlgType;

static int g_fail = 0;

static void expect(const char* what, int got, int want) {
	printf("  %-48s cycles=%d (expected %d) %s\n", what, got, want, got == want ? "ok" : "FAIL");
	if (got != want) g_fail++;
}

static INDEX_TYPE cell(MeshType* m, int x, int y, int z) {
	return m->coords2Cellid(Vec3l(x, y, z));
}

static void pair_cells(GradType* g, INDEX_TYPE a, INDEX_TYPE b) {
	g->setPair(a, b);
	g->setPair(b, a);
	g->setAssigned(a, true);
	g->setAssigned(b, true);
}

static void all_critical(MeshType* m, GradType* g) {
	g->ClearAllGradient();
	MeshType::AllCellsIterator it(m);
	for (it.begin(); it.valid(); it.advance()) {
		g->setCritical(it.value(), true);
		g->setAssigned(it.value(), true);
	}
}

int main() {
	printf("loopcheck_smoke\n");

	// 3x3x3 vertices -> doubled cells 5x5x5, centre vertex at (2,2,2).
	GridType* grid = new GridType(Vec3l(3, 3, 3), Vec3b(0, 0, 0));
	MeshType* mesh = new MeshType(grid);
	GradType* grad = new GradType(mesh);
	AlgType alg(nullptr, mesh, grad);

	// Edges radiating from the centre in the z=2 plane, and the quads between them.
	INDEX_TYPE E[4] = { cell(mesh, 3, 2, 2), cell(mesh, 2, 3, 2), cell(mesh, 1, 2, 2), cell(mesh, 2, 1, 2) };
	INDEX_TYPE Q[4] = { cell(mesh, 3, 3, 2), cell(mesh, 1, 3, 2), cell(mesh, 1, 1, 2), cell(mesh, 3, 1, 2) };
	// Qi contains Ei and E(i+1): pairing Ei~Qi makes E1->Q1->E2->Q2->E3->Q3->E4->Q4->E1.
	all_critical(mesh, grad);
	for (int i = 0; i < 4; i++) pair_cells(grad, E[i], Q[i]);
	expect("closed V-path of four quads", alg.CheckGradientForLoops(), 1);

	// Drop the E0~Q0 pair: E1->Q1->E2->Q2->E3->Q3 then reaches the critical E0 and stops.
	all_critical(mesh, grad);
	pair_cells(grad, E[1], Q[1]);
	pair_cells(grad, E[2], Q[2]);
	pair_cells(grad, E[3], Q[3]);
	expect("open V-path (no cycle)", alg.CheckGradientForLoops(), 0);

	// An all-critical field has no arrows at all.
	all_critical(mesh, grad);
	expect("all cells critical", alg.CheckGradientForLoops(), 0);

	delete grad;
	delete mesh;
	delete grid;

	// MSCEER's own 3D lower-star gradient on random distinct fields.
	int dims[][3] = { {3, 3, 4}, {5, 6, 7}, {9, 8, 7}, {12, 10, 11} };
	for (auto& d : dims) {
		GridType* g = new GridType(Vec3l(d[0], d[1], d[2]), Vec3b(0, 0, 0));
		GridFuncType* f = new GridFuncType(g);
		const INDEX_TYPE n = (INDEX_TYPE)d[0] * d[1] * d[2];
		std::vector<float> vals(n);
		for (INDEX_TYPE i = 0; i < n; i++) vals[i] = (float)i;
		std::mt19937 rng(1234u + (unsigned)n);
		std::shuffle(vals.begin(), vals.end(), rng);
		f->DeepCopyImage(vals.data());
		MeshType* m = new MeshType(g);
		MaxVLType* maxv = new MaxVLType(m, f);
		maxv->ComputeOutput();
		GradType* gr = new GradType(m);
		gr->ClearAllGradient();
		RobinsType robins(m, maxv, gr);
		robins.ComputePairing();
		AlgType a(nullptr, m, gr);
		char what[64];
		snprintf(what, sizeof(what), "MyRobinsNoalloc %dx%dx%d", d[0], d[1], d[2]);
		expect(what, a.CheckGradientForLoops(), 0);
		delete gr;
		delete maxv;
		delete m;
		delete f;
		delete g;
	}

	printf("loopcheck_smoke: %s\n", g_fail ? "FAIL" : "ok");
	return g_fail ? 1 : 0;
}
