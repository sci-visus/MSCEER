// The star tables built from TopologicalRegularGrid3D's own iterators must
// describe the same topology as the canonical tables the device-safe core uses.
//
// Exit code: 0 iff they are equivalent and internally consistent.

#include "gi_basic_types.h"
#include "gi_regular_grid.h"
#include "gi_topological_regular_grid.h"
#include "dmt/gi_dmt_star27_mesh_tables.h"

#include <cstdio>

using namespace GInt;

int main() {
	int bad = 0;
	::dmt::Star27Tables C;
	::dmt::BuildStar27TablesCanonical(C);
	for (int s = 0; s < 27; s++) {
		const int d = C.dim[s];
		if (C.nfac[s] != d || C.ncof[s] != 2 * (3 - d) || C.nvert[s] != (1 << d)) {
			printf("  canonical slot %d inconsistent\n", s);
			bad++;
		}
	}
	if (C.dim[::dmt::kStarCentre] != 0) bad++;

	// Several mesh sizes, every interior vertex as centre.
	int dims[][3] = { {3, 3, 3}, {4, 5, 6}, {7, 3, 4} };
	for (auto& dd : dims) {
		RegularGrid3D grid(Vec3l(dd[0], dd[1], dd[2]), Vec3b(0, 0, 0));
		TopologicalRegularGrid3D mesh(&grid);
		int checked = 0;
		for (int z = 1; z < dd[2] - 1; z++)
			for (int y = 1; y < dd[1] - 1; y++)
				for (int x = 1; x < dd[0] - 1; x++) {
					::dmt::Star27Tables M;
					const INDEX_TYPE centre = mesh.coords2Cellid(Vec3l(2 * x, 2 * y, 2 * z));
					if (!GInt::dmt::BuildStar27TablesFromMesh(&mesh, centre, M) || !GInt::dmt::TablesEquivalent(M, C)) {
						printf("  mesh %dx%dx%d centre (%d,%d,%d): tables differ\n", dd[0], dd[1], dd[2], x, y, z);
						bad++;
					}
					checked++;
				}
		printf("  mesh %dx%dx%d: %d interior centres checked\n", dd[0], dd[1], dd[2], checked);
	}
	printf("dmt_tables_validate: %s\n", bad ? "FAIL" : "ok");
	return bad ? 1 : 0;
}
