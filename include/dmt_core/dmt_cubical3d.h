#ifndef DMT_CUBICAL3D_H
#define DMT_CUBICAL3D_H

// The cubical complex of an X x Y x Z vertex grid in doubled coordinates.
//
// A cell is (cx, cy, cz) with each coordinate in [0, 2n-2]; its dimension is the
// number of odd coordinates, and a vertex (x, y, z) is the cell (2x, 2y, 2z).
// Cell ids are the MSCEER linear ids of TopologicalRegularGrid3D,
//     id = cx + Cx * (cy + Cy * cz),   C* = 2 n* - 1,
// and vertex ids are the RegularGrid3D ids, v = x + X * (y + Y * z).
//
// Facets and cofacets are enumerated in a fixed order: axes ascending, and on
// each axis the -1 neighbour before the +1 neighbour. Direction codes match
// TopologicalRegularGrid3D::Compress6NeighborOffsetToByte:
//     0:+x 1:-x 2:+y 3:-y 4:+z 5:-z.

#include "dmt_hd.h"

namespace dmt {

	struct Cubical3D {
		int64_t n[3];     // vertices per axis
		int64_t c[3];     // cells per axis, 2n-1
		int64_t stride[3];// cell-id strides: 1, Cx, Cx*Cy

		DMT_HD static Cubical3D Make(int64_t X, int64_t Y, int64_t Z) {
			Cubical3D g;
			g.n[0] = X; g.n[1] = Y; g.n[2] = Z;
			for (int i = 0; i < 3; i++) g.c[i] = 2 * g.n[i] - 1;
			g.stride[0] = 1;
			g.stride[1] = g.c[0];
			g.stride[2] = g.c[0] * g.c[1];
			return g;
		}

		DMT_HD int64_t NumCells() const { return c[0] * c[1] * c[2]; }
		DMT_HD int64_t NumVertices() const { return n[0] * n[1] * n[2]; }
		DMT_HD int64_t Top(int axis) const { return c[axis] - 1; }

		DMT_HD int64_t Id(int64_t x, int64_t y, int64_t z) const {
			return x + c[0] * (y + c[1] * z);
		}
		DMT_HD void Coords(int64_t id, int64_t out[3]) const {
			out[0] = id % c[0];
			int64_t r = id / c[0];
			out[1] = r % c[1];
			out[2] = r / c[1];
		}
		DMT_HD bool InDomain(int64_t x, int64_t y, int64_t z) const {
			return x >= 0 && y >= 0 && z >= 0 && x < c[0] && y < c[1] && z < c[2];
		}

		DMT_HD int Dim(int64_t id) const {
			int64_t q[3];
			Coords(id, q);
			return (int)((q[0] & 1) + (q[1] & 1) + (q[2] & 1));
		}

		// On the boundary of the domain: some coordinate is 0 or its maximum.
		// Only even coordinates can be extreme, so this is the subcomplex of
		// cells lying in the boundary surface.
		DMT_HD bool OnBoundary(int64_t id) const {
			int64_t q[3];
			Coords(id, q);
			for (int i = 0; i < 3; i++)
				if (q[i] == 0 || q[i] == Top(i)) return true;
			return false;
		}

		// Number of axes on which the cell sits at an extreme (0..3); this is
		// TopologicalRegularGrid3D::boundaryValue for a non-periodic grid.
		DMT_HD int BoundaryValue(int64_t id) const {
			int64_t q[3];
			Coords(id, q);
			int b = 0;
			for (int i = 0; i < 3; i++)
				if (q[i] == 0 || q[i] == Top(i)) b++;
			return b;
		}

		// Steps from the domain boundary in doubled units; 0 means on it.
		DMT_HD int64_t Bdist(int64_t id) const {
			int64_t q[3];
			Coords(id, q);
			int64_t best = q[0] < Top(0) - q[0] ? q[0] : Top(0) - q[0];
			for (int i = 1; i < 3; i++) {
				int64_t d = q[i] < Top(i) - q[i] ? q[i] : Top(i) - q[i];
				if (d < best) best = d;
			}
			return best;
		}

		// Facets, axes ascending, -1 before +1. Returns the count (= dim).
		DMT_HD int Facets(int64_t id, int64_t out[6]) const {
			int64_t q[3];
			Coords(id, q);
			int k = 0;
			for (int i = 0; i < 3; i++) {
				if (q[i] & 1) {
					out[k++] = id - stride[i];
					out[k++] = id + stride[i];
				}
			}
			return k;
		}

		// Cofacets inside the domain, axes ascending, -1 before +1.
		DMT_HD int Cofacets(int64_t id, int64_t out[6]) const {
			int64_t q[3];
			Coords(id, q);
			int k = 0;
			for (int i = 0; i < 3; i++) {
				if ((q[i] & 1) == 0) {
					if (q[i] - 1 >= 0) out[k++] = id - stride[i];
					if (q[i] + 1 <= Top(i)) out[k++] = id + stride[i];
				}
			}
			return k;
		}

		// Vertex ids (x + X*(y + Y*z)) of the cell, in lexicographic order of
		// (x, y, z) with x slowest. Returns the count (2^dim).
		DMT_HD int Vertices(int64_t id, int64_t out[8]) const {
			int64_t q[3];
			Coords(id, q);
			int64_t lo[3], hi[3];
			for (int i = 0; i < 3; i++) {
				lo[i] = (q[i] & 1) ? (q[i] - 1) / 2 : q[i] / 2;
				hi[i] = (q[i] & 1) ? (q[i] + 1) / 2 : q[i] / 2;
			}
			int k = 0;
			for (int64_t x = lo[0]; x <= hi[0]; x++)
				for (int64_t y = lo[1]; y <= hi[1]; y++)
					for (int64_t z = lo[2]; z <= hi[2]; z++)
						out[k++] = x + n[0] * (y + n[1] * z);
			return k;
		}

		DMT_HD int64_t VertexCell(int64_t vid) const {
			int64_t x = vid % n[0];
			int64_t r = vid / n[0];
			int64_t y = r % n[1];
			int64_t z = r / n[1];
			return Id(2 * x, 2 * y, 2 * z);
		}
		DMT_HD int64_t CellVertex(int64_t vertex_cell) const {
			int64_t q[3];
			Coords(vertex_cell, q);
			return q[0] / 2 + n[0] * (q[1] / 2 + n[1] * (q[2] / 2));
		}

		// Cell reached from `id` along direction code 0..5.
		DMT_HD int64_t Neighbor(int64_t id, int dir) const {
			int axis = dir >> 1;
			return (dir & 1) ? id - stride[axis] : id + stride[axis];
		}

		// Direction code from a to an adjacent cell b (one axis, +-1); 7 if
		// they are not adjacent that way.
		DMT_HD int DirBetween(int64_t a, int64_t b) const {
			int64_t qa[3], qb[3];
			Coords(a, qa);
			Coords(b, qb);
			int axis = -1;
			for (int i = 0; i < 3; i++) {
				if (qa[i] != qb[i]) {
					if (axis >= 0) return 7;
					axis = i;
				}
			}
			if (axis < 0) return 7;
			int64_t d = qb[axis] - qa[axis];
			if (d == 1) return 2 * axis;
			if (d == -1) return 2 * axis + 1;
			return 7;
		}

		// Forman's cubical incidence <d tau, sigma>: nonzero only when sigma
		// differs from tau on exactly one axis a, on which tau is odd and the
		// difference is +-1. Then it is (-1)^k * (sigma_a - tau_a), where k is
		// the number of odd coordinates of tau before axis a.
		DMT_HD int Incidence(int64_t tau, int64_t sigma) const {
			int64_t qt[3], qs[3];
			Coords(tau, qt);
			Coords(sigma, qs);
			int a = -1;
			for (int i = 0; i < 3; i++) {
				if (qt[i] != qs[i]) {
					if (a >= 0) return 0;
					a = i;
				}
			}
			if (a < 0 || (qt[a] & 1) == 0) return 0;
			int64_t d = qs[a] - qt[a];
			if (d != 1 && d != -1) return 0;
			int k = 0;
			for (int i = 0; i < a; i++)
				if (qt[i] & 1) k++;
			return (k & 1) ? (int)-d : (int)d;
		}
	};

} // namespace dmt

#endif
