#ifndef DMT_STAR27_TABLES_H
#define DMT_STAR27_TABLES_H

// Topology of the star of a vertex in the cubical complex, as fixed tables.
//
// The star of vertex v is every cell containing v; in doubled coordinates
// these are the 27 cells v + (dx, dy, dz), d in {-1, 0, 1}. Cell slots and
// vertex slots use the same numbering,
//     s = (dx + 1) + 3 (dy + 1) + 9 (dz + 1),   centre = 13,
// cell slot offsets in doubled units, vertex slot offsets in vertex units.
// This is the 27-slot convention planned for the 3D GPU kernel.
//
// Within the star, a cell's facets that still contain v are obtained by
// zeroing one of its nonzero offsets; its cofacets in the star by setting one
// zero offset to -1 or +1. Facets that do not contain v lie outside the star
// and are never needed by a lower-star rule.

#include "dmt_hd.h"

namespace dmt {

	constexpr int kStarSlots = 27;
	constexpr int kStarCentre = 13;

	struct Star27Tables {
		int8_t off[27][3];      // doubled-coordinate offset of each cell slot
		uint8_t dim[27];        // dimension of the cell in each slot
		uint8_t nfac[27];       // facets inside the star (= dim)
		uint8_t fac[27][3];
		uint8_t ncof[27];       // cofacets inside the star (= 2 * (3 - dim))
		uint8_t cof[27][6];
		uint8_t nvert[27];      // vertices of the cell (= 2^dim)
		uint8_t vert[27][8];    // as vertex-block slots
		uint32_t facmask[27];   // fac[] as a bitmask over slots
		uint32_t cofmask[27];   // cof[] as a bitmask over slots
		uint32_t vertmask[27];  // vert[] as a bitmask over vertex slots
		uint32_t dimmask[4];    // slots of each dimension
	};

	DMT_HD int StarSlot(int dx, int dy, int dz) { return (dx + 1) + 3 * (dy + 1) + 9 * (dz + 1); }

	// Fills the masks from the list fields; shared by every builder.
	DMT_HD void FinishStar27Masks(Star27Tables& T) {
		for (int d = 0; d < 4; d++) T.dimmask[d] = 0;
		for (int s = 0; s < 27; s++) {
			T.facmask[s] = 0;
			T.cofmask[s] = 0;
			T.vertmask[s] = 0;
			for (int i = 0; i < T.nfac[s]; i++) T.facmask[s] |= 1u << T.fac[s][i];
			for (int i = 0; i < T.ncof[s]; i++) T.cofmask[s] |= 1u << T.cof[s][i];
			for (int i = 0; i < T.nvert[s]; i++) T.vertmask[s] |= 1u << T.vert[s][i];
			T.dimmask[T.dim[s]] |= 1u << s;
		}
	}

	// The tables by arithmetic, in a fixed order (axes ascending, -1 before +1).
	DMT_HD void BuildStar27TablesCanonical(Star27Tables& T) {
		for (int dz = -1; dz <= 1; dz++)
			for (int dy = -1; dy <= 1; dy++)
				for (int dx = -1; dx <= 1; dx++) {
					const int s = StarSlot(dx, dy, dz);
					const int o[3] = { dx, dy, dz };
					T.off[s][0] = (int8_t)dx;
					T.off[s][1] = (int8_t)dy;
					T.off[s][2] = (int8_t)dz;
					int dim = 0;
					for (int i = 0; i < 3; i++)
						if (o[i] != 0) dim++;
					T.dim[s] = (uint8_t)dim;

					T.nfac[s] = 0;
					for (int i = 0; i < 3; i++) {
						if (o[i] == 0) continue;
						int p[3] = { o[0], o[1], o[2] };
						p[i] = 0;
						T.fac[s][T.nfac[s]++] = (uint8_t)StarSlot(p[0], p[1], p[2]);
					}

					T.ncof[s] = 0;
					for (int i = 0; i < 3; i++) {
						if (o[i] != 0) continue;
						for (int sg = -1; sg <= 1; sg += 2) {
							int p[3] = { o[0], o[1], o[2] };
							p[i] = sg;
							T.cof[s][T.ncof[s]++] = (uint8_t)StarSlot(p[0], p[1], p[2]);
						}
					}

					// vertices: each nonzero offset spans vertex offsets {0, o}
					T.nvert[s] = 0;
					for (int vx = 0; vx <= 1; vx++) {
						if (vx == 1 && o[0] == 0) continue;
						for (int vy = 0; vy <= 1; vy++) {
							if (vy == 1 && o[1] == 0) continue;
							for (int vz = 0; vz <= 1; vz++) {
								if (vz == 1 && o[2] == 0) continue;
								T.vert[s][T.nvert[s]++] = (uint8_t)StarSlot(vx * o[0], vy * o[1], vz * o[2]);
							}
						}
					}
				}
		FinishStar27Masks(T);
	}

	// Direction code from cell slot a to an adjacent cell slot b (0:+x 1:-x
	// 2:+y 3:-y 4:+z 5:-z), or 7 if they are not adjacent.
	DMT_HD int StarDir(const Star27Tables& T, int a, int b) {
		int axis = -1;
		for (int i = 0; i < 3; i++) {
			if (T.off[a][i] != T.off[b][i]) {
				if (axis >= 0) return 7;
				axis = i;
			}
		}
		if (axis < 0) return 7;
		int d = T.off[b][axis] - T.off[a][axis];
		if (d == 1) return 2 * axis;
		if (d == -1) return 2 * axis + 1;
		return 7;
	}

} // namespace dmt

#endif
