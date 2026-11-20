#ifndef DMT_DMF_KEYCACHE_H
#define DMT_DMF_KEYCACHE_H

// DMF keys read back from the top-facet table instead of recomputed.
//
// The G0 of a cell's key (dmt_dmf_cell.h) is the facet with the largest key,
// which is exactly the g0 direction DMFTopFacets stores in top[] (low nibble),
// and G1 is the facet on the other side of the same axis (direction g0 ^ 1).
// So once top[] is filled for the cells of dimension d-1, the key of a d-cell
// is a walk of d levels down the g0 / g0^1 tree: 2^d - 1 byte reads of top[]
// and 2^d reads of ranks[], with no comparisons. Filling top[] one dimension
// at a time (DMFTopFacetsFromTop over dimension 1, then 2, then 3) yields the
// same table as DMFTopFacets, and the selection below the same cofacet as
// DMFLowestSelectingCofacet; the extra state is the one byte per cell top[]
// already holds.

#include "dmt_hd.h"
#include "dmt_cubical3d.h"
#include "dmt_dmf_cell.h"

namespace dmt {

	// Flattened F_w of `cell` from top[] (filled for every face of the cell of
	// dimension >= 1); returns its length 2^dim. Same result as DMFKey.
	DMT_HD int DMFKeyFromTop(const Cubical3D& G, const uint32_t* ranks, const uint8_t* top,
		int64_t cell, uint32_t out[8]) {
		int64_t q[3];
		G.Coords(cell, q);
		const int dim = (int)((q[0] & 1) + (q[1] & 1) + (q[2] & 1));
		// Cells of the current level as coordinate triples, in key order.
		int64_t lv[8][3];
		lv[0][0] = q[0]; lv[0][1] = q[1]; lv[0][2] = q[2];
		int n = 1;
		for (int l = 0; l < dim; l++) {
			// Expand in place from the back so entry i becomes 2i (G0), 2i+1 (G1).
			for (int i = n - 1; i >= 0; i--) {
				const int64_t x = lv[i][0], y = lv[i][1], z = lv[i][2];
				const int g0 = top[G.Id(x, y, z)] & 15;
				const int axis = g0 >> 1;
				const int64_t s = (g0 & 1) ? -1 : 1;
				lv[2 * i][0] = x; lv[2 * i][1] = y; lv[2 * i][2] = z;
				lv[2 * i + 1][0] = x; lv[2 * i + 1][1] = y; lv[2 * i + 1][2] = z;
				lv[2 * i][axis] += s;
				lv[2 * i + 1][axis] -= s;
			}
			n *= 2;
		}
		for (int i = 0; i < n; i++)
			out[i] = ranks[lv[i][0] / 2 + G.n[0] * (lv[i][1] / 2 + G.n[1] * (lv[i][2] / 2))];
		return n;
	}

	// Component `comp` (0 .. 2^dim - 1) of the flattened key: bit dim-1-l of
	// comp picks G0 (0) or G1 (1) at level l. Same as DMFKey(...)[comp].
	DMT_HD uint32_t DMFKeyComponentFromTop(const Cubical3D& G, const uint32_t* ranks, const uint8_t* top,
		int64_t cell, int comp) {
		int64_t q[3];
		G.Coords(cell, q);
		const int dim = (int)((q[0] & 1) + (q[1] & 1) + (q[2] & 1));
		for (int l = dim - 1; l >= 0; l--) {
			const int g0 = top[G.Id(q[0], q[1], q[2])] & 15;
			const int axis = g0 >> 1;
			const int64_t s = (g0 & 1) ? -1 : 1;
			q[axis] += ((comp >> l) & 1) ? -s : s;
		}
		return ranks[q[0] / 2 + G.n[0] * (q[1] / 2 + G.n[1] * (q[2] / 2))];
	}

	// DMFTopFacets with facet keys read from top[] (filled for the cells of
	// dimension dim(cell) - 1). Same g0 / hi2 directions.
	DMT_HD bool DMFTopFacetsFromTop(const Cubical3D& G, const uint32_t* ranks, const uint8_t* top,
		int64_t cell, uint8_t& g0_dir, uint8_t& hi2_dir) {
		int64_t f[6];
		const int nf = G.Facets(cell, f);
		if (nf == 0) return false;
		uint32_t fk[6][8];
		int flen = 1;
		for (int i = 0; i < nf; i++) flen = DMFKeyFromTop(G, ranks, top, f[i], fk[i]);
		int a = 0, b = -1;
		for (int i = 1; i < nf; i++) {
			if (DMFKeyCmp(fk[i], fk[a], flen) > 0) { b = a; a = i; }
			else if (b < 0 || DMFKeyCmp(fk[i], fk[b], flen) > 0) b = i;
		}
		g0_dir = (uint8_t)G.DirBetween(cell, f[a]);
		hi2_dir = b < 0 ? (uint8_t)7 : (uint8_t)G.DirBetween(cell, f[b]);
		return true;
	}

	// DMFLowestSelectingCofacet with cofacet keys read from top[] (filled for
	// every cell). Same cofacet.
	DMT_HD int64_t DMFLowestSelectingCofacetFromTop(const Cubical3D& G, const uint32_t* ranks,
		const uint8_t* top, int64_t a, int which) {
		int64_t c[6];
		const int nc = G.Cofacets(a, c);
		int64_t best = -1;
		uint32_t bk[8], k[8];
		int len = 0;
		for (int i = 0; i < nc; i++) {
			const uint8_t t = top[c[i]];
			const uint8_t dir = which == 0 ? (uint8_t)(t & 15u) : (uint8_t)(t >> 4);
			if (dir > 5 || G.Neighbor(c[i], dir) != a) continue;
			len = DMFKeyFromTop(G, ranks, top, c[i], k);
			if (best < 0 || DMFKeyCmp(k, bk, len) < 0) {
				best = c[i];
				for (int j = 0; j < len; j++) bk[j] = k[j];
			}
		}
		return best;
	}

	// ---------------------------------------------------------------------
	// Coordinate variants for any index width (Cubical3DT<I>): the same walks
	// and selections, but the cell is given and returned as doubled
	// coordinates, so no id -> coordinate division is ever needed. The
	// tiled GPU DMF (gpu_dgrad/dmt3d_tiled.cu) runs these with int32 on
	// tile-local grids; dmf3d_tiled_validate checks them byte for byte
	// against the id-based functions above (through the in-core GPU driver)
	// and the CPU driver.
	// ---------------------------------------------------------------------

	template <class I>
	DMT_HD int DMFKeyFromTopQ(const Cubical3DT<I>& G, const uint32_t* ranks, const uint8_t* top,
		I x, I y, I z, uint32_t out[8]) {
		const int dim = (int)((x & 1) + (y & 1) + (z & 1));
		I lv[8][3];
		lv[0][0] = x; lv[0][1] = y; lv[0][2] = z;
		int n = 1;
		for (int l = 0; l < dim; l++) {
			for (int i = n - 1; i >= 0; i--) {
				const I cx = lv[i][0], cy = lv[i][1], cz = lv[i][2];
				const int g0 = top[G.Id(cx, cy, cz)] & 15;
				const int axis = g0 >> 1;
				const I s = (g0 & 1) ? (I)-1 : (I)1;
				lv[2 * i][0] = cx; lv[2 * i][1] = cy; lv[2 * i][2] = cz;
				lv[2 * i + 1][0] = cx; lv[2 * i + 1][1] = cy; lv[2 * i + 1][2] = cz;
				lv[2 * i][axis] += s;
				lv[2 * i + 1][axis] -= s;
			}
			n *= 2;
		}
		for (int i = 0; i < n; i++)
			out[i] = ranks[lv[i][0] / 2 + G.n[0] * (lv[i][1] / 2 + G.n[1] * (lv[i][2] / 2))];
		return n;
	}

	// DMFTopFacetsFromTop by coordinates. Facets in the same order as
	// Cubical3D::Facets (axes ascending, -1 before +1).
	template <class I>
	DMT_HD bool DMFTopFacetsFromTopQ(const Cubical3DT<I>& G, const uint32_t* ranks, const uint8_t* top,
		I x, I y, I z, uint8_t& g0_dir, uint8_t& hi2_dir) {
		const I q[3] = { x, y, z };
		I f[6][3];
		uint8_t fd[6];
		int nf = 0;
		for (int i = 0; i < 3; i++) {
			if (q[i] & 1) {
				for (int side = 0; side < 2; side++) {
					f[nf][0] = q[0]; f[nf][1] = q[1]; f[nf][2] = q[2];
					f[nf][i] += side == 0 ? (I)-1 : (I)1;
					fd[nf] = (uint8_t)(side == 0 ? 2 * i + 1 : 2 * i);
					nf++;
				}
			}
		}
		if (nf == 0) return false;
		uint32_t fk[6][8];
		int flen = 1;
		for (int i = 0; i < nf; i++) flen = DMFKeyFromTopQ(G, ranks, top, f[i][0], f[i][1], f[i][2], fk[i]);
		int a = 0, b = -1;
		for (int i = 1; i < nf; i++) {
			if (DMFKeyCmp(fk[i], fk[a], flen) > 0) { b = a; a = i; }
			else if (b < 0 || DMFKeyCmp(fk[i], fk[b], flen) > 0) b = i;
		}
		g0_dir = fd[a];
		hi2_dir = b < 0 ? (uint8_t)7 : fd[b];
		return true;
	}

	// DMFLowestSelectingCofacetFromTop by coordinates: returns the direction
	// code from the cell to the selected cofacet, or -1. Cofacets in the same
	// order as Cubical3D::Cofacets. A cofacet selects the cell iff its stored
	// direction is the step back to the cell (the id test Neighbor(c, dir) ==
	// a of the id variant; the two agree because stored directions only ever
	// point along axes with odd coordinates, so no two codes alias).
	template <class I>
	DMT_HD int DMFLowestSelectingCofacetFromTopQ(const Cubical3DT<I>& G, const uint32_t* ranks,
		const uint8_t* top, I x, I y, I z, int which) {
		const I q[3] = { x, y, z };
		int best = -1;
		uint32_t bk[8], k[8];
		for (int i = 0; i < 3; i++) {
			if (q[i] & 1) continue;
			for (int side = 0; side < 2; side++) {
				I c[3] = { q[0], q[1], q[2] };
				if (side == 0) { if (q[i] - 1 < 0) continue; c[i] -= 1; }
				else { if (q[i] + 1 > G.Top(i)) continue; c[i] += 1; }
				const int to_c = side == 0 ? 2 * i + 1 : 2 * i;   // a -> c
				const uint8_t t = top[G.Id(c[0], c[1], c[2])];
				const int dir = which == 0 ? (int)(t & 15u) : (int)(t >> 4);
				if (dir != (to_c ^ 1)) continue;                   // c -> a
				const int len = DMFKeyFromTopQ(G, ranks, top, c[0], c[1], c[2], k);
				if (best < 0 || DMFKeyCmp(k, bk, len) < 0) {
					best = to_c;
					for (int j = 0; j < len; j++) bk[j] = k[j];
				}
			}
		}
		return best;
	}

} // namespace dmt

#endif
