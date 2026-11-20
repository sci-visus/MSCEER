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
	// Pattern-templated variants for the tiled GPU kernels (any index width,
	// Cubical3DT<I>). P is the parity pattern of the cell (bit i set = odd
	// coordinate on axis i), so its dimension and its odd axes are known at
	// compile time and every walk below unrolls into straight-line code with
	// static indices: on the GPU the walks, keys and comparisons stay in
	// registers. Cells are carried as (cell id, corner-vertex id) pairs, the
	// corner vertex being floor(q / 2) per axis, so no id -> coordinate
	// division is ever needed: a step along axis i moves the cell id by the
	// cell stride and, when it lands on the upper vertex, the vertex id by
	// the vertex stride. Same results as the id-based functions above
	// (dmt3d_tiled_validate checks the tiled kernels built on these against
	// the in-core GPU driver and the CPU driver).
	// ---------------------------------------------------------------------

	constexpr int DMFPop3(int p) { return (p & 1) + ((p >> 1) & 1) + ((p >> 2) & 1); }

	template <class I>
	struct DMFStrides {
		I cs[3];   // cell-id strides
		I vs[3];   // vertex-id strides
	};

	template <class I>
	DMT_HD DMFStrides<I> MakeDMFStrides(const Cubical3DT<I>& G) {
		DMFStrides<I> S;
		S.cs[0] = 1; S.cs[1] = G.c[0]; S.cs[2] = G.c[0] * G.c[1];
		S.vs[0] = 1; S.vs[1] = G.n[0]; S.vs[2] = G.n[0] * G.n[1];
		return S;
	}

	template <class I>
	DMT_HD I DMFSelAxis(const I a[3], int axis) { return axis == 0 ? a[0] : (axis == 1 ? a[1] : a[2]); }

	// Flattened F_w of a D-cell (DMFKeyFromTop), as 2^D values in out[].
	template <int D, class I>
	DMT_HD void DMFKeyWalk(const DMFStrides<I>& S, const uint32_t* ranks, const uint8_t* top,
		I cid, I vb, uint32_t out[1 << D]) {
		I id[1 << D], vo[1 << D];
		id[0] = cid;
		vo[0] = vb;
		DMT_UNROLL
		for (int l = 0; l < D; l++) {
			// expand from the back so entry i becomes 2i (G0) and 2i+1 (G1)
			DMT_UNROLL
			for (int i = (1 << l) - 1; i >= 0; i--) {
				const int g0 = top[id[i]] & 15;
				const int axis = g0 >> 1;
				const bool neg = (g0 & 1) != 0;
				const I cst = DMFSelAxis(S.cs, axis), vst = DMFSelAxis(S.vs, axis);
				const I base = id[i], vbase = vo[i];
				id[2 * i] = neg ? base - cst : base + cst;
				id[2 * i + 1] = neg ? base + cst : base - cst;
				vo[2 * i] = neg ? vbase : vbase + vst;
				vo[2 * i + 1] = neg ? vbase + vst : vbase;
			}
		}
		DMT_UNROLL
		for (int k = 0; k < (1 << D); k++) out[k] = ranks[vo[k]];
	}

	template <int N>
	DMT_HD int DMFKeyCmpN(const uint32_t* a, const uint32_t* b) {
		DMT_UNROLL
		for (int i = 0; i < N; i++) {
			if (a[i] < b[i]) return -1;
			if (a[i] > b[i]) return 1;
		}
		return 0;
	}

	// DMFTopFacetsFromTop for a cell of pattern P (P != 0). Facets in the order
	// of Cubical3D::Facets: odd axes ascending, -1 before +1.
	template <int P, class I>
	DMT_HD void DMFTopFacetsP(const DMFStrides<I>& S, const uint32_t* ranks, const uint8_t* top,
		I cid, I vb, uint8_t& g0_dir, uint8_t& hi2_dir) {
		constexpr int D = DMFPop3(P);
		constexpr int FL = 1 << (D - 1);
		// running largest (ka, da) and second largest (kb, db): the same
		// a / b updates as DMFTopFacetsFromTop, with static copies instead of
		// indices into a key array
		uint32_t ka[FL], kb[FL], kf[FL];
		uint8_t da = 7, db = 7;
		bool first = true, hasb = false;
		DMT_UNROLL
		for (int i = 0; i < 3; i++) {
			if (!((P >> i) & 1)) continue;
			DMT_UNROLL
			for (int side = 0; side < 2; side++) {
				if (side == 0) DMFKeyWalk<D - 1>(S, ranks, top, cid - S.cs[i], vb, kf);
				else DMFKeyWalk<D - 1>(S, ranks, top, cid + S.cs[i], vb + S.vs[i], kf);
				const uint8_t df = (uint8_t)(side == 0 ? 2 * i + 1 : 2 * i);
				if (first) {
					DMT_UNROLL
					for (int j = 0; j < FL; j++) ka[j] = kf[j];
					da = df;
					first = false;
				} else if (DMFKeyCmpN<FL>(kf, ka) > 0) {
					DMT_UNROLL
					for (int j = 0; j < FL; j++) { kb[j] = ka[j]; ka[j] = kf[j]; }
					db = da;
					da = df;
					hasb = true;
				} else if (!hasb || DMFKeyCmpN<FL>(kf, kb) > 0) {
					DMT_UNROLL
					for (int j = 0; j < FL; j++) kb[j] = kf[j];
					db = df;
					hasb = true;
				}
			}
		}
		g0_dir = da;
		hi2_dir = hasb ? db : (uint8_t)7;
	}

	// DMFLowestSelectingCofacetFromTop for a cell of pattern P (dimension < 3)
	// at doubled coordinates q, selector W (0: g0, 1: hi2). Returns the
	// direction code from the cell to the selected cofacet, or -1. Cofacets in
	// the order of Cubical3D::Cofacets (even axes ascending, -1 before +1).
	// For W = 0 every candidate c has G0(c) = the cell itself, so the first
	// half of every candidate key is the same and the comparison only needs
	// the second halves, F(G1(c)), G1(c) being the facet of c opposite the
	// cell: keys of dimension D instead of D + 1, half the walk.
	template <int P, int W, class I>
	DMT_HD int DMFLowestSelectingP(const Cubical3DT<I>& G, const DMFStrides<I>& S, const uint32_t* ranks,
		const uint8_t* top, I cid, I vb, const I q[3]) {
		constexpr int D = DMFPop3(P);
		constexpr int KL = W == 0 ? (1 << D) : (1 << (D + 1));
		int best = -1;
		uint32_t bk[KL], k[KL];
		DMT_UNROLL
		for (int i = 0; i < 3; i++) {
			if ((P >> i) & 1) continue;
			DMT_UNROLL
			for (int side = 0; side < 2; side++) {
				I c, cvb;
				if (side == 0) {
					if (q[i] - 1 < 0) continue;
					c = cid - S.cs[i];
					cvb = vb - S.vs[i];
				} else {
					if (q[i] + 1 > G.Top(i)) continue;
					c = cid + S.cs[i];
					cvb = vb;
				}
				const int to_c = side == 0 ? 2 * i + 1 : 2 * i;   // cell -> c
				const uint8_t t = top[c];
				const int dir = W == 0 ? (int)(t & 15u) : (int)(t >> 4);
				if (dir != (to_c ^ 1)) continue;                  // c selects the cell
				if (W == 0) {
					const I opp = side == 0 ? c - S.cs[i] : c + S.cs[i];
					const I ovb = side == 0 ? cvb : cvb + S.vs[i];
					DMFKeyWalk<(W == 0 ? D : D + 1)>(S, ranks, top, opp, ovb, k);
				} else {
					DMFKeyWalk<(W == 0 ? D : D + 1)>(S, ranks, top, c, cvb, k);
				}
				if (best < 0 || DMFKeyCmpN<KL>(k, bk) < 0) {
					best = to_c;
					DMT_UNROLL
					for (int j = 0; j < KL; j++) bk[j] = k[j];
				}
			}
		}
		return best;
	}

} // namespace dmt

#endif
