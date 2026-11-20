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

	template <int N>
	DMT_HD int DMFKeyCmpN(const uint32_t* a, const uint32_t* b) {
		DMT_UNROLL
		for (int i = 0; i < N; i++) {
			if (a[i] < b[i]) return -1;
			if (a[i] > b[i]) return 1;
		}
		return 0;
	}

	// Keys from the cell's own corner ranks. A cell's key depends only on the
	// values of its vertices (F_w is defined recursively: G0 = the facet with
	// the largest key, G1 = the facet opposite it; DMFKey in dmt_dmf_cell.h),
	// so instead of walking top[] -- a chain of dependent byte reads -- the key
	// is computed from the 2^D corner ranks, loaded independently, by the same
	// recursion unrolled at compile time. v[c] is the rank of corner c (bit i
	// set = vertex on the upper side of odd axis i); B holds the corner bits
	// already fixed by the enclosing cells. Facets go in Cubical3D::Facets
	// order (odd axes ascending, lower side first), so the facet opposite
	// facet f is f ^ 1, and the first strict maximum wins as in
	// DMFKeyFromFacets (keys of one dimension are distinct anyway).
	template <int P, int B>
	struct DMFCornerKey {
		static constexpr int D = DMFPop3(P);
		static constexpr int FL = D > 0 ? (1 << (D - 1)) : 1;

		// the facets across axis I: lower (corner bit I = 0) and upper
		template <int I>
		static DMT_HD void Facets(const uint32_t v[8], uint32_t lo[FL], uint32_t hi[FL]) {
			constexpr bool has = ((P >> I) & 1) != 0;
			constexpr int PF = has ? (P & ~(1 << I)) : 0;   // a never-used leaf otherwise
			constexpr int BH = has ? (B | (1 << I)) : B;
			DMFCornerKey<PF, B>::Get(v, lo);
			DMFCornerKey<PF, BH>::Get(v, hi);
		}

		template <int I>
		static DMT_HD void Step(const uint32_t v[8], bool& first, uint32_t g0[FL], uint32_t g1[FL]) {
			if (!((P >> I) & 1)) return;
			uint32_t lo[FL], hi[FL];
			Facets<I>(v, lo, hi);
			if (first || DMFKeyCmpN<FL>(lo, g0) > 0) {
				DMT_UNROLL
				for (int j = 0; j < FL; j++) { g0[j] = lo[j]; g1[j] = hi[j]; }
			}
			first = false;
			if (DMFKeyCmpN<FL>(hi, g0) > 0) {
				DMT_UNROLL
				for (int j = 0; j < FL; j++) { g0[j] = hi[j]; g1[j] = lo[j]; }
			}
		}

		static DMT_HD void Get(const uint32_t v[8], uint32_t out[2 * FL]) {
			uint32_t g0[FL], g1[FL];
			bool first = true;
			Step<0>(v, first, g0, g1);
			Step<1>(v, first, g0, g1);
			Step<2>(v, first, g0, g1);
			DMT_UNROLL
			for (int j = 0; j < FL; j++) { out[j] = g0[j]; out[FL + j] = g1[j]; }
		}
	};

	template <int B>
	struct DMFCornerKey<0, B> {
		static constexpr int FL = 1;
		static DMT_HD void Get(const uint32_t v[8], uint32_t* out) { out[0] = v[B]; }
	};

	// Ranks of the corners of a cell of pattern P whose corner vertex
	// (floor(q / 2) per axis) has id vb: v[c] for every c whose bits lie in P.
	template <int P, class I>
	DMT_HD void DMFLoadCorners(const DMFStrides<I>& S, const uint32_t* ranks, I vb, uint32_t v[8]) {
		DMT_UNROLL
		for (int c = 0; c < 8; c++) {
			if (c & ~P) continue;
			const I off = ((c & 1) ? S.vs[0] : (I)0) + ((c & 2) ? S.vs[1] : (I)0) + ((c & 4) ? S.vs[2] : (I)0);
			v[c] = ranks[vb + off];
		}
	}

	// Key of the cell of pattern P with corner vertex vb.
	template <int P, class I>
	DMT_HD void DMFKeyFromCorners(const DMFStrides<I>& S, const uint32_t* ranks, I vb,
		uint32_t out[1 << DMFPop3(P)]) {
		uint32_t v[8];
		DMFLoadCorners<P>(S, ranks, vb, v);
		DMFCornerKey<P, 0>::Get(v, out);
	}

	// DMFTopFacetsFromTop for a cell of pattern P (P != 0): the same g0 / hi2
	// directions, with the facet keys from the cell's corner ranks (top[] and
	// the cell id are not needed; the signature matches the walk version).
	template <int P, class I>
	DMT_HD void DMFTopFacetsP(const DMFStrides<I>& S, const uint32_t* ranks, const uint8_t* top,
		I cid, I vb, uint8_t& g0_dir, uint8_t& hi2_dir) {
		(void)top;
		(void)cid;
		typedef DMFCornerKey<P, 0> CK;
		constexpr int FL = CK::FL;
		uint32_t v[8];
		DMFLoadCorners<P>(S, ranks, vb, v);
		// running largest (ka, da) and second largest (kb, db), the a / b
		// updates of DMFTopFacetsFromTop over the facets in Facets order
		uint32_t ka[FL], kb[FL], lo[FL], hi[FL];
		uint8_t da = 7, db = 7;
		bool first = true, hasb = false;
		auto consider = [&](const uint32_t* kf, uint8_t df) {
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
		};
		if (P & 1) { CK::template Facets<0>(v, lo, hi); consider(lo, 1); consider(hi, 0); }
		if (P & 2) { CK::template Facets<1>(v, lo, hi); consider(lo, 3); consider(hi, 2); }
		if (P & 4) { CK::template Facets<2>(v, lo, hi); consider(lo, 5); consider(hi, 4); }
		g0_dir = da;
		hi2_dir = hasb ? db : (uint8_t)7;
	}

	// One candidate cofacet of DMFLowestSelectingP, across even axis I.
	template <int P, int W, int I, class Idx>
	DMT_HD void DMFSelectStep(const DMFStrides<Idx>& S, const uint32_t* ranks, Idx vb, unsigned sel,
		int& best, uint32_t* bk) {
		constexpr int D = DMFPop3(P);
		constexpr int KL = W == 0 ? (1 << D) : (1 << (D + 1));
		constexpr int PC = ((P >> I) & 1) ? P : (P | (1 << I));   // c's pattern (P never has I here)
		if ((P >> I) & 1) return;
		DMT_UNROLL
		for (int side = 0; side < 2; side++) {
			const int to_c = side == 0 ? 2 * I + 1 : 2 * I;   // cell -> c
			if (!((sel >> to_c) & 1u)) continue;
			const Idx cvb = side == 0 ? vb - S.vs[I] : vb;     // corner vertex of c
			uint32_t k[KL];
			if (W == 0) {
				// G0(c) is the cell, so only F(G1(c)) differs between candidates:
				// G1(c) is the facet of c opposite the cell, of the cell's pattern
				DMFKeyFromCorners<P>(S, ranks, side == 0 ? cvb : cvb + S.vs[I], k);
			} else {
				DMFKeyFromCorners<PC>(S, ranks, cvb, k);
			}
			if (best < 0 || DMFKeyCmpN<KL>(k, bk) < 0) {
				best = to_c;
				DMT_UNROLL
				for (int j = 0; j < KL; j++) bk[j] = k[j];
			}
		}
	}

	// DMFLowestSelectingCofacetFromTop for a cell of pattern P (dimension < 3)
	// at doubled coordinates q, selector W (0: g0, 1: hi2). Returns the
	// direction code from the cell to the selected cofacet, or -1. Candidates
	// in the order of Cubical3D::Cofacets (even axes ascending, -1 before +1).
	template <int P, int W, class I>
	DMT_HD int DMFLowestSelectingP(const Cubical3DT<I>& G, const DMFStrides<I>& S, const uint32_t* ranks,
		const uint8_t* top, I cid, I vb, const I q[3]) {
		constexpr int D = DMFPop3(P);
		constexpr int KL = W == 0 ? (1 << D) : (1 << (D + 1));
		// First the selecting cofacets, as a mask of direction codes (cell -> c):
		// one top[] byte each. Every (d+1)-cell selects exactly one facet, so a
		// cell has about one candidate on average, and one candidate needs no key.
		unsigned sel = 0;
		DMT_UNROLL
		for (int i = 0; i < 3; i++) {
			if ((P >> i) & 1) continue;
			DMT_UNROLL
			for (int side = 0; side < 2; side++) {
				I c;
				if (side == 0) {
					if (q[i] - 1 < 0) continue;
					c = cid - S.cs[i];
				} else {
					if (q[i] + 1 > G.Top(i)) continue;
					c = cid + S.cs[i];
				}
				const int to_c = side == 0 ? 2 * i + 1 : 2 * i;
				const uint8_t t = top[c];
				const int dir = W == 0 ? (int)(t & 15u) : (int)(t >> 4);
				if (dir == (to_c ^ 1)) sel |= 1u << to_c;          // c selects the cell
			}
		}
		if (sel == 0) return -1;
		if ((sel & (sel - 1)) == 0) {
			int b = 0;
			while (!((sel >> b) & 1u)) b++;
			return b;
		}
		// two or more: the lowest key, keys from the corner ranks
		int best = -1;
		uint32_t bk[KL];
		DMFSelectStep<P, W, 0>(S, ranks, vb, sel, best, bk);
		DMFSelectStep<P, W, 1>(S, ranks, vb, sel, best, bk);
		DMFSelectStep<P, W, 2>(S, ranks, vb, sel, best, bk);
		return best;
	}

} // namespace dmt

#endif
