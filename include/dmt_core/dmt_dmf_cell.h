#ifndef DMT_DMF_CELL_H
#define DMT_DMF_CELL_H

// Per-cell pieces of the Shivashankar-Natarajan discrete gradient
// (N. Shivashankar, V. Natarajan, "Parallel computation of 3D Morse-Smale
// complexes", Computer Graphics Forum 31(3), 2012: ASSIGNGRADIENT and
// ASSIGNGRADIENT2).
//
// Every cell gets the recursive key
//     F_w(v)     = (f(v))                       for a vertex,
//     F_w(sigma) = (F_w(G0), F_w(G1))           otherwise,
// G0 the facet of sigma with the largest key and G1 the largest facet sharing
// no vertex with G0 -- on a cube, the facet opposite G0. Flattened, the key of
// a d-cell is a sequence of its 2^d vertex values; cells of one dimension are
// compared lexicographically on it, which is the order of the nested tuples.
// With pairwise distinct vertex values no two cells of one dimension have
// equal keys, so "largest" and "lowest" are always unique.
//
// These functions read the cell's own vertices only (keys are computed on the
// fly, no table of keys is stored), so they can run per cell on a GPU.

#include "dmt_hd.h"
#include "dmt_cubical3d.h"

namespace dmt {

	// Lexicographic comparison of two keys of length n: -1, 0 or 1.
	DMT_HD int DMFKeyCmp(const uint32_t* a, const uint32_t* b, int n) {
		for (int i = 0; i < n; i++) {
			if (a[i] < b[i]) return -1;
			if (a[i] > b[i]) return 1;
		}
		return 0;
	}

	DMT_HD int DMFKey0(const Cubical3D& G, const uint32_t* ranks, int64_t cell, uint32_t out[8]) {
		out[0] = ranks[G.CellVertex(cell)];
		return 1;
	}

	// Key of a cell of dimension >= 1 from the keys of its facets. Facets come
	// in opposite pairs per odd axis (lower, upper), so the facet opposite
	// facet i is facet i ^ 1.
	DMT_HD int DMFKeyFromFacets(const uint32_t fk[6][8], int nf, int flen, uint32_t out[8]) {
		int g0 = 0;
		for (int i = 1; i < nf; i++)
			if (DMFKeyCmp(fk[i], fk[g0], flen) > 0) g0 = i;
		const int g1 = g0 ^ 1;
		for (int i = 0; i < flen; i++) out[i] = fk[g0][i];
		for (int i = 0; i < flen; i++) out[flen + i] = fk[g1][i];
		return 2 * flen;
	}

	DMT_HD int DMFKey1(const Cubical3D& G, const uint32_t* ranks, int64_t cell, uint32_t out[8]) {
		int64_t f[6];
		const int nf = G.Facets(cell, f);
		uint32_t fk[6][8];
		for (int i = 0; i < nf; i++) DMFKey0(G, ranks, f[i], fk[i]);
		return DMFKeyFromFacets(fk, nf, 1, out);
	}

	DMT_HD int DMFKey2(const Cubical3D& G, const uint32_t* ranks, int64_t cell, uint32_t out[8]) {
		int64_t f[6];
		const int nf = G.Facets(cell, f);
		uint32_t fk[6][8];
		for (int i = 0; i < nf; i++) DMFKey1(G, ranks, f[i], fk[i]);
		return DMFKeyFromFacets(fk, nf, 2, out);
	}

	DMT_HD int DMFKey3(const Cubical3D& G, const uint32_t* ranks, int64_t cell, uint32_t out[8]) {
		int64_t f[6];
		const int nf = G.Facets(cell, f);
		uint32_t fk[6][8];
		for (int i = 0; i < nf; i++) DMFKey2(G, ranks, f[i], fk[i]);
		return DMFKeyFromFacets(fk, nf, 4, out);
	}

	// Flattened F_w of any cell; returns its length 2^dim.
	DMT_HD int DMFKey(const Cubical3D& G, const uint32_t* ranks, int64_t cell, uint32_t out[8]) {
		switch (G.Dim(cell)) {
		case 0: return DMFKey0(G, ranks, cell, out);
		case 1: return DMFKey1(G, ranks, cell, out);
		case 2: return DMFKey2(G, ranks, cell, out);
		default: return DMFKey3(G, ranks, cell, out);
		}
	}

	// The facets of `cell` with the largest and second-largest keys, as
	// direction codes from the cell (0..5); hi2 is 7 when the cell has a
	// single facet. Returns false for a vertex.
	DMT_HD bool DMFTopFacets(const Cubical3D& G, const uint32_t* ranks, int64_t cell,
		uint8_t& g0_dir, uint8_t& hi2_dir) {
		int64_t f[6];
		const int nf = G.Facets(cell, f);
		if (nf == 0) return false;
		uint32_t fk[6][8];
		int flen = 1;
		for (int i = 0; i < nf; i++) flen = DMFKey(G, ranks, f[i], fk[i]);
		int a = 0, b = -1;
		for (int i = 1; i < nf; i++) {
			if (DMFKeyCmp(fk[i], fk[a], flen) > 0) { b = a; a = i; }
			else if (b < 0 || DMFKeyCmp(fk[i], fk[b], flen) > 0) b = i;
		}
		g0_dir = (uint8_t)G.DirBetween(cell, f[a]);
		hi2_dir = b < 0 ? (uint8_t)7 : (uint8_t)G.DirBetween(cell, f[b]);
		return true;
	}

	// Among the cofacets b of `a` with sel(b) == a -- sel being the stored g0
	// or hi2 direction of b -- the one with the lowest key, or -1 if none.
	// top[] holds, per cell, (g0_dir | hi2_dir << 4); which = 0 selects g0,
	// 1 selects hi2.
	DMT_HD int64_t DMFLowestSelectingCofacet(const Cubical3D& G, const uint32_t* ranks,
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
			len = DMFKey(G, ranks, c[i], k);
			if (best < 0 || DMFKeyCmp(k, bk, len) < 0) {
				best = c[i];
				for (int j = 0; j < len; j++) bk[j] = k[j];
			}
		}
		return best;
	}

} // namespace dmt

#endif
