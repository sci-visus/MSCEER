#ifndef DMT_LOWERSTAR_PQ_H
#define DMT_LOWERSTAR_PQ_H

// Priority-queue lower-star gradient (Robins, Wood, Sheppard, "Theory and
// algorithms for constructing discrete Morse complexes from grayscale digital
// images", IEEE TPAMI 2011, ProcessLowerStars), on the 27-slot star of one
// vertex, with no heap and no recursion.
//
// The lower star of v is the set of cells whose highest vertex is v. Cells are
// ordered by their Robins key: the values of their vertices sorted in
// decreasing order, compared lexicographically, a proper prefix comparing
// smaller. Values must be pairwise distinct; then every choice below is a
// unique minimum and the result does not depend on iteration order.
//
//   seed:  pair v with its steepest (lowest-key) edge, or make v critical
//          when the star has no edge;
//   pq1:   cells with exactly one unclassified facet in the part;
//   pq0:   cells with none.
//   Repeatedly pop the lowest-key cell of pq1 and pair it with its one
//   unclassified facet; when pq1 is empty, make the lowest-key cell of pq0
//   critical. After each step, cofacets that now have exactly one
//   unclassified facet join pq1.
//
// Two rules are provided:
//   PLS     one part: the whole lower star, always seeded.
//   PLS-BR  boundary-restricted: the lower star is split into the cells on the
//           domain boundary and the interior cells, and each part is collapsed
//           on its own, seeding v only in the part that contains it. No pair
//           ever joins a boundary cell to an interior cell.
//
// Keys are local: within one star only the relative order of the (at most
// 27) vertices matters, so each is replaced by a 5-bit local rank and a
// cell's key packs its sorted local ranks (+1, 0 for absent) into one
// integer, most significant first. Integer comparison of packed keys is then
// exactly lexicographic comparison with the shorter-prefix-first rule.

#include "dmt_hd.h"
#include "dmt_cubical3d.h"
#include "dmt_star27_tables.h"

namespace dmt {

	// Result for one star. pair[s] is the partner slot, -1 for a critical
	// cell and -2 for a slot not in the lower star. step[s] is the index of
	// the decision that classified the cell within its part (a pair is one
	// decision; both cells share it), kOrphanStep for a cell no rule reached.
	// part[s] is 0 for PLS and for the boundary part of PLS-BR, 1 for the
	// interior part.
	constexpr uint8_t kOrphanStep = 255;

	struct StarResult {
		int8_t pair[27];
		uint8_t step[27];
		uint8_t part[27];
		uint8_t n_orphans;
	};

	DMT_HD int Popc27(uint32_t m) {
		int c = 0;
		while (m) { m &= m - 1; c++; }
		return c;
	}
	DMT_HD int LowBit(uint32_t m) {
		int i = 0;
		while (((m >> i) & 1u) == 0) i++;
		return i;
	}

	DMT_HD void ClearStarResult(StarResult& R) {
		for (int s = 0; s < 27; s++) {
			R.pair[s] = -2;
			R.step[s] = 0;
			R.part[s] = 0;
		}
		R.n_orphans = 0;
	}

	// Ranks of the 27 vertices around vertex (x, y, z); vvalid has bit j set
	// for each vertex slot inside the domain.
	DMT_HD void GatherStarKeys(const Cubical3D& G, const uint32_t* ranks,
		int64_t x, int64_t y, int64_t z, uint32_t key[27], uint32_t& vvalid) {
		vvalid = 0;
		for (int dz = -1; dz <= 1; dz++)
			for (int dy = -1; dy <= 1; dy++)
				for (int dx = -1; dx <= 1; dx++) {
					const int j = StarSlot(dx, dy, dz);
					const int64_t X = x + dx, Y = y + dy, Z = z + dz;
					if (X < 0 || Y < 0 || Z < 0 || X >= G.n[0] || Y >= G.n[1] || Z >= G.n[2]) {
						key[j] = 0;
						continue;
					}
					key[j] = ranks[X + G.n[0] * (Y + G.n[1] * Z)];
					vvalid |= 1u << j;
				}
	}

	// Cell slots of the star of vertex (x, y, z) that lie on the domain boundary.
	DMT_HD uint32_t StarBoundaryMask(const Star27Tables& T, const Cubical3D& G,
		int64_t x, int64_t y, int64_t z) {
		uint32_t m = 0;
		const int64_t v[3] = { 2 * x, 2 * y, 2 * z };
		for (int s = 0; s < 27; s++) {
			for (int i = 0; i < 3; i++) {
				const int64_t q = v[i] + T.off[s][i];
				if (q == 0 || q == G.Top(i)) { m |= 1u << s; break; }
			}
		}
		return m;
	}

	// lowerv: valid vertex slots below the centre; L: cell slots of the lower
	// star (inside the domain, every vertex either the centre or below it).
	DMT_HD void StarLower(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		uint32_t& lowerv, uint32_t& L) {
		lowerv = 0;
		const uint32_t kc = key[kStarCentre];
		for (int j = 0; j < 27; j++)
			if (j != kStarCentre && ((vvalid >> j) & 1u) && key[j] < kc) lowerv |= 1u << j;
		const uint32_t allowed = lowerv | (1u << kStarCentre);
		L = 0;
		for (int s = 0; s < 27; s++)
			if ((T.vertmask[s] & ~allowed) == 0) L |= 1u << s;
	}

	// Packed Robins keys of the cells in L.
	DMT_HD void RobinsKeys(const Star27Tables& T, const uint32_t key[27], uint32_t lowerv,
		uint32_t L, uint64_t rkey[27]) {
		uint8_t lrank[27];
		for (int j = 0; j < 27; j++) lrank[j] = 0;
		for (int j = 0; j < 27; j++) {
			if (!((lowerv >> j) & 1u)) continue;
			int r = 0;
			for (int i = 0; i < 27; i++)
				if (((lowerv >> i) & 1u) && key[i] < key[j]) r++;
			lrank[j] = (uint8_t)r;
		}
		lrank[kStarCentre] = (uint8_t)Popc27(lowerv);
		for (int s = 0; s < 27; s++) {
			rkey[s] = 0;
			if (!((L >> s) & 1u)) continue;
			uint8_t r[8];
			const int nv = T.nvert[s];
			for (int i = 0; i < nv; i++) r[i] = lrank[T.vert[s][i]];
			for (int i = 1; i < nv; i++) {          // insertion sort, descending
				uint8_t x = r[i];
				int k = i - 1;
				while (k >= 0 && r[k] < x) { r[k + 1] = r[k]; k--; }
				r[k + 1] = x;
			}
			uint64_t p = 0;
			for (int i = 0; i < 8; i++) p = (p << 5) | (uint64_t)(i < nv ? r[i] + 1 : 0);
			rkey[s] = p;
		}
	}

	DMT_HD int ArgMinKey(const uint64_t rkey[27], uint32_t m) {
		int best = -1;
		uint64_t bk = 0;
		while (m) {
			const int s = LowBit(m);
			m &= m - 1;
			if (best < 0 || rkey[s] < bk) { best = s; bk = rkey[s]; }
		}
		return best;
	}

	// Collapse one part L of a lower star. A port of the reference
	// ProcessLowerStars kernel, decision for decision.
	DMT_HD void CollapsePart(const Star27Tables& T, const uint64_t rkey[27], uint32_t L,
		bool seed_vertex, uint8_t part_id, StarResult& R) {
		uint32_t cls = 0;
		uint8_t step = 0;

#define DMT_NUN(a) Popc27(T.facmask[(a)] & L & ~cls)

		if (seed_vertex) {
			const uint32_t edges = L & T.dimmask[1];
			if (!edges) {
				R.pair[kStarCentre] = -1;
				R.step[kStarCentre] = step++;
				R.part[kStarCentre] = part_id;
				return;
			}
			const int delta = ArgMinKey(rkey, edges);
			R.pair[kStarCentre] = (int8_t)delta;
			R.pair[delta] = (int8_t)kStarCentre;
			R.step[kStarCentre] = R.step[delta] = step++;
			R.part[kStarCentre] = R.part[delta] = part_id;
			cls |= (1u << kStarCentre) | (1u << delta);
		}

		uint32_t pq1 = 0, pq0 = 0;
		{
			uint32_t m = L & ~cls;
			while (m) {
				const int s = LowBit(m);
				m &= m - 1;
				const int k = DMT_NUN(s);
				if (k == 1) pq1 |= 1u << s;
				else if (k == 0) pq0 |= 1u << s;
			}
		}

		while (pq1 || pq0) {
			while (pq1) {
				const int a = ArgMinKey(rkey, pq1);
				pq1 &= ~(1u << a);
				if ((cls >> a) & 1u) continue;
				if (DMT_NUN(a) == 0) { pq0 |= 1u << a; continue; }
				const int t = LowBit(T.facmask[a] & L & ~cls);
				R.pair[t] = (int8_t)a;
				R.pair[a] = (int8_t)t;
				R.step[t] = R.step[a] = step++;
				R.part[t] = R.part[a] = part_id;
				cls |= (1u << t) | (1u << a);
				pq0 &= ~(1u << t);
				const int both[2] = { a, t };
				for (int w = 0; w < 2; w++) {
					uint32_t m = T.cofmask[both[w]] & L & ~cls;
					while (m) {
						const int b = LowBit(m);
						m &= m - 1;
						if (DMT_NUN(b) == 1) pq1 |= 1u << b;
					}
				}
			}
			pq0 &= ~cls;
			if (pq0) {
				const int g = ArgMinKey(rkey, pq0);
				pq0 &= ~(1u << g);
				R.pair[g] = -1;
				R.step[g] = step++;
				R.part[g] = part_id;
				cls |= 1u << g;
				uint32_t m = T.cofmask[g] & L & ~cls;
				while (m) {
					const int b = LowBit(m);
					m &= m - 1;
					if (DMT_NUN(b) == 1) pq1 |= 1u << b;
				}
			}
		}
#undef DMT_NUN

		// Cells no step reached are made critical, as the reference does in
		// its final sweep. On a lower star this does not happen.
		uint32_t left = L & ~cls;
		while (left) {
			const int s = LowBit(left);
			left &= left - 1;
			R.pair[s] = -1;
			R.step[s] = kOrphanStep;
			R.part[s] = part_id;
			R.n_orphans++;
		}
	}

	// PLS on the star whose vertex keys are key[] (vvalid: slots in the domain).
	DMT_HD void PLSStar(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		StarResult& R) {
		ClearStarResult(R);
		uint32_t lowerv, L;
		StarLower(T, key, vvalid, lowerv, L);
		uint64_t rkey[27];
		RobinsKeys(T, key, lowerv, L, rkey);
		CollapsePart(T, rkey, L, true, 0, R);
	}

	// PLS-BR: boundary part first, then the interior part; bnd marks the cell
	// slots on the domain boundary.
	DMT_HD void PLSBRStar(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		uint32_t bnd, StarResult& R) {
		ClearStarResult(R);
		uint32_t lowerv, L;
		StarLower(T, key, vvalid, lowerv, L);
		uint64_t rkey[27];
		RobinsKeys(T, key, lowerv, L, rkey);
		const uint32_t Lb = L & bnd, Li = L & ~bnd;
		if (Lb) CollapsePart(T, rkey, Lb, ((Lb >> kStarCentre) & 1u) != 0, 0, R);
		if (Li) CollapsePart(T, rkey, Li, ((Li >> kStarCentre) & 1u) != 0, 1, R);
	}

} // namespace dmt

#endif
