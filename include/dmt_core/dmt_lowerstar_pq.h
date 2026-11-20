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
#if defined(__CUDA_ARCH__)
		return __popc(m);
#else
		int c = 0;
		while (m) { m &= m - 1; c++; }
		return c;
#endif
	}
	DMT_HD int LowBit(uint32_t m) {
#if defined(__CUDA_ARCH__)
		return __ffs((int)m) - 1;
#else
		int i = 0;
		while (((m >> i) & 1u) == 0) i++;
		return i;
#endif
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
		DMT_UNROLL
		for (int dz = -1; dz <= 1; dz++)
			DMT_UNROLL
			for (int dy = -1; dy <= 1; dy++)
				DMT_UNROLL
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
		DMT_UNROLL
		for (int j = 0; j < 27; j++)
			if (j != kStarCentre && ((vvalid >> j) & 1u) && key[j] < kc) lowerv |= 1u << j;
		const uint32_t allowed = lowerv | (1u << kStarCentre);
		L = 0;
		DMT_UNROLL
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

	// Rank-mask keys: the same order as RobinsKeys in 32 bits, with no sort.
	// Give each vertex of the lower star its local rank r (0 = lowest; the
	// centre is above all of them) and let a cell's key be the set of its
	// vertices' ranks as a bitmask. For two distinct vertex sets A and B the
	// Robins comparison (descending sorted values, lexicographic, a proper
	// prefix smaller) is decided by the largest element of A xor B: the set
	// holding it is the larger one. That element is the highest set bit of
	// mask(A) ^ mask(B), so comparing the masks as integers gives exactly the
	// Robins order. Every index below is static once the loops are unrolled,
	// so on the GPU the keys and ranks live in registers.
	//
	// The vertex set of a cell is its far corner (vertex slot s for cell slot
	// s) united with the vertex sets of its facets in the star, so the masks
	// are built in increasing dimension from the canonical star topology
	// (BuildStar27TablesCanonical; gi_dmt_star27_mesh_tables.h checks mesh-
	// built tables against it).
	DMT_HD void RobinsMaskKeys(const uint32_t key[27], uint32_t lowerv, uint32_t rk[27]) {
		uint32_t less[27];
		DMT_UNROLL
		for (int j = 0; j < 27; j++) less[j] = 0;
		DMT_UNROLL
		for (int j = 0; j < 27; j++) {
			DMT_UNROLL
			for (int i = j + 1; i < 27; i++) {
				// the centre's rank is fixed below; its pairs are not needed
				if (i == kStarCentre || j == kStarCentre) continue;
				// only read for pairs of lower vertices, whose keys are distinct
				const bool c = key[i] < key[j];
				less[j] |= (uint32_t)c << i;
				less[i] |= (uint32_t)!c << j;
			}
		}
		uint32_t rb[27];
		DMT_UNROLL
		for (int j = 0; j < 27; j++) rb[j] = 1u << (Popc27(less[j] & lowerv) & 31);
		rb[kStarCentre] = 1u << Popc27(lowerv);
		DMT_UNROLL
		for (int s = 0; s < 27; s++) rk[s] = 0;
		rk[kStarCentre] = rb[kStarCentre];
		DMT_UNROLL
		for (int d = 1; d <= 3; d++) {
			DMT_UNROLL
			for (int s = 0; s < 27; s++) {
				const int dx = s % 3 - 1, dy = (s / 3) % 3 - 1, dz = s / 9 - 1;
				if ((dx != 0) + (dy != 0) + (dz != 0) != d) continue;
				uint32_t m = rb[s];
				if (dx != 0) m |= rk[s - dx];
				if (dy != 0) m |= rk[s - 3 * dy];
				if (dz != 0) m |= rk[s - 9 * dz];
				rk[s] = m;
			}
		}
	}

	// Argmin over the slots of m with every index static (keys stay in
	// registers on the GPU). Keys of distinct cells are distinct.
	template <class K>
	DMT_HD int ArgMinKeyT(const K rkey[27], uint32_t m) {
		if ((m & (m - 1)) == 0) return m ? LowBit(m) : -1;
		int best = -1;
		K bk = 0;
		DMT_UNROLL
		for (int s = 0; s < 27; s++)
			if (((m >> s) & 1u) && (best < 0 || rkey[s] < bk)) { best = s; bk = rkey[s]; }
		return best;
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

	// Where the decisions of a collapse go. StarResultRecorder fills a
	// StarResult (partner, step, part: the CPU drivers and the journals);
	// StarCodeRecorder keeps only the pair codes, as three bit planes over the
	// slots, so on the GPU the whole collapse stays in registers.
	struct StarResultRecorder {
		StarResult& R;
		DMT_HD explicit StarResultRecorder(StarResult& r) : R(r) {}
		DMT_HD void Pair(int a, int b, uint8_t step, uint8_t part) {
			R.pair[a] = (int8_t)b;
			R.pair[b] = (int8_t)a;
			R.step[a] = R.step[b] = step;
			R.part[a] = R.part[b] = part;
		}
		DMT_HD void Critical(int g, uint8_t step, uint8_t part) {
			R.pair[g] = -1;
			R.step[g] = step;
			R.part[g] = part;
		}
		DMT_HD void Orphan(int g, uint8_t part) {
			R.pair[g] = -1;
			R.step[g] = kOrphanStep;
			R.part[g] = part;
			R.n_orphans++;
		}
	};

	// Pair codes as StarDir gives them (0:+x 1:-x 2:+y 3:-y 4:+z 5:-z, 7 =
	// critical): slots of adjacent cells differ by +-1, +-3 or +-9 along x, y, z.
	struct StarCodeRecorder {
		uint32_t c0 = 0, c1 = 0, c2 = 0;
		uint32_t done = 0;   // slots a decision wrote (StarResult: pair != -2)
		int n_orphans = 0;
		DMT_HD void Set(int s, uint32_t code) {
			done |= 1u << s;
			c0 |= (code & 1u) << s;
			c1 |= ((code >> 1) & 1u) << s;
			c2 |= ((code >> 2) & 1u) << s;
		}
		DMT_HD void Pair(int a, int b, uint8_t, uint8_t) {
			const int d = b - a;
			const int ad = d < 0 ? -d : d;
			const uint32_t axis = ad == 1 ? 0u : (ad == 3 ? 1u : 2u);
			Set(a, 2u * axis + (d < 0 ? 1u : 0u));
			Set(b, 2u * axis + (d > 0 ? 1u : 0u));
		}
		DMT_HD void Critical(int g, uint8_t, uint8_t) { Set(g, 7u); }
		DMT_HD void Orphan(int g, uint8_t) { Set(g, 7u); n_orphans++; }
		DMT_HD uint32_t Code(int s) const {
			return ((c0 >> s) & 1u) | (((c1 >> s) & 1u) << 1) | (((c2 >> s) & 1u) << 2);
		}
	};

	// Collapse one part L of a lower star. A port of the reference
	// ProcessLowerStars kernel, decision for decision.
	template <class K, class Rec>
	DMT_HD void CollapsePartT(const Star27Tables& T, const K rkey[27], uint32_t L,
		bool seed_vertex, uint8_t part_id, Rec& rec) {
		uint32_t cls = 0;
		uint8_t step = 0;

#define DMT_NUN(a) Popc27(T.facmask[(a)] & L & ~cls)

		if (seed_vertex) {
			const uint32_t edges = L & T.dimmask[1];
			if (!edges) {
				rec.Critical(kStarCentre, step++, part_id);
				return;
			}
			const int delta = ArgMinKeyT(rkey, edges);
			rec.Pair(kStarCentre, delta, step++, part_id);
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
				const int a = ArgMinKeyT(rkey, pq1);
				pq1 &= ~(1u << a);
				if ((cls >> a) & 1u) continue;
				if (DMT_NUN(a) == 0) { pq0 |= 1u << a; continue; }
				const int t = LowBit(T.facmask[a] & L & ~cls);
				rec.Pair(t, a, step++, part_id);
				cls |= (1u << t) | (1u << a);
				pq0 &= ~(1u << t);
				// cofacets of a, then of t, that now miss exactly one facet
				uint32_t m = T.cofmask[a] & L & ~cls;
				while (m) {
					const int c = LowBit(m);
					m &= m - 1;
					if (DMT_NUN(c) == 1) pq1 |= 1u << c;
				}
				m = T.cofmask[t] & L & ~cls;
				while (m) {
					const int c = LowBit(m);
					m &= m - 1;
					if (DMT_NUN(c) == 1) pq1 |= 1u << c;
				}
			}
			pq0 &= ~cls;
			if (pq0) {
				const int g = ArgMinKeyT(rkey, pq0);
				pq0 &= ~(1u << g);
				rec.Critical(g, step++, part_id);
				cls |= 1u << g;
				uint32_t m = T.cofmask[g] & L & ~cls;
				while (m) {
					const int c = LowBit(m);
					m &= m - 1;
					if (DMT_NUN(c) == 1) pq1 |= 1u << c;
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
			rec.Orphan(s, part_id);
		}
	}

	DMT_HD void CollapsePart(const Star27Tables& T, const uint64_t rkey[27], uint32_t L,
		bool seed_vertex, uint8_t part_id, StarResult& R) {
		StarResultRecorder rec(R);
		CollapsePartT<uint64_t>(T, rkey, L, seed_vertex, part_id, rec);
	}

	// PLS on the star whose vertex keys are key[] (vvalid: slots in the domain).
	DMT_HD void PLSStar(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		StarResult& R) {
		ClearStarResult(R);
		uint32_t lowerv, L;
		StarLower(T, key, vvalid, lowerv, L);
		uint32_t rkey[27];
		RobinsMaskKeys(key, lowerv, rkey);
		StarResultRecorder rec(R);
		CollapsePartT<uint32_t>(T, rkey, L, true, 0, rec);
	}

	// PLS-BR: boundary part first, then the interior part; bnd marks the cell
	// slots on the domain boundary.
	DMT_HD void PLSBRStar(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		uint32_t bnd, StarResult& R) {
		ClearStarResult(R);
		uint32_t lowerv, L;
		StarLower(T, key, vvalid, lowerv, L);
		uint32_t rkey[27];
		RobinsMaskKeys(key, lowerv, rkey);
		StarResultRecorder rec(R);
		const uint32_t Lb = L & bnd, Li = L & ~bnd;
		if (Lb) CollapsePartT<uint32_t>(T, rkey, Lb, ((Lb >> kStarCentre) & 1u) != 0, 0, rec);
		if (Li) CollapsePartT<uint32_t>(T, rkey, Li, ((Li >> kStarCentre) & 1u) != 0, 1, rec);
	}

	// Pair codes only (no step / part journal): the same decisions as PLSStar
	// and PLSBRStar (bnd = 0 and boundary_restricted = false gives PLS).
	// L returns the lower-star slots; rec.done the slots a decision wrote
	// (all of L unless a part holding the vertex has no edge, which a true
	// boundary mask never produces); rec.Code(s) is the pair code of slot s.
	DMT_HD void LowerStarCodes(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		bool boundary_restricted, uint32_t bnd, uint32_t& L, StarCodeRecorder& rec) {
		uint32_t lowerv;
		StarLower(T, key, vvalid, lowerv, L);
		uint32_t rkey[27];
		RobinsMaskKeys(key, lowerv, rkey);
		if (!boundary_restricted) {
			CollapsePartT<uint32_t>(T, rkey, L, true, 0, rec);
			return;
		}
		const uint32_t Lb = L & bnd, Li = L & ~bnd;
		if (Lb) CollapsePartT<uint32_t>(T, rkey, Lb, ((Lb >> kStarCentre) & 1u) != 0, 0, rec);
		if (Li) CollapsePartT<uint32_t>(T, rkey, Li, ((Li >> kStarCentre) & 1u) != 0, 1, rec);
	}

	// The same two rules with the packed 64-bit keys of RobinsKeys (the
	// original formulation), kept as the reference the rank-mask path is
	// checked against star by star (dmt3d_tiled_validate).
	DMT_HD void PLSStarPacked(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		StarResult& R) {
		ClearStarResult(R);
		uint32_t lowerv, L;
		StarLower(T, key, vvalid, lowerv, L);
		uint64_t rkey[27];
		RobinsKeys(T, key, lowerv, L, rkey);
		CollapsePart(T, rkey, L, true, 0, R);
	}

	DMT_HD void PLSBRStarPacked(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
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
