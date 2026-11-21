#ifndef DMT_PLS_BR_H
#define DMT_PLS_BR_H

// PLS-BR: boundary-restricted PLS. The lower star is split into the cells on
// the domain boundary and the interior cells, and each part is collapsed on
// its own, seeding the vertex only in the part that contains it. No pair ever
// joins a boundary cell to an interior cell.
//
// The scaffold (lower-star masks, Robins keys, the collapse itself) is in
// dmt_lowerstar.h.
//
// aligned_first selects the aligned-first variant (see CollapsePartT). Only
// the part that holds the vertex is seeded, so only that part has a first
// generation to align.

#include "dmt_lowerstar.h"

namespace dmt {

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

	// PLS-BR: boundary part first, then the interior part; bnd marks the cell
	// slots on the domain boundary.
	DMT_HD void PLSBRStar(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		uint32_t bnd, StarResult& R, bool aligned_first = false) {
		ClearStarResult(R);
		uint32_t lowerv, L;
		StarLower(T, key, vvalid, lowerv, L);
		uint32_t rkey[27];
		RobinsMaskKeys(key, lowerv, rkey);
		StarResultRecorder rec(R);
		const uint32_t Lb = L & bnd, Li = L & ~bnd;
		if (Lb) CollapsePartT<uint32_t>(T, rkey, Lb, ((Lb >> kStarCentre) & 1u) != 0, 0, rec, aligned_first);
		if (Li) CollapsePartT<uint32_t>(T, rkey, Li, ((Li >> kStarCentre) & 1u) != 0, 1, rec, aligned_first);
	}

	// Pair codes only (no step / part journal): the same decisions as
	// PLSBRStar. L returns the lower-star slots; rec.done the slots a decision
	// wrote (all of L unless a part holding the vertex has no edge, which a
	// true boundary mask never produces); rec.Code(s) is the pair code of slot s.
	DMT_HD void PLSBRStarCodes(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		uint32_t bnd, uint32_t& L, StarCodeRecorder& rec, bool aligned_first = false) {
		uint32_t lowerv;
		StarLower(T, key, vvalid, lowerv, L);
		uint32_t rkey[27];
		RobinsMaskKeys(key, lowerv, rkey);
		const uint32_t Lb = L & bnd, Li = L & ~bnd;
		if (Lb) CollapsePartT<uint32_t>(T, rkey, Lb, ((Lb >> kStarCentre) & 1u) != 0, 0, rec, aligned_first);
		if (Li) CollapsePartT<uint32_t>(T, rkey, Li, ((Li >> kStarCentre) & 1u) != 0, 1, rec, aligned_first);
	}

	// The same rule with the packed 64-bit keys of RobinsKeys (the original
	// formulation), kept as the reference the rank-mask path is checked
	// against star by star (dmt3d_tiled_validate).
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
