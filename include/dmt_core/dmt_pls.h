#ifndef DMT_PLS_H
#define DMT_PLS_H

// PLS: the priority-queue lower-star gradient of Robins, Wood and Sheppard
// ("Theory and algorithms for constructing discrete Morse complexes from
// grayscale digital images", IEEE TPAMI 2011, ProcessLowerStars), on the
// 27-slot star of one vertex.
//
// One part: the whole lower star, always seeded. The scaffold (lower-star
// masks, Robins keys, the collapse itself) is in dmt_lowerstar.h.
//
// aligned_first selects the aligned-first variant: the quads containing the
// seed edge are paired before the queues take over (see CollapsePartT).

#include "dmt_lowerstar.h"

namespace dmt {

	// PLS on the star whose vertex keys are key[] (vvalid: slots in the domain).
	DMT_HD void PLSStar(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		StarResult& R, bool aligned_first = false) {
		ClearStarResult(R);
		uint32_t lowerv, L;
		StarLower(T, key, vvalid, lowerv, L);
		uint32_t rkey[27];
		RobinsMaskKeys(key, lowerv, rkey);
		StarResultRecorder rec(R);
		CollapsePartT<uint32_t>(T, rkey, L, true, 0, rec, aligned_first);
	}

	// Pair codes only (no step / part journal): the same decisions as PLSStar.
	// L returns the lower-star slots; rec.done the slots a decision wrote;
	// rec.Code(s) is the pair code of slot s.
	DMT_HD void PLSStarCodes(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		uint32_t& L, StarCodeRecorder& rec, bool aligned_first = false) {
		uint32_t lowerv;
		StarLower(T, key, vvalid, lowerv, L);
		uint32_t rkey[27];
		RobinsMaskKeys(key, lowerv, rkey);
		CollapsePartT<uint32_t>(T, rkey, L, true, 0, rec, aligned_first);
	}

	// The same rule with the packed 64-bit keys of RobinsKeys (the original
	// formulation), kept as the reference the rank-mask path is checked
	// against star by star (dmt3d_tiled_validate).
	DMT_HD void PLSStarPacked(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		StarResult& R) {
		ClearStarResult(R);
		uint32_t lowerv, L;
		StarLower(T, key, vvalid, lowerv, L);
		uint64_t rkey[27];
		RobinsKeys(T, key, lowerv, L, rkey);
		CollapsePart(T, rkey, L, true, 0, R);
	}

} // namespace dmt

#endif
