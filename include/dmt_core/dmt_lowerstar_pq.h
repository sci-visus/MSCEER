// This header was split into dmt_lowerstar.h (star scaffold), dmt_pls.h (PLS)
// and dmt_pls_br.h (PLS-BR). The old path and the names below are kept so that
// existing code keeps compiling; new code includes the header of the rule it
// uses.
#ifndef GI_FORWARD_DMT_CORE_DMT_LOWERSTAR_PQ_H
#define GI_FORWARD_DMT_CORE_DMT_LOWERSTAR_PQ_H

#include "dmt_pls.h"
#include "dmt_pls_br.h"

namespace dmt {

	// PLSStar / PLSBRStar with aligned_first = true.
	DMT_HD void PLSAlignedStar(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		StarResult& R) {
		PLSStar(T, key, vvalid, R, true);
	}

	DMT_HD void PLSBRAlignedStar(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		uint32_t bnd, StarResult& R) {
		PLSBRStar(T, key, vvalid, bnd, R, true);
	}

	// PLSStarCodes or PLSBRStarCodes, chosen by boundary_restricted.
	DMT_HD void LowerStarCodes(const Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
		bool boundary_restricted, uint32_t bnd, uint32_t& L, StarCodeRecorder& rec) {
		if (boundary_restricted) PLSBRStarCodes(T, key, vvalid, bnd, L, rec);
		else PLSStarCodes(T, key, vvalid, L, rec);
	}

} // namespace dmt

#endif
