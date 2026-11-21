#ifndef GI_PLS_GRADIENT_H
#define GI_PLS_GRADIENT_H

// PLS: host drivers for the priority-queue lower-star gradient of Robins,
// Wood and Sheppard (dmt_core/dmt_pls.h).
//
// Two entry points with identical results:
//   ComputePLSGradient3D   flat uint32 ranks in, GradBitfield bytes out. The
//                          GPU entry points have the same contract; no MSCEER
//                          mesh is built.
//   PLSGradient<...>       the same rule on an MSCEER mesh, max/min vertex
//                          labeling and gradient labeling (ComputeLowerStar per
//                          vertex and an OpenMP ComputePairing).
//
// aligned_first selects the aligned-first variant in both.

#include "gi_lowerstar_driver.h"
#include "dmt_core/dmt_pls.h"

namespace GInt {
	namespace dmt {

		struct PLSRule {
			static uint32_t Boundary(const ::dmt::Star27Tables&, const ::dmt::Cubical3D&,
				int64_t, int64_t, int64_t) {
				return 0u;
			}
			static void Star(const ::dmt::Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
				uint32_t /*bnd*/, ::dmt::StarResult& R, bool aligned_first) {
				::dmt::PLSStar(T, key, vvalid, R, aligned_first);
			}
		};

		template <class M, class V, class G> using PLSGradient = LowerStarPQGradient<M, V, G, PLSRule>;

		// ranks: one per vertex, x fastest, pairwise distinct (only their order
		// matters). out_grad: (2X-1)(2Y-1)(2Z-1) bytes, x fastest. out_step and
		// out_part are optional (nullptr) per-cell journal arrays. Returns true
		// iff every cell was assigned.
		inline bool ComputePLSGradient3D(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
			uint8_t* out_grad, uint8_t* out_step, uint8_t* out_part,
			LowerStarStats* stats, int num_threads = 0, bool aligned_first = false) {
			return LowerStarGradient3D<PLSRule>(ranks, X, Y, Z, aligned_first,
				out_grad, out_step, out_part, stats, num_threads);
		}

	} // namespace dmt
} // namespace GInt

#endif
