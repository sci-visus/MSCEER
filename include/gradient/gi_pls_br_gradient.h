#ifndef GI_PLS_BR_GRADIENT_H
#define GI_PLS_BR_GRADIENT_H

// PLS-BR: host drivers for the boundary-restricted priority-queue lower-star
// gradient (dmt_core/dmt_pls_br.h). Cells on the domain boundary and interior
// cells are collapsed separately, so no pair joins the two.
//
// Two entry points with identical results:
//   ComputePLSBRGradient3D   flat uint32 ranks in, GradBitfield bytes out. The
//                            GPU entry points have the same contract; no
//                            MSCEER mesh is built.
//   PLSBRGradient<...>       the same rule on an MSCEER mesh, max/min vertex
//                            labeling and gradient labeling.
//
// aligned_first selects the aligned-first variant in both.

#include "gi_lowerstar_driver.h"
#include "dmt_core/dmt_pls_br.h"

namespace GInt {
	namespace dmt {

		struct PLSBRRule {
			static uint32_t Boundary(const ::dmt::Star27Tables& T, const ::dmt::Cubical3D& G,
				int64_t x, int64_t y, int64_t z) {
				return ::dmt::StarBoundaryMask(T, G, x, y, z);
			}
			static void Star(const ::dmt::Star27Tables& T, const uint32_t key[27], uint32_t vvalid,
				uint32_t bnd, ::dmt::StarResult& R, bool aligned_first) {
				::dmt::PLSBRStar(T, key, vvalid, bnd, R, aligned_first);
			}
		};

		template <class M, class V, class G> using PLSBRGradient = LowerStarPQGradient<M, V, G, PLSBRRule>;

		// Same contract as ComputePLSGradient3D (gi_pls_gradient.h).
		inline bool ComputePLSBRGradient3D(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
			uint8_t* out_grad, uint8_t* out_step, uint8_t* out_part,
			LowerStarStats* stats, int num_threads = 0, bool aligned_first = false) {
			return LowerStarGradient3D<PLSBRRule>(ranks, X, Y, Z, aligned_first,
				out_grad, out_step, out_part, stats, num_threads);
		}

	} // namespace dmt
} // namespace GInt

#endif
