// This header was split into gradient/gi_lowerstar_driver.h (shared drivers),
// gradient/gi_pls_gradient.h (PLS) and gradient/gi_pls_br_gradient.h (PLS-BR).
// The old path and the names below are kept so that existing code keeps
// compiling; new code includes the header of the method it uses.
#ifndef GI_FORWARD_GI_DMT_LOWERSTAR_GRADIENT_H
#define GI_FORWARD_GI_DMT_LOWERSTAR_GRADIENT_H

#include "gradient/gi_pls_gradient.h"
#include "gradient/gi_pls_br_gradient.h"

namespace GInt {
	namespace dmt {

		// PLS_ALIGNED and PLS_BR_ALIGNED are the aligned-first variants, now the
		// aligned_first argument of the PLS and PLS-BR entry points.
		enum class LowerStarRule { PLS = 0, PLS_BR = 1, PLS_ALIGNED = 2, PLS_BR_ALIGNED = 3 };

		inline bool RuleIsBoundaryRestricted(LowerStarRule r) {
			return r == LowerStarRule::PLS_BR || r == LowerStarRule::PLS_BR_ALIGNED;
		}
		inline bool RuleIsAligned(LowerStarRule r) {
			return r == LowerStarRule::PLS_ALIGNED || r == LowerStarRule::PLS_BR_ALIGNED;
		}

		// ComputePLSGradient3D or ComputePLSBRGradient3D, chosen by rule.
		inline bool ComputeLowerStarGradient3D(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
			LowerStarRule rule, uint8_t* out_grad, uint8_t* out_step, uint8_t* out_part,
			LowerStarStats* stats, int num_threads = 0) {
			return RuleIsBoundaryRestricted(rule)
				? ComputePLSBRGradient3D(ranks, X, Y, Z, out_grad, out_step, out_part, stats, num_threads, RuleIsAligned(rule))
				: ComputePLSGradient3D(ranks, X, Y, Z, out_grad, out_step, out_part, stats, num_threads, RuleIsAligned(rule));
		}

		// PLSGradient / PLSBRGradient constructed with aligned_first = true.
		template <class M, class V, class G>
		class PLSAlignedGradient : public PLSGradient<M, V, G> {
		public:
			PLSAlignedGradient(M* mesh, V* maxvl, G* grad, uint8_t* step_out = nullptr, uint8_t* part_out = nullptr)
				: PLSGradient<M, V, G>(mesh, maxvl, grad, step_out, part_out, true) {}
		};

		template <class M, class V, class G>
		class PLSBRAlignedGradient : public PLSBRGradient<M, V, G> {
		public:
			PLSBRAlignedGradient(M* mesh, V* maxvl, G* grad, uint8_t* step_out = nullptr, uint8_t* part_out = nullptr)
				: PLSBRGradient<M, V, G>(mesh, maxvl, grad, step_out, part_out, true) {}
		};

	} // namespace dmt
} // namespace GInt

#endif
