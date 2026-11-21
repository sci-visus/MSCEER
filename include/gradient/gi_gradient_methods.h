#ifndef GI_GRADIENT_METHODS_H
#define GI_GRADIENT_METHODS_H

// The discrete gradient methods of this library, by ID.
//
// An ID names a pairing rule, that is, whatever determines the output bytes.
// It has the form <CORE>[_<BOUNDARY>][_<VARIANT>]:
//
//   PLS   ProcessLowerStars: the priority-queue lower-star rule of V. Robins,
//         P. J. Wood and A. P. Sheppard, "Theory and algorithms for constructing
//         discrete Morse complexes from grayscale digital images", IEEE TPAMI
//         2011. gradient/gi_pls_gradient.h, dmt_core/dmt_pls.h.
//   CLS   Conforming lower stars (A. Gyulassy): a lower-star gradient that
//         conforms to a cell labeling. It uses the lower-star decomposition of
//         Robins, Wood and Sheppard with a different pairing order. The
//         labeling is an input, not part of the ID: the boundary class of each
//         cell by default, restriction labels when they are given.
//         gradient/gi_cls_gradient.h, gradient/gi_cls_sliding_gradient.h.
//   DMF   The discrete Morse function assignment of N. Shivashankar and
//         V. Natarajan, "Parallel computation of 3D Morse-Smale complexes",
//         Computer Graphics Forum 2012. gradient/gi_dmf_gradient.h,
//         dmt_core/dmt_dmf.h. DMF1 is its first pass alone.
//
// Boundary tokens, on PLS only:
//   BR    binary split: cells on the domain boundary and interior cells are
//         collapsed separately. gradient/gi_pls_br_gradient.h.
//   BC    boundary classes: a pair must stay inside one of the mesh's boundary
//         classes (interior, face, box edge, corner).
//         gradient/gi_pls_bc_gradient.h.
//
// Variant token:
//   AF    aligned-first: the quads containing the vertex's edge are paired
//         first, then the usual order. In every API this is an aligned_first
//         argument of the base method, not a separate function.
//
// Words that describe an implementation (sliding, tiled, streamed, GPU) are
// never part of an ID.
//
// This header is self-contained (<cstdint> only) so device code may include it.

#include <cstdint>

namespace GInt {

	enum class GradientMethod : uint8_t {
		PLS = 0,
		PLS_BR,
		PLS_BC,
		PLS_AF,
		PLS_BR_AF,
		CLS,
		CLS_AF,
		DMF,
		DMF1
	};

	constexpr int kNumGradientMethods = 9;

	inline const char* MethodName(GradientMethod m) {
		switch (m) {
		case GradientMethod::PLS: return "PLS";
		case GradientMethod::PLS_BR: return "PLS_BR";
		case GradientMethod::PLS_BC: return "PLS_BC";
		case GradientMethod::PLS_AF: return "PLS_AF";
		case GradientMethod::PLS_BR_AF: return "PLS_BR_AF";
		case GradientMethod::CLS: return "CLS";
		case GradientMethod::CLS_AF: return "CLS_AF";
		case GradientMethod::DMF: return "DMF";
		case GradientMethod::DMF1: return "DMF1";
		}
		return "";
	}

	// Exact, case-sensitive match of an ID. Returns false and leaves out
	// unchanged when name is not an ID.
	inline bool ParseMethod(const char* name, GradientMethod& out) {
		if (!name) return false;
		for (int i = 0; i < kNumGradientMethods; i++) {
			const GradientMethod m = (GradientMethod)i;
			const char* a = MethodName(m);
			const char* b = name;
			while (*a && *a == *b) { a++; b++; }
			if (*a == 0 && *b == 0) { out = m; return true; }
		}
		return false;
	}

	inline bool IsAlignedFirst(GradientMethod m) {
		return m == GradientMethod::PLS_AF || m == GradientMethod::PLS_BR_AF || m == GradientMethod::CLS_AF;
	}

	// The method an aligned-first ID varies; every other ID maps to itself.
	inline GradientMethod BaseMethod(GradientMethod m) {
		switch (m) {
		case GradientMethod::PLS_AF: return GradientMethod::PLS;
		case GradientMethod::PLS_BR_AF: return GradientMethod::PLS_BR;
		case GradientMethod::CLS_AF: return GradientMethod::CLS;
		default: return m;
		}
	}

	// The ID of base with aligned_first applied. Returns false when
	// aligned_first is set and base has no aligned-first form.
	inline bool WithAlignedFirst(GradientMethod base, bool aligned_first, GradientMethod& out) {
		if (!aligned_first) { out = base; return true; }
		switch (base) {
		case GradientMethod::PLS: out = GradientMethod::PLS_AF; return true;
		case GradientMethod::PLS_BR: out = GradientMethod::PLS_BR_AF; return true;
		case GradientMethod::CLS: out = GradientMethod::CLS_AF; return true;
		default: return false;
		}
	}

} // namespace GInt

#endif
