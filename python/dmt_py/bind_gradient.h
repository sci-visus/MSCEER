#ifndef DMT_PY_BIND_GRADIENT_H
#define DMT_PY_BIND_GRADIENT_H

// The C++ object behind msceer_dmt.Gradient3D, shared by every binding that
// reads or modifies a gradient.

#include "bind_common.h"
#include "dmt/gi_dmt_gradient_view.h"
#include "dmt/gi_dmt_vpath_counter.h"

#include <cstdint>
#include <memory>
#include <vector>

namespace dmt_py {

	struct PyGradient {
		Cubical3D G;
		std::vector<uint8_t> bytes;   // GradBitfield per cell, x fastest
		std::vector<uint8_t> step;    // optional journal: decision index within the lower star
		std::vector<uint8_t> part;    // optional journal: part of the lower star
		std::vector<uint8_t> stage;   // DMF: 1 = Algorithm 1 pair, 2 = Algorithm 2 pair, 0 = critical
		std::vector<int64_t> order;   // DMF (optional): position in the (dimension, key) order
		uint64_t version = 0;         // bumped by every change to the pairing
		std::shared_ptr<GInt::dmt::PairScorer> scorer;   // reused scratch for score()

		GInt::dmt::GradientView3D View() const {
			return GInt::dmt::GradientView3D(G, const_cast<uint8_t*>(bytes.data()));
		}
	};

} // namespace dmt_py

#endif
