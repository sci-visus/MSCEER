#ifndef NEIGHBORHOOD_EXTRACTOR_H	
#define NEIGHBORHOOD_EXTRACTOR_H

#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_topological_regular_grid.h"
#include <vector>

namespace GInt {

	class VertexNeighborhoodState {
	public:
		class LabelGradState {
			BYTE_TYPE grad[27];
			char label[27];
		};
		INDEX_TYPE baseid;

		float vals[9];
		char lstar[27];
		std::vector<LabelGradState> states;

	
	};


};

#endif
