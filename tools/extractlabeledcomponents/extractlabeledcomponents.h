// extract all labeled components in the 'filename' file as subregion volumes, and do likewise for the optional additional fields;
//	values without label are set to zero in the mask file, but the additional fields are kept intact;
//	all volumes must have the same dimensions
//	assumes label 0 is the background, and other labels are integers

#pragma once

#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif

#include <algorithm>
#include <assert.h>
#include <stdio.h>
#include <unordered_map>
#include <vector>

#include "gi_basic_types.h"
#include "gi_timing.h"
#include "gi_topological_explicit_mesh_function.h"

#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"

#include "gi_topological_utility_functions.h"
#include "gi_numeric_integrator_expanding_region_stop.h" // not where comparer should be
#include "gi_topological_regular_masked_grid.h"
#include "gi_topological_regular_masked_restricted_grid.h"

#include "gi_max_vertex_labeling.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_modified_robins.h"
#include "gi_morse_smale_complex_basic.h"

typedef GInt::RegularGrid3D GridType;

using namespace GInt;

struct AABB {
	Vec3l low, high;
};


struct LabeledComponent {
	int64_t label;
	float* mask_field;
	float* additional_field;
	int64_t offset[3];
	int64_t dims[3];
};


std::vector<struct LabeledComponent> extract_labeled_components(const float* mask_field, const float* additional_field, const int64_t* dims, float background_value)
{
	// set up structures to navigate grid, and load the 3d image
	GridType* underlying_grid = new GridType(Vec3l(dims[0], dims[1], dims[2]), Vec3b(0, 0, 0));

	// compute axis-aligned bounding boxes for each component
	std::unordered_map<INDEX_TYPE, AABB> aabbs;
	for (INDEX_TYPE i = 0; i < underlying_grid->NumElements(); i++) {
		if (mask_field[i] == background_value) {
			continue;
		}

		INDEX_TYPE id = (INDEX_TYPE)mask_field[i];
		assert((float)id == mask_field[i]); // no fractional labels allowed
		auto aabb = aabbs.find(id);
		auto position = underlying_grid->Coords(i);
		if (aabb == aabbs.end()) {
			aabbs[id] = {position, position};
		} else {
			aabb->second = {
				{std::min(position[0], aabb->second.low[0]), std::min(position[1], aabb->second.low[1]), std::min(position[2], aabb->second.low[2])},
				{std::max(position[0], aabb->second.high[0]), std::max(position[1], aabb->second.high[1]), std::max(position[2], aabb->second.high[2])},
			};
		}
	}

	std::vector<struct LabeledComponent> labeled_components;

	for (auto item : aabbs) {
		auto id = item.first;
		auto aabb = item.second;
#if !defined(NDEBUG)
		printf("AABB: id %lld, low (%lld, %lld, %lld), high (%lld, %lld, %lld)\n", item.first, aabb.low[0], aabb.low[1], aabb.low[2], aabb.high[0], aabb.high[1], aabb.high[2]);
#endif

		struct LabeledComponent component = {0};

		auto subgrid = GridType(aabb.high + Vec3l(1, 1, 1) - aabb.low, Vec3b(0, 0, 0));
		component.label = id;
		component.offset[0] = aabb.low[0];
		component.offset[1] = aabb.low[1];
		component.offset[2] = aabb.low[2];
		component.dims[0] = subgrid.XYZ()[0];
		component.dims[1] = subgrid.XYZ()[1];
		component.dims[2] = subgrid.XYZ()[2];

		component.mask_field = new float[subgrid.NumElements()];

		for (INDEX_TYPE z = aabb.low[2]; z <= aabb.high[2]; z++) {
			for (INDEX_TYPE y = aabb.low[1]; y <= aabb.high[1]; y++) {
				for (INDEX_TYPE x = aabb.low[0]; x <= aabb.high[0]; x++) {
					auto idx = underlying_grid->Index3d({x, y, z});
					if ((INDEX_TYPE)mask_field[idx] == id) {
						component.mask_field[subgrid.Index3d({x - aabb.low[0], y - aabb.low[1], z - aabb.low[2]})] = id;
					} else {
						component.mask_field[subgrid.Index3d({x - aabb.low[0], y - aabb.low[1], z - aabb.low[2]})] = background_value;
					}
				}
			}
		}

		if (additional_field != NULL) {
			auto subgrid = GridType(aabb.high + Vec3l(1, 1, 1) - aabb.low, Vec3b(0, 0, 0));
			component.additional_field = new float[subgrid.NumElements()];

			for (INDEX_TYPE z = aabb.low[2]; z <= aabb.high[2]; z++) {
				for (INDEX_TYPE y = aabb.low[1]; y <= aabb.high[1]; y++) {
					for (INDEX_TYPE x = aabb.low[0]; x <= aabb.high[0]; x++) {
						component.additional_field[subgrid.Index3d({x - aabb.low[0], y - aabb.low[1], z - aabb.low[2]})] = additional_field[underlying_grid->Index3d({x, y, z})];
					}
				}
			}
		}

		labeled_components.push_back(component);
	}
  
	return labeled_components;
}
