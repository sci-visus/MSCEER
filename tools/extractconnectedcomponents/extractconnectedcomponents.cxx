// extract all connected components above the specified 'threshold' in the 'filename' file as a subregion volumes, and do likewise for the optional additional fields;
//	values below the threshold are set to zero in the 'filename' file, but the additional fields are kept intact;
//	all volumes must have the same dimensions


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


int main(int argc, char** argv) {

	// READ IN THE COMMAND LINE ARGUMENTS
	INDEX_TYPE X, Y, Z;
	if (argc < 6) { printf("Usage: X Y Z filename threshold [additional_fields]\n"); return 0; }
	sscanf(argv[1], "%lld", &X);
	sscanf(argv[2], "%lld", &Y);
	sscanf(argv[3], "%lld", &Z);
	char const* filename = argv[4];
	float threshold;
	sscanf(argv[5], "%f", &threshold);
  
	float* mask_field = new float[X*Y*Z];
	FILE* fin = fopen(filename, "rb");
	assert(fin != NULL);
	fread(mask_field, sizeof(float), X*Y*Z, fin);
	fclose(fin);
  
	// set up structures to navigate grid, and load the 3d image
	GridType* underlying_grid = new GridType(Vec3l(X, Y, Z), Vec3b(0, 0, 0));
	
	printf("loaded field function\n");
  
	VolumeConnectedComponents cc(underlying_grid);
	DenseLabeling<char> *maskvol = new DenseLabeling<char>(underlying_grid->NumElements());
  
	for (INDEX_TYPE i=0; i< underlying_grid->NumElements();i++) {
		maskvol->SetLabel(i, mask_field[i] > threshold);
	}

	cc.PerformUnionFind(maskvol);

	// compute axis-aligned bounding boxes for each component
	std::unordered_map<INDEX_TYPE, AABB> aabbs;
	for (INDEX_TYPE i = 0; i < underlying_grid->NumElements(); i++) {
		if (maskvol->GetLabel(i) == 0) {
			continue;
		}

		INDEX_TYPE id = cc.Find(i);
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


	// expand the axis-aligned bounding box by 1 voxel on each axis if not on domain boundary to get correct isosurface
	//	for each aabb in the subsequent analysis
	for (auto& item : aabbs) {
		auto& aabb = item.second;
		for (int i = 0; i < 3; i++) {
		  INDEX_TYPE ind_type = 0;
			aabb.low[i] = std::max(aabb.low[i] - 1,ind_type);
			aabb.high[i] = std::min(aabb.high[i] + 1, underlying_grid->XYZ()[i] - 1);
		}
	}


	// save mask field
	for (auto item : aabbs) {
		auto id = item.first;
		auto aabb = item.second;
		printf("AABB: id %lld, low (%lld, %lld, %lld), high (%lld, %lld, %lld)\n", item.first, aabb.low[0], aabb.low[1], aabb.low[2], aabb.high[0], aabb.high[1], aabb.high[2]);

		// set to zero all values in the aabb that do not have aabb's label
		auto subgrid = GridType(aabb.high + Vec3l(1, 1, 1) - aabb.low, Vec3b(0, 0, 0));
		float* tmp = new float[subgrid.NumElements()];

		for (INDEX_TYPE z = aabb.low[2]; z <= aabb.high[2]; z++) {
			for (INDEX_TYPE y = aabb.low[1]; y <= aabb.high[1]; y++) {
				for (INDEX_TYPE x = aabb.low[0]; x <= aabb.high[0]; x++) {
					auto idx = underlying_grid->Index3d({x, y, z});
					if (maskvol->GetLabel(idx) != 0 && cc.Find(idx) == id) {
						tmp[subgrid.Index3d({x - aabb.low[0], y - aabb.low[1], z - aabb.low[2]})] = mask_field[idx];
					} else {
						tmp[subgrid.Index3d({x - aabb.low[0], y - aabb.low[1], z - aabb.low[2]})] = 0;
					}					
				}
			}
		}

		char filename[2048];
		auto ret = snprintf(filename, sizeof filename, "label%d_mask_position%lldx%lldx%lld_%lldx%lldx%lld_float32.raw", (int)mask_field[id], aabb.low[0], aabb.low[1], aabb.low[2], subgrid.XYZ()[0], subgrid.XYZ()[1], subgrid.XYZ()[2]);
		assert(ret < sizeof filename);
		FILE *fp = fopen(filename, "wb");
		fwrite(tmp, sizeof *tmp, subgrid.NumElements(), fp);
		fclose(fp);

		delete[] tmp;
	}


	// save additional fields
	float* additional_field = new float[X*Y*Z];
	for (int i = 6; i < argc; i++) {
		FILE* fp = fopen(argv[i], "rb");
		assert(fp != NULL);
		fread(additional_field, sizeof *additional_field, underlying_grid->NumElements(), fp);
		fclose(fp);

		for (auto item : aabbs) {
			auto aabb = item.second;			

			auto subgrid = GridType(aabb.high + Vec3l(1, 1, 1) - aabb.low, Vec3b(0, 0, 0));
			float* tmp = new float[subgrid.NumElements()];

			for (INDEX_TYPE z = aabb.low[2]; z <= aabb.high[2]; z++) {
				for (INDEX_TYPE y = aabb.low[1]; y <= aabb.high[1]; y++) {
					for (INDEX_TYPE x = aabb.low[0]; x <= aabb.high[0]; x++) {
						tmp[subgrid.Index3d({x - aabb.low[0], y - aabb.low[1], z - aabb.low[2]})] = additional_field[underlying_grid->Index3d({x, y, z})];
					}
				}
			}

			char filename[2048];
			auto ret = snprintf(filename, sizeof filename, "label%d_field%d_position%lldx%lldx%lld_%lldx%lldx%lld_float32.raw", (int)mask_field[item.first], i - 6, aabb.low[0], aabb.low[1], aabb.low[2], subgrid.XYZ()[0], subgrid.XYZ()[1], subgrid.XYZ()[2]);
			assert(ret < sizeof filename);
			FILE *fp = fopen(filename, "wb");
			fwrite(tmp, sizeof *tmp, subgrid.NumElements(), fp);
			fclose(fp);

			delete[] tmp;
		}
	}
	delete[] additional_field;
	delete[] mask_field;
  
	return 0;
}


