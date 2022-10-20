// TODO(2/24/2022): do not compute gradient on boundary for internal tiles
// TODO(2/23/2022): measure memory overhead during construction and check if we could reduce it
// TODO(2/23/2022): provide API to build the localized morse complex directly from a file (so the global data does not need to be in memory)


#if !defined(LOCALIZED_MORSE_COMPLEX)
#define LOCALIZED_MORSE_COMPLEX

#include <stdint.h>

#include <queue>
#include <stack>
#include <vector>


// avoid any printing inside the headers
#define printf

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
#include "gi_bifiltration_pairing.h"
#include "gi_fast_robins_noalloc.h"

#undef printf
#define USEMAXV

using namespace GInt;
typedef RegularGrid3D GridType;
typedef RegularGridTrilinearFunction GridFuncType;
typedef TopologicalRegularGrid3D MeshType;
typedef DiscreteGradientLabeling<MeshType> GradType;
typedef RegularGridMaxMinVertexLabeling3D<MeshType, GridFuncType> MaxVLType;
typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType, 5, 4> RobinsType;
typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> TopoFuncType;
typedef SlidingWindowRobinsNoalloc<GridType, GridFuncType, MeshType, MaxVLType, GradType> NewRobinsType;


int64_t const tile_x = 64;
int64_t const tile_y = 64;
int64_t const tile_z = 64;


typedef float data_t;
typedef int32_t local_t;


struct cell {
	local_t local_index;
	local_t tile_index;
};


struct tile {
	int64_t position[3]; // in terms of the global index space
	int64_t size[3];

	int64_t offset[3];
	int64_t bridge_set_offset[3];

	data_t *data;

	GridType* grid;
	GridFuncType* func;
	GradType* labeling;
	MeshType* tgrid;

	double temporary_allocations_time_seconds;
	double max_labeling_time_seconds;
	double robins_time_seconds;
	double total_time_seconds;
};

struct localized_morse_complex {
	struct tile *tiles;
	int64_t tile_count;
	int64_t dims[3];
};


// TODO(2/14/2022): data argument is not const due to the GridFuncType constructor argument not being const
void
compute_localized_morse_complex_from_grid(data_t *data, int64_t const *dims,
	struct localized_morse_complex *out_morse_complex)
{
	// TODO(1/21/2022): handle nontile shape multiple grids
	assert(dims[0]%tile_x == 0);
	assert(dims[1]%tile_y == 0);
	assert(dims[2]%tile_z == 0);

	// guarantees we can always construct 2 vertex wide ghost layer on +x,+y,+z side
	assert(tile_x > 1);
	assert(tile_y > 1);
	assert(tile_z > 1);

	out_morse_complex->dims[0] = dims[0]/tile_x;
	out_morse_complex->dims[1] = dims[1]/tile_y;
	out_morse_complex->dims[2] = dims[2]/tile_z;

	GridType *m_grid = new GridType(Vec3l(dims[0], dims[1], dims[2]), Vec3b(0, 0, 0));
	GridFuncType *m_func2 = new GridFuncType(m_grid, data);

	std::vector<struct tile> tiles;
	for (int64_t k = 0; k < dims[2]; k += tile_z) {
		for (int64_t j = 0; j < dims[1]; j += tile_y) {
			for (int64_t i = 0; i < dims[0]; i+= tile_x) {
				struct tile tile = {
					.position = {i, j, k},
					.size = {tile_x, tile_y, tile_z},
					.offset = {(i == 0) ? 0 : 1, (j == 0) ? 0 : 1, (k == 0) ? 0 : 1},
					.bridge_set_offset = {(i + tile_x >= dims[0]) ? 0 : 1, (j + tile_y >= dims[1]) ? 0 : 1, (k + tile_z >= dims[2]) ? 0 : 1},
				};

				int64_t ghost_low[3] = {i - tile.offset[0], j - tile.offset[1], k - tile.offset[2]}; // inclusive
				int64_t ghost_high[3] = {i + tile_x + 2*tile.bridge_set_offset[0], j + tile_y + 2*tile.bridge_set_offset[1], k + tile_z + 2*tile.bridge_set_offset[2]}; // noninclusive, +1 for bridge set, +1 for ghost

				tile.data = new data_t[(ghost_high[0] - ghost_low[0])*(ghost_high[1] - ghost_low[1])*(ghost_high[2] - ghost_low[2])];
				assert(tile.data != NULL);
				// read tile function from global function
				{
					int64_t offset = 0;
					for (int64_t k = ghost_low[2]; k < ghost_high[2]; k++) {
						for (int64_t j = ghost_low[1]; j < ghost_high[1]; j++) {
							for (int64_t i = ghost_low[0]; i < ghost_high[0]; i++) {
								tile.data[offset++] = m_func2->SampleImage({i, j, k});
							}
						}
					}
				}

				tile.grid = new GridType({ghost_high[0] - ghost_low[0], ghost_high[1] - ghost_low[1], ghost_high[2] - ghost_low[2]}, {0, 0, 0});
				tile.func = new GridFuncType(tile.grid, tile.data);
				tile.tgrid = new MeshType(tile.grid);

				tiles.push_back(tile);
			}
		}
	}

	delete m_func2;
	delete m_grid;

	auto total_t0 = std::chrono::high_resolution_clock::now();
	int tile_i;
#pragma omp parallel for schedule(dynamic)
	for (tile_i = 0; tile_i < (int)tiles.size(); tile_i++) {
		auto &tile = tiles[tile_i];

		int tid = omp_get_thread_num();

		// memory allocations not part of the Robins algorithm implementation
		{
			tile.labeling = new GradType(tile.tgrid);
			tile.labeling->ClearAllGradient();
		}

		// Robins algorithm implementation
		{
			auto total_t0 = std::chrono::high_resolution_clock::now();

			MaxVLType* maxv = new MaxVLType(tile.tgrid, tile.func);

			// TODO(1/21/2022): reuse memory allocation
			auto allocations_t0 = std::chrono::high_resolution_clock::now();		
			maxv->allocate();
			auto allocations_t1 = std::chrono::high_resolution_clock::now();

			auto max_labeling_t0 = std::chrono::high_resolution_clock::now();
			maxv->ComputeOutputSerial();
			auto max_labeling_t1 = std::chrono::high_resolution_clock::now();

			auto robins_t0 = std::chrono::high_resolution_clock::now();
			NewRobinsType* trobins =  new NewRobinsType(tile.grid, tile.func, tile.tgrid, maxv, tile.labeling);
			trobins->ComputePairing_slidingSerial();
			delete trobins;
			auto robins_t1 = std::chrono::high_resolution_clock::now();

			delete maxv;

			auto total_t1 = std::chrono::high_resolution_clock::now();

			tile.temporary_allocations_time_seconds = std::chrono::duration<double>(allocations_t1 - allocations_t0).count();
			tile.max_labeling_time_seconds = std::chrono::duration<double>(max_labeling_t1 - max_labeling_t0).count();
			tile.robins_time_seconds = std::chrono::duration<double>(robins_t1 - robins_t0).count();
			tile.total_time_seconds = std::chrono::duration<double>(total_t1 - total_t0).count();
		}
	}
	auto total_t1 = std::chrono::high_resolution_clock::now();

	// copy tiles to C array
	out_morse_complex->tiles = (struct tile *)malloc(tiles.size()*sizeof tiles[0]);
	out_morse_complex->tile_count = 0;
	for (auto tile : tiles) {
		out_morse_complex->tiles[out_morse_complex->tile_count++] = tile;
	}
}


// supports nonpower-of-two tiles
// NOTE(2/6/2022): buggy if the last tile is of width or height or depth of 1 (assertion on ghost_high catches the case)
void
compute_localized_morse_complex_from_grid_notilemultiple(data_t *data, int64_t const *dims,
	struct localized_morse_complex *out_morse_complex)
{
	// round up the number of tiles
	// TODO: maybe do remainder and if statement for better readability
	out_morse_complex->dims[0] = (dims[0] + tile_x - 1)/tile_x;
	out_morse_complex->dims[1] = (dims[1] + tile_y - 1)/tile_y;
	out_morse_complex->dims[2] = (dims[2] + tile_z - 1)/tile_z;

	GridType *m_grid = new GridType(Vec3l(dims[0], dims[1], dims[2]), Vec3b(0, 0, 0));
	GridFuncType *m_func2 = new GridFuncType(m_grid, data);

	std::vector<struct tile> tiles;
	for (int64_t k = 0; k < dims[2]; k += tile_z) {
		for (int64_t j = 0; j < dims[1]; j += tile_y) {
			for (int64_t i = 0; i < dims[0]; i+= tile_x) {
				struct tile tile = {
					.position = {i, j, k},
					.size = {(i + tile_x > dims[0]) ? (dims[0] - i) : tile_x, (j + tile_y > dims[1]) ? (dims[1] - j) : tile_y, (k + tile_z > dims[2]) ? (dims[2] - k) : tile_z},
					.offset = {(i == 0) ? 0 : 1, (j == 0) ? 0 : 1, (k == 0) ? 0 : 1},
					.bridge_set_offset = {(i + tile_x >= dims[0]) ? 0 : 1, (j + tile_y >= dims[1]) ? 0 : 1, (k + tile_z >= dims[2]) ? 0 : 1},
				};

				int64_t ghost_low[3] = {i - tile.offset[0], j - tile.offset[1], k - tile.offset[2]}; // inclusive
				int64_t ghost_high[3] = {i + tile.size[0] + 2*tile.bridge_set_offset[0], j + tile.size[1] + 2*tile.bridge_set_offset[1], k + tile.size[2] + 2*tile.bridge_set_offset[2]}; // noninclusive, +1 for bridge set, +1 for ghost

				assert(ghost_high[0] <= dims[0] && ghost_high[1] <= dims[1] && ghost_high[2] <= dims[2]);

				tile.data = new data_t[(ghost_high[0] - ghost_low[0])*(ghost_high[1] - ghost_low[1])*(ghost_high[2] - ghost_low[2])];
				assert(tile.data != NULL);
				// read tile function from global function
				{
					int64_t offset = 0;
					for (int64_t k = ghost_low[2]; k < ghost_high[2]; k++) {
						for (int64_t j = ghost_low[1]; j < ghost_high[1]; j++) {
							for (int64_t i = ghost_low[0]; i < ghost_high[0]; i++) {
								tile.data[offset++] = m_func2->SampleImage({i, j, k});
							}
						}
					}
				}

				tile.grid = new GridType({ghost_high[0] - ghost_low[0], ghost_high[1] - ghost_low[1], ghost_high[2] - ghost_low[2]}, {0, 0, 0});
				tile.func = new GridFuncType(tile.grid, tile.data);
				tile.tgrid = new MeshType(tile.grid);

				tiles.push_back(tile);
			}
		}
	}

	delete m_func2;
	delete m_grid;

	auto total_t0 = std::chrono::high_resolution_clock::now();
	int tile_i;
#pragma omp parallel for schedule(dynamic)
	for (tile_i = 0; tile_i < (int)tiles.size(); tile_i++) {
		auto &tile = tiles[tile_i];

		int tid = omp_get_thread_num();

		// memory allocations not part of the Robins algorithm implementation
		{
			tile.labeling = new GradType(tile.tgrid);
			tile.labeling->ClearAllGradient();
		}

		// Robins algorithm implementation
		{
			auto total_t0 = std::chrono::high_resolution_clock::now();

			MaxVLType* maxv = new MaxVLType(tile.tgrid, tile.func);

			// TODO(1/21/2022): reuse memory allocation
			auto allocations_t0 = std::chrono::high_resolution_clock::now();		
			maxv->allocate();
			auto allocations_t1 = std::chrono::high_resolution_clock::now();

			auto max_labeling_t0 = std::chrono::high_resolution_clock::now();
			maxv->ComputeOutputSerial();
			auto max_labeling_t1 = std::chrono::high_resolution_clock::now();

			auto robins_t0 = std::chrono::high_resolution_clock::now();
			NewRobinsType* trobins =  new NewRobinsType(tile.grid, tile.func, tile.tgrid, maxv, tile.labeling);
			trobins->ComputePairing_slidingSerial();
			delete trobins;
			auto robins_t1 = std::chrono::high_resolution_clock::now();

			delete maxv;

			auto total_t1 = std::chrono::high_resolution_clock::now();

			tile.temporary_allocations_time_seconds = std::chrono::duration<double>(allocations_t1 - allocations_t0).count();
			tile.max_labeling_time_seconds = std::chrono::duration<double>(max_labeling_t1 - max_labeling_t0).count();
			tile.robins_time_seconds = std::chrono::duration<double>(robins_t1 - robins_t0).count();
			tile.total_time_seconds = std::chrono::duration<double>(total_t1 - total_t0).count();
		}
	}
	auto total_t1 = std::chrono::high_resolution_clock::now();

	// copy tiles to C array
	out_morse_complex->tiles = (struct tile *)malloc(tiles.size()*sizeof tiles[0]);
	out_morse_complex->tile_count = 0;
	for (auto tile : tiles) {
		out_morse_complex->tiles[out_morse_complex->tile_count++] = tile;
	}
}


// naive version that reads each tile separately
void
compute_localized_morse_complex_from_file(char const *file, int64_t const *dims,
	struct localized_morse_complex *out_morse_complex)
{
	// round up the number of tiles
	// TODO: maybe do remainder and if statement for better readability
	out_morse_complex->dims[0] = (dims[0] + tile_x - 1)/tile_x;
	out_morse_complex->dims[1] = (dims[1] + tile_y - 1)/tile_y;
	out_morse_complex->dims[2] = (dims[2] + tile_z - 1)/tile_z;

	std::vector<struct tile> tiles;
	{
		FILE *fp = fopen(file, "rb");
		assert(fp != NULL);

		for (int64_t k = 0; k < dims[2]; k += tile_z) {
			for (int64_t j = 0; j < dims[1]; j += tile_y) {
				for (int64_t i = 0; i < dims[0]; i+= tile_x) {
					struct tile tile = {
						.position = {i, j, k},
						.size = {(i + tile_x > dims[0]) ? (dims[0] - i) : tile_x, (j + tile_y > dims[1]) ? (dims[1] - j) : tile_y, (k + tile_z > dims[2]) ? (dims[2] - k) : tile_z},
						.offset = {(i == 0) ? 0 : 1, (j == 0) ? 0 : 1, (k == 0) ? 0 : 1},
						.bridge_set_offset = {(i + tile_x >= dims[0]) ? 0 : 1, (j + tile_y >= dims[1]) ? 0 : 1, (k + tile_z >= dims[2]) ? 0 : 1},
					};

					int64_t ghost_low[3] = {i - tile.offset[0], j - tile.offset[1], k - tile.offset[2]}; // inclusive
					int64_t ghost_high[3] = {i + tile.size[0] + 2*tile.bridge_set_offset[0], j + tile.size[1] + 2*tile.bridge_set_offset[1], k + tile.size[2] + 2*tile.bridge_set_offset[2]}; // noninclusive, +1 for bridge set, +1 for ghost

					tile.data = new data_t[(ghost_high[0] - ghost_low[0])*(ghost_high[1] - ghost_low[1])*(ghost_high[2] - ghost_low[2])];
					assert(tile.data != NULL);
					// read tile function from global function
					{
						int64_t offset = 0;
						for (int64_t k = ghost_low[2]; k < ghost_high[2]; k++) {
							for (int64_t j = ghost_low[1]; j < ghost_high[1]; j++) {
								int64_t file_offset = sizeof *tile.data*(ghost_low[0] + j*dims[0] + k*dims[0]*dims[1]);
								fseek(fp, file_offset, SEEK_SET);
								fread(&tile.data[offset], sizeof *tile.data, ghost_high[0] - ghost_low[0], fp);
								offset += ghost_high[0] - ghost_low[0];
							}
						}
					}

					tile.grid = new GridType({ghost_high[0] - ghost_low[0], ghost_high[1] - ghost_low[1], ghost_high[2] - ghost_low[2]}, {0, 0, 0});
					tile.func = new GridFuncType(tile.grid, tile.data);
					tile.tgrid = new MeshType(tile.grid);

					tiles.push_back(tile);
				}
			}
		}

		fclose(fp);
	}

	auto total_t0 = std::chrono::high_resolution_clock::now();
	int tile_i;
#pragma omp parallel for schedule(dynamic)
	for (tile_i = 0; tile_i < (int)tiles.size(); tile_i++) {
		auto &tile = tiles[tile_i];

		int tid = omp_get_thread_num();

		// memory allocations not part of the Robins algorithm implementation
		{
			tile.labeling = new GradType(tile.tgrid);
			tile.labeling->ClearAllGradient();
		}

		// Robins algorithm implementation
		{
			auto total_t0 = std::chrono::high_resolution_clock::now();

			MaxVLType* maxv = new MaxVLType(tile.tgrid, tile.func);

			// TODO(1/21/2022): reuse memory allocation
			auto allocations_t0 = std::chrono::high_resolution_clock::now();		
			maxv->allocate();
			auto allocations_t1 = std::chrono::high_resolution_clock::now();

			auto max_labeling_t0 = std::chrono::high_resolution_clock::now();
			maxv->ComputeOutputSerial();
			auto max_labeling_t1 = std::chrono::high_resolution_clock::now();

			auto robins_t0 = std::chrono::high_resolution_clock::now();
			NewRobinsType* trobins =  new NewRobinsType(tile.grid, tile.func, tile.tgrid, maxv, tile.labeling);
			trobins->ComputePairing_slidingSerial();
			delete trobins;
			auto robins_t1 = std::chrono::high_resolution_clock::now();

			delete maxv;

			auto total_t1 = std::chrono::high_resolution_clock::now();

			tile.temporary_allocations_time_seconds = std::chrono::duration<double>(allocations_t1 - allocations_t0).count();
			tile.max_labeling_time_seconds = std::chrono::duration<double>(max_labeling_t1 - max_labeling_t0).count();
			tile.robins_time_seconds = std::chrono::duration<double>(robins_t1 - robins_t0).count();
			tile.total_time_seconds = std::chrono::duration<double>(total_t1 - total_t0).count();
		}
	}
	auto total_t1 = std::chrono::high_resolution_clock::now();

	// copy tiles to C array
	out_morse_complex->tiles = (struct tile *)malloc(tiles.size()*sizeof tiles[0]);
	out_morse_complex->tile_count = 0;
	for (auto tile : tiles) {
		out_morse_complex->tiles[out_morse_complex->tile_count++] = tile;
	}
}


// global discrete gradient is reconstructed from the tiles and saved in row-major order
// NOTE(2/5/2022): assumes row-major tile layout and tiles be in row-major order too
void
save_discrete_gradient_to_file(struct localized_morse_complex const *morse_complex, char const *file)
{
	FILE *fp = fopen(file, "wb");
	assert(fp != NULL);

	struct tile *last_tile = &morse_complex->tiles[morse_complex->tile_count - 1];
	int64_t dims[3] = {
		last_tile->position[0] + last_tile->size[0],
		last_tile->position[1] + last_tile->size[1],
		last_tile->position[2] + last_tile->size[2],
	};

	printf("dims %d %d %d\n", dims[0], dims[1], dims[2]);
	printf("tile grid dims %d %d %d\n", morse_complex->dims[0], morse_complex->dims[1], morse_complex->dims[2]);

	int64_t in_tile_k = 0;
	for (int64_t tile_j = 0; tile_j < morse_complex->dims[1]; tile_j++) {

		// write plane of a tile row
		struct tile *first_tile_in_row = &morse_complex->tiles[tile_j*morse_complex->dims[0]];
	
		printf("limit j %d\n", 2*first_tile_in_row->size[1] - 1 + first_tile_in_row->bridge_set_offset[1]);
		printf("tile %d, offset[1] %d, bridge_set_offset[1] %d\n", tile_j*morse_complex->dims[0], first_tile_in_row->offset[1], first_tile_in_row->bridge_set_offset[1]);
		printf("tile dims %d %d %d\n", first_tile_in_row->size[0], first_tile_in_row->size[1], first_tile_in_row->size[2]);

		for (int64_t in_tile_j = 0; in_tile_j < 2*first_tile_in_row->size[1] - 1 + first_tile_in_row->bridge_set_offset[1]; in_tile_j++) {

			for (int64_t tile_i = 0; tile_i < morse_complex->dims[0]; tile_i++) {
				struct tile *tile = &morse_complex->tiles[tile_i + tile_j*morse_complex->dims[0]];
				printf("tile size %d %d %d\n", tile->offset[0] + tile->size[0] + 2*tile->bridge_set_offset[0], tile->offset[1] + tile->size[1] + 2*tile->bridge_set_offset[1], tile->offset[2] + tile->size[2] + 2*tile->bridge_set_offset[2]);

				// TODO(3/5/2022): set dimension of ascending manifold
				//	topo_algs = new TopologicalGradientUsingAlgorithms<MeshType, TopoFuncType, GradType>(m_topofunc, m_tgrid, labeling);
				//	topo_algs->setAscendingManifoldDimensions();
				int64_t row_count = 2*(tile->offset[0] + tile->size[0] + 2*tile->bridge_set_offset[0]) - 1;
				int64_t offset = 2*tile->offset[0] + (2*tile->offset[1] + in_tile_j)*row_count;
				int64_t count = 2*tile->size[0] - 1 + tile->bridge_set_offset[0];

				// set boundary flag
				{
					// x-axis
					if (tile->offset[0] == 0) {
						tile->labeling->m_dgrad->LabelArray()[offset].flag = true;
					}
					if (tile->bridge_set_offset[0] == 0) {
						tile->labeling->m_dgrad->LabelArray()[offset + count - 1].flag = true;	
					}
					// y-axis or z-axis
					for (int64_t i = 0; i < count; i++) {
						if ((tile->offset[1] == 0 && in_tile_j == 0) || (tile->bridge_set_offset[1] == 0 && in_tile_j + 1 == 2*first_tile_in_row->size[1] - 1)) {
							tile->labeling->m_dgrad->LabelArray()[offset + i].flag = true;
						}
						if ((tile->offset[2] == 0 && in_tile_k == 0) || (tile->bridge_set_offset[2] == 0 && in_tile_k + 1 == 2*first_tile_in_row->size[2] - 1)) {
							tile->labeling->m_dgrad->LabelArray()[offset + i].flag = true;
						}
					}
				}

				fwrite(&tile->labeling->m_dgrad->LabelArray()[offset], sizeof (GradBitfield), count, fp);
			}
		}
	}

	fclose(fp);
}


void
localized_morse_complex_free(struct localized_morse_complex *morse_complex)
{
	for (int64_t i = 0; i < morse_complex->tile_count; i++) {
		struct tile *tile = &morse_complex->tiles[i];

		delete[] tile->data;
		delete tile->grid;
		delete tile->func;
		delete tile->labeling;
		delete tile->tgrid;
	}
	free(morse_complex->tiles);
}



void
cell_to_global_coords(struct localized_morse_complex const *morse_complex, struct cell cell, int64_t *out_coords)
{
	struct tile *tile = &morse_complex->tiles[cell.tile_index];

	Vec3l coords;
	tile->tgrid->cellid2Coords(cell.local_index, coords);
	for (int64_t i = 0; i < 3; i++) {
		if (tile->offset[i] == 1) {
			coords[i] -= 2;
		}
	}

	int64_t tile_coords[3] = {
		cell.tile_index%morse_complex->dims[0],
		(cell.tile_index/morse_complex->dims[0])%morse_complex->dims[1],
		cell.tile_index/(morse_complex->dims[0]*morse_complex->dims[1]),
	};

	// TODO(3/5/2022): use tile->position instead of tile_coords...
	out_coords[0] = coords[0] + tile_coords[0]*2*tile_x;
	out_coords[1] = coords[1] + tile_coords[1]*2*tile_y;
	out_coords[2] = coords[2] + tile_coords[2]*2*tile_z;
}









std::vector<struct cell>
query_critical_cells(struct localized_morse_complex const *morse_complex, int64_t dimension)
{
	std::vector<struct cell> critical_cells;

	for (int64_t tile_i = 0; tile_i < morse_complex->tile_count; tile_i++) {
		struct tile tile = morse_complex->tiles[tile_i];

		// TODO(3/5/2022): start i,j,k at 2*tile.offset[0]...
		for (int64_t k = 2*tile.position[2]; k <= 2*(tile.position[2] + tile.size[2] - 1) + tile.bridge_set_offset[2]; k++) {
			for (int64_t j = 2*tile.position[1]; j <= 2*(tile.position[1] + tile.size[1] - 1) + tile.bridge_set_offset[1]; j++) {
				for (int64_t i = 2*tile.position[0]; i <= 2*(tile.position[0] + tile.size[0] - 1) + tile.bridge_set_offset[0]; i++) {
					int64_t i_tile = i - 2*tile.position[0] + 2*tile.offset[0];
					int64_t j_tile = j - 2*tile.position[1] + 2*tile.offset[1];
					int64_t k_tile = k - 2*tile.position[2] + 2*tile.offset[2];
					auto index_tile = tile.tgrid->coords2Cellid({i_tile, j_tile, k_tile});

					if (tile.labeling->getCritical(index_tile) && tile.tgrid->dimension(index_tile) == dimension) {
						struct cell c = {.local_index = (local_t)index_tile, .tile_index = (local_t)tile_i};
						critical_cells.push_back(c);
					}
				}
			}
		}
	}

	return critical_cells;
}




// TODO(2/14/2022): return as C array
std::vector<struct cell>
query_descending_manifold(struct localized_morse_complex const *morse_complex, struct cell maximum)
{
	// NOTE(2/14/2022): descending 3-manifolds can only split, thus we do not use a visited set and check that
	//	indeed the queried critical cell is a maximum
	assert(morse_complex->tiles[maximum.tile_index].labeling->getCritical(maximum.local_index));
	assert(morse_complex->tiles[maximum.tile_index].tgrid->dimension(maximum.local_index) == 3);

	std::stack<struct cell> tile_queue;
	tile_queue.push(maximum);

	std::vector<struct cell> cells;
	while (!tile_queue.empty()) {
		struct cell tile_cell = tile_queue.top();
		tile_queue.pop();

		struct tile *tile = &morse_complex->tiles[tile_cell.tile_index];

		// TODO(2/14/2022): we could cache the tile_index here and thus make the queue take less memory because
		//	only cells in the same tile are processed
		std::stack<struct cell> queue;
		queue.push(tile_cell);

		while (!queue.empty()) {
			struct cell cell = queue.top();
			queue.pop();

			cells.push_back(cell);

			MeshType::FacetsIterator facets(tile->tgrid);
			for (facets.begin(cell.local_index); facets.valid(); facets.advance()) {
				local_t c = facets.value();

				if (tile->labeling->getCritical(c)) {
					continue;
				}

				local_t index = tile->labeling->getPair(c);

				// one of the facets points to the current cell so we skip it, also we skip traversal to cells of different dimension
				if (index == cell.local_index || tile->tgrid->dimension(index) != tile->tgrid->dimension(cell.local_index)) {
					continue;
				}

				// determine tile coordinates
				bool outside_tile = false;
				{
					Vec3l coords;
					tile->tgrid->cellid2Coords(index, coords);

					outside_tile = coords[0] < 2*tile->offset[0] || coords[0] >= (2*(tile->offset[0] + tile->size[0]) - 1 + tile->bridge_set_offset[0])
						|| coords[1] < 2*tile->offset[1] || coords[1] >= (2*(tile->offset[1] + tile->size[1]) - 1 + tile->bridge_set_offset[1])
						|| coords[2] < 2*tile->offset[2] || coords[2] >= (2*(tile->offset[2] + tile->size[2]) - 1 + tile->bridge_set_offset[2]);
				}


				if (outside_tile) {
					Vec3l coords;
					tile->tgrid->cellid2Coords(index, coords);
					int64_t directions[3] = {0, 0, 0};
					for (int64_t i = 0; i < 3; i++) {
						if (coords[i] < 2*tile->offset[i]) {
							directions[i] = -1;
						} else if (coords[i] >= (2*(tile->offset[i] + tile->size[i]) - 1 + tile->bridge_set_offset[i])) {
							directions[i] = 1;
						}
					};

					int64_t tile_coords[3] = {
						tile_cell.tile_index%morse_complex->dims[0],
						(tile_cell.tile_index/morse_complex->dims[0])%morse_complex->dims[1],
						tile_cell.tile_index/(morse_complex->dims[0]*morse_complex->dims[1]),
					};

					int64_t next_tile_index = (tile_coords[0] + directions[0])
						+ (tile_coords[1] + directions[1])*morse_complex->dims[0]
						+ (tile_coords[2] + directions[2])*morse_complex->dims[0]*morse_complex->dims[2];
					struct tile *next_tile = &morse_complex->tiles[next_tile_index];

					// compute the first or last cell coordinate alongside each dimension and add/subtract 3 cells
					int64_t local_coords[3] = {coords[0], coords[1], coords[2]};
					for (int64_t i = 0; i < 3; i++) {
						if (directions[i] == -1) {
							local_coords[i] = (2*(next_tile->offset[i] + next_tile->size[i] + 2*next_tile->bridge_set_offset[i]) - 1 - 1) - 3;
						} else if (directions[i] == 1) {
							local_coords[i] = 3;
						}
					}
					int64_t next_local_index = next_tile->tgrid->coords2Cellid({local_coords[0], local_coords[1], local_coords[2]});
					
					struct cell c = {.local_index = (local_t)next_local_index, .tile_index = (local_t)next_tile_index};
					tile_queue.push(c);
				} else {
					struct cell c = {.local_index = index, .tile_index = cell.tile_index};
					queue.push(c);
				}
			}
		}
	}

	return cells;
}









// localized Morse complex version that does not store ghost layer after tile's gradient is computed
struct tile1 {
	int64_t position[3]; // in terms of the global index space
	int64_t size[3];

	int64_t gradient_size[3]; // 2*size - 1 (+ 1 if has bridge set)
	GInt::GradBitfield *gradient;
};

struct localized_morse_complex1 {
	struct tile1 *tiles;
	int64_t tile_count;
	int64_t dims[3];
};


int64_t
min(int64_t a, int64_t b)
{
	if (a < b) {
		return a;
	}
	return b;
}

int64_t
max(int64_t a, int64_t b)
{
	if (a > b) {
		return a;
	}
	return b;
}

void
compute_localized_morse_complex_from_grid1(data_t *data, int64_t const *dims,
	struct localized_morse_complex1 *out_morse_complex)
{
	// round up the number of tiles
	// TODO: maybe do remainder and if statement for better readability
	out_morse_complex->dims[0] = (dims[0] + tile_x - 1)/tile_x;
	out_morse_complex->dims[1] = (dims[1] + tile_y - 1)/tile_y;
	out_morse_complex->dims[2] = (dims[2] + tile_z - 1)/tile_z;

	GridType *m_grid = new GridType(Vec3l(dims[0], dims[1], dims[2]), Vec3b(0, 0, 0));
	GridFuncType *m_func2 = new GridFuncType(m_grid, data);

	out_morse_complex->tile_count = out_morse_complex->dims[0]*out_morse_complex->dims[1]*out_morse_complex->dims[2];
	out_morse_complex->tiles = (struct tile1 *)malloc(out_morse_complex->tile_count*sizeof out_morse_complex->tiles[0]);

	// TODO(3/22/2022): initialize these tiles the the parallel loop below?
	int64_t tile_offset = 0;
	for (int64_t k = 0; k < dims[2]; k += tile_z) {
		for (int64_t j = 0; j < dims[1]; j += tile_y) {
			for (int64_t i = 0; i < dims[0]; i+= tile_x) {
				int64_t size[3] = {
					(i + tile_x > dims[0]) ? (dims[0] - i) : tile_x,
					(j + tile_y > dims[1]) ? (dims[1] - j) : tile_y,
					(k + tile_z > dims[2]) ? (dims[2] - k) : tile_z,
				};
				struct tile1 tile = {
					.position = {i, j, k},
					.size = {size[0], size[1], size[2]},
					.gradient_size = {
						// TODO: exploit the fact that 2*size[0] is power of two if tile_x is power of two
						(i + size[0] == dims[0]) ? 2*size[0] - 1 : 2*size[0],
						(j + size[1] == dims[1]) ? 2*size[1] - 1 : 2*size[1],
						(k + size[2] == dims[2]) ? 2*size[2] - 1 : 2*size[2],
					},
				};
				out_morse_complex->tiles[tile_offset++] = tile;
			}
		}
	}

#pragma omp parallel for schedule(dynamic)
	for (int64_t tile_offset = 0; tile_offset < out_morse_complex->tile_count; tile_offset++) {
		struct tile1 *tile = &out_morse_complex->tiles[tile_offset];

		// inclusive
		int64_t low[3] = {
			(tile->position[0] == 0) ? tile->position[0] : tile->position[0] - 1,
			(tile->position[1] == 0) ? tile->position[1] : tile->position[1] - 1,
			(tile->position[2] == 0) ? tile->position[2] : tile->position[2] - 1,
		};
		// exclusive, +1 for bridge set, +1 for ghost
		int64_t high[3] = {
			min(tile->position[0] + tile->size[0] + 2, dims[0]),
			min(tile->position[1] + tile->size[1] + 2, dims[1]),
			min(tile->position[2] + tile->size[2] + 2, dims[2]),
		};

		// TODO: temporary allocations should use per thread allocator
		data_t *data = (data_t *)malloc((high[0] - low[0])*(high[1] - low[1])*(high[2] - low[2])*sizeof *data);
		assert(data != NULL);
		// read tile function from global function
		{
			int64_t offset = 0;
			for (int64_t k = low[2]; k < high[2]; k++) {
				for (int64_t j = low[1]; j < high[1]; j++) {
					for (int64_t i = low[0]; i < high[0]; i++) {
						data[offset++] = m_func2->SampleImage({i, j, k});
					}
				}
			}
		}

		GradType *labeling = NULL;
		{
			GridType *grid = new GridType({high[0] - low[0], high[1] - low[1], high[2] - low[2]}, {0, 0, 0});
			GridFuncType *func = new GridFuncType;
			func->RegularGridTrilinearFunctionSerial(grid, data);
			MeshType *tgrid = new MeshType(grid);

			// memory allocations not part of the Robins algorithm implementation
			labeling = new GradType(tgrid);
			labeling->ClearAllGradient();

			// Robins algorithm implementation
			{
				MaxVLType* maxv = new MaxVLType(tgrid, func);

				// TODO(1/21/2022): reuse memory allocation
				maxv->allocate();
				maxv->ComputeOutputSerial();

				NewRobinsType* trobins =  new NewRobinsType(grid, func, tgrid, maxv, labeling);
				trobins->ComputePairing_slidingSerial();

				delete trobins;
				delete maxv;
			}

			delete tgrid;
			delete func;
			delete grid;
		}

		free(data);

		assert(sizeof *tile->gradient == sizeof *labeling->m_dgrad->LabelArray());
		tile->gradient = (GInt::GradBitfield *)malloc(tile->gradient_size[0]*tile->gradient_size[1]*tile->gradient_size[2]*sizeof *tile->gradient);
		assert(tile->gradient != NULL);

		// copy discrete gradient to tile (we remove the ghost layer discrete gradient)
		{
			// TODO(3/6/2022): optimize away the indexing operations
			int64_t offset = 0;
			for (int64_t k = 0; k < tile->gradient_size[2]; k++) {
				for (int64_t j = 0; j < tile->gradient_size[1]; j++) {
					for (int64_t i = 0; i < tile->gradient_size[0]; i++) {
						int64_t index = (i + ((tile->position[0] == 0) ? 0 : 2))
							+ (j + ((tile->position[1] == 0) ? 0 : 2))*(2*(high[0] - low[0]) - 1)
							+ (k + ((tile->position[2] == 0) ? 0 : 2))*(2*(high[0] - low[0]) - 1)*(2*(high[1] - low[1]) - 1);
						tile->gradient[offset++] = labeling->m_dgrad->LabelArray()[index];
					}
				}
			}
		}
	}
	delete m_func2;
	delete m_grid;
}

void
localized_morse_complex_free1(struct localized_morse_complex1 *morse_complex)
{
	for (int64_t i = 0; i < morse_complex->tile_count; i++) {
		struct tile1 *tile = &morse_complex->tiles[i];

		free(tile->gradient);
	}
	free(morse_complex->tiles);
}


bool
is_critical(GInt::GradBitfield gradient)
{
	return gradient.pair == 7;
}

int64_t
dimension(int64_t i, int64_t j, int64_t k)
{
	return i%2 + j%2 + k%2;
}


std::vector<struct cell>
query_critical_cells1(struct localized_morse_complex1 const *morse_complex, int64_t dim)
{
	std::vector<struct cell> critical_cells;

	for (int64_t tile_i = 0; tile_i < morse_complex->tile_count; tile_i++) {
		struct tile1 tile = morse_complex->tiles[tile_i];

		// TODO(3/6/2022): we could iterate only through the correct dimension cells
		for (int64_t k = 0; k < tile.gradient_size[2]; k++) {
			for (int64_t j = 0; j < tile.gradient_size[1]; j++) {
				for (int64_t i = 0; i < tile.gradient_size[0]; i++) {
					int64_t index = i + j*tile.gradient_size[0] + k*tile.gradient_size[0]*tile.gradient_size[1];

					if (is_critical(tile.gradient[index]) && dimension(i, j, k) == dim) {
						struct cell c = {.local_index = (local_t)index, .tile_index = (local_t)tile_i};
						critical_cells.push_back(c);
					}
				}
			}
		}
	}

	return critical_cells;
}



/*
std::vector<struct cell>
query_descending_manifold1(struct localized_morse_complex1 const *morse_complex, struct cell maximum)
{
	// NOTE(2/14/2022): descending 3-manifolds can only split, thus we do not use a visited set and check that
	//	indeed the queried critical cell is a maximum
	assert(morse_complex->tiles[maximum.tile_index].labeling->getCritical(maximum.local_index));
	assert(morse_complex->tiles[maximum.tile_index].tgrid->dimension(maximum.local_index) == 3);

	std::stack<struct cell> tile_queue;
	tile_queue.push(maximum);

	std::vector<struct cell> cells;
	while (!tile_queue.empty()) {
		struct cell tile_cell = tile_queue.top();
		tile_queue.pop();

		struct tile *tile = &morse_complex->tiles[tile_cell.tile_index];

		// TODO(2/14/2022): we could cache the tile_index here and thus make the queue take less memory because
		//	only cells in the same tile are processed
		std::stack<struct cell> queue;
		queue.push(tile_cell);

		while (!queue.empty()) {
			struct cell cell = queue.top();
			queue.pop();

			cells.push_back(cell);

			MeshType::FacetsIterator facets(tile->tgrid);
			for (facets.begin(cell.local_index); facets.valid(); facets.advance()) {
				local_t c = facets.value();

				struct cell1 facet = ;

				// facet can be outside if it is of bridge set cell
				if (facet is outside) {
					tile_queue.push(facet);
					continue;
				}

				if (is_critical(facet)) {
					continue;
				}

				struct cell1 pair = pair(tile, facet);

				// one of the facets points to the current cell so we skip it, also we skip traversal to cells of different dimension
				if (pair == cell || dimension(pair) != dimension(cell)) {
					continue;
				}

				// determine tile coordinates
				bool outside_tile = false;
				{
					Vec3l coords;
					tile->tgrid->cellid2Coords(index, coords);

					outside_tile = coords[0] < 0 || coords[0] >= tile->gradient_size[0]
						|| coords[1] < 0 || coords[1] >= tile->gradient_size[1]
						|| coords[2] < 0 || coords[2] >= tile->gradient_size[2];
				}


				if (outside_tile) {
					Vec3l coords;
					tile->tgrid->cellid2Coords(index, coords);
					int64_t directions[3] = {0, 0, 0};
					for (int64_t i = 0; i < 3; i++) {
						if (coords[i] < 2*tile->offset[i]) {
							directions[i] = -1;
						} else if (coords[i] >= (2*(tile->offset[i] + tile->size[i]) - 1 + tile->bridge_set_offset[i])) {
							directions[i] = 1;
						}
					};

					int64_t tile_coords[3] = {
						tile_cell.tile_index%morse_complex->dims[0],
						(tile_cell.tile_index/morse_complex->dims[0])%morse_complex->dims[1],
						tile_cell.tile_index/(morse_complex->dims[0]*morse_complex->dims[1]),
					};

					int64_t next_tile_index = (tile_coords[0] + directions[0])
						+ (tile_coords[1] + directions[1])*morse_complex->dims[0]
						+ (tile_coords[2] + directions[2])*morse_complex->dims[0]*morse_complex->dims[2];
					struct tile *next_tile = &morse_complex->tiles[next_tile_index];

					// compute the first or last cell coordinate alongside each dimension and add/subtract 3 cells
					int64_t local_coords[3] = {coords[0], coords[1], coords[2]};
					for (int64_t i = 0; i < 3; i++) {
						if (directions[i] == -1) {
							local_coords[i] = (2*(next_tile->offset[i] + next_tile->size[i] + 2*next_tile->bridge_set_offset[i]) - 1 - 1) - 3;
						} else if (directions[i] == 1) {
							local_coords[i] = 3;
						}
					}
					int64_t next_local_index = next_tile->tgrid->coords2Cellid({local_coords[0], local_coords[1], local_coords[2]});
					
					struct cell c = {.local_index = (local_t)next_local_index, .tile_index = (local_t)next_tile_index};
					tile_queue.push(c);
				} else {
					struct cell c = {.local_index = index, .tile_index = cell.tile_index};
					queue.push(c);
				}
			}
		}
	}

	return cells;
}
*/



int64_t
cell_to_global_index1(struct localized_morse_complex1 const *morse_complex, struct cell cell)
{
	struct tile1 *last_tile = &morse_complex->tiles[morse_complex->tile_count - 1];
	int64_t dims[3] = {
		last_tile->position[0] + last_tile->size[0],
		last_tile->position[1] + last_tile->size[1],
		last_tile->position[2] + last_tile->size[2],
	};

	struct tile1 *tile = &morse_complex->tiles[cell.tile_index];

	int64_t local_coords[3] = {
		cell.local_index%tile->gradient_size[0],
		(cell.local_index/tile->gradient_size[0])%tile->gradient_size[1],
		cell.local_index/(tile->gradient_size[0]*tile->gradient_size[1]),
	};
	int64_t coords[3] = {
		local_coords[0] + 2*tile->position[0],
		local_coords[1] + 2*tile->position[1],
		local_coords[2] + 2*tile->position[2],
	};

	return coords[0] + coords[1]*(2*dims[0] - 1) + coords[2]*(2*dims[0] - 1)*(2*dims[1] - 1);
}



// global discrete gradient is reconstructed from the tiles and saved in row-major order
// TODO(3/8/2022): modifies the flag bit of discrete gradient, maybe we should allocate temporary tile
//	and keep the input gradient intact
// NOTE(3/5/2022): assumes row-major tile layout and tiles be in row-major order too
void
save_discrete_gradient_to_file1(struct localized_morse_complex1 *morse_complex, char const *file)
{
	FILE *fp = fopen(file, "wb");
	assert(fp != NULL);

	struct tile1 *last_tile = &morse_complex->tiles[morse_complex->tile_count - 1];
	int64_t dims[3] = {
		last_tile->position[0] + last_tile->size[0],
		last_tile->position[1] + last_tile->size[1],
		last_tile->position[2] + last_tile->size[2],
	};

	//printf("dims %d %d %d\n", dims[0], dims[1], dims[2]);
	//printf("tile grid dims %d %d %d\n\n", morse_complex->dims[0], morse_complex->dims[1], morse_complex->dims[2]);

	for (int64_t tile_k = 0; tile_k < morse_complex->dims[2]; tile_k++) {

		struct tile1 *first_tile_in_column = &morse_complex->tiles[tile_k*morse_complex->dims[0]*morse_complex->dims[1]];
		for (int64_t in_tile_k = 0; in_tile_k < first_tile_in_column->gradient_size[2]; in_tile_k++) {

			for (int64_t tile_j = 0; tile_j < morse_complex->dims[1]; tile_j++) {

				// write planes of a tile row
				struct tile1 *first_tile_in_row = &morse_complex->tiles[tile_j*morse_complex->dims[0] + tile_k*morse_complex->dims[0]*morse_complex->dims[1]];
		
				//printf("limit j %d\n", first_tile_in_row->gradient_size[1]);
				//printf("tile %d\n", tile_j*morse_complex->dims[0]);
				//printf("tile dims %d %d %d\n", first_tile_in_row->size[0], first_tile_in_row->size[1], first_tile_in_row->size[2]);

				for (int64_t in_tile_j = 0; in_tile_j < first_tile_in_row->gradient_size[1]; in_tile_j++) {
					for (int64_t tile_i = 0; tile_i < morse_complex->dims[0]; tile_i++) {
						struct tile1 *tile = &morse_complex->tiles[tile_i + tile_j*morse_complex->dims[0] + tile_k*morse_complex->dims[0]*morse_complex->dims[1]];

						// TODO(3/5/2022): set dimension of ascending manifold (this is used to restrict traversals of descending manifolds to
						//	only the necessary traversal, so instead of 3-manifold only 1-manifold is traversed for maxima for example)
						//	topo_algs = new TopologicalGradientUsingAlgorithms<MeshType, TopoFuncType, GradType>(m_topofunc, m_tgrid, labeling);
						//	topo_algs->setAscendingManifoldDimensions();
						int64_t offset = in_tile_j*tile->gradient_size[0] + in_tile_k*tile->gradient_size[0]*tile->gradient_size[1];
						int64_t count = tile->gradient_size[0];

						// set boundary flag
						// TODO(3/8/2022): why is the boundary flag not already set after running Robins algorithm?
						//	(it is likely because the gradient file is precomputed and stored on disk and never directly used after Robins algorithm)
						{
							// x-axis
							if (tile->position[0] == 0) {
								tile->gradient[offset].flag = true;
							}
							if (tile->gradient_size[0] == 2*tile->size[0] - 1) { // last tile in a row
								tile->gradient[offset + count - 1].flag = true;	
							}
							// y-axis or z-axis
							for (int64_t i = 0; i < count; i++) {
								if ((tile->position[1] == 0 && in_tile_j == 0) || (tile->gradient_size[1] == 2*tile->size[1] - 1 && in_tile_j + 1 == first_tile_in_row->gradient_size[1])) {
									tile->gradient[offset + i].flag = true;
								}
								if ((tile->position[2] == 0 && in_tile_k == 0) || (tile->gradient_size[2] == 2*tile->size[2] - 1 && in_tile_k + 1 == first_tile_in_row->gradient_size[2])) {
									tile->gradient[offset + i].flag = true;
								}
							}
						}

						fwrite(&tile->gradient[offset], sizeof (GradBitfield), count, fp);
					}
				}
				//printf("\n");
			}
		}
	}

	fclose(fp);
}


#endif