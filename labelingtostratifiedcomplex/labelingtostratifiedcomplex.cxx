
#include <vector>
#include <set>
#include <queue>
#include <time.h>
#include "gi_numeric_integrator_region_stop.h"
#include "gi_numeric_integrator_expanding_region_stop.h"
#include "gi_timing.h"
#include "gi_labeling_to_bounary_labeling.h"
#include "gi_topological_explicit_mesh_function.h"
#include "gi_robin_labeling.h"
#include "gi_adaptive_in_quad_euler_advector.h"
#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_isolated_region_remover.h"
#include "gi_isolated_region_remover2.h"
#include "gi_topological_utility_functions.h"
#include "gi_numeric_streamline_integrator.h"
#include "gi_numeric_integrator_expanding_region_stop_filtered2.h"
#include "gi_numeric_integrator_expanding_region_stop_filtered.h"
#include <vector>
#include <set>
#include <queue>
#include <time.h>
#include "gi_strictly_numeric_integrator.h"
#include "gi_numeric_integrator_region_stop.h"
#include "gi_numeric_integrator_expanding_region_stop.h"
#include "gi_timing.h"
#include "gi_labeling_to_bounary_labeling.h"
#include "gi_topological_explicit_mesh_function.h"
#include "gi_robin_labeling.h"
#include "gi_adaptive_in_quad_euler_advector.h"
#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_isolated_region_remover.h"
#include "gi_topological_utility_functions.h"
#include "gi_numeric_streamline_integrator.h"
#include "gi_bifiltration_pairing.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_extrema_region_builder.h"
#include "gi_numeric_integrator_path_compressing.h"
#include "gi_fast_robins_noalloc.h"
#include "gi_graphs.h"
#include "gi_morse_smale_complex_basic.h"


using namespace GInt;

typedef RegularGrid3D GridType;
typedef TopologicalRegularGrid3D MeshType;
typedef MeshCellsGraphBuilder<MeshType> SkeletonType;

// annoyingly, we need to define these types even though we never use MaxVLType...
// to be fixed in the future
typedef RegularGridTrilinearFunction GridFuncType;
typedef RegularGridMaxMinVertexLabeling3D<MeshType, GridFuncType> MaxVLType;
// end annoyingly

int X, Y, Z;
int per_x, per_y, per_z;
std::string filename;
int parallelism = -1;
int maxdimofcomplex = 1;

bool GetOptions(int argc, char** argv) {
	if (argc < 6) { printf("Usage: X Y Z filename maxdimofcomplex [parallelism=ompmaxnumthreads]\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
	sscanf(argv[5], "%f", &maxdimofcomplex);

	if (argc >= 7)
		sscanf(argv[6], "%d", &parallelism);

	// set remaining values
	if (parallelism != -1) {
		omp_set_num_threads(parallelism);
	}

}


RegularGrid3D* g_grid;
MeshType *g_topo_grid;

int main(int argc, char** argv) {

	// read command line options
	GetOptions(argc, argv);

	// START IO ---------------------------

	g_grid = new RegularGrid3D(Vec3i(X, Y, Z), Vec3b(per_x, per_y, per_z));
	DenseLabeling<INT_TYPE>* region_labels = new DenseLabeling<INT_TYPE>(g_grid->NumElements());
	region_labels->ReadFromFile(filename.c_str());

	// END IO -----------------------------

	g_topo_grid = new MeshType(g_grid);

	VertexLabelingToBoundaryLabeling<int, MaxVLType>* g_edge_map = new
		VertexLabelingToBoundaryLabeling<int, MaxVLType>(g_topo_grid);
	g_edge_map->ComputeRegionBoundaryKind(region_labels);

	char label_name[2048];
	sprintf(label_name, "%s._edges_%dx%dx%d.raw", filename.c_str(), g_topo_grid->XYZ()[0], g_topo_grid->XYZ()[1], g_topo_grid->XYZ()[2]);
	printf("outputting labels:\n%s\n", label_name);
	g_edge_map->GetOutputLabels()->OutputToFile(label_name);

	// here is an msc that will be hand constructed
	GInt::NoFuncNoGradMSC<float, MeshType>* test_msc =
		new GInt::NoFuncNoGradMSC<float, MeshType>(g_topo_grid);

	// gather the cells forming the stratified complex
	std::vector< INDEX_TYPE > maxima_cells;
	std::vector<INDEX_TYPE> topo_index_partition;
	int num_threads;
#pragma omp parallel
	{
#pragma omp single
		{
			num_threads = omp_get_num_threads();
			ArrayIndexPartitioner::EvenChunkSplit(g_topo_grid->numCells(), num_threads, topo_index_partition);
		}

		int thread_num = omp_get_thread_num();
		// in parallel go through and find all 2-saddles
		std::vector<INDEX_TYPE> local_maxima;
		MeshType::DCellsIterator hexes( g_topo_grid, 3,topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
		for (hexes.begin(); hexes.valid(); hexes.advance()) {
			INDEX_TYPE cell_id = hexes.value();
			auto tval = g_edge_map->GetOutputLabels()->GetLabel(cell_id);
			if (tval == 3) {
				local_maxima.push_back(cell_id);
			}
		}
#pragma omp critical
		{
			maxima_cells.insert(maxima_cells.end(), local_maxima.begin(), local_maxima.end());
		}
	}

	// first value is id of the lowest cell
	typedef pair<int, vector<INDEX_TYPE>> SadLigPair;
	INT_TYPE num_maxima = maxima_cells.size();
#pragma omp parallel for schedule(dynamic)
	for (int max_id = 0; max_id < num_maxima; max_id++) {

		MeshType::FacetsIterator facets(g_topo_grid);

		for (facets.begin(max_id); facets.valid(); facets.advance()) {
			INDEX_TYPE quad_id = facets.value();
			auto quad_kind = g_edge_map->GetOutputLabels()->GetLabel(quad_id);
			if (quad_kind != 2) continue; // this is not part of a  linear structure

			// start a line HERE - now the next max is going to be the other end of the line
			vector<INDEX_TYPE> local_line;
			local_line.push_back(max_id); // start with maximum on the line
			INDEX_TYPE lowest_quad_lineid = 1; // initialize this to the first quad
			local_line.push_back(quad_id);

			// []| --> now we want to get next []
			while (true) {
				bool has_next_hex = false;
				INDEX_TYPE next_hex;
				BYTE_TYPE next_hex_kind;
				MeshType::CofacetsIterator hexes(g_topo_grid);
				hexes.begin(quad_id); 
				INDEX_TYPE possible_next_hex = hexes.value();
				// don't allow going back on the line, so check
				if (possible_next_hex == local_line[local_line.size() - 2]) hexes.advance();
					
				// check if we hit a dead end
				if (!hexes.valid()) break;
				next_hex = hexes.value();

				// this means that the next one MUST be the only possible other hex
				next_hex_kind = g_edge_map->GetOutputLabels()->GetLabel(possible_next_hex);
				if (next_hex_kind < 2) break; // this is not part of a  linear structure, so dead end

				local_line.push_back(next_hex); // so push on the end of the line []|[]|[]

				// do we keep looking or are we done?
				if (next_hex_kind == 3) {
					// reached another critical point!!! []|[]|[]
					// local_line is the geometric segment, with lowest_quad_lineid the lowest quad
					//SadLigPair
				}
				
				

			}

			// 


		}

		MeshType::CofacetsIterator cofacets(g_topo_grid);

	}



	// gather the cells forming the stratified complex
	std::vector<std::pair<int, INDEX_TYPE> > skeleton_cells;
	std::vector<INDEX_TYPE> topo_index_partition;
	int num_threads;
#pragma omp parallel
	{
#pragma omp single
		{
			num_threads = omp_get_num_threads();
			ArrayIndexPartitioner::EvenChunkSplit(g_topo_grid->numCells(), num_threads, topo_index_partition);
		}

		int thread_num = omp_get_thread_num();
		// in parallel go through and find all 2-saddles
		std::vector<std::pair<int, INDEX_TYPE> > local_skeleton_cells;
		MeshType::AllCellsIterator all_cells(g_topo_grid, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
		for (all_cells.begin(); all_cells.valid(); all_cells.advance()) {
			INDEX_TYPE cell_id = all_cells.value();
			auto tval = g_edge_map->GetOutputLabels()->GetLabel(cell_id);
			if (tval > 2) {
				std::pair<int, INDEX_TYPE> p(tval, cell_id);
				local_skeleton_cells.push_back(p);
			}
		}
#pragma omp critical
		{
			skeleton_cells.insert(skeleton_cells.end(), local_skeleton_cells.begin(), local_skeleton_cells.end());
		}
	}


















	// gather the cells forming the stratified complex
	std::vector<std::pair<int, INDEX_TYPE> > skeleton_cells; 
	std::vector<INDEX_TYPE> topo_index_partition;
	int num_threads;
#pragma omp parallel
	{
#pragma omp single
		{
			num_threads = omp_get_num_threads();
			ArrayIndexPartitioner::EvenChunkSplit(g_topo_grid->numCells(), num_threads, topo_index_partition);
		}

		int thread_num = omp_get_thread_num();
		// in parallel go through and find all 2-saddles
		std::vector<std::pair<int, INDEX_TYPE> > local_skeleton_cells;
		MeshType::AllCellsIterator all_cells(g_topo_grid, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
		for (all_cells.begin(); all_cells.valid(); all_cells.advance()) {
			INDEX_TYPE cell_id = all_cells.value();
			auto tval = g_edge_map->GetOutputLabels()->GetLabel(cell_id);
			if (tval > 2) {
				std::pair<int, INDEX_TYPE> p(tval, cell_id);
				local_skeleton_cells.push_back(p);
			}
		}
#pragma omp critical
		{
			skeleton_cells.insert(skeleton_cells.end(), local_skeleton_cells.begin(), local_skeleton_cells.end());
		}
	}


	// now build the complex
	MeshCellsGraphBuilder<MeshType>* skeleton_graph =
		new MeshCellsGraphBuilder<MeshType>(g_topo_grid);
	for (auto& p : skeleton_cells) {
		skeleton_graph->AddCellIndex(p.second);
		if (p.first == 4) skeleton_graph->AddToForcedVertices(p.second);
	}
	auto* result_graph = skeleton_graph->ComputeStuff();

	

};


