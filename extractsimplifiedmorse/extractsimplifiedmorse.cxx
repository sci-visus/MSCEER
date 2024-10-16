#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif
#include "extractsimplifiedmorse.h"

using namespace GInt;



int main(int argc, char** argv) {

	// READ IN THE COMMAND LINE ARGUMENTS
	int X, Y, Z;
	std::string filename;
	if (argc < 6) { printf("Usage: X Y Z filename persistence mode[0=basins,1=mountains,2=both] \n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
	float persistence;
	sscanf(argv[5], "%f", &persistence);
  

	int extract_mode = 0;
	sscanf(argv[6], "%d", &extract_mode);
  std::string mask_filename;
  
  int data_dims[3]={X,Y,Z};
  printf("extractmsc persistence is %f\n", persistence);

	GridType* underlying_grid;
	GridFuncType* grid_function;
	MeshType *topological_grid;
	TopoFuncType* topological_grid_function;
	GradType *discrete_gradient;

	// ------ Load a discrete gradient ------------ 

	// set up structures to navigate grid, and load the 3d image
	underlying_grid = new GridType(GInt::Vec3i(X, Y, Z), GInt::Vec3b(0, 0, 0));
	grid_function = new GridFuncType(underlying_grid);
	grid_function->LoadImageFromFile(filename.c_str());

	printf("loaded cont function\n");

	// now set up an indexing scheme to use in a topological interpretation of the 
	// regular grid
	topological_grid = new MeshType(underlying_grid);

	// we will use a lazy-evaluation max vertex mesh function - i.e. the value of a 
	// cell in the topological grid is the maximum value of its vertices in the input
	// image
	MaxVLType* maxv = new MaxVLType(topological_grid, grid_function);
	maxv->ComputeOutput();
	topological_grid_function = new TopoFuncType();
	topological_grid_function->setMeshAndFuncAndMaxVLabeling(topological_grid, grid_function, maxv);

	// read the discrete gradient from disk
	discrete_gradient = new GradType(topological_grid);
	printf("created dgrad struct\n");
	char gradname[2048];
	sprintf(gradname, "%s.grad", argv[4]);
	discrete_gradient->load_from_file(gradname);


	// ----------- method 1: label based on simplified extremum graph ---------------------
	// make volumes of indices based on the discrete gradient
	
	// find all critical cells of dim 0 or 3
	std::vector<INDEX_TYPE> criticals[4];
	std::vector<INDEX_TYPE> topo_index_partition;
	int num_threads;
#pragma omp parallel
	{
#pragma omp single
		{
			num_threads = omp_get_num_threads();
			ArrayIndexPartitioner::EvenChunkSplit(topological_grid->numCells(), num_threads, topo_index_partition);
		}
		int thread_num = omp_get_thread_num();
		typename MeshType::AllCellsIterator all_cells_iterator(topological_grid, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
		std::vector<INDEX_TYPE> thread_criticals[4];
		for (all_cells_iterator.begin(); all_cells_iterator.valid(); all_cells_iterator.advance()) {
			INDEX_TYPE cell_id = all_cells_iterator.value();
			if (discrete_gradient->getCritical(cell_id)) {
				DIM_TYPE tdim = topological_grid->dimension(cell_id);
				thread_criticals[tdim].push_back(cell_id);
			}
		}
#pragma omp critical
		{
			criticals[0].insert(criticals[0].end(), thread_criticals[0].begin(), thread_criticals[0].end());
		}
	}

	// figure out index remapping
	SimplifiedExtremumGraph<MeshType, TopoFuncType, GradType>* simplified_ext_graph =
		new SimplifiedExtremumGraph<MeshType, TopoFuncType, GradType>(topological_grid, topological_grid_function, discrete_gradient);

	switch (extract_mode) {
	case 0:
		simplified_ext_graph->SetMode(SimplifiedExtremumGraph<MeshType, TopoFuncType, GradType>::EXTGRAPHMODE::MINS);
		break;
	case 1:
		simplified_ext_graph->SetMode(SimplifiedExtremumGraph<MeshType, TopoFuncType, GradType>::EXTGRAPHMODE::MAXS);
		break;
	default:
		simplified_ext_graph->SetMode(SimplifiedExtremumGraph<MeshType, TopoFuncType, GradType>::EXTGRAPHMODE::BOTH);
	}

	simplified_ext_graph->ComputeMinMapFromGradient(persistence);

	simplified_ext_graph->mMaxGraph.

	// create empty msc just to get access to tracing manifolds
	MscType* msc = new MscType(discrete_gradient, topological_grid, topological_grid_function);
	if (extract_mode == 0 || extract_mode == 2) {
		// extract basins

		// make space to store ids
		DenseLabeling<int>* basins = new DenseLabeling<int>(topological_grid->numCells(0));
		
		int num_mins = criticals[0].size();
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < num_mins; i++) {

			INDEX_TYPE min_id = criticals[0][i];
			std::set<INDEX_TYPE> manifold;
			msc->rec_man_trace_up(min_id, manifold);

			for (auto cell_id : manifold) {
				Vec3l coords;
				topological_grid->cellid2Coords(cell_id, coords);
				if (topological_grid->dimension(coords) != 0) continue; // skip non-vertex cells in asc 3-manifold
				INDEX_TYPE data_id = topological_grid->VertexNumberFromCellID(min_id);
				basins[data_id] = i;
			}
		}
		basins->OutputToIntFile("output_basins.raw");

	}



	
	
	
	
	// now compute the Morse-Smale complex
		MscType* msc = new MscType(discrete_gradient, topological_grid, topological_grid_function);
		msc->SetBuildArcGeometry(Vec3b(false, false, false)); // we only need geometric realizations of 2saddle-max arcs
		msc->ComputeFromGrad();

		// simplify to persistence and dump record to file
		char persname[2048];
		sprintf(persname, "%s.pers", argv[4]);
		msc->set_output_cancellation_records(persname);
		
		msc->ComputeHierarchy(persistence);
		// get persistence to get to the number of maxima requested
		msc->SetSelectPersAbs(persistence);

		//msc->WriteComplex1Skeleton("testcomplex.bin");
		//msc->print_complex_info(true);
		//
		////sanity check arcs
		//unordered_set<INT_TYPE> living_arcs;
		//MscType::LivingArcsIterator liv_arcs_it(msc);
		//for (liv_arcs_it.begin(); liv_arcs_it.valid(); liv_arcs_it.advance()) {
		//	living_arcs.insert(liv_arcs_it.value());
		//}

		//char pairsname[1024];
		//sprintf(pairsname, "%s.pairs.txt", argv[4]);
		//FILE* f_pairs = fopen(pairsname, "w");
		//for (auto arcid : living_arcs) {
		//	auto& a = msc->getArc(arcid);
		//	if (a.dim != 2) continue;
		//	if (a.boundary) continue;
		//	float low_val = msc->getNode(a.lower).value;
		//	float high_val = msc->getNode(a.upper).value;
		//	fprintf(f_pairs, "%f %f\n", low_val, high_val);
		//}
		//fclose(f_pairs);

	



	return 0;

}


