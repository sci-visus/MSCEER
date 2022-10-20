#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif
#include "extractfilteredlinesoup.h"

using namespace GInt;


std::string filename;
std::string filter_filename;
float filter_value;

int main(int argc, char** argv) {

	// READ IN THE COMMAND LINE ARGUMENTS
	int X, Y;
	if (argc < 6) { printf("Usage: X Y topo_filename persistence filter_value [filter_filename] \n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	filename = std::string(argv[3]);
	float persistence;
	sscanf(argv[4], "%f", &persistence);


	sscanf(argv[5], "%f", &filter_value);


	std::string mask_filename;

	int data_dims[3] = { X,Y,1 };


	GridType* underlying_grid;
	GridFuncType* grid_function;
	MeshType *topological_grid;
	TopoFuncType* topological_grid_function;
	GradType *discrete_gradient;

	GridFuncType* filter_function;
	TopoFuncType* topological_grid_filter_function;

	// set up structures to navigate grid, and load the 3d image
	underlying_grid = new GridType(GInt::Vec2i(X, Y), GInt::Vec2b(0, 0));
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

	if (argc == 7) {
		filter_filename = std::string(argv[6]);
		filter_function = new GridFuncType(underlying_grid);
		filter_function->LoadImageFromFile(filter_filename.c_str());
		topological_grid_filter_function = new TopoFuncType();
		topological_grid_filter_function->setMeshAndFuncAndMaxVLabeling(topological_grid, filter_function, new MaxVLType(topological_grid, filter_function));

	}
	else {
		filter_function = grid_function;
		topological_grid_filter_function = topological_grid_function;
	}

	// read the discrete gradient from disk
	discrete_gradient = new GradType(topological_grid);
	printf("created dgrad struct\n");
	char gradname[2048];
	sprintf(gradname, "%s.randomized", filename.c_str());
	discrete_gradient->load_from_file(gradname);

	// now compute the Morse-Smale complex
	MscType* msc = new MscType(discrete_gradient, topological_grid, topological_grid_function);
	msc->SetBuildArcGeometry(Vec3b(false, true, true)); // we only need geometric realizations of 2d saddle-max arcs
	msc->ComputeFromGrad();

	// simplify to persistence and dump record to file
	char persname[2048];
	sprintf(persname, "%s.pers", argv[3]);
	msc->set_output_cancellation_records(persname);

	// get persistence to simplify to:
	float maxval = grid_function->GetMaxValue();
	float minval = grid_function->GetMinValue();

	float pers_limit = 0.01 * persistence * (maxval - minval);
	printf("funciton range: [%f, %f], target persistence = %f\n", minval, maxval, pers_limit);

	msc->ComputeHierarchy(pers_limit);
	// get persistence to get to the number of maxima requested
	msc->SetSelectPersAbs(pers_limit);

	printf("filtering arcs:\n");
	// now we need to compute average values
	printf(" -- identifying living arcs\n");
	MscType::LivingArcsIterator ait(msc);
	vector<int> living_arcs;
	for (ait.begin(); ait.valid(); ait.advance()) living_arcs.push_back(ait.value());
	int num_arcs = living_arcs.size();
	printf(" -- parallel computing averages on %d arcs\n", num_arcs);

	vector<INDEX_TYPE> global_soup;
#pragma omp parallel
	{
		std::vector<INDEX_TYPE> index_partition;
		int num_threads = omp_get_num_threads();
		ArrayIndexPartitioner::EvenChunkSplit(num_arcs, num_threads, index_partition);
		int thread_num = omp_get_thread_num();

		// iterate over all vertices
		unordered_set<INDEX_TYPE> local_soup;
		for (int arc_pos = index_partition[thread_num]; arc_pos < index_partition[thread_num + 1]; arc_pos++) {
			
			// fill geometry
			int arc_id = living_arcs[arc_pos];
			vector<INDEX_TYPE> arc_geom;
			msc->fillArcGeometry(arc_id, arc_geom);

			// compute average;
			float sum = 0;
			for (auto id : arc_geom) {
				sum += topological_grid_filter_function->cellValue(id);
			}
			float ave = sum / arc_geom.size();

			// filter
			if (ave > filter_value) {
				local_soup.insert(arc_geom.begin(), arc_geom.end());
			}
		}


#pragma omp critical
		{
			printf("   -- thread %d adding %d cells\n", thread_num, local_soup.size());
			std::copy(local_soup.begin(), local_soup.end(), std::back_inserter(global_soup));
		}

		
	}

	// end parallel section, write out results

	char soupname[2048];
	sprintf(soupname, "%s.soup", argv[3]);
	FILE* fout = fopen(soupname, "wb");
	int count = global_soup.size();
	fwrite(&count, sizeof(int), 1, fout);
	fwrite(&(global_soup[0]), sizeof(INDEX_TYPE), global_soup.size(), fout);
	fclose(fout);
	printf("wrote %d elements\n", global_soup.size());

	return 0;

}


