#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif
#include "filteredlinesoup.h"

using namespace GInt;



void compute_sampled_soups(std::string basename, float persistence, float filter_value,
	float* fdata, GridFuncType* raw_function, GridFuncType* gauss_function, long long mX, long long mY);

	
	
template<class MSCType, class FuncType, typename dtype>
void mscToCellSoup(MSCType* msc, FuncType* func, dtype filter_value, vector<INDEX_TYPE>& global_soup) {
	printf("filtering arcs:\n");
	// now we need to compute average values
	printf(" -- identifying living arcs\n");
	MSCType::LivingArcsIterator ait(msc);
	vector<int> living_arcs;
	for (ait.begin(); ait.valid(); ait.advance()) living_arcs.push_back(ait.value());
	int num_arcs = living_arcs.size();
	printf(" -- parallel computing averages on %d arcs\n", num_arcs);


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
				sum += func->cellValue(id);
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
		printf("all threads joined\n");
}



int main(int argc, char** argv) {





	// READ IN THE COMMAND LINE ARGUMENTS
	int X, Y;
	if (argc < 6) { printf("Usage: X Y topo_filename persistence filter_value [filter_filename] \n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	auto filename = std::string(argv[3]);
	float persistence;
	sscanf(argv[4], "%f", &persistence);

	float filter_value;
	sscanf(argv[5], "%f", &filter_value);
	auto gauss_filename = std::string(argv[6]);
	auto raw_filename = std::string(argv[7]);
	auto outputstring = std::string(argv[8]);






	std::string mask_filename;

	int data_dims[3] = { X,Y,1 };


	GridType* underlying_grid;
	GridFuncType* grid_function;
	GridFuncType* raw_function;
	GridFuncType* gauss_function;


	// set up structures to navigate grid, and load the 3d image
	underlying_grid = new GridType(GInt::Vec2i(X, Y), GInt::Vec2b(0, 0));
	grid_function = new GridFuncType(underlying_grid);
	raw_function = new GridFuncType(underlying_grid);
	gauss_function = new GridFuncType(underlying_grid);

//#define DEBUG_ALL
#ifdef DEBUG_ALL
	// READ IN THE COMMAND LINE ARGUMENTS
	filename = "testoutput_" + std::to_string(X) + "x" + std::to_string(Y) + ".raw";
	float* newdat = new float[underlying_grid->NumElements()];
	persistence = 5.0;
	filter_value = 0.3;
	while (true) {
		
		// make random data
		for (int i = 0; i < underlying_grid->NumElements(); i++) {
			newdat[i] = float(rand() % 223) / 223.0;
		}
		FILE* fout = fopen(filename.c_str(), "wb");
		fwrite(newdat, sizeof(float), underlying_grid->NumElements(), fout);
		fclose(fout);
#endif

	grid_function->LoadImageFromFile(filename.c_str());
	raw_function->LoadImageFromFile(raw_filename.c_str());
	gauss_function->LoadImageFromFile(gauss_filename.c_str());

	compute_sampled_soups(outputstring, persistence, filter_value, grid_function->GetImage(), raw_function, gauss_function, X, Y);
	printf("ALL DONE\n");
#ifdef DEBUG_ALL
}
#endif

}


// User supplies 3 numpy 2d arrays -- laplacian, gaussian, raw
// computes msc on first, samples all and outputs soups and values files
void compute_sampled_soups(std::string basename, float persistence, float filter_value,
	float* fdata, GridFuncType* raw_function, GridFuncType* gauss_function, long long mX, long long mY) {

	//GridFuncType* grid_function = new GridFuncType(underlying_grid, fdata);
	
	Accurate2D::DiscreteGradientBuilder* dgb = new Accurate2D::DiscreteGradientBuilder();
	dgb->SetFloadArrayAndDims(mX, mY, fdata);
	dgb->SetNeededAccuracy(1, 0);
	dgb->SetParallelism(1);
	dgb->ComputeDiscreteGradient();

	auto grid = dgb->GetGrid();
	auto gridfunc = dgb->GetGridFunc();
	auto mesh = dgb->GetTopoMesh();
	auto meshfunc = dgb->GetMeshFunc();
	auto grad = dgb->GetGrad();
	typedef GInt::MorseSmaleComplexBasic<float, Accurate2D::MeshType, Accurate2D::MeshFuncType, Accurate2D::GradType> MyMscType;

	MyMscType* msc = new MyMscType(grad, mesh, meshfunc);
	msc->SetBuildArcGeometry(Vec3b(false, true, true)); // we only need geometric realizations of 2d saddle-max arcs
	msc->ComputeFromGrad();

	// get persistence to simplify to:
	float maxval = gridfunc->GetMaxValue();
	float minval = gridfunc->GetMinValue();

	float pers_limit = 0.01 * persistence * (maxval - minval);
	printf("funciton range: [%f, %f], target persistence = %f\n", minval, maxval, pers_limit);

	msc->ComputeHierarchy(pers_limit);
	// get persistence to get to the number of maxima requested
	msc->SetSelectPersAbs(pers_limit);

	vector<INDEX_TYPE> global_soup;
	mscToCellSoup<MyMscType, Accurate2D::MeshFuncType, float>(msc, meshfunc, filter_value, global_soup);

	if (global_soup.size() == 0) return;
	char soupname[2048];
	sprintf(soupname, "%s.soup", basename.c_str());
	FILE* fout = fopen(soupname, "wb");
	int count = global_soup.size();
	fwrite(&count, sizeof(int), 1, fout);
	fwrite(&(global_soup[0]), sizeof(INDEX_TYPE), global_soup.size(), fout);
	fclose(fout);
	printf("wrote %d soup elements\n", global_soup.size());
	
	MaxVLType* raw_maxv = new MaxVLType(mesh, raw_function);
	raw_maxv->ComputeOutput();
	MaxVLType* gauss_maxv = new MaxVLType(mesh, gauss_function);
	gauss_maxv->ComputeOutput();


	TopoFuncType* raw_mesh_func = new TopoFuncType();
	raw_mesh_func->setMeshAndFuncAndMaxVLabeling(mesh, raw_function, raw_maxv);
	TopoFuncType* gauss_mesh_func = new TopoFuncType();
	gauss_mesh_func->setMeshAndFuncAndMaxVLabeling(mesh, gauss_function, gauss_maxv);

	sprintf(soupname, "%s.soupvals", basename.c_str());
	fout = fopen(soupname, "wb");

	fwrite(&count, sizeof(int), 1, fout);

	for (int i = 0; i < count; i++) {
		float vals[3];
		auto id = global_soup[i];
		vals[0] = meshfunc->cellValue(id);
		vals[1] = raw_mesh_func->cellValue(id);
		vals[2] = gauss_mesh_func->cellValue(id);
		fwrite(vals, sizeof(float), 3, fout);
	}
	fclose(fout);
	printf("wrote %d soup x 3 values\n", global_soup.size());


}


