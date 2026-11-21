#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif
#include "filteredlinesoup_py.h"

using namespace GInt;





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
}



//int main(int argc, char** argv) {
//
//	// READ IN THE COMMAND LINE ARGUMENTS
//	int X, Y;
//	if (argc < 6) { printf("Usage: X Y topo_filename persistence filter_value [filter_filename] \n"); return 0; }
//	sscanf(argv[1], "%d", &X);
//	sscanf(argv[2], "%d", &Y);
//	filename = std::string(argv[3]);
//	float persistence;
//	sscanf(argv[4], "%f", &persistence);
//
//
//	sscanf(argv[5], "%f", &filter_value);
//
//
//	std::string mask_filename;
//
//	int data_dims[3] = { X,Y,1 };
//
//
//	GridType* underlying_grid;
//	GridFuncType* grid_function;
//	MeshType *topological_grid;
//	TopoFuncType* topological_grid_function;
//	GradType *discrete_gradient;
//
//	GridFuncType* filter_function;
//	TopoFuncType* topological_grid_filter_function;
//
//	// set up structures to navigate grid, and load the 3d image
//	underlying_grid = new GridType(GInt::Vec2i(X, Y), GInt::Vec2b(0, 0));
//	grid_function = new GridFuncType(underlying_grid);
//	grid_function->LoadImageFromFile(filename.c_str());
//
//
//	//py::scoped_interpreter guard{}; // start the interpreter and keep it alive
//	//py::print("Hello, World!"); // use the Python API
//
//	//py::module_ imstuff = py::module_::import("process_image");
//
//	//// make some data Y, X order
//	//py::array_t<float> arr({ 300, 400 });
//	//float* ptr = arr.mutable_data();
//	//// make some image
//	//for (int i = 0; i < 300; i++)
//	//	for (int j = 0; j < 400; j++)
//	//	{
//	//		float val = 0;
//	//		if (i == j) val = 1;
//	//		if ((i + 50) % 100 == 0) val = 1.1f;
//	//		ptr[i * 400 + j] = val; // column major order
//	//	}
//
//	//py::object result = imstuff.attr("ExistingProcess")(arr, true, 16, 3, 30, 2, 1.25, "test");
//
//	//auto result_t = result.cast<py::tuple>();
//	//auto r1 = result_t[0].cast<py::array_t<float>>();
//	//auto r2 = result_t[1].cast<py::array_t<float>>();
//	//auto r3 = result_t[2].cast<py::array_t<float>>();
//
//
//
//	printf("loaded cont function\n");
//
//	// now set up an indexing scheme to use in a topological interpretation of the 
//	// regular grid
//	topological_grid = new MeshType(underlying_grid);
//
//	// we will use a lazy-evaluation max vertex mesh function - i.e. the value of a 
//	// cell in the topological grid is the maximum value of its vertices in the input
//	// image
//	MaxVLType* maxv = new MaxVLType(topological_grid, grid_function);
//	maxv->ComputeOutput();
//	topological_grid_function = new TopoFuncType();
//	topological_grid_function->setMeshAndFuncAndMaxVLabeling(topological_grid, grid_function, maxv);
//
//	if (argc == 7) {
//		filter_filename = std::string(argv[6]);
//		filter_function = new GridFuncType(underlying_grid);
//		filter_function->LoadImageFromFile(filter_filename.c_str());
//		topological_grid_filter_function = new TopoFuncType();
//		topological_grid_filter_function->setMeshAndFuncAndMaxVLabeling(topological_grid, filter_function, new MaxVLType(topological_grid, filter_function));
//
//	}
//	else {
//		filter_function = grid_function;
//		topological_grid_filter_function = topological_grid_function;
//	}
//
//	// read the discrete gradient from disk
//	discrete_gradient = new GradType(topological_grid);
//	printf("created dgrad struct\n");
//	char gradname[2048];
//	sprintf(gradname, "%s.randomized", filename.c_str());
//	discrete_gradient->load_from_file(gradname);
//
//	// now compute the Morse-Smale complex
//	MscType* msc = new MscType(discrete_gradient, topological_grid, topological_grid_function);
//	msc->SetBuildArcGeometry(Vec3b(false, true, true)); // we only need geometric realizations of 2d saddle-max arcs
//	msc->ComputeFromGrad();
//
//	// simplify to persistence and dump record to file
//	char persname[2048];
//	sprintf(persname, "%s.pers", argv[3]);
//	msc->set_output_cancellation_records(persname);
//
//	// get persistence to simplify to:
//	float maxval = grid_function->GetMaxValue();
//	float minval = grid_function->GetMinValue();
//
//	float pers_limit = 0.01 * persistence * (maxval - minval);
//	printf("funciton range: [%f, %f], target persistence = %f\n", minval, maxval, pers_limit);
//
//	msc->ComputeHierarchy(pers_limit);
//	// get persistence to get to the number of maxima requested
//	msc->SetSelectPersAbs(pers_limit);
//
//	vector<INDEX_TYPE> global_soup;
//	mscToCellSoup<MscType, TopoFuncType, float>(msc, topological_grid_filter_function, filter_value, global_soup);
//
//	// end parallel section, write out results
//
//	char soupname[2048];
//	sprintf(soupname, "%s.soup", argv[3]);
//	FILE* fout = fopen(soupname, "wb");
//	int count = global_soup.size();
//	fwrite(&count, sizeof(int), 1, fout);
//	fwrite(&(global_soup[0]), sizeof(INDEX_TYPE), global_soup.size(), fout);
//	fclose(fout);
//	printf("wrote %d elements\n", global_soup.size());
//
//	return 0;
//
//}


// User supplies 3 numpy 2d arrays -- laplacian, gaussian, raw
// computes msc on first, samples all and outputs soups and values files
void compute_sampled_soups(std::string basename, float persistence, float filter_value,
	py::array_t<float> laplace,
	py::array_t<float> gaussian,
	py::array_t<float> raw) {

	Py_BEGIN_ALLOW_THREADS

	auto la = laplace.unchecked<2>(); // make fast access
	auto ga = gaussian.unchecked<2>();
	auto ra = raw.unchecked<2>();
	// do dumb copy
	// in python y moves fastest
	auto pX = la.shape(0);
	auto pY = la.shape(1);
	auto mX = pY;
	auto mY = pX;	
	
	//GridType* underlying_grid = new GridType({ la.shape(0), la.shape(1) }, { false, false });
	float* fdata = new float[pX * pY];
	float* frawdata = new float[pX * pY];
	float* fgaussdata = new float[pX * pY];


	printf("FROM PYTHON ARRAY SHAPE (%d, %d)\n", pX , pY);
	for (py::ssize_t i = 0; i < pX; i++)
		for (py::ssize_t j = 0; j <pY; j++) {
			fdata[j + i * mX] = la(i, j);
			frawdata[j + i * mX] = ra(i, j);
			fgaussdata[j + i * mX] = ga(i, j);
		}
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

	char soupname[2048];
	sprintf(soupname, "%s.soup", basename.c_str());
	FILE* fout = fopen(soupname, "wb");
	int count = global_soup.size();
	fwrite(&count, sizeof(int), 1, fout);
	fwrite(&(global_soup[0]), sizeof(INDEX_TYPE), global_soup.size(), fout);
	fclose(fout);
	printf("wrote %d soup elements\n", global_soup.size());
	
	printf("writing soupvals file...\n");

	auto X = grid->XY()[0];
	auto Y = grid->XY()[1];

	GridFuncType* raw_function;
	GridFuncType* gauss_function;

	raw_function = new GridFuncType(grid, frawdata);
	gauss_function = new GridFuncType(grid, fgaussdata);

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

	delete[] fdata;
	delete[] frawdata;
	delete[] fgaussdata;
	
	delete dgb;
	delete msc;
	delete raw_function;
	delete gauss_function;

	delete raw_maxv;
	delete gauss_maxv;

	delete raw_mesh_func;
	delete gauss_mesh_func;

	delete grid;
	delete gridfunc;
	delete mesh;
	delete meshfunc;
	delete grad;

	Py_END_ALLOW_THREADS

}

void print_values(std::string str) {
	py::print("str =", str);
}

PYBIND11_MODULE(filteredlinesoup_py, m) {
	m.def("compute_sampled_soups", &compute_sampled_soups, "A function that returns a 2D NumPy array");
	m.def("print_values", &print_values, "A function that prints an integer and a string");
}


