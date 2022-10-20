//#include "integrate3.hpp"
//#include "vector2.hpp"
//#include <stdio.h>
#include <vector>
#include <algorithm>

#include <set>
#include <queue>
#include <time.h>
#include "gi_strictly_numeric_integrator.h"
#include "gi_numeric_integrator_region_stop.h"
#include "gi_numeric_integrator_expanding_region_stop.h"
#include "gi_timing.h"
#include "gi_labeling_to_bounary_labeling.h"
#include "gi_topological_explicit_mesh_function.h"
//#include "gi_topological_region_growing_simple_gradient_builder.h"
//#include "gi_topological_convergent_gradient_builder.h"
#include "gi_robin_labeling.h"
#include "gi_adaptive_in_quad_euler_advector.h"
//#include "gi_numeric_integrator_2d_restricted_expanding_region.h"
//#include "gi_topological_2d_restricted_expanding_regions.h"
#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_isolated_region_remover.h"
#include "gi_isolated_region_remover2.h"
#include "gi_topological_utility_functions.h"
//#include "gi_morse_smale_complex_basic.h"
//#include "gi_morse_smale_complex_restricted.h"
//#include "gi_kdtree.h"
//#include "gi_msc_selectors.h"
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
//#include "gi_topological_region_growing_simple_gradient_builder.h"
//#include "gi_topological_convergent_gradient_builder.h"
#include "gi_robin_labeling.h"
#include "gi_adaptive_in_quad_euler_advector.h"
//#include "gi_numeric_integrator_2d_restricted_expanding_region.h"
//#include "gi_topological_2d_restricted_expanding_regions.h"
#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_isolated_region_remover.h"
#include "gi_topological_utility_functions.h"
//#include "gi_morse_smale_complex_basic.h"
//#include "gi_morse_smale_complex_restricted.h"
//#include "gi_kdtree.h"
//#include "gi_msc_selectors.h"
#include "gi_numeric_streamline_integrator.h"
//#include "gi_experimental.h"
//#include "gi_experimental2.h"
//#include "gi_experimental3.h"
//#include "gi_convergent_2dmanifold_gradient_builder.h"
//#include "gi_experimental5.h"
#include "gi_bifiltration_pairing.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_extrema_region_builder.h"
#include "gi_numeric_integrator_path_compressing.h"
#include "gi_fast_robins_noalloc.h"

using namespace GInt;


//typedef IsolatedRegionRemover<ComparerASC> RegionRemoverTypeASC;

//#include "integrate3.hpp"
//#include "vector2.hpp"
//#include <stdio.h>


//
//
//int *certains;
//int *dests;
//
//
//
//
//
//RegularGridTrilinearFunction* GlobalGrid;
//
//INT_TYPE highest_neighbor(Vec3i xyz) {
//	INT_TYPE tid = GlobalGrid->Index3d(xyz);
//	INT_TYPE thighest = tid;
//	Vec3i neighbors[6];
//	int nn = GlobalGrid->GatherExistingNeighborsSameBdry6(xyz, neighbors);
//	for (int i = 0; i < nn; i++) {
//		INT_TYPE tneg = GlobalGrid->Index3d(neighbors[i]);
//		if (GlobalGrid->is_greater(tneg, thighest)) thighest = tneg;
//	}
//	return thighest;
//}
//
//// does path tracing up with path compression
//int find(int s, int* a) {
//	if (a[s] == s) {
//		return s;
//	}
//	else if (a[a[s]] == s) {
//		return s;
//	}
//	a[s] = find(a[s], a);
//	return a[s];
//}

using namespace GInt;


typedef RegularGrid2D GridType;
typedef RegularGridBilinearFunction GridFuncType;

//typedef MorseSmaleComplexBasic<FLOATTYPE, MeshType, MeshFuncType, GradType> MSCType;
//typedef NumericIntegratorExpandingRegionStopWithCutoff<AdaptiveEulerAdvector3D<-1>, ComparerASC> IntegratorTypeWC;
//typedef NumericIntegratorExpandingRegionStop<AdaptiveEulerAdvector3D<-1>, ComparerASC> IntegratorTypeASC;
typedef NumericStreamlineIntegrator2D< AdaptiveEulerAdvector2D<GridFuncType, 1> > StreamlineIntegratorTypeASC;
typedef NumericStreamlineIntegrator2D<AdaptiveEulerAdvector2D<GridFuncType, -1> > StreamlineIntegratorTypeDSC;

int X, Y;
int iteration_limit;
int per_x, per_y;
float error_threshold, gradient_threshold;
std::string filename;
int parallelism = -1;
int outputdebug = 0;
int hacknum = 1000;
int stride = 1;

bool GetOptions(int argc, char** argv) {
	if (argc < 8) { printf("Usage: X Y filename error_threshold grad_threshold maxnumiter  [parallelism=ompmaxnumthreads] [outputdebug=0] [integrationinteraltimer=0]\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	filename = std::string(argv[3]);
	sscanf(argv[4], "%f", &error_threshold);
	sscanf(argv[5], "%f", &gradient_threshold);
	sscanf(argv[6], "%d", &iteration_limit);
	sscanf(argv[7], "%d", &stride);
	
	printf("dims=(%d,%d)\nfile=%s\nintegration parameters: e=%f, gt=%f, il=%d, stride=%d\n",
		X, Y, argv[3], error_threshold, gradient_threshold, iteration_limit,  stride);

}


GridType* g_grid;
GridFuncType* g_rgt_func;

//RobinsLabelingAlgorithm<MeshType, MeshFuncType> *g_robin_alg;
//TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>* g_topo_alg;
StreamlineIntegratorTypeASC* g_digitizing_streamline_integrator_asc;
StreamlineIntegratorTypeDSC* g_digitizing_streamline_integrator_dsc;



int main(int argc, char** argv) {

	// read command line options
	GetOptions(argc, argv);


	// start timing overall algorithm


	char linesname[2048];
	sprintf(linesname, "%s.slines", argv[3]);


	// will write timing to this file

	// START IO ---------------------------
	//StartTask();
	g_grid = new GridType(Vec2i(X, Y), Vec2b(per_x, per_y));
	g_rgt_func = new GridFuncType(g_grid);
	g_rgt_func->LoadImageFromFile(filename.c_str());
	//EndTaskAndRecord(Read Data");

	printf("here\n");
	g_rgt_func->ComputeGradFromImage(1);
	
	g_digitizing_streamline_integrator_asc = new StreamlineIntegratorTypeASC(g_grid, g_rgt_func, error_threshold, gradient_threshold, iteration_limit);
	g_digitizing_streamline_integrator_dsc = new StreamlineIntegratorTypeDSC(g_grid, g_rgt_func, error_threshold, gradient_threshold, iteration_limit);

	
	//g_rgt_func->Negate();
	printf("here %s\n", linesname);
	FILE* fout = fopen(linesname, "wb");
	printf("here 2\n");
#pragma omp parallel for
	for (int i = stride; i < X; i += stride) {
		for (int j = stride; j < Y; j += stride) {
			
			Vec2d seed(i, j);
			float wigglex = (rand() % 101) / 103.f - .51f;
			float wiggley = (rand() % 101) / 103.f - .51f;
			seed += Vec2d(wigglex* stride, wiggley * stride);
			seed.PrintFloat();
			std::vector<Vec2d> pointsA;
			std::vector<Vec2d> pointsD;
		//	printf("%d %d\n", i, j);
			g_digitizing_streamline_integrator_asc->IntegrateStreamline(seed, pointsA);
			g_digitizing_streamline_integrator_dsc->IntegrateStreamline(seed, pointsD);
			std::reverse(pointsD.begin(), pointsD.end());

#pragma omp critical
			{
				int num = pointsA.size() + pointsD.size() - 1 ;
				if (num > 1) {
					fwrite(&num, sizeof(int), 1, fout);
					fwrite(pointsD.data(), sizeof(Vec2d), pointsD.size() - 1, fout);
					fwrite(pointsA.data(), sizeof(Vec2d), pointsA.size(), fout);
				}
			}

		}

	}
	fclose(fout);


	return 1;
};


