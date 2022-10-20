
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


using namespace GInt;


typedef RegularGrid3D GridType;

#ifndef LOW_MEMORY_MODE
typedef RegularGridTrilinearFunction GridFuncType;
#else
typedef UncachedRegularGridTrilinearFunction GridFuncType;
#endif

typedef TopologicalRegularGrid3D MeshType;
typedef IndexCompareLessThan<GridFuncType> ComparerASC;
typedef NumericIntegratorPathCompressingToTerminal<AdaptiveEulerAdvector3D<GridFuncType, -1>, GridFuncType > IntegratorTypeASC;
typedef IsolatedRegionRemoverMasked<ComparerASC> RegionRemoverTypeASC;
typedef DiscreteGradientLabeling<MeshType> GradType;
#ifndef LOW_MEMORY_MODE
typedef RegularGridMaxMinVertexLabeling3D<MeshType, GridFuncType> MaxVLType;
#else
typedef UncachedMaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
#endif
typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType, 4, 6> PointwiseSteepestDescentAlgType;
typedef SlidingWindowRobinsNoalloc < RegularGrid3D, RegularGridTrilinearFunction, MeshType, MaxVLType, GradType> FastSteepestDescentAlgType;
typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> MeshFuncType;

int X, Y, Z;
int iteration_limit;
int per_x, per_y, per_z;
float error_threshold, gradient_threshold;
std::string filename;
int parallelism = -1;
int outputdebug = 0;
int hacknum = 1000;
bool usecutoff = false;
float g_pre_simp_threshold = 0.0f;
int asc_or_dsc = 0;

bool GetOptions(int argc, char** argv) {
	if (argc < 11) { printf("Usage: X Y Z filename error_threshold grad_threshold maxnumiter PresimpThesh AscOrDsc(0,1) [parallelism=ompmaxnumthreads] [outputdebug=0] [integrationinteraltimer=0]\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
	sscanf(argv[5], "%f", &error_threshold);
	sscanf(argv[6], "%f", &gradient_threshold);
	sscanf(argv[7], "%d", &iteration_limit);
	sscanf(argv[8], "%f", &g_pre_simp_threshold);
	sscanf(argv[9], "%d", &asc_or_dsc);
	if (argc >= 11)
		sscanf(argv[10], "%d", &parallelism);

	// set remaining values
	if (parallelism != -1) {
		omp_set_num_threads(parallelism);
	}

	printf("dims=(%d,%d,%d)\nfile=%s\nintegration parameters: e=%f, gt=%f, il=%d, ps=%f, ascdsc=%d,\npar=%d\n",
		X, Y, Z, argv[4], error_threshold, gradient_threshold, iteration_limit, g_pre_simp_threshold, asc_or_dsc, parallelism);

}


std::chrono::steady_clock::time_point g_start_time;
std::chrono::steady_clock::time_point task_start_time;
std::chrono::steady_clock::time_point prior_time;
std::chrono::steady_clock::time_point now_time;
FILE* gtiming;

void StartTask() {
	task_start_time = std::chrono::steady_clock::now();
	printf("starting task\n");
}

std::vector<int> timings_sequence;

void EndTaskAndRecord(const char* s) {
	now_time = std::chrono::steady_clock::now();
	// format: [global activity name] [task] [start] [end] [dration]
	timings_sequence.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count());
	fprintf(gtiming, "%s %d %d %d\n", s,
		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	printf("TIMING: %s %d %d %d\n", s,
		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
}

void bigEnd(int parallelism) {
	for (auto i : timings_sequence) {
		printf("%d ", i);
	}
	printf("\n");
}


RegularGrid3D* g_grid;
GridFuncType* g_rgt_func;
MeshType *g_topo_grid;
IntegratorTypeASC* g_num_integrator;
//IntegratorTypeWC* g_num_integrator_with_cutoff;
RegionRemoverTypeASC* g_region_remover;
#ifdef USE_REGION_CLEANER
VertexLabelingToBoundaryLabeling<INDEX_TYPE>* g_edge_map;
#else
VertexLabelingToBoundaryLabeling<int, MaxVLType>* g_edge_map;
#endif

MeshFuncType* g_topo_func;
GradType *base_grad;
//RobinsLabelingAlgorithm<MeshType, MeshFuncType> *g_robin_alg;
TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>* g_topo_alg;
MaxVLType* g_maxv_labeling;


void RecordGrad(GradType* grad, const char* gradname) {
	printf("setting dim asc man\n");
	g_topo_alg->setAscendingManifoldDimensions();
	printf("outputting to file %s\n", gradname);
	grad->outputToFile(gradname);
	//return 1;
}


int main(int argc, char** argv) {

	// read command line options
	GetOptions(argc, argv);


	// start timing overall algorithm
	ThreadedTimer timer(1);
	timer.StartGlobal();


	char gradname[2048];
	sprintf(gradname, "%s.grad", argv[4]);


	// will write timing to this file
	char timingname[2048];
	sprintf(timingname, "%s.%03d.gtime.txt", argv[4], parallelism);
	gtiming = fopen(timingname, "w");

	// START IO ---------------------------
	//StartTask();
	g_grid = new RegularGrid3D(Vec3i(X, Y, Z), Vec3b(per_x, per_y, per_z));
	g_rgt_func = new GridFuncType(g_grid);
	g_rgt_func->LoadImageFromFile(filename.c_str());
	// if want dsc man, then flip
	if (asc_or_dsc == 1)
		g_rgt_func->Negate();
	//EndTaskAndRecord(Read Data");


	g_start_time = std::chrono::steady_clock::now();
	// CREATE MAX VERTEX LABELING
	StartTask();
	g_topo_grid = new MeshType(g_grid);

	g_maxv_labeling = new MaxVLType(g_topo_grid, g_rgt_func);
	g_maxv_labeling->ComputeOutput();
	EndTaskAndRecord("Topological MaxVlabel");
	// create a topology function
	StartTask();
	printf("topo function...\n");
	g_topo_func = new MeshFuncType();
	g_topo_func->setMeshAndFuncAndMaxVLabeling(g_topo_grid, g_rgt_func, g_maxv_labeling);
	EndTaskAndRecord("Topological Function");

	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	// DO FIRST DISCRETE GRADIENT COMPUTATION WITH NO RESTRICTION
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------

	//StartTask();
	//base_grad = new GradType(g_topo_grid);
	//base_grad->ClearAllGradient();
	//PointwiseSteepestDescentAlgType* first_robins = new PointwiseSteepestDescentAlgType(g_topo_grid, g_maxv_labeling, base_grad);
	//first_robins->ComputePairing();
	//EndTaskAndRecord("Topological Robins0");

	StartTask();
	base_grad = new GradType(g_topo_grid);
	base_grad->ClearAllGradient();
	FastSteepestDescentAlgType* compare_robins = new FastSteepestDescentAlgType(g_grid, g_rgt_func, g_topo_grid, g_maxv_labeling, base_grad);
	compare_robins->ComputePairing_sliding();
	EndTaskAndRecord("Topological SlidingRobins");

	//return 1;


	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	// DO VOLUME ACCURATE GRADIENT COMPUTATION WITH NO RESTRICTION
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------

	//-------------------------------------------------------------
	// IF WE WANT NUMERIC ACCURACY, WE WILL NEED NUMEIC GRADIENT
	//-------------------------------------------------------------

	// Do gradient vectors computation from raw image data
	printf("computing gradient\n");
	StartTask();
	g_rgt_func->ComputeGradFromImage(1);
	//g_rgt_func->Negate();
	EndTaskAndRecord("NumericalTracing GradCompute");

	// CREATE RESTRICTION MAP - WILL NEED
	DenseLabeling<char>* restriction_labels = new DenseLabeling<char>(g_topo_grid->numCells());
	restriction_labels->SetAll(0);
	// we will always need a constrained robins alg
	PointwiseSteepestDescentAlgType* constrained_robins = new PointwiseSteepestDescentAlgType(g_topo_grid, g_maxv_labeling, restriction_labels, base_grad);

	// simplify the extremum graph so we can get bigger integration targets
	StartTask();
	SimplifiedExtremumGraph<MeshType, MeshFuncType, GradType>* simplified_ext_graph =
		new SimplifiedExtremumGraph<MeshType, MeshFuncType, GradType>(g_topo_grid, g_topo_func, base_grad);

	simplified_ext_graph->SetMode(SimplifiedExtremumGraph<MeshType, MeshFuncType, GradType>::EXTGRAPHMODE::MINS);

	simplified_ext_graph->ComputeMinMapFromGradient(g_pre_simp_threshold);
	EndTaskAndRecord("Topological SimplifiedGraph");
	printf("done creating simplified extremum graph\n");

	// a simplified map means that we only care about accurate boundaries for 
	// extrema that persist above a given threshold. 
	// per Julien's observation, we can compute the extremal simplification graphs
	// without spending the cost of doing a full MS Complex (even though we have the gradient)


	StartTask();
	std::unordered_map<INDEX_TYPE, INT_TYPE> extmapASC;
	for (auto p : simplified_ext_graph->mMinGraph->mCellIndexToListIdMap) {
		extmapASC[p.first] = simplified_ext_graph->mMinGraph->Representative(p.second);
	}
	GridSimplifiedExtremalRegionBuilder<ComparerASC, GridFuncType, MeshType>* test_simp_reg_builder_asc =
		new GridSimplifiedExtremalRegionBuilder<ComparerASC, GridFuncType, MeshType>(g_rgt_func, g_grid, g_topo_grid);
	test_simp_reg_builder_asc->BeginIntegration(extmapASC);
	EndTaskAndRecord("Topological SimExtRBASC");

#ifdef EXTRA_DEBUG_OUTPUT
	char label_name1[2048];
	sprintf(label_name1, "%s._terminal_%dx%dx%d.raw", filename.c_str(), g_grid->XYZ()[0], g_grid->XYZ()[1], g_grid->XYZ()[2]);
	printf("outputting labels:\n%s\n", label_name1);
	test_simp_reg_builder_asc->GetOutputLabels()->OutputToIntFile(label_name1);
#endif

	//test_simp_reg_builder_asc->GetOutputLabels()->OutputToIntFile("extremal_asc.raw");

	// now do actual numeric integration using simplified extremum regions map as target
	StartTask();
	IntegratorTypeASC* newintegrator_asc =
		new IntegratorTypeASC(g_rgt_func, g_grid, error_threshold,
			gradient_threshold, iteration_limit);
	newintegrator_asc->BeginIntegration(test_simp_reg_builder_asc->GetIdMap(),
		test_simp_reg_builder_asc->GetOutputLabels(), true);
	//newintegrator_asc->GetOutputLabels()->OutputToIntFile("newintegrator_asc.raw");
	EndTaskAndRecord("NumericalTracing NewIntegrationASC");
	//test_simp_reg_builder_asc->GetOutputLabels()->OutputToIntFile("integrated_asc.raw");
	
#ifdef EXTRA_DEBUG_OUTPUT
	char label_name2[2048];
	sprintf(label_name2, "%s._unclean_%dx%dx%d.raw", filename.c_str(), g_grid->XYZ()[0], g_grid->XYZ()[1], g_grid->XYZ()[2]);
	printf("outputting labels:\n%s\n", label_name2);
	test_simp_reg_builder_asc->GetOutputLabels()->OutputToIntFile(label_name2);
#endif

	// REMOVE DISCONNECTED COMPONENTS
	StartTask();
	IsolatedCCRegionRemoverNEW<ComparerASC, GridFuncType>* cleaner1_asc =
		new IsolatedCCRegionRemoverNEW<ComparerASC, GridFuncType>(g_rgt_func, newintegrator_asc->GetOutputLabels());
	printf("Removing Disconnected Regions\n");
	cleaner1_asc->ComputeOutput();
	printf("here2\n");
	EndTaskAndRecord("Removed CleanerASC");


	char label_name[2048];
	sprintf(label_name, "%s._int_%dx%dx%d.raw", filename.c_str(), g_grid->XYZ()[0], g_grid->XYZ()[1], g_grid->XYZ()[2]);
	printf("outputting labels:\n%s\n", label_name);
	newintegrator_asc->GetOutputLabels()->OutputToIntFile(label_name);

	VertexLabelingToBoundaryLabeling<int, MaxVLType>* g_edge_map = new
		VertexLabelingToBoundaryLabeling<int, MaxVLType>(g_topo_grid);
	g_edge_map->ComputeRegionBoundaryKind(newintegrator_asc->GetOutputLabels());

	sprintf(label_name, "%s._edges_%dx%dx%d.raw", filename.c_str(), g_topo_grid->XYZ()[0], g_topo_grid->XYZ()[1], g_topo_grid->XYZ()[2]);
	printf("outputting labels:\n%s\n", label_name);
	g_edge_map->GetOutputLabels()->OutputToFile(label_name);


};


