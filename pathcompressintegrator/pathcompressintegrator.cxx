//#include "integrate3.hpp"
//#include "vector2.hpp"
//#include <stdio.h>
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
#include "gi_topological_region_growing_simple_gradient_builder.h"
#include "gi_topological_convergent_gradient_builder.h"
#include "gi_robin_labeling.h"
#include "gi_adaptive_in_quad_euler_advector.h"
#include "gi_numeric_integrator_2d_restricted_expanding_region.h"
#include "gi_topological_2d_restricted_expanding_regions.h"
#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_isolated_region_remover2.h"
#include "gi_topological_utility_functions.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_morse_smale_complex_restricted.h"
#include "gi_kdtree.h"
#include "gi_msc_selectors.h"
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
#include "gi_topological_region_growing_simple_gradient_builder.h"
#include "gi_topological_convergent_gradient_builder.h"
#include "gi_robin_labeling.h"
#include "gi_adaptive_in_quad_euler_advector.h"
#include "gi_numeric_integrator_2d_restricted_expanding_region.h"
#include "gi_topological_2d_restricted_expanding_regions.h"
#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_isolated_region_remover.h"
#include "gi_topological_utility_functions.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_morse_smale_complex_restricted.h"
#include "gi_kdtree.h"
#include "gi_msc_selectors.h"
#include "gi_numeric_streamline_integrator.h"
#include "gi_experimental.h"
#include "gi_experimental2.h"
#include "gi_experimental3.h"
#include "gi_convergent_2dmanifold_gradient_builder_grad_align.h"
#include "gi_experimental5.h"
#include "gi_max_vertex_labeling.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_bifiltration_pairing.h"

using namespace GInt;
//typedef IndexCompareLessThan Comparer;
//#ifdef BADER_CONSTRAINT
//typedef TopologicalRegularGridRestricted GridType;
//#else
//typedef TopologicalRegularGridRestricted GridType;
//#endif
//typedef MorseSmaleComplexBasic<FLOATTYPE, GridType, TopologicalExplicitDenseMeshFunction, DiscreteGradientLabeling> MSCType;
//typedef NumericIntegratorExpandingRegionStopFiltered2<AdaptiveEulerAdvector<-1>, Comparer> IntegratorType;
//typedef IsolatedRegionRemover<Comparer> RegionRemoverType;

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
//	int nn = GlobalGrid->GatherExistingNeighbors6(xyz, neighbors);
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

// rename discretegradientlabeling to regulargriddiscretegradientlabeling

using namespace GInt;
typedef IndexCompareLessThan Comparer;
typedef RegularGrid GridType;
typedef TopologicalRegularGridRestricted MeshType;
typedef RegularGridTrilinearFunction GridFuncType;
typedef DiscreteGradientLabeling<MeshType> GradType;

//typedef NumericIntegratorExpandingRegionStopWithCutoff<AdaptiveEulerAdvector<-1>, Comparer> IntegratorTypeWC;
typedef NumericIntegratorExpandingRegionStopFiltered2<AdaptiveEulerAdvector<-1>, Comparer> IntegratorType;
//typedef NumericIntegratorRegionStop<AdaptiveEulerAdvector<-1>, Comparer> IntegratorType2;
//typedef NumericStreamlineIntegrator<AdaptiveEulerAdvector<-1> > StreamlineIntegratorType;
typedef IsolatedRegionRemoverMasked<Comparer> RegionRemoverType;
//typedef NumericStreamlineIntegrator<AdaptiveEulerAdvector<-1> > StreamlineIntegratorType;
typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLabelingType;
typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLabelingType, GridFuncType, float> TopoFuncType;
typedef BoundaryEdgeGraphSparse<MeshType, TopoFuncType> BoundaryEdgeGraphType;
typedef Convergent2ManifoldGradientBuilderGradAlign<BoundaryEdgeGraphType, MeshType, TopoFuncType, GridFuncType, MaxVLabelingType> ManifoldLabelerType;

int HashingWaveFront::count = 0;
int HashingWaveFrontTry2::count = 0;
DenseLabeling<char>* HashingWaveFrontTry2::mProcessed = NULL;
int ExplicitWaveFrontTry2::count = 0;
DenseLabeling<char>* ExplicitWaveFrontTry2::mProcessed = NULL;


int X, Y, Z;
int iteration_limit;
int per_x, per_y, per_z;
float error_threshold, gradient_threshold;
std::string filename;
int parallelism = -1;
int outputdebug = 0;
int hacknum = 1000;
//bool usecutoff = false;
FLOATTYPE cutoffvalue;
int internaltimer = 0;

bool GetOptions(int argc, char** argv) {
	if (argc < 11) { printf("Usage: X Y Z filename grad_threshold error_threshold maxnumiter per_x per_y per_z [parallelism=ompmaxnumthreads] [outputdebug=0]\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
	sscanf(argv[5], "%f", &error_threshold);
	sscanf(argv[6], "%f", &gradient_threshold);
	sscanf(argv[7], "%d", &iteration_limit);
	sscanf(argv[8], "%d", &per_x);
	sscanf(argv[9], "%d", &per_y);
	sscanf(argv[10], "%d", &per_z);
	if (argc >= 12)
		sscanf(argv[11], "%d", &parallelism);
	if (argc >= 13)
		sscanf(argv[12], "%d", &outputdebug);
	if (parallelism != -1) {
		omp_set_num_threads(parallelism);
	}

	sscanf(argv[13], "%d", &hacknum);
	sscanf(argv[14], "%d", &internaltimer);
}
const char* ITimingStrings[] = {
	"FindMinima",
	"ExpandCertain",
	"EarlyTerminationIntegration",
	"ReintegrateBoundaries",
	"Done"
};

const char* TimingStrings[] = {
	"GlobalIO",
	"IExpand",
	"Integration",
	"Wait",
	"Cleanup",
	"ExpandRegion",
	"MergeRegions",
	"CopyGlobal",
	"Robins"
};

int main(int argc, char** argv) {

	GetOptions(argc, argv);





	ThreadedTimer timer(1);
	timer.StartGlobal();





	GridType* m_grid;
	//RegularGridTrilinearFunction* m_func;
	GridFuncType* g_grid_function;
	MeshType *m_tgrid;
	IntegratorType* i2;
	//IntegratorTypeWC* i2wc;
	RegionRemoverType* i2clean;
	VertexLabelingToBoundaryLabeling<INDEX_TYPE>* edgemap;
	TopoFuncType* m_topofunc;
	GradType *labeling;
	MyRobinsNoalloc<MeshType, MaxVLabelingType, GradType> *robin;
	TopologicalGradientUsingAlgorithms<MeshType, TopoFuncType, GradType>* topo_algs;

	NumericIntegratorNew* sinteegrator;
	//StreamlineIntegratorType* sintegrator;
	TerminateNearExtrema* extremumtermination;


	// added by Harsh

	char timingname[2048];
	sprintf(timingname, "%s.%03d.%04d.gtime.txt", argv[4], parallelism, hacknum);
	FILE* gtiming = fopen(timingname, "w");
	// format: [global activity name] [task] [start] [end] [dration]

	// START IO ---------------------------
	// -- start timing IO
	//prior_time = std::chrono::steady_clock::now();
	m_grid = new GridType(Vec3i(X, Y, Z), Vec3b(per_x, per_y, per_z));
	g_grid_function = new GridFuncType(m_grid);
	g_grid_function->LoadImageFromFile(filename.c_str());

	// START TIMING
	std::chrono::steady_clock::time_point g_start_time = std::chrono::steady_clock::now();;
	std::chrono::steady_clock::time_point task_start_time;
	std::chrono::steady_clock::time_point prior_time;
	std::chrono::steady_clock::time_point now_time;
	// -- end timing IO
	//now_time = std::chrono::steady_clock::now();
	//fprintf(gtiming, "IO Read %d %d %d\n", 
	//	std::chrono::duration_cast<std::chrono::milliseconds>(prior_time - g_start_time).count(),
	//	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
	//	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - prior_time).count());
	// END IO ------------------------------

	// START GRAD INTEGRATION --------------
	prior_time = std::chrono::steady_clock::now();
	task_start_time = std::chrono::steady_clock::now();


	g_grid_function->ComputeGradFromImage(1);
	now_time = std::chrono::steady_clock::now();
	fprintf(gtiming, "NumericalTracing GradCompute %d %d %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	task_start_time = std::chrono::steady_clock::now();

	//g_grid_function->Negate();
	m_tgrid = new MeshType(m_grid);

	DenseLabeling<INDEX_TYPE>* outlabels;
	DenseLabeling<INT_TYPE>* otherlabels;
	//if (usecutoff) {
	//	i2wc = new IntegratorTypeWC(g_grid_function, m_grid, error_threshold, gradient_threshold, iteration_limit, cutoffvalue); //
	//	i2wc->BeginIntegration();
	//	outlabels = i2wc->GetOutputLabels();
	//}
	//else {
		i2 = new IntegratorType(g_grid_function, m_grid, error_threshold, gradient_threshold, iteration_limit); //
		if (internaltimer != 0) {
			ThreadedTimer expand_timer(parallelism);
			expand_timer.StartGlobal();

			i2->BeginIntegration(&expand_timer);

			expand_timer.EndGlobal();
			char tname[1024];
			sprintf(tname, "%s.%d.itime.txt", argv[4], parallelism);
			expand_timer.PrintAllToFile(tname, ITimingStrings);
		}
		else {
			i2->BeginIntegration();
		}
		otherlabels = i2->GetOutputLabels();
	//}

	i2clean = new RegionRemoverType(g_grid_function, otherlabels, &(i2->GetExtrema()));
	i2clean->ComputeOutput();

	////if (outputdebug) {
	////	outlabels->OutputToIntFile("classes.raw");
	////}
	//printf("done with numeric integration!\ncreating vertex to boundary labeling...\n");
	//now_time = std::chrono::steady_clock::now();
	//fprintf(gtiming, "NumericalTracing Integration %d %d %d\n",
	//	std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
	//	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
	//	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	//task_start_time = std::chrono::steady_clock::now();

	//printf("removing regions...\n");
	//i2clean = new RegionRemoverType(g_grid_function, outlabels);
	//i2clean->ComputeOutput();

	//if (usecutoff) {
	//	delete i2wc;
	//}
	//else {
	//	delete i2;
	//}
	//if (outputdebug) {
	//	i2clean->GetOutputLabels()->OutputToIntFile("regions_clean.raw");
	//}
	// END GRAD INTEGRATION ----------------
	now_time = std::chrono::steady_clock::now();
	//fprintf(gtiming, "NumericalTracing Cleaning %d %d %d\n",
	//	std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
	//	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
	//	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	task_start_time = std::chrono::steady_clock::now();
	fprintf(gtiming, "NumericalTracing Overall %d %d %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(prior_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - prior_time).count());

	// start topology here ---------------
	prior_time = std::chrono::steady_clock::now();
	task_start_time = std::chrono::steady_clock::now();

	printf("edge labeling\n");
	edgemap = new VertexLabelingToBoundaryLabeling<INDEX_TYPE>(i2clean->GetOutputLabels(), m_tgrid);
	edgemap->ComputeBoundary();


	//delete i2clean;


	m_tgrid->set_restriction(edgemap->GetOutputLabels());
	now_time = std::chrono::steady_clock::now();
	fprintf(gtiming, "RegionGrow EdgeMap %d %d %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	task_start_time = std::chrono::steady_clock::now();

	if (outputdebug) {
		edgemap->GetOutputLabels()->OutputToFile("boundary_labels.raw");
	}
	//printf("done with boundary map creation!\ncreating topology function and copying values...\n");
	//printf("\n --------------------------------------- \n");

	//printf(" Creating TopologicalExplicitDenseMeshFunction ...");
	//fflush(stdout);
	MaxVLabelingType* maxVLabeling = new MaxVLabelingType(m_tgrid, g_grid_function);
	maxVLabeling->ComputeOutput();


	printf("topo function...\n");
	m_topofunc = new TopoFuncType();
	m_topofunc->setMeshAndFuncAndMaxVLabeling(m_tgrid, g_grid_function, maxVLabeling);
	//m_topofunc->setMeshAndAllocate(m_tgrid);
	//m_topofunc->copyVertexValuesFromGridFunction(g_grid_function);
	//m_topofunc->setCellValuesMaxOfVerts();







	now_time = std::chrono::steady_clock::now();
	fprintf(gtiming, "RegionGrow TopoFunc %d %d %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	task_start_time = std::chrono::steady_clock::now();

	printf("expansion...\n");
	BoundaryEdgeGraphType* d = new BoundaryEdgeGraphType(m_tgrid, m_topofunc);
	d->Init();

	edgemap->GetOutputLabels()->OutputToFile("before_expand.raw");
	ManifoldLabelerType* texpand =
		new ManifoldLabelerType(d, m_tgrid, m_topofunc, maxVLabeling, g_grid_function, edgemap->GetOutputLabels());
	texpand->hacknum = hacknum;
	//Convergent2ManifoldGradientBuilder/*< BoundaryEdgeGraph>*/* texpand = new Convergent2ManifoldGradientBuilder/*< BoundaryEdgeGraph>*/(d, m_tgrid, m_topofunc, labeling);

	if (internaltimer == 2) {
		ThreadedTimer expand_timer(parallelism);
		expand_timer.StartGlobal();

		//texpand->BeginIntegration(&expand_timer);

		expand_timer.EndGlobal();
		expand_timer.PrintAllToFile("test_timings.txt", TimingStrings);
	}
	else {
		//texpand->BeginIntegration();

	}



	if (outputdebug) {
		edgemap->GetOutputLabels()->OutputToFile("after_expand.raw");
	}

	//delete d;
	//delete texpand;











	now_time = std::chrono::steady_clock::now();
	fprintf(gtiming, "RegionGrow Assign %d %d %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	task_start_time = std::chrono::steady_clock::now();

	printf("robins...\n");
	labeling = new GradType(m_tgrid);
	labeling->ClearAllGradient();
	robin = new MyRobinsNoalloc<MeshType, MaxVLabelingType, GradType>(m_tgrid, maxVLabeling, edgemap->GetOutputLabels(),labeling);
	robin->ComputePairing();

	now_time = std::chrono::steady_clock::now();
	fprintf(gtiming, "RegionGrow Robins %d %d %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	task_start_time = std::chrono::steady_clock::now();

	fprintf(gtiming, "RegionGrow Overall %d %d %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(prior_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - prior_time).count());
	//printf("final write\n");


	//topo_algs = new TopologicalGradientUsingAlgorithms(m_topofunc, m_tgrid, labeling);
	//topo_algs->setAscendingManifoldDimensions();
	//topo_algs->CheckGradientConsistency();

	topo_algs = new TopologicalGradientUsingAlgorithms<MeshType, TopoFuncType, GradType>(m_topofunc, m_tgrid, labeling);
	topo_algs->setAscendingManifoldDimensions();
	char gradname[1024];
	sprintf(gradname, "%s.grad", argv[4]);
	labeling->outputToFile(gradname);
	//fprintf(gtiming, "IO Write %d %d %d\n",
	//	std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
	//	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
	//	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	task_start_time = std::chrono::steady_clock::now();

	fprintf(gtiming, "Total Overall %d %d %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(g_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count());


	fclose(gtiming);


	return 1;
};




