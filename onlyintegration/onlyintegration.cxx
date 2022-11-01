/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

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
#include "gi_isolated_region_remover.h"
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
#include "gi_convergent_2dmanifold_gradient_builder.h"
#include "gi_experimental5.h"
#include "gi_bifiltration_pairing.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_extrema_region_builder.h"
#include "gi_numeric_integrator_path_compressing.h"

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

typedef RegularGridTrilinearFunction GridFuncType;
typedef IndexCompareLessThan<GridFuncType> ComparerASC;
typedef IndexCompareGreaterThan<GridFuncType> ComparerDSC;
typedef TopologicalRegularGrid3D GridType;
//typedef MorseSmaleComplexBasic<FLOATTYPE, GridType, FuncType, GradType> MSCType;
//typedef NumericIntegratorExpandingRegionStopWithCutoff<AdaptiveEulerAdvector3D<-1>, ComparerASC> IntegratorTypeWC;
//typedef NumericIntegratorExpandingRegionStop<AdaptiveEulerAdvector3D<-1>, ComparerASC> IntegratorTypeASC;
typedef NumericIntegratorPathCompressingToTerminal<AdaptiveEulerAdvector3D<GridFuncType, -1>, GridFuncType > IntegratorTypeASC;
typedef NumericIntegratorPathCompressingToTerminal<AdaptiveEulerAdvector3D<GridFuncType, 1>, GridFuncType > IntegratorTypeDSC;
typedef IsolatedRegionRemoverMasked<ComparerASC> RegionRemoverTypeASC;
typedef IsolatedRegionRemoverMasked<ComparerDSC> RegionRemoverTypeDSC;
typedef DigitizingNumericStreamlineIntegrator3dASC<GridType, GridFuncType, AdaptiveEulerAdvector3D<GridFuncType, 1>  > StreamlineIntegratorTypeASC;
typedef DigitizingNumericStreamlineIntegrator3dDSC<GridType, GridFuncType  , AdaptiveEulerAdvector3D<GridFuncType, -1>> StreamlineIntegratorTypeDSC;
typedef DiscreteGradientLabeling<GridType> GradType;
typedef RegularGridMaxMinVertexLabeling3D<GridType, GridFuncType> MaxVLType;
typedef MyRobinsNoalloc<GridType, MaxVLType, GradType, 4, 6> RobinsType;
typedef TopologicalMaxVertexMeshFunction<GridType, MaxVLType, GridFuncType, float> FuncType;

int X, Y, Z;
int iteration_limit;
int per_x, per_y, per_z;
float error_threshold, gradient_threshold;
std::string filename;
int parallelism = -1;
int outputdebug = 0;
int hacknum = 1000;
bool usecutoff = false;
FLOATTYPE cutoffvalue;
int internaltimer = 0;

bool GetOptions(int argc, char** argv) {
	if (argc < 11) { printf("Usage: X Y Z filename grad_threshold error_threshold maxnumiter per_x per_y per_z [parallelism=ompmaxnumthreads] [outputdebug=0] [integrationinteraltimer=0]\n"); return 0; }
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
	if (argc >= 14) {
		sscanf(argv[14], "%d", &internaltimer);
	}
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


std::chrono::steady_clock::time_point g_start_time;
std::chrono::steady_clock::time_point task_start_time;
std::chrono::steady_clock::time_point prior_time;
std::chrono::steady_clock::time_point now_time;
FILE* gtiming;

void StartTask() {
	task_start_time = std::chrono::steady_clock::now();
	printf("starting task\n");
}

void EndTaskAndRecord(const char* s) {
	now_time = std::chrono::steady_clock::now();
	// format: [global activity name] [task] [start] [end] [dration]
	fprintf(gtiming, "%s %d %d %d\n", s,
		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	printf("%s %d %d %d\n", s,
		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
}

RegularGrid3D* g_grid;
GridFuncType* g_rgt_func;
GridType *g_topo_grid;
IntegratorTypeASC* g_num_integrator;
//IntegratorTypeWC* g_num_integrator_with_cutoff;
RegionRemoverTypeASC* g_region_remover;
#ifdef USE_REGION_CLEANER
VertexLabelingToBoundaryLabeling<INDEX_TYPE>* g_edge_map;
#else
VertexLabelingToBoundaryLabeling<int, MaxVLType>* g_edge_map;
#endif

FuncType* g_topo_func;
GradType *base_grad;
//RobinsLabelingAlgorithm<GridType, FuncType> *g_robin_alg;
TopologicalGradientUsingAlgorithms<GridType, FuncType, GradType>* g_topo_alg;
StreamlineIntegratorTypeASC* g_digitizing_streamline_integrator_asc;
StreamlineIntegratorTypeDSC* g_digitizing_streamline_integrator_dsc;


//void CombineLabels() {
//	INDEX_TYPE num = g_topo_grid->numCells();
//
//#pragma omp parallel for
//	for (INDEX_TYPE i = 0; i < num; i++) {
//		char v1 = g_digitizing_streamline_integrator_asc->get_output()->GetLabel(i);
//		char v2 = g_edge_map->GetOutputLabels()->GetLabel(i);
//
//		g_edge_map->GetOutputLabels()->SetLabel(i, max(v1, v2));
//	}
//}

class NeighborState {

public:
	GradType::GradBitfield grad[5 * 5 * 5];
	char labs[5 * 5 * 5];

	NeighborState() {
		for (int i = 0; i < 5 * 5 * 5; i++) {
			labs[i] = 0;
		}
	}

};

class VertexNeighborhood {
public:
	std::vector<NeighborState> states;
	float vert_v[27];
	Vec3d vec_v[27];
	char lstar[5 * 5 * 5];
	VertexNeighborhood() {
		for (int i = 0; i < 5 * 5 * 5; i++) {
			lstar[i] = 0;
		}
	}

};


void ReIntegrateUpFrom2Saddles() {

	// FIRST gather all the critical 2-saddles from the discrete gradient
	StartTask();
	std::vector<INDEX_TYPE> topo_index_partition;
	int num_threads;
	std::vector<std::pair<float, INDEX_TYPE> > criticals;
	printf("gothere 2\n");
#pragma omp parallel
	{
#pragma omp single
		{
			num_threads = omp_get_num_threads();
			ArrayIndexPartitioner::EvenChunkSplit(g_topo_grid->numCells(), num_threads, topo_index_partition);
		}

		int thread_num = omp_get_thread_num();
		// in parallel go through and find all 2-saddles
		std::vector<std::pair<float, INDEX_TYPE> > lcriticals;
		GridType::DCellsIterator face_iterator(g_topo_grid, 2, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
		for (face_iterator.begin(); face_iterator.valid(); face_iterator.advance()) {
			INDEX_TYPE cell_id = face_iterator.value();

			if (base_grad->getCritical(cell_id)) {
				std::pair<float, INDEX_TYPE> p(g_topo_func->cellValue(cell_id), cell_id);

				lcriticals.push_back(p);

			}
		}
#pragma omp critical
		{
			criticals.insert(criticals.end(), lcriticals.begin(), lcriticals.end());
		}
	}
	printf("gothere 3\n");
	EndTaskAndRecord("Topological GatherCritical2Saddles");

	StartTask();
	std::sort(criticals.begin(), criticals.end());
	EndTaskAndRecord("Topological Sort2Saddles");
	//for (auto p : criticals) printf("cp: %f %llu\n", p.first, p.second);

	// so criticals is now sorted list, highest value last

	INDEX_TYPE total_left = criticals.size() - 1;
#pragma omp parallel 
	{
		while (true) {
			INDEX_TYPE local_id;
			//int tots;
#pragma omp critical
			{
				local_id = total_left;
				total_left -= 1;
				//tots = rand() % 10000;
				//printf("%d doing %llu %llu\n", omp_get_thread_num(), local_id, total_left);
			}

			if (local_id < 0) break;

			//std::vector<int> fff;
			//fff.clear();
			//for (int i = 0; i < tots; i++) {
			//	fff.push_back(i * i);
			//}
			//printf("fff size %d\n", fff.size());

			INDEX_TYPE sad_id = criticals[local_id].second;
			// 
			//// now find each arc
			std::vector<INDEX_TYPE> result;
			std::queue<INDEX_TYPE> cell_queue;
			cell_queue.push(sad_id);

			// THIS IS A SUPER FAST WAY OF filling in gometry of arc... only 2 in each direction
			// gather the 4 hexes along paths on either side of the critical saddle
			std::set<INDEX_TYPE> cell_visited;
			int counter = 4;
			while (!cell_queue.empty() && counter >= 0) {
				INDEX_TYPE current = cell_queue.front();
				cell_queue.pop();

				cell_visited.insert(current);
				//result.push_back(current);

				GridType::CofacetsIterator cofacets(g_topo_grid);
				for (cofacets.begin(current); cofacets.valid(); cofacets.advance()) {
					INDEX_TYPE temp_id = cofacets.value();

					if (base_grad->getCritical(temp_id) &&
						cell_visited.count(temp_id) == 0) {
						result.push_back(temp_id);
						cell_visited.insert(temp_id);
					}
					else if (cell_visited.count(temp_id) == 0) {
						INDEX_TYPE pair = base_grad->getPair(temp_id);
						result.push_back(temp_id);
						//result.push_back(pair);
						cell_visited.insert(temp_id);
						cell_visited.insert(pair);
						cell_queue.push(pair);
					}
				}
				counter--;
			}

			//printf("result size %d\n", result.size());
			for (auto arc_hex_id : result) {
				//printf("%llu ", arc_hex_id);
				if (g_topo_grid->dimension(arc_hex_id) != 3) continue;
				std::vector<Vec3d> points;
				std::vector<INDEX_TYPE> dline;
				Vec3d seed; g_topo_grid->centroid(arc_hex_id, seed);
				seed = seed * 0.5; // back to grid coordinates
				g_digitizing_streamline_integrator_asc->IntegrateStreamline(seed, points, dline);
			}
			//printf("\n");
			g_digitizing_streamline_integrator_asc->set_label(sad_id);



		}

	}


}

void ReIntegrateDownFrom1Saddles() {

	// FIRST gather all the critical 2-saddles from the discrete gradient
	StartTask();
	std::vector<INDEX_TYPE> topo_index_partition;
	int num_threads;
	std::vector<std::pair<float, INDEX_TYPE> > criticals;
	printf("gothere 2\n");
#pragma omp parallel
	{
#pragma omp single
		{
			num_threads = omp_get_num_threads();
			ArrayIndexPartitioner::EvenChunkSplit(g_topo_grid->numCells(), num_threads, topo_index_partition);
		}

		int thread_num = omp_get_thread_num();
		// in parallel go through and find all 2-saddles
		std::vector<std::pair<float, INDEX_TYPE> > lcriticals;
		GridType::DCellsIterator face_iterator(g_topo_grid, 1, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
		for (face_iterator.begin(); face_iterator.valid(); face_iterator.advance()) {
			INDEX_TYPE cell_id = face_iterator.value();

			if (base_grad->getCritical(cell_id)) {
				std::pair<float, INDEX_TYPE> p(-1 * g_topo_func->cellValue(cell_id), cell_id);

				lcriticals.push_back(p);

			}
		}
#pragma omp critical
		{
			criticals.insert(criticals.end(), lcriticals.begin(), lcriticals.end());
		}
	}
	printf("gothere 3\n");
	EndTaskAndRecord("Topological GatherCritical1Saddles");

	StartTask();
	std::sort(criticals.begin(), criticals.end());
	EndTaskAndRecord("Topological Sort2Saddles");
	//for (auto p : criticals) printf("cp: %f %llu\n", p.first, p.second);

	// so criticals is now sorted list, highest value last

	INDEX_TYPE total_left = criticals.size() - 1;
#pragma omp parallel 
	{
		while (true) {
			INDEX_TYPE local_id;
			//int tots;
#pragma omp critical
			{
				local_id = total_left;
				total_left -= 1;
				//tots = rand() % 10000;
				//printf("%d doing %llu %llu\n", omp_get_thread_num(), local_id, total_left);
			}

			if (local_id < 0) break;

			//std::vector<int> fff;
			//fff.clear();
			//for (int i = 0; i < tots; i++) {
			//	fff.push_back(i * i);
			//}
			//printf("fff size %d\n", fff.size());

			INDEX_TYPE sad_id = criticals[local_id].second;
			// 
			//// now find each arc
			std::vector<INDEX_TYPE> result;
			std::queue<INDEX_TYPE> cell_queue;
			cell_queue.push(sad_id);

			// THIS IS A SUPER FAST WAY OF filling in gometry of arc... only 2 in each direction
			// gather the 4 vertices along paths on either side of the critical saddle
			std::set<INDEX_TYPE> cell_visited;
			int counter = 4;
			while (!cell_queue.empty() && counter >= 0) {
				INDEX_TYPE current = cell_queue.front();
				cell_queue.pop();

				cell_visited.insert(current);
				//result.push_back(current);

				GridType::FacetsIterator facets(g_topo_grid);
				for (facets.begin(current); facets.valid(); facets.advance()) {
					INDEX_TYPE temp_id = facets.value();

					if (base_grad->getCritical(temp_id) &&
						cell_visited.count(temp_id) == 0) {
						result.push_back(temp_id);
						cell_visited.insert(temp_id);
					}
					else if (cell_visited.count(temp_id) == 0) {
						INDEX_TYPE pair = base_grad->getPair(temp_id);
						result.push_back(temp_id);
						//result.push_back(pair);
						cell_visited.insert(temp_id);
						cell_visited.insert(pair);
						cell_queue.push(pair);
					}
				}
				counter--;
			}

			//printf("result size %d\n", result.size());
			for (auto arc_vert_id : result) {
				//printf("%llu ", arc_hex_id);
				if (g_topo_grid->dimension(arc_vert_id) != 0) continue;
				std::vector<Vec3d> points;
				std::vector<INDEX_TYPE> dline;
				Vec3d seed; g_topo_grid->centroid(arc_vert_id, seed);
				seed = seed * 0.5; // back to grid coordinates
				g_digitizing_streamline_integrator_dsc->IntegrateStreamline(seed, points, dline);
			}
			//printf("\n");
			g_digitizing_streamline_integrator_dsc->set_label(sad_id);



		}

	}


}

//// we need asc man 3 dsc man 3?
//bool needasc3 = false;
//bool needdsc3 = false;
//bool needasc1 = false;
//bool needdsc1 = false;
//bool needsad = false;
//


GradType* stageGrad[4];
DenseLabeling<char>* stageLabel[4];

int main(int argc, char** argv) {

	// read command line options
	GetOptions(argc, argv);

	// start timing overall algorithm
	ThreadedTimer timer(1);
	timer.StartGlobal();


	char gradname[1024];
	sprintf(gradname, "%s.grad", argv[4]);


	// will write timing to this file
	char timingname[2048];
	sprintf(timingname, "%s.%03d.gtime.txt", argv[4], parallelism);
	gtiming = fopen(timingname, "w");

	// START IO ---------------------------
	StartTask();
	g_grid = new RegularGrid3D(Vec3i(X, Y, Z), Vec3b(per_x, per_y, per_z));
	g_rgt_func = new RegularGridTrilinearFunction(g_grid);
	g_rgt_func->LoadImageFromFile(filename.c_str());
	EndTaskAndRecord("I/O - Read Data");

	// CREATE MAX VERTEX LABELING
	StartTask();
	g_topo_grid = new GridType(g_grid);
	DenseLabeling<char>* nolabel = new 	DenseLabeling<char>(g_topo_grid->numCells());
	nolabel->SetAll(0);
	MaxVLType* g_maxv_labeling = new MaxVLType(g_topo_grid, g_rgt_func);
	g_maxv_labeling->ComputeOutput();
	EndTaskAndRecord("Topological MaxVlabel");
	// create a topology function
	StartTask();
	printf("topo function...\n");
	g_topo_func = new FuncType();
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

	StartTask();
	base_grad = new GradType(g_topo_grid);
	base_grad->ClearAllGradient();
	RobinsType* first_robins = new RobinsType(g_topo_grid, g_maxv_labeling, base_grad);
	first_robins->ComputePairing();
	EndTaskAndRecord("Topological Robins0");

	// so i guess the max is not always on the +1,+1,+1 coordinate of the vertex whose lower star
	// it is part of...
	//typename GridType::DCellsIterator hexes(g_topo_grid, 3);
	//for (hexes.begin(); hexes.valid(); hexes.advance()) {
	//	INDEX_TYPE hid = hexes.value();
	//	if (base_grad->getCritical(hid)) {
	//		Vec3l hc; g_topo_grid->cellid2Coords(hid, hc);
	//		INDEX_TYPE vid = g_maxv_labeling->Cell2HighestVertex(hid);
	//		Vec3l vc; g_topo_grid->cellid2Coords(vid, vc);
	//		hc.PrintInt();
	//		vc.PrintInt();
	//		printf("\n");
	//	}
	//}

	//return 1;

	g_topo_alg = new TopologicalGradientUsingAlgorithms<GridType, FuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
	printf("after base first robins:\n");
	g_topo_alg->count_critical_points(4);

	//g_topo_alg->setAscendingManifoldDimensions();
	//base_grad->outputToFile(gradname);

	// === some debugging here ===
	std::set<INDEX_TYPE> original_cp;
	GridType::AllCellsIterator allit1(g_topo_grid);
	for (allit1.begin(); allit1.valid(); allit1.advance()) {
		INDEX_TYPE cid = allit1.value();
		if (base_grad->getCritical(cid)) original_cp.insert(cid);
	}


	stageGrad[0] = new GradType(g_topo_grid);
	stageGrad[0]->copyValues(base_grad);
	stageLabel[0] = nolabel;


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
	g_rgt_func->ComputeGradFromImage(2);
	//g_rgt_func->Negate();
	EndTaskAndRecord("NumericalTracing GradCompute");

	//-------------------------------------------------------------
	// IF WE WANT ACCURATE 3-MANIFOLDS, COMPUTE SIMPLIFIED MAPS
	//-------------------------------------------------------------
	{
		// a simplified map means that we only care about accurate boundaries for 
		// extrema that persist above a given threshold. 
		// per Julien's observation, we can compute the extremal simplification graphs
		// without spending the cost of doing a full MS Complex (even though we have the gradient)

		// create a map that takes extremum to its simplified representative.
		StartTask();
		printf("testing....\n");
		SimplifiedExtremumGraph<GridType, FuncType, GradType>* test =
			new SimplifiedExtremumGraph<GridType, FuncType, GradType>(g_topo_grid, g_topo_func, base_grad);
		test->SetMode(SimplifiedExtremumGraph<GridType, FuncType, GradType>::EXTGRAPHMODE::BOTH);
		printf("about to create thingy\n");

		test->ComputeMinMapFromGradient(0.001003);
		EndTaskAndRecord("Simplified Min/Max Merge Graph");
		printf("done creating thingy\n");
		FILE* ftest = fopen("TESTMERGE.txt", "w");
		for (int i = 0; i < test->mMinGraph->NumExtrema(); i++) {
			fprintf(ftest, "%llu %llu\n", test->mMinGraph->TopoIndex(i), test->mMinGraph->TopoIndex(test->mMinGraph->Representative(i)));
		}
		fclose(ftest);
		printf("done testing\n");
		//UFMergeGraph<FuncType>* asdf = new UFMergeGraph<FuncType>();

		// NOW BUILD A TERMINAL MAP FROM THE SIMPLIFIED MIN/MAX HIERARCHIES
		// DO ASCENDING MANIFOLDS
		StartTask();
		unordered_map<INDEX_TYPE, INT_TYPE> extmapASC;
		for (auto p : test->mMinGraph->mCellIndexToListIdMap) {
			extmapASC[p.first] = test->mMinGraph->Representative(p.second);
		}
		GridSimplifiedExtremalRegionBuilder<ComparerASC, GridFuncType, GridType>* test_simp_reg_builder_asc =
			new GridSimplifiedExtremalRegionBuilder<ComparerASC, GridFuncType, GridType>(g_rgt_func, g_grid, g_topo_grid);
		test_simp_reg_builder_asc->BeginIntegration(extmapASC);
		EndTaskAndRecord("Simplified Min/Max Termination Labeling - ASC");

		//test_simp_reg_builder_asc->GetOutputLabels()->OutputToIntFile("simp_terminal_regions_asc.raw");
		//// test region identifier here
		//GridExtremalRegionBuilder<ComparerASC>* terminal_region_builder = new GridExtremalRegionBuilder<ComparerASC>(g_rgt_func, g_grid);
		//terminal_region_builder->BeginIntegration();
		//terminal_region_builder->GetOutputLabels()->OutputToIntFile("terminal_regions.raw");

		// DO DESCENDING MANIFOLDS
		StartTask();
		unordered_map<INDEX_TYPE, INT_TYPE> extmapDSC;
		for (auto p : test->mMaxGraph->mCellIndexToListIdMap) {
			extmapDSC[g_maxv_labeling->Cell2HighestVertex(p.first)] = test->mMaxGraph->Representative(p.second);
		}
		GridSimplifiedExtremalRegionBuilder<ComparerDSC, GridFuncType, GridType>* test_simp_reg_builder_dsc =
			new GridSimplifiedExtremalRegionBuilder<ComparerDSC, GridFuncType, GridType>(g_rgt_func, g_grid, g_topo_grid);
		test_simp_reg_builder_dsc->BeginIntegration(extmapDSC);
		EndTaskAndRecord("Simplified Min/Max Termination Labeling - DSC");

		test_simp_reg_builder_dsc->GetOutputLabels()->OutputToIntFile("simp_terminal_regions_dsc.raw");

		//

		StartTask();
		//return 1;

		IntegratorTypeASC* newintegrator_asc =
			new IntegratorTypeASC(g_rgt_func, g_grid, error_threshold,
				gradient_threshold, iteration_limit);
		newintegrator_asc->BeginIntegration(test_simp_reg_builder_asc->GetIdMap(),
			test_simp_reg_builder_asc->GetOutputLabels(), true);

		newintegrator_asc->GetOutputLabels()->OutputToIntFile("newintegrator_asc.raw");
		EndTaskAndRecord("NumericalTracing NewIntegration - ASC");

		//return 1;

		StartTask();

		IntegratorTypeDSC* newintegrator_dsc =
			new IntegratorTypeDSC(g_rgt_func, g_grid, error_threshold,
				gradient_threshold, iteration_limit);
		newintegrator_dsc->BeginIntegration(test_simp_reg_builder_dsc->GetIdMap(),
			test_simp_reg_builder_dsc->GetOutputLabels(), true);

		newintegrator_dsc->GetOutputLabels()->OutputToIntFile("newintegrator_dsc.raw");
		EndTaskAndRecord("NumericalTracing NewIntegration - DSC");

		//return 1;

		// REMOVE DISCONNECTED COMPONENTS
		StartTask();
		IsolatedCCRegionRemoverNEW<ComparerASC, GridFuncType >* cleaner1_asc =
			new IsolatedCCRegionRemoverNEW<ComparerASC, GridFuncType>(g_rgt_func, newintegrator_asc->GetOutputLabels());
		printf("Removing Disconnected Regions\n");
		cleaner1_asc->ComputeOutput();
		printf("here2\n");
		EndTaskAndRecord("Removed Disconnected Components - ASC");
		newintegrator_asc->GetOutputLabels()->OutputToIntFile("newintegrator_cleaned_asc.raw");

		// REMOVE DISCONNECTED COMPONENTS
		StartTask();
		IsolatedCCRegionRemoverNEW<ComparerDSC, GridFuncType>* cleaner1_dsc =
			new IsolatedCCRegionRemoverNEW<ComparerDSC, GridFuncType>(g_rgt_func, newintegrator_dsc->GetOutputLabels());
		printf("Removing Disconnected Regions\n");
		cleaner1_dsc->ComputeOutput();
		printf("here2\n");
		EndTaskAndRecord("Removed Disconnected Components - DSC");
		newintegrator_dsc->GetOutputLabels()->OutputToIntFile("newintegrator_cleaned_dsc.raw");

		//return 1;

		// NOW MAKE EDGEMAP
		StartTask();
		DenseLabeling<int>* otherlabels = newintegrator_asc->GetOutputLabels();

		DenseLabeling<INDEX_TYPE>* outlabels;

		g_edge_map = new VertexLabelingToBoundaryLabeling<int, MaxVLType>( g_topo_grid);

		g_edge_map->ComputeMINBoundary(otherlabels);
		g_edge_map->ComputeMAXBoundary(newintegrator_dsc->GetOutputLabels(), g_maxv_labeling);
		
		//g_topo_grid->set_restriction(g_edge_map->GetOutputLabels());
		EndTaskAndRecord("Topological EdgeMap");

	}

	

	// REBUILD GRADIENT WITH BOUNDARY MAP
	RobinsType* constrained_robins = new RobinsType(g_topo_grid, g_maxv_labeling, g_edge_map->GetOutputLabels(), base_grad);
	if (outputdebug) {
		g_edge_map->GetOutputLabels()->OutputToFile("boundary_labels.raw");
	}
	StartTask();
	printf("redoing discrete gradient in changed boundarids\n");
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
		INDEX_TYPE threadfixcont = 0;
		// iterate over all vertices
		GridType::DCellsIterator verts(g_topo_grid, 0, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
		for (verts.begin(); verts.valid(); verts.advance()) {
			INDEX_TYPE vert_GI = verts.value();
			bool hasdiff = false;
			if (g_edge_map->GetOutputLabels()->GetLabel(vert_GI) > 0) hasdiff = true;
			if (!hasdiff) {
				GridType::CofacetsIterator edgeit(g_topo_grid);
				//bool hasdiff = false;

				for (edgeit.begin(vert_GI); edgeit.valid(); edgeit.advance()) {
					INDEX_TYPE edge_GI = edgeit.value();
					if (g_edge_map->GetOutputLabels()->GetLabel(edge_GI) > 0 &&
						g_maxv_labeling->Cell2HighestVertex(edge_GI) == vert_GI) {
						hasdiff = true;
						break;
					}
				}
			}
			if (hasdiff) {
				constrained_robins->ComputeLowerStar(vert_GI);
				threadfixcont++;
			}
		}
#pragma omp critical
		{
			printf("thread %d did %llu fixed vertices\n", thread_num, threadfixcont);
		}
	}
	EndTaskAndRecord("Topological Robins round 2");

	//g_topo_alg->setAscendingManifoldDimensions();

	//base_grad->outputToFile(gradname);
	//return 1;


	stageGrad[1] = new GradType(g_topo_grid);
	stageGrad[1]->copyValues(base_grad);
	stageLabel[1] = new DenseLabeling<char>(g_topo_grid->numCells());;
	stageLabel[1]->CopyValues(g_edge_map->GetOutputLabels());


	g_topo_alg = new TopologicalGradientUsingAlgorithms<GridType, FuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
	printf("after Second robins:\n");
	g_topo_alg->count_critical_points(4);

	//g_topo_alg->setAscendingManifoldDimensions();
	//base_grad->outputToFile(gradname);
	//return 1;

	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	// DO LINES ACCURATE GRADIENT COMPUTATION 
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------


	StartTask();
	int labeltarget = 4;
	g_digitizing_streamline_integrator_asc = new StreamlineIntegratorTypeASC(g_grid, g_rgt_func, g_topo_grid, error_threshold, gradient_threshold, iteration_limit);
	g_digitizing_streamline_integrator_asc->SetDigitizingTarget(g_edge_map->GetOutputLabels(), labeltarget);

	ReIntegrateUpFrom2Saddles();
	EndTaskAndRecord("Numerical SaddleInt");

	g_digitizing_streamline_integrator_asc->get_output()->OutputToFile("digitized_lines.raw");

	//StartTask();
	//CombineLabels();
	//EndTaskAndRecord("Topological CombineLabels");

	g_edge_map->GetOutputLabels()->OutputToFile("boundary_labels2.raw");

	//g_discrete_grad = new GradType(g_topo_grid);
	//g_discrete_grad->ClearAllGradient();
	//g_robin_alg = new RobinsType(g_topo_grid, g_maxv_labeling, g_edge_map->GetOutputLabels(), g_discrete_grad);
	//g_robin_alg->ComputePairing();


	StartTask();
	printf("robins...\n");
	set<INDEX_TYPE> reverts;
	GridType::AllCellsIterator allit(g_topo_grid);
	for (allit.begin(); allit.valid(); allit.advance()) {
		INDEX_TYPE cid = allit.value();
		if (g_edge_map->GetOutputLabels()->GetLabel(cid) == labeltarget) {
			reverts.insert(g_maxv_labeling->Cell2HighestVertex(cid));
		}
	}
	for (auto id : reverts) {
		constrained_robins->ComputeLowerStar(id);
	}

	EndTaskAndRecord("Topological Robins3");

	g_digitizing_streamline_integrator_asc->get_output()->OutputToFile("digitized_lines.raw");

	g_topo_alg = new TopologicalGradientUsingAlgorithms<GridType, FuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
	printf("after Second robins:\n");
	g_topo_alg->count_critical_points(4);


	//g_topo_alg->setAscendingManifoldDimensions();


	////// record stage
	//stageGrad[2] = new GradType(g_topo_grid);
	//stageGrad[2]->copyValues(base_grad);
	//stageLabel[2] = new DenseLabeling<char>(g_topo_grid->numCells());;
	//stageLabel[2]->CopyValues(g_edge_map->GetOutputLabels());

	//std::set<INDEX_TYPE> new_cp;
	//GridType::AllCellsIterator allit2(g_topo_grid);
	//int countboundary = 0;
	//for (allit2.begin(); allit2.valid(); allit2.advance()) {
	//	INDEX_TYPE cid = allit2.value();
	//	if (stageGrad[1]->getCritical(cid) && original_cp.count(cid) == 0) new_cp.insert(cid);
	//	if (stageGrad[1]->getCritical(cid) && g_topo_grid->boundaryValue(cid) != 0) countboundary++;
	//}


	//printf("have %d new cps \n", new_cp.size());
	//printf("there were %d boundary cps\n", countboundary);
	//// pick which ones to write
	//std::set<INDEX_TYPE> crit_v;
	//for (auto id : new_cp) {
	//	crit_v.insert(g_maxv_labeling->Cell2HighestVertex(id));
	//}


	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	// DO LINES ACCURATE GRADIENT COMPUTATION 
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------


	StartTask();
	labeltarget = 5;
	g_digitizing_streamline_integrator_dsc = new StreamlineIntegratorTypeDSC(g_grid, g_rgt_func, g_topo_grid, error_threshold, gradient_threshold, iteration_limit);
	g_digitizing_streamline_integrator_dsc->SetDigitizingTarget(g_edge_map->GetOutputLabels(), labeltarget);

	ReIntegrateDownFrom1Saddles();
	EndTaskAndRecord("Numerical SaddleInt");

	g_digitizing_streamline_integrator_dsc->get_output()->OutputToFile("digitized_lines_dsc.raw");

	//StartTask();
	//CombineLabels();
	//EndTaskAndRecord("Topological CombineLabels");

	g_edge_map->GetOutputLabels()->OutputToFile("boundary_labels3.raw");

	//g_discrete_grad = new GradType(g_topo_grid);
	//g_discrete_grad->ClearAllGradient();
	//g_robin_alg = new RobinsType(g_topo_grid, g_maxv_labeling, g_edge_map->GetOutputLabels(), g_discrete_grad);
	//g_robin_alg->ComputePairing();


	StartTask();
	printf("robins...\n");
	//set<INDEX_TYPE> reverts;
	reverts.clear();
	//GridType::AllCellsIterator allit(g_topo_grid);
	for (allit.begin(); allit.valid(); allit.advance()) {
		INDEX_TYPE cid = allit.value();
		if (g_edge_map->GetOutputLabels()->GetLabel(cid) == labeltarget) {
			reverts.insert(g_maxv_labeling->Cell2HighestVertex(cid));
		}
	}
	for (auto id : reverts) {
		constrained_robins->ComputeLowerStar(id);
	}

	EndTaskAndRecord("Topological Robins4");

	g_digitizing_streamline_integrator_dsc->get_output()->OutputToFile("digitized_lines4.raw");

	g_topo_alg = new TopologicalGradientUsingAlgorithms<GridType, FuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
	printf("after third robins:\n");
	g_topo_alg->count_critical_points(4);


	//g_topo_alg->setAscendingManifoldDimensions();


	//// record stage
	stageGrad[2] = new GradType(g_topo_grid);
	stageGrad[2]->copyValues(base_grad);
	stageLabel[2] = new DenseLabeling<char>(g_topo_grid->numCells());;
	stageLabel[2]->CopyValues(g_edge_map->GetOutputLabels());

	std::set<INDEX_TYPE> new_cp;
	GridType::AllCellsIterator allit2(g_topo_grid);
	int countboundary = 0;
	for (allit2.begin(); allit2.valid(); allit2.advance()) {
		INDEX_TYPE cid = allit2.value();
		if (stageGrad[1]->getCritical(cid) && original_cp.count(cid) == 0) new_cp.insert(cid);
		if (stageGrad[1]->getCritical(cid) && g_topo_grid->boundaryValue(cid) != 0) countboundary++;
	}


	printf("have %d new cps \n", new_cp.size());
	printf("there were %d boundary cps\n", countboundary);
	// pick which ones to write
	std::set<INDEX_TYPE> crit_v;
	for (auto id : new_cp) {
		crit_v.insert(g_maxv_labeling->Cell2HighestVertex(id));
	}



	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	// DO LOCAL SIMPLIFICATION
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------

	printf("before local cancellations:\n");
	g_topo_alg->count_critical_points(4);
	printf("checking for loops:\n");
	g_topo_alg->CheckGradientForLoops();

	int counts[3];

	for (int i = 0; i < 3; i++) counts[i] = 0;
	for (int k = 0; k < 2; k++) {
		GridType::AllCellsIterator eit(g_topo_grid);
		for (eit.begin(); eit.valid(); eit.advance()) {
			INDEX_TYPE eid = eit.value();
			INDEX_TYPE v1 = g_maxv_labeling->Cell2HighestVertex(eid);
			if (base_grad->getCritical(eid)) {
				char lab1 = g_edge_map->GetOutputLabels()->GetLabel(eid);
				std::vector<INDEX_TYPE> candidates;
				GridType::CofacetsIterator cfit(g_topo_grid);
				for (cfit.begin(eid); cfit.valid(); cfit.advance()) {
					INDEX_TYPE fid = cfit.value();



					if (base_grad->getCritical(fid)) {
						INDEX_TYPE v2 = g_maxv_labeling->Cell2HighestVertex(fid);

						//if (v1 == v2)
						candidates.push_back(fid);


					}

				}
				//if (candidates.size() == 1) {
				for (auto fid : candidates) {
					//INDEX_TYPE fid = candidates[0];
					char lab2 = g_edge_map->GetOutputLabels()->GetLabel(fid);
					//
					if (g_topo_grid->dimension(eid) == 1) {
						if (g_topo_grid->boundaryValue(eid) == g_topo_grid->boundaryValue(fid)
							//lab1 != lab2 &&
							/*g_maxv_labeling->Before(g_topo_grid->VertexNumberFjkromCellID(g_robin_alg->lowest_vertex(fid)),
							g_topo_grid->VertexNumberFromCellID(g_robin_alg->lowest_vertex(eid)))*/) {
							base_grad->setPair(eid, fid);
							base_grad->setPair(fid, eid);
							break;
							//printf("%d: %llu -> %llu : %d - %d\n",k, eid, fid, lab1, lab2);
						}
					}
					else if (lab1 != lab2) {
						counts[g_topo_grid->dimension(eid)]++;
						continue;
					}
					else {
						if (g_topo_grid->boundaryValue(eid) == g_topo_grid->boundaryValue(fid)
							//lab1 != lab2 &&
							/*g_maxv_labeling->Before(g_topo_grid->VertexNumberFjkromCellID(g_robin_alg->lowest_vertex(fid)),
							g_topo_grid->VertexNumberFromCellID(g_robin_alg->lowest_vertex(eid)))*/) {
							base_grad->setPair(eid, fid);
							base_grad->setPair(fid, eid);
							break;
							//printf("%d: %llu -> %llu : %d - %d\n",k, eid, fid, lab1, lab2);
						}
					}
				}
				//}


			}



		}

	}
	for (int i = 0; i < 3; i++) printf("dim-%d crits from mismatch\n", counts[i]);

	printf("after local cancellations:\n");
	g_topo_alg->count_critical_points(4);

	printf("chekcing for loops:\n");
	g_topo_alg->CheckGradientForLoops();


	stageGrad[3] = new GradType(g_topo_grid);
	stageGrad[3]->copyValues(base_grad);
	stageLabel[3] = new DenseLabeling<char>(g_topo_grid->numCells());;
	stageLabel[3]->CopyValues(g_edge_map->GetOutputLabels());


	// do the writing here!
	std::vector<VertexNeighborhood*> vns;
	for (auto id : crit_v) {
		if (g_topo_grid->boundaryValue(id) != 0) continue;
		Vec3l coords; g_topo_grid->cellid2Coords(id, coords);

		VertexNeighborhood* vn = new VertexNeighborhood();
		vns.push_back(vn);
		for (int i = 0; i < 4; i++) {
			vn->states.push_back(NeighborState());
		}
		// function stuff;
		int lid = 0;// local index in 5 cube
		INDEX_TYPE fid = g_topo_grid->VertexNumberFromCellID(id);
		Vec3l fc = g_grid->XYZ3d(fid);
		for (int z = -1; z <= 1; z++) {
			for (int y = -1; y <= 1; y++) {
				for (int x = -1; x <= 1; x++) {
					Vec3l nfc = fc + Vec3l(x, y, z);
					vn->vec_v[lid] = g_rgt_func->SampleGrad(nfc);
					vn->vert_v[lid] = g_rgt_func->SampleImage(nfc);
					lid++;
				}
			}
		}

		lid = 0;
		for (int z = -2; z <= 2; z++) {
			for (int y = -2; y <= 2; y++) {
				for (int x = -2; x <= 2; x++) {
					Vec3l offset(x, y, z);
					Vec3l ncoords = coords + offset;
					INDEX_TYPE nid = g_topo_grid->coords2Cellid(ncoords);
					if (g_maxv_labeling->Cell2HighestVertex(nid) == id) {
						vn->lstar[lid] = 1;


					}
					for (int i = 0; i < 4; i++) {
						vn->states[i].grad[lid] = stageGrad[i]->getBitfield(nid);
						vn->states[i].labs[lid] = stageLabel[i]->GetLabel(nid);
					}


					lid++;
				}
			}
		}





	}

	FILE* fvn = fopen("VertexNeighborhoods.bin", "wb");
	int numverts = vns.size();
	printf("sizes: char=%d, Vec3d=%d, float=%d, GradBitField=%d\n", sizeof(char), sizeof(GInt::Vec3d), sizeof(float), sizeof(GradType::GradBitfield));

	printf("writing %d vertex neighborhoods\n", numverts);
	fwrite(&numverts, sizeof(int), 1, fvn);
	int numstages = 4;
	fwrite(&numstages, sizeof(int), 1, fvn);
	for (auto vn : vns) {
		fwrite(vn->lstar, sizeof(char), 5 * 5 * 5, fvn);
		fwrite(vn->vec_v, sizeof(Vec3d), 3 * 3 * 3, fvn);
		fwrite(vn->vert_v, sizeof(float), 3 * 3 * 3, fvn);
		for (int i = 0; i < numstages; i++) {
			fwrite(vn->states[i].grad, sizeof(GradType::GradBitfield), 5 * 5 * 5, fvn);
			fwrite(vn->states[i].labs, sizeof(char), 5 * 5 * 5, fvn);
		}
	}
	fclose(fvn);
	//


	g_topo_alg->setAscendingManifoldDimensions();
	base_grad->outputToFile(gradname);



	//g_topo_alg->setAscendingManifoldDimensions();



	//char gradname[1024];
	//sprintf(gradname, "%s.grad", argv[4]);
	//g_discrete_grad->outputToFile(gradname);

	//g_topo_alg->setAscendingManifoldDimensions();
	//char gradname[1024];
	//sprintf(gradname, "%s.grad", argv[4]);
	//g_discrete_grad->outputToFile(gradname);
	////fprintf(gtiming, "IO Write %d %d %d\n",
	////	std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
	////	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
	////	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
	//task_start_time = std::chrono::steady_clock::now();

	//fprintf(gtiming, "Total Overall %d %d %d\n",
	//	std::chrono::duration_cast<std::chrono::milliseconds>(g_start_time - g_start_time).count(),
	//	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
	//	std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count());


	fclose(gtiming);


	return 1;
};


