//#include "integrate3.hpp"
//#include "vector2.hpp"
//#include <stdio.h>
#include <vector>
#include <set>
#include <queue>
#include <time.h>
#include "gi_timing.h"
#include "gi_topological_explicit_mesh_function.h"
#include "gi_topological_region_growing_simple_gradient_builder.h"
#include "gi_topological_convergent_gradient_builder.h"
#include "gi_robin_labeling.h"
#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_isolated_region_remover.h"
#include "gi_topological_utility_functions.h"
#include "gi_numeric_integrator_expanding_region_stop.h" // not where comparer should be
#include "gi_morse_smale_complex_basic.h"
#include "gi_morse_smale_complex_restricted.h"
#include "gi_kdtree.h"
#include "gi_msc_selectors.h"
#include "gi_numeric_streamline_integrator.h"
//#include "gi_dataflow_objects.h"
//#include "gi_dfo_morse_smale_selectors.h"
#include "gi_fast_robins_noalloc.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_graphs.h"
#include "gi_ms_complex_to_graph.h"


using namespace GInt;
typedef RegularGrid2D GridType;
typedef RegularGridBilinearFunction GridFuncType;

//typedef UncachedRegularGridTrilinearFunction GridFuncType;
typedef TopologicalRegularGrid2D MeshType;
typedef IndexCompareLessThan<GridFuncType> ComparerASC;
typedef IndexCompareGreaterThan<GridFuncType> ComparerDSC;
//typedef MorseSmaleComplexBasic<FLOATTYPE, MeshType, MeshFuncType, GradType> MSCType;
//typedef NumericIntegratorExpandingRegionStopWithCutoff<AdaptiveEulerAdvector3D<-1>, ComparerASC> IntegratorTypeWC;
//typedef NumericIntegratorExpandingRegionStop<AdaptiveEulerAdvector3D<-1>, ComparerASC> IntegratorTypeASC;
typedef DigitizingNumericStreamlineIntegrator2dASC<MeshType, GridFuncType, AdaptiveEulerAdvector2D<GridFuncType, 1> > StreamlineIntegratorTypeASC;
typedef DigitizingNumericStreamlineIntegrator2dDSC<MeshType, GridFuncType, AdaptiveEulerAdvector2D<GridFuncType, -1> > StreamlineIntegratorTypeDSC;
typedef DiscreteGradientLabeling<MeshType> GradType;
//typedef UncachedMaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
//typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
typedef RegularGridMaxMinVertexLabeling2D<MeshType, GridFuncType> MaxVLType;
typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType, 4, 6> RobinsType;
typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> MeshFuncType;
//typedef SlidingWindowRobinsNoalloc < RegularGrid2D, RegularGridTrilinearFunction, MeshType, MaxVLType, GradType> NewRobinsType;

typedef MorseSmaleComplexBasic<FLOATTYPE, MeshType, MeshFuncType, GradType> MSCType;

std::chrono::steady_clock::time_point g_start_time;
std::chrono::steady_clock::time_point task_start_time;
std::chrono::steady_clock::time_point prior_time;
std::chrono::steady_clock::time_point now_time;
FILE* gtiming;

void StartTask(const char* str = " ") {
	task_start_time = std::chrono::steady_clock::now();
	printf("starting task %s\n", str);
}
std::vector<int> timings_sequence;

void EndTaskAndRecord(const char* s) {
	now_time = std::chrono::steady_clock::now();
	// format: [global activity name] [task] [start] [end] [dration]
	timings_sequence.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count());

	printf("TIMING: %s %d %d %d\n", s,
		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
}
StreamlineIntegratorTypeASC* g_digitizing_streamline_integrator_asc;
StreamlineIntegratorTypeDSC* g_digitizing_streamline_integrator_dsc;
GridType* g_grid;
GridFuncType* g_rgt_func;
MeshType* g_topo_grid;

#ifdef USE_REGION_CLEANER
VertexLabelingToBoundaryLabeling<INDEX_TYPE>* g_edge_map;
#else
VertexLabelingToBoundaryLabeling<int, MaxVLType>* g_edge_map;
#endif

GradType* base_grad;
//RobinsLabelingAlgorithm<MeshType, MeshFuncType> *g_robin_alg;
TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>* g_topo_alg;

MaxVLType* g_maxv_labeling;
MeshFuncType* g_topo_func;
int iteration_limit = 1000;
int per_x, per_y;
float error_threshold = 0.1;
float gradient_threshold = 0.01;

void ReIntegrateUpFromSaddles(GradType* base_grad) {


#ifdef OUTPUT_DEBUG_2LINES
	FILE* fout = fopen("test_uplines.bin", "wb");
	int countwrite = 0;
	int writelength = 0;
#endif

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
		MeshType::DCellsIterator face_iterator(g_topo_grid, 1, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
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
			int counter = 10;
			while (!cell_queue.empty() && counter >= 0) {
				INDEX_TYPE current = cell_queue.front();
				cell_queue.pop();

				cell_visited.insert(current);
				//result.push_back(current);

				MeshType::CofacetsIterator cofacets(g_topo_grid);
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
				if (g_topo_grid->dimension(arc_hex_id) != 2) continue;
				std::vector<Vec2d> points;
				std::vector<INDEX_TYPE> dline;
				Vec2d seed; g_topo_grid->centroid(arc_hex_id, seed);
				seed = seed * 0.5; // back to grid coordinates
				g_digitizing_streamline_integrator_asc->IntegrateStreamline(seed, points, dline);
#ifdef OUTPUT_DEBUG_2LINES
#pragma omp critical
				{
					int num = points.size();
					fwrite(&num, sizeof(int), 1, fout);
					fwrite(points.data(), sizeof(Vec2d), points.size(), fout);
					countwrite++; writelength += points.size();
				}
#endif

			}
			//printf("\n");
			g_digitizing_streamline_integrator_asc->set_label(sad_id);



		}

	}

#ifdef OUTPUT_DEBUG_2LINES
	printf("countwrite = %d, length = %d\n", countwrite, writelength);
	fclose(fout);
#endif

}

void ReIntegrateDownFromSaddles(GradType* base_grad) {
#ifdef OUTPUT_DEBUG_2LINES
	FILE* fout = fopen("test_dnlines.bin", "wb");
#endif
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
		MeshType::DCellsIterator face_iterator(g_topo_grid, 1, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
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
			int counter = 10;
			while (!cell_queue.empty() && counter >= 0) {
				INDEX_TYPE current = cell_queue.front();
				cell_queue.pop();

				cell_visited.insert(current);
				//result.push_back(current);

				MeshType::FacetsIterator facets(g_topo_grid);
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
				std::vector<Vec2d> points;
				std::vector<INDEX_TYPE> dline;
				Vec2d seed; g_topo_grid->centroid(arc_vert_id, seed);
				seed = seed * 0.5; // back to grid coordinates
				g_digitizing_streamline_integrator_dsc->IntegrateStreamline(seed, points, dline);
#ifdef OUTPUT_DEBUG_2LINES
#pragma omp critical
				{
					int num = points.size();
					fwrite(&num, sizeof(int), 1, fout);
					fwrite(points.data(), sizeof(Vec2d), points.size(), fout);
				}
#endif
			}
			//printf("\n");
			g_digitizing_streamline_integrator_dsc->set_label(sad_id);



		}
#ifdef OUTPUT_DEBUG_2LINES



		for (auto p : criticals) {
			INDEX_TYPE id = p.second;
			g_digitizing_streamline_integrator_dsc->get_output()->SetLabel(id, 4);
		}

		fclose(fout);
#endif
	}


}



int use_valleys = 0;
GInt::Geometric2DGraph* ComputeFromFloatSlice(float* buffer, int x, int y, float persistence, std::vector<float>& cut_points) {
	g_start_time = std::chrono::steady_clock::now();


	// create structures for input and output
	printf("making grid...\n");
	GridType* grid = new GridType({ x,y }, { 0,0 });
	printf("making grid function...\n");
	GridFuncType* grid_func = new GridFuncType(grid, buffer);
	printf("making mesh...\n");
	MeshType* mesh = new MeshType(grid);
	printf("making discrete grad structure...\n");
	GradType* grad = new GradType(mesh);
	printf("making max vertex labeling...\n");
	MaxVLType* max_vl = new MaxVLType(mesh, grid_func);
	max_vl->ComputeOutput();
	printf("making mesh function...\n");
	MeshFuncType* mesh_func = new MeshFuncType();
	mesh_func->setMeshAndFuncAndMaxVLabeling(mesh, grid_func, max_vl);

	// build discrete gradient
	printf("computing discrete gradient...\n");
	//if (! grad->load_from_file())
	RobinsType* compute_grad = new RobinsType(mesh, max_vl, grad);
	compute_grad->ComputePairing();


	g_grid = grid;
	g_rgt_func = grid_func;
	g_topo_grid = mesh;
	g_maxv_labeling = max_vl;
	g_topo_func = mesh_func;
	base_grad = grad;
	printf("computing gradient\n");
	StartTask();
	grid_func->ComputeGradFromImage(1);
	//g_rgt_func->Negate();
	EndTaskAndRecord("NumericalTracing GradCompute");

	// CREATE RESTRICTION MAP - WILL NEED
	DenseLabeling<char>* restriction_labels = new DenseLabeling<char>(mesh->numCells());
	restriction_labels->SetAll(0);
	// we will always need a constrained robins alg
	RobinsType* constrained_robins = new RobinsType(mesh, max_vl, restriction_labels, grad);

	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	// DO LINES ACCURATE GRADIENT COMPUTATION 
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------
	//-------------------------------------------------------------

	if (true) {
		StartTask("asc 1");
		int labeltarget = 1;
		g_digitizing_streamline_integrator_asc = new StreamlineIntegratorTypeASC(g_grid, g_rgt_func, g_topo_grid, error_threshold, gradient_threshold, iteration_limit);
		g_digitizing_streamline_integrator_asc->SetDigitizingTarget(restriction_labels, labeltarget);
		ReIntegrateUpFromSaddles(base_grad);
		EndTaskAndRecord("Numerical SaddleIntASC");
	}

	if (true) {
		StartTask("dsc 1");
		int labeltarget = 2;
		g_digitizing_streamline_integrator_dsc = new StreamlineIntegratorTypeDSC(g_grid, g_rgt_func, g_topo_grid, error_threshold, gradient_threshold, iteration_limit);
		g_digitizing_streamline_integrator_dsc->SetDigitizingTarget(restriction_labels, labeltarget);
		ReIntegrateDownFromSaddles(base_grad);
		EndTaskAndRecord("Numerical SaddleIntDSC");
	}
#ifdef OUTPUT_DEBUG_2LINES
	restriction_labels->OutputToFile("labs.bin");
#endif

	//-------------------------------------------------------------
	// NOW REDO DISCRETE GRADIENT IN NEIGHBORHOOD
	//-------------------------------------------------------------
	StartTask();
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
		MeshType::DCellsIterator verts(g_topo_grid, 0, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
		for (verts.begin(); verts.valid(); verts.advance()) {
			INDEX_TYPE vert_GI = verts.value();
			bool hasdiff = false;
			if (restriction_labels->GetLabel(vert_GI) > 0) hasdiff = true;
			if (!hasdiff) {
				MeshType::AdjacentCellsIterator cocells(g_topo_grid);
				//bool hasdiff = false;

				for (cocells.begin(vert_GI); cocells.valid(); cocells.advance()) {
					INDEX_TYPE cocell_GI = cocells.value();
					if (restriction_labels->GetLabel(cocell_GI) > 0 &&
						g_maxv_labeling->Cell2HighestVertex(cocell_GI) == vert_GI) {
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
	EndTaskAndRecord("Topological ConformingGrad");













	g_topo_alg = 
		new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(mesh_func, mesh, grad);
	g_topo_alg->setAscendingManifoldDimensions();
	// build ms complex
	printf("computing ms complex...\n");
	MSCType* msc = new MSCType(grad, mesh, mesh_func);
	msc->SetBuildArcGeometry({ 1,1,1 });
	msc->ComputeFromGrad();
	printf("computing msc hierarchy...\n");
	msc->ComputeHierarchy(persistence);
	msc->SetSelectPersAbs(persistence);
	// build ridge graph

	printf("computing mesh cells graph...\n");
	printf(" -- cut points size: %d\n", cut_points.size());
	MeshCellsGraph* graph;
	if (! use_valleys) {
		if (cut_points.size() == 0) {
			graph = GInt::BuildMeshCellsGraphFromMSCRidges<MSCType, MeshType>(msc, mesh);
		}
		else {
			printf(" -- using cut poitns\n");
			graph = GInt::BuildMeshCellsGraphFromMSCRidgesWithCuts<MSCType, MeshType, MeshFuncType>(msc, mesh, mesh_func, cut_points);
		}
	}
	else {
		if (cut_points.size() == 0) {
			graph = GInt::BuildMeshCellsGraphFromMSCValleys<MSCType, MeshType>(msc, mesh);
		}
		else {
			printf(" -- using cut poitns\n");
			graph = GInt::BuildMeshCellsGraphFromMSCValleysWithCuts<MSCType, MeshType, MeshFuncType>(msc, mesh, mesh_func, cut_points);
		}
	}

	
	// make a geometric graph
	printf("computing geometric graph...\n");
    auto geom_graph = GInt::BuildGeometricGraphFromMeshGraph<MeshType>(graph, mesh, 10);
    printf("finished compute graph.");
	return geom_graph;



}

void one_step_neighbors(INT_TYPE arcid, MSCType* msc, set<INT_TYPE>& arcs) {
	MSCType::SurroundingArcsIterator sarcs(msc);

	auto lower = msc->getArc(arcid).lower;
	for (sarcs.begin(lower); sarcs.valid(); sarcs.advance()) {
		INT_TYPE sarcid = sarcs.value();
		arcs.insert(sarcid);
	}

	auto upper = msc->getArc(arcid).upper;
	for (sarcs.begin(upper); sarcs.valid(); sarcs.advance()) {
		INT_TYPE sarcid = sarcs.value();
		arcs.insert(sarcid);
	}

}

bool arc_less_than(INT_TYPE arc1, INT_TYPE arc2, MSCType* msc) {

	auto& a1 = msc->getArc(arc1);
	auto& a2 = msc->getArc(arc2);

	if (a1.persistence < a2.persistence) return true;
	if (a1.persistence > a2.persistence) return false;

	auto& n1l = msc->getNode(a1.lower);
	auto& n2l = msc->getNode(a2.lower);

	if (n1l.value < n2l.value) return true;
	if (n1l.value > n2l.value) return false;

	if (a1.lower < a2.lower) return true;
	if (a1.lower > a2.lower) return false;
	if (a1.upper < a2.upper) return true;
	if (a1.upper > a2.upper) return false;
	return false;
}

bool IsIndependentArc(INT_TYPE arcid, MSCType* msc) {
	
	set<INT_TYPE> surrounding_arcs;
	one_step_neighbors(arcid, msc, surrounding_arcs);
	//printf("negsied = %d\n", surrounding_arcs.size());

	for (auto id : surrounding_arcs) {
		if (arcid != id && !arc_less_than(arcid, id, msc)) return false;
	}
	return true;
}


void DoAccuarte() {

}

GInt::Geometric2DGraph* GetIndependentArcs(float* buffer, int x, int y, float persistence) {

	// create structures for input and output
	printf("making grid...\n");
	GridType* grid = new GridType({ x,y }, { 0,0 });
	printf("making grid function...\n");
	GridFuncType* grid_func = new GridFuncType(grid, buffer);
	printf("making mesh...\n");
	MeshType* mesh = new MeshType(grid);
	printf("making discrete grad structure...\n");
	GradType* grad = new GradType(mesh);
	printf("making max vertex labeling...\n");
	MaxVLType* max_vl = new MaxVLType(mesh, grid_func);
	max_vl->ComputeOutput();
	printf("making mesh function...\n");
	MeshFuncType* mesh_func = new MeshFuncType();
	mesh_func->setMeshAndFuncAndMaxVLabeling(mesh, grid_func, max_vl);

	// build discrete gradient
	printf("computing discrete gradient...\n");
	//if (! grad->load_from_file())
	RobinsType* compute_grad = new RobinsType(mesh, max_vl, grad);
	compute_grad->ComputePairing();
	TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>* g_topo_alg =
		new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(mesh_func, mesh, grad);
	g_topo_alg->setAscendingManifoldDimensions();
	// build ms complex
	printf("computing ms complex...\n");
	MSCType* msc = new MSCType(grad, mesh, mesh_func);
	msc->ComputeFromGrad();

	set<INT_TYPE> independent_arcs;
	int modcount[8];
	for (int i = 0; i < 8; i++) modcount[i] = 0;
	MSCType::ArcsIterator ait(msc);
	for (ait.begin(); ait.valid(); ait.advance()) {
		INT_TYPE aid = ait.value();

		if (IsIndependentArc(aid, msc)) {
			independent_arcs.insert(aid);
			modcount[(aid * 8) / msc->numArcs()]++;
		}
	}
	printf(" ##-- "); for (int i = 0; i < 8; i++) printf("%d ", modcount[i]); printf("\n");

	printf("##-- num msc base arcs %d\n", msc->numArcs());

	printf("computing msc hierarchy...\n");
	msc->ComputeHierarchy(persistence);

	int num_matching_indep = 0;
	int num_not_matching = 0;
	for (auto cr : msc->GetCancellationRecords()) {
		//if (cr.persistence > 0) break;

		auto& a = msc->getArc(cr.arcp);
		if (a.created > 0) continue;  /// this is not a base arc
		
		if (independent_arcs.count(cr.arcp) != 0) {
			num_matching_indep++;
		}
		else {
			num_not_matching++;
		}
	}
	
	printf("##-- Num Indep = %d\n", independent_arcs.size());
	printf("##-- Num matching = %d\n", num_matching_indep);
	printf("##-- Num not matching = %d\n", num_not_matching);
	return NULL;
	
	msc->SetSelectPersAbs(persistence);
	// build ridge graph

	printf("computing mesh cells graph...\n");
	auto graph = GInt::BuildMeshCellsGraphFromMSCRidges<MSCType, MeshType>(msc, mesh);

	// make a geometric graph
	printf("computing geometric graph...\n");
	auto geom_graph = GInt::BuildGeometricGraphFromMeshGraph<MeshType>(graph, mesh, 10);
	return geom_graph;



}


int main(int argc, char** argv) {

	ThreadedTimer timer(1);
	timer.StartGlobal();

	int X, Y;
	float pers;
	std::string filename;

    if (argc < 5) { printf("Usage: filename X Y persistence use_valleys\n"); return 0; }
	filename = std::string(argv[1]);
	sscanf(argv[2], "%d", &X);
	sscanf(argv[3], "%d", &Y);
	sscanf(argv[4], "%f", &pers);
	
	if (argc >= 6) {
		sscanf(argv[5], "%d", &use_valleys);
	}

	vector<float> cut_values;
	if (argc >= 7) {
		for (int i = 6; i < argc; i++) {
			float val;
			sscanf(argv[i], "%f", &val);
			cut_values.push_back(val);
			printf("cut value = %f\n", val);
			
		}
	}

	FILE* fin = fopen(filename.c_str(), "rb");
	float* buff = new float[X*Y];
    printf("read %d bytes\n", fread(buff, sizeof(float), X*Y, fin));
	fclose(fin);

	//GetIndependentArcs(buff, X, Y, pers);
	//return 1;

	auto geometric_graph = ComputeFromFloatSlice(buff, X, Y, pers, cut_values);
	printf("sanity checking1 ... \n");
	geometric_graph->CheckConnections();
	printf("done\n");
	// the output graph for ML
	GraphForMLLabeling outgraph;
	outgraph.CreateFrom2DPolylineGraph(geometric_graph);
	outgraph.Write(filename);
	
	if (true) {
		GraphForMLLabeling testgraph;
		testgraph.Read(filename);
		auto* ng = testgraph.Create2DPolylineGraph();
		printf("sanity checking2 ... \n");
		ng->CheckConnections();
		printf("done\n");
		printf("tested\n");
	}

	printf("Done!!\n");

}

////#include "integrate3.hpp"
////#include "vector2.hpp"
////#include <stdio.h>
//#include <vector>
//#include <set>
//#include <queue>
//#include <time.h>
//#include "gi_strictly_numeric_integrator.h"
//#include "gi_numeric_integrator_region_stop.h"
//#include "gi_numeric_integrator_expanding_region_stop.h"
//#include "gi_timing.h"
//#include "gi_labeling_to_bounary_labeling.h"
//#include "gi_topological_explicit_mesh_function.h"
//#include "gi_topological_region_growing_simple_gradient_builder.h"
//#include "gi_topological_convergent_gradient_builder.h"
//#include "gi_robin_labeling.h"
//#include "gi_adaptive_in_quad_euler_advector.h"
//#include "gi_numeric_integrator_2d_restricted_expanding_region.h"
//#include "gi_topological_2d_restricted_expanding_regions.h"
//#include "gi_topological_gradient_using_algorithms.h"
//#include "gi_topological_regular_grid_restricted.h"
//#include "gi_isolated_region_remover.h"
//#include "gi_topological_utility_functions.h"
//#include "gi_morse_smale_complex_basic.h"
//#include "gi_morse_smale_complex_restricted.h"
//#include "gi_kdtree.h"
//#include "gi_msc_selectors.h"
//#include "gi_numeric_streamline_integrator.h"
//
//#include "gi_experimental.h"
//
//#include "concurrentqueue.h"
//#include <thread>
//#include <mutex>
//#include <map>
//#include <chrono>
//#include <unordered_map>
//
////
////
////int *certains;
////int *dests;
////
////
////
////
////
////RegularGridTrilinearFunction* GlobalGrid;
////
////INT_TYPE highest_neighbor(Vec3i xyz) {
////	INT_TYPE tid = GlobalGrid->Index3d(xyz);
////	INT_TYPE thighest = tid;
////	Vec3i neighbors[6];
////	int nn = GlobalGrid->GatherExistingNeighbors6(xyz, neighbors);
////	for (int i = 0; i < nn; i++) {
////		INT_TYPE tneg = GlobalGrid->Index3d(neighbors[i]);
////		if (GlobalGrid->is_greater(tneg, thighest)) thighest = tneg;
////	}
////	return thighest;
////}
////
////// does path tracing up with path compression
////int find(int s, int* a) {
////	if (a[s] == s) {
////		return s;
////	}
////	else if (a[a[s]] == s) {
////		return s;
////	}
////	a[s] = find(a[s], a);
////	return a[s];
////}
//
////using namespace GInt;
////typedef IndexCompareLessThan Comparer;
////typedef TopologicalRegularGridRestricted GridType;
////typedef MorseSmaleComplexBasic<FLOATTYPE, TopologicalRegularGrid, TopologicalExplicitDenseMeshFunction, DiscreteGradientLabeling> MSCType;
////typedef NumericIntegratorExpandingRegionStop<AdaptiveEulerAdvector<-1>, Comparer> IntegratorType;
////typedef NumericIntegratorRegionStop<AdaptiveEulerAdvector<-1>, Comparer> IntegratorType2;
////typedef NumericStreamlineIntegrator<AdaptiveEulerAdvector<-1> > StreamlineIntegratorType;
////typedef IsolatedRegionRemover<Comparer> RegionRemoverType;
////typedef NumericStreamlineIntegrator<AdaptiveEulerAdvector<-1> > StreamlineIntegratorType;
//
//typedef std::chrono::time_point < std::chrono::high_resolution_clock> TIMEP;
//#define  maptype std::map
//
//struct barf {
//	barf(){}
//	barf(int a, int b, bool c) : a(a), b(b), c(c) {}
//	int a; 
//	int b;
//	bool c; 
//};
//
//
//int me[1000];
//
//class testm {
//public:
//	int id;
//	std::mutex m;
//	testm(int id) :id(id), counter(0){}
//
//
//	void work() {
//		
//		for (int i = 0; i < 1000; i++) {
//			me[id]++;
//		}
//	}
//
//	inline void increment() {
//		counter++;
//	}
//	inline void increment_TS() {
//		work();
//			{
//				lock_guard<std::mutex> guard(m);
//				counter++;
//			}
//
//	}
//	int counter;
//};
//
//std::mutex m;
//int uniqueid = 0;
//void addstuff(moodycamel::ConcurrentQueue<barf>* q, int i) {
//	int myid;
//		{
//			std::lock_guard<std::mutex> guard(m);
//			myid = uniqueid++;
//		}
//	q->enqueue(barf(i, myid, false));
//
//}
//
//int gnumiter = 1000000;
//int gnumbins = 10;
//void teststuff(maptype<int, testm*>* map, int tid) {
//	for (int i = 0; i < gnumiter; i++) {
//		for (int j = 0; j < gnumbins; j++) {
//			map->operator[](j)->increment_TS();
//		}
//	}
//
//}
//
//using namespace GInt;
//typedef IndexCompareLessThan Comparer;
//typedef TopologicalRegularGrid GridType;
//typedef MorseSmaleComplexBasic<FLOATTYPE, TopologicalRegularGrid, TopologicalExplicitDenseMeshFunction, DiscreteGradientLabeling> MSCType;
//typedef NumericIntegratorExpandingRegionStop<AdaptiveEulerAdvector<-1>, Comparer> IntegratorType;
//typedef NumericIntegratorRegionStop<AdaptiveEulerAdvector<-1>, Comparer> IntegratorType2;
//typedef NumericStreamlineIntegrator<AdaptiveEulerAdvector<-1> > StreamlineIntegratorType;
//typedef IsolatedRegionRemover<Comparer> RegionRemoverType;
//typedef NumericStreamlineIntegrator<AdaptiveEulerAdvector<-1> > StreamlineIntegratorType;
//
//
//int main(int argc, char** argv) {
//
//	int X, Y, Z;
//	int per_x, per_y, per_z;
//	std::string filename;
//
//	if (argc < 8) { printf("Usage: X Y Z filename per_x per_y per_z\n"); return 0; }
//	sscanf(argv[1], "%d", &X);
//	sscanf(argv[2], "%d", &Y);
//	sscanf(argv[3], "%d", &Z);
//	filename = std::string(argv[4]);
//	sscanf(argv[5], "%d", &per_x);
//	sscanf(argv[6], "%d", &per_y);
//	sscanf(argv[7], "%d", &per_z);
//	
//	
//	RegularGrid* m_grid;
//	//RegularGridTrilinearFunction* m_func;
//	RegularGridTrilinearFunction* m_func2;
//	GridType *m_tgrid;
//
//	m_grid = new RegularGrid(Vec3i(X, Y, Z), Vec3b(per_x, per_y, per_z));
//
//	m_func2 = new RegularGridTrilinearFunction(m_grid);
//	m_func2->LoadImageFromFile(filename.c_str());
//	m_tgrid = new GridType(m_grid);
//
//	TopologicalExplicitDenseMeshFunction* m_topofunc;
//
//	m_topofunc = new TopologicalExplicitDenseMeshFunction();
//	m_topofunc->setMeshAndAllocate(m_tgrid);
//	m_topofunc->copyVertexValuesFromGridFunction(m_func2);
//	m_topofunc->setCellValuesMaxOfVerts();
//
//	DEPGRAPH_TYPE* tgraph = new DEPGRAPH_TYPE(m_tgrid, m_topofunc);
//
//	TestingDependencyGraphIndependentTraversal* traversal = new TestingDependencyGraphIndependentTraversal(tgraph);
//
//	auto starttimer = chrono::steady_clock::now();
//	traversal->doStuff();
//
//	auto endtimer = chrono::steady_clock::now();
//
//	printf("took %d ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(endtimer - starttimer).count());
//
//
//	
//	printf("Done\n");
//	return 1;
//	//int numthreads; 
//	//sscanf(argv[1], "%d", &numthreads);
//	//sscanf(argv[2], "%d", &gnumbins);
//	//sscanf(argv[3], "%d", &gnumiter);
//
//	//maptype<int, testm*> mymap;
//
//	//for (int i = 0; i < gnumbins; i++) {
//	//	mymap[i] = new testm(i);
//	//}
//	//TIMEP m_global_start = std::chrono::high_resolution_clock::now();
//
//	////int numthreads =  std::thread::hardware_concurrency();
//
//
//
//	//std::thread* t = new std::thread[numthreads];
//
//	//for (int i = 0; i < numthreads; i++) {
//	//	t[i] = std::thread(teststuff, &mymap, i);
//	//}
//
//	//for (int i = 0; i < numthreads; i++) {
//	//	t[i].join();
//	//}
//
//	//TIMEP m_global_end = std::chrono::high_resolution_clock::now();
//
//
//	////for (int i = 0; i < 10; i++) {
//	////	printf("%d %d\n", i, mymap[i]->counter);
//	////}
//	//printf("%d %d %d %d\n", std::chrono::duration_cast<std::chrono::milliseconds>(m_global_end - m_global_start).count(), numthreads, gnumiter, gnumbins);
//
//	//moodycamel::ConcurrentQueue<barf> testqueue;
//
//	//{
//	//	std::lock_guard<std::mutex> guard(m);
//	//	testqueue.enqueue(barf(99, uniqueid++, false));
//	//}
//
//	//int numthreads = std::thread::hardware_concurrency();
//	//std::thread* t = new std::thread[numthreads];
//
//	//for (int i = 0; i < numthreads; i++) {
//	//	t[i] = std::thread(addstuff, &testqueue, i);
//	//}
//
//	//for (int i = 0; i < numthreads; i++) {
//	//	t[i].join();
//	//}
//
//	//barf placeholder;
//	//while (testqueue.try_dequeue(placeholder)) {
//	//	printf("%d %d\n", placeholder.a, placeholder.b);
//	//}
//
//	return 0;
//
//} 
//
//
//
////
////
////	ThreadedTimer timer(1);
////	timer.StartGlobal();
////
////    int X, Y, Z;
////    int iteration_limit;
////    int per_x, per_y, per_z;
////    float error_threshold, gradient_threshold;
////    std::string filename;
////
////	if (argc < 11) { printf("Usage: X Y Z filename grad_threshold error_threshold maxnumiter per_x per_y per_z\n"); return 0; }
////    sscanf(argv[1], "%d", &X);
////	sscanf(argv[2], "%d", &Y);
////	sscanf(argv[3], "%d", &Z);
////    filename = std::string(argv[4]);
////	sscanf(argv[5], "%f", &error_threshold);
////	sscanf(argv[6], "%f", &gradient_threshold);
////	sscanf(argv[7], "%d", &iteration_limit);
////	sscanf(argv[8], "%d", &per_x);
////	sscanf(argv[9], "%d", &per_y);
////	sscanf(argv[10], "%d", &per_z);
////
////#if 1
////	RegularGrid* m_grid;
////	//RegularGridTrilinearFunction* m_func;
////	RegularGridTrilinearFunction* m_func2;
////	GridType *m_tgrid;
////	IntegratorType* i2;
////	RegionRemoverType* i2clean;
////	VertexLabelingToBoundaryLabeling<INDEX_TYPE>* edgemap;
////	TopologicalExplicitDenseMeshFunction* m_topofunc;
////	DiscreteGradientLabeling *labeling;
////	RobinsLabelingAlgorithm *robin;
////	TopologicalGradientUsingAlgorithms* topo_algs;
////	MSCType* msc;
////	MSCSelectorLivingNodes<MSCType>* s_nodes;
////	MSCSelectorNodeIndex<MSCType>* s_1saddles;
////	MSCSelectorRepresentative1Saddle<MSCType>* s_r1saddles;
////
////	NumericIntegratorNew* sinteegrator;
////	StreamlineIntegratorType* sintegrator;
////	TerminateNearExtrema* extremumtermination;
////
////
////	
////	// added by Harsh
////    m_grid = new RegularGrid(Vec3i(X, Y, Z), Vec3b(per_x, per_y, per_z));
////
////    m_func2 = new RegularGridTrilinearFunction(m_grid);
////    m_func2->LoadImageFromFile(filename.c_str());
////    m_func2->ComputeGradFromImage(1);
////
////
////#ifdef BADER_CONSTRAINT
////    m_func2->Negate();
////
////    m_tgrid = new GridType(m_grid);
////
////	//RunMeshConsistencyChecks(m_tgrid);
////
////
////
////
////
////
////    printf("\n --------------------------------------- \n");
////
////	i2 = new IntegratorType(m_func2, m_grid, error_threshold, gradient_threshold, iteration_limit);
////   
////	i2->BeginIntegration();
////	i2->GetOutputLabels()->OutputToIntFile("classes.raw");
////    printf("done with numeric integration!\ncreating vertex to boundary labeling...\n");
////
////	i2clean = new RegionRemoverType(m_func2, i2->GetOutputLabels());
////	i2clean->ComputeOutput();
////	i2clean->GetOutputLabels()->OutputToIntFile("classes_clean.raw");
////
////	edgemap = new VertexLabelingToBoundaryLabeling<INDEX_TYPE>(i2clean->GetOutputLabels(), m_tgrid);
////    edgemap->ComputeBoundary();
////    edgemap->GetOutputLabels()->OutputToFile("boundary_labels.raw");
////    printf("done with boundary map creation!\ncreating topology function and copying values...\n");
////    //edgemap->OutputEdgesToFile("surfin.raw");
////
////    printf("\n --------------------------------------- \n");
////    m_tgrid->set_restriction (edgemap->GetOutputLabels ());
////#else
////    TopologicalRegularGrid *m_tgrid = new TopologicalRegularGrid(m_grid);
////#endif
////
////    printf(" Creating TopologicalExplicitDenseMeshFunction ...");
////    fflush(stdout);
////
////    ThreadedTimer timer2(1);
////    timer2.StartGlobal();
////
////   m_topofunc = new TopologicalExplicitDenseMeshFunction();
////    m_topofunc->setMeshAndAllocate(m_tgrid);
////    m_topofunc->copyVertexValuesFromGridFunction(m_func2);
////    m_topofunc->setCellValuesMaxOfVerts();
////
////    timer2.EndGlobal ();
////    printf(" Done! ");
////    timer2.PrintAll ();
////
////    labeling = new DiscreteGradientLabeling(m_tgrid);
////    labeling->ClearAllGradient();
////
////    robin = new RobinsLabelingAlgorithm(m_topofunc, m_tgrid, labeling);
////    robin->compute_output();
////
////
////
////
////
////
////    //TIMEP start = std::chrono::high_resolution_clock::now();
////
////    //DiscreteGradientLabeling* topo_grad2 = new DiscreteGradientLabeling(m_tgrid);
////    //TopologicalRegionGrowingSimpleGradientBuilder* topo_cg_builder = new TopologicalRegionGrowingSimpleGradientBuilder(m_topofunc, m_tgrid, topo_grad2);
////    //topo_cg_builder->computeGradient();
////
////    //TIMEP end = std::chrono::high_resolution_clock::now();
////
////    //printf("serial alg took %d ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
////
////	//FILE* TEST_GRAD = fopen("TESTGRAD.raw", "wb");
////
////	//fclose(TEST_GRAD);
////	//printf("about to set dim asc man\n");
////	//TopologicalGradientUsingAlgorithms* topo_algs = new TopologicalGradientUsingAlgorithms(m_topofunc, m_tgrid, topo_grad2);
////	//topo_algs->setAscendingManifoldDimensions();
////
////	printf("about to set dim asc man\n");
////	topo_algs = new TopologicalGradientUsingAlgorithms(m_topofunc, m_tgrid, labeling);
////	topo_algs->setAscendingManifoldDimensions();
////	topo_algs->CheckGradientConsistency();
////
////    char gradname2[2048];
////    sprintf(gradname2, "%s.grad", filename.c_str ());
////    labeling->outputToFile(gradname2);
////    printf("output to %s\n", gradname2);
////
////
////	msc =
////		new MSCType(labeling, m_tgrid, m_topofunc);
////
////	msc->ComputeFromGrad();
////
////	// ATOM LOCALIZATION BEGINS HERE
////
////	// add in selector shit here
////	kdtree* mkt = kd_create(3);
////
////	// add all my minima to the kdtree
////	for (int i = 0; i < msc->numNodes(); i++) {
////		GInt::node<FLOATTYPE>& n = msc->getNode(i);
////		if (n.dim != 0) continue;
////		Vec3l ncoords;
////		m_tgrid->cellid2Coords(n.cellindex, ncoords);
////		Vec3d ndcoords = ncoords;
////		ndcoords *= 0.5;
////
////		// add in periodic versions (assums periodic = 1,1,1)
////		int* data = new int[1];
////		data[0] = i;
////		for (int px = -1; px <= 1; px++)
////			for (int py = -1; py <= 1; py++)
////				for (int pz = -1; pz <= 1; pz++) {
////			double dd[3];
////			dd[0] = ndcoords[0] + X * px;
////			dd[1] = ndcoords[1] + Y * py;
////			dd[2] = ndcoords[2] + Z * pz;
////			kd_insert(mkt, dd, data);
////				}
////
////	}
////	// now i added all the "atoms" 
////	// for each node find closest "atom"
////	// because i dont have atoms file, simply use a perturbed version of each node
////	
////	vector<pair<Vec3d, INT_TYPE> > restrictnodes; // i save the atom positions as well as the index so that i can render later the atom position.
////	for (int i = 0; i < msc->numNodes(); i++) {
////		GInt::node<FLOATTYPE>& n = msc->getNode(i);
////		if (n.dim != 0) continue;
////		Vec3l ncoords;
////		m_tgrid->cellid2Coords(n.cellindex, ncoords);
////		Vec3d ndcoords = ncoords;
////		ndcoords *= 0.5;
////		Vec3d rando((rand() % 100) / 100.0, (rand() % 100) / 100.0, (rand() % 100) / 100.0);
////		ndcoords += rando;
////		double p[3]; p[0] = ndcoords[0]; p[1] = ndcoords[1]; p[2] = ndcoords[2];
////		// now find closest to ndcoords
////		kdres* res = kd_nearest(mkt, p);
////		if (res == NULL) {
////			printf("ERROR: null point found in nearest\n");
////		}
////		else {
////			int* cc = (int*)kd_res_item(res, 0);
////			restrictnodes.push_back(pair<Vec3d, INT_TYPE>(ndcoords, cc[0]));
////		}
////	}
////	// now this code below actually forces the MSC to avoid cancelling the right atoms
////	for (int i = 0; i < restrictnodes.size(); i++) {
////		INT_TYPE nodeid = restrictnodes[i].second;
////		msc->getNode(nodeid).boundary = 128; // make boundary type
////	}
////
////
////	// ATOM LOCALIZATION ENDS HERE
////
////	msc->ComputeHierarchy(1.0);
////	s_nodes = new MSCSelectorLivingNodes<MSCType>(msc);
////	s_1saddles = new MSCSelectorNodeIndex<MSCType>(msc, 1);
////	s_1saddles->add_parent(s_nodes);
////	s_r1saddles = new MSCSelectorRepresentative1Saddle<MSCType>(msc);
////	s_r1saddles->add_parent(s_1saddles);
////	
////	// here is how to actually extract arcs
////	sintegrator = new StreamlineIntegratorType(m_grid, m_func2, 0.01, 0, 5000);
////	extremumtermination = new TerminateNearExtrema(i2->GetExtrema(), m_grid);
////	sintegrator->SetAdvectionChecker(extremumtermination);//sinteegrator->
////
////	s_r1saddles->compute_output();
////	for (auto it = s_r1saddles->output.begin();
////		it != s_r1saddles->output.end(); it++) {
////		GInt::node<FLOATTYPE>& n = msc->getNode(*it);
////
////		// iterate over arcs around minimum
////		MSCType::SurroundingArcsIterator sit(msc);
////		for (sit.begin(*it); sit.valid(); sit.advance()) {
////			//printf("asdf %d\n", sit.value());
////			INT_TYPE aid = sit.value();
////
////			GInt::arc<FLOATTYPE>& a = msc->getArc(aid);
////			// only look at saddle-minimum arcs
////			if (a.upper == *it) {
////				vector<INDEX_TYPE> g;
////				msc->fillArcGeometry(aid, g); // combinatorial arc geometry
////
////
////				vector<Vec3d> numgeom; // structure to hold numerical path
////				INDEX_TYPE lowerid = msc->getNode(a.lower).cellindex; // combinatorial says this is where wer should end up
////				INDEX_TYPE stopid = m_tgrid->VertexNumberFromCellID(lowerid); // get vertex id in grid coordinates from combinatorial min from topological coordinates
////
////				// now try integrating down starting from cominatorial geomery
////				for (int j = 0; j < g.size(); j++) {
////
////					Vec3l c;
////					m_tgrid->cellid2Coords(g[j], c);
////					vector<Vec3d> gg;
////					Vec3d cf = c;
////					cf *= 0.5;
////					numgeom.push_back(cf); // record combinatorial step as part of geometric path
////
////					ADVECTION_EVENT event = sintegrator->IntegrateStreamline(cf, gg); // do numerical streamline integration
////					if (event == ADVECTION_EVENT::LOW_GRADIENT) {
////						//glColor3f(0.5, .5, .9);
////					}
////					else if (event == ADVECTION_EVENT::HIT_PREASSIGNED) {
////						//glColor3f(0.5, .9, .2);
////						// THIS IS THE ONLY ONE THAT I'VE SEEN HAPPEN IN PRACTICE FOR BES DATA
////					}
////					else if (event == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
////						//glColor3f(0.9, 0.5, 0.2);
////					}
////					else {
////						//glColor3f(.2, .2, .2);
////					}
////
////
////					Vec3d cf2 = gg.back(); // last point on numerically integrated line
////					//printf("endpoint found %d -> (%f, %f, %f)\n", extremumtermination->WhichExtremum(gg.back()), cf2[0], cf2[1], cf2[2]);
////					// test if this is located near a critical vertex that was used as a stopping criteria in region integration
////					if (extremumtermination->WhichExtremum(gg.back()) == stopid) {
////						numgeom.insert(numgeom.end(), gg.begin(), gg.end()); // this numerical path is good! so paste it onto rest of geometry and we're done with this arc
////						break;
////					}
////
////				}
////
////				// This code here that is commented is what renders the geometry - notice restarting line strip when distance is too big- likey periodic boundary crossing
////				//glColor3f(1, .5, .1);
////				//glBegin(GL_LINE_STRIP);
////				//for (int j = 0; j < numgeom.size(); j++) {
////				//	if (j > 0 && (numgeom[j] - numgeom[j - 1]).Mag() > 4.0) {
////				//		glEnd();
////				//		glBegin(GL_LINE_STRIP);
////				//	}
////				//	glVertex3f(numgeom[j][0], numgeom[j][1], numgeom[j][2]);
////				//}
////				//glEnd();
////
////
////
////			}
////		}
////
////
////	}
////
////
////#else
////	//StrictlyNumericIntegrator* i = new StrictlyNumericIntegrator(Vec3i(X, Y, Z), Vec3b(true, true, true), argv[4], error_threshold, gradient_threshold);
////    //i->BeginIntegration();
////    RegularGrid* m_grid = new RegularGrid(Vec3i(X, Y, Z), Vec3b(per_x, per_y, per_z));
////	RegularGridTrilinearFunction* m_func = new RegularGridTrilinearFunction(m_grid);
////    m_func->LoadImageFromFile(argv[4]);
////    m_func->ComputeGradFromImage(1);
////
////	NumericIntegratorExpandingRegionStop* i2 = new NumericIntegratorExpandingRegionStop(m_func, m_grid, error_threshold, gradient_threshold, iteration_limit);
////	//NumericIntegratorRegionStop* i2 = new NumericIntegratorRegionStop(Vec3i(X, Y, Z), Vec3b(true, true, true), argv[4], error_threshold, gradient_threshold);
////	i2->BeginIntegration();
////
////	//i2->GetOutputLabels()->OutputToIntFile("vertex_labels.raw");
////	printf("done with numeric integration!\ncreating vertex to boundary labeling...\n");
////	TopologicalRegularGrid* topo_grid = new TopologicalRegularGrid(i2->GetGrid());
////	
////	VertexLabelingToBoundaryLabeling<INDEX_TYPE>* edgemap = new VertexLabelingToBoundaryLabeling<INDEX_TYPE>(i2->GetOutputLabels(), topo_grid);
////	edgemap->ComputeBoundary();
////	//edgemap->GetOutputLabels()->OutputToFile("boundary_labels.raw");
////	printf("done with boundary map creation!\ncreating topology function and copying values...\n");
////	edgemap->OutputEdgesToFile("surfin.raw");
////
////	TopologicalExplicitDenseMeshFunction* topo_func = new TopologicalExplicitDenseMeshFunction();
////	topo_func->setMeshAndAllocate(topo_grid);
////	topo_func->copyVertexValuesFromGridFunction(i2->GetFunction());
////	topo_func->setCellValuesMaxOfVerts();	
////
////	Topological2DRestrictedExpandingRegions* topo_2d_expand =
////		new Topological2DRestrictedExpandingRegions(topo_func, topo_grid, edgemap->GetOutputLabels());
////	topo_2d_expand->BeginIntegration();
////	//topo_2d_expand->TEST_OUTPUT->OutputToFile("edgemapexpansiontest.raw");
////	int countassigned = 0;
////	int countunassigned = 0;
////	
////	FILE* fcertains = fopen("certainout.raw", "wb");
////	for (INDEX_TYPE i = 0; i < topo_grid->numCells(); i++) {
////
////		if (topo_grid->dimension(i) == 1 && edgemap->GetOutputLabels()->GetLabel(i) == 1) {
////			if (topo_2d_expand->TEST_OUTPUT->GetLabel(i) == 0) {
////				countunassigned++;
////				//printf("ERROR: edge that is part of boundar is not certain!\n");
////			}
////			else {
////				countassigned++;
////			}
////		}
////
////		if (topo_grid->dimension(i) == 1 && topo_2d_expand->TEST_OUTPUT->GetLabel(i) > 0){
////			fwrite(&i, sizeof(INDEX_TYPE), 1, fcertains);
////			int tmpclass = topo_2d_expand->TEST_OUTPUT->GetLabel(i);
////			fwrite(&tmpclass, sizeof(int), 1, fcertains);
////		}
////	}
////	fclose(fcertains);
////
////	printf("%d assigned, %d unassigned\n", countassigned, countunassigned);
////
////	FILE* fdelayed = fopen("delayedout.raw", "wb");
////	for (int i = 0; i < topo_2d_expand->delayed_cells.size(); i++) {
////		if (topo_2d_expand->TEST_OUTPUT->GetLabel(topo_2d_expand->delayed_cells[i]) == 0) {
////			fwrite(&topo_2d_expand->delayed_cells[i], sizeof(INDEX_TYPE), 1, fdelayed);
////		} 
////	}
////	fclose(fdelayed);
////	fdelayed = fopen("crititcal1cells.raw", "wb");
////	for (int i = 0; i < topo_2d_expand->critical_primal_edges.size(); i++) {
////		fwrite(&topo_2d_expand->critical_primal_edges[i], sizeof(INDEX_TYPE), 1, fdelayed);
////	}
////	fclose(fdelayed);
////	//fdelayed = fopen("STARTPOINTS.raw", "wb");
////	//for (auto it = topo_2d_expand->bottleneckcount.begin(); it != topo_2d_expand->bottleneckcount.end(); it++) {
////	//	if ((*it).second > 1) 
////	//		fwrite(&((*it).first), sizeof(INDEX_TYPE), 1, fdelayed);
////	//}
////	//fclose(fdelayed);
////	//fdelayed = fopen("STARTPOINTS1.raw", "wb");
////	//for (auto it = topo_2d_expand->bottleneckcount.begin(); it != topo_2d_expand->bottleneckcount.end(); it++) {
////	//	if ((*it).second == 1)
////	//		fwrite(&((*it).first), sizeof(INDEX_TYPE), 1, fdelayed);
////	//}
////	//fclose(fdelayed);	
////	//NumericIntegrator2DRestrictedExpandingRegion* edgeintegrator =
////	//	new NumericIntegrator2DRestrictedExpandingRegion(m_func, m_grid, topo_grid,
////	//	edgemap->GetOutputLabels(), error_threshold, gradient_threshold, iteration_limit);
////
////	//edgeintegrator->BeginIntegration();
////	//edgeintegrator->TEST_OUTPUT->OutputToFile("edgemapexpansiontest.raw");
////	//edgeintegrator->TEST_LINEOUT->OutputToFile("edgemapquadintegration.raw"); 
////	//// great sin here- atually negating function!
////	//m_func->Negate();
////	//char negname[2048];
////	//sprintf(negname, "%s.neg.raw", argv[4]);
////
////	//}
////	
////	
////	
////
//////	// teh following just tests partition
//////	std::vector<INDEX_TYPE> topo_index_partition;
//////	DenseLabeling<char>* TEST_PARTITION = new DenseLabeling<char>(topo_grid->numCells());
//////	TEST_PARTITION->SetAll(0);
//////
//////#pragma omp parallel
//////	{
//////#pragma omp single
//////				{
//////					int num_threads = omp_get_num_threads();
//////					ArrayIndexPartitioner::EvenChunkSplit(topo_grid->numCells(), num_threads, topo_index_partition);
//////				}
//////
//////		int thread_num = omp_get_thread_num();
//////		TopologicalRegularGrid::AllCellsIterator all_cell_iter(topo_grid, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
//////		for (all_cell_iter.begin(); all_cell_iter.valid(); all_cell_iter.advance()) {
//////			INDEX_TYPE cell_id = all_cell_iter.value();
//////			TEST_PARTITION->SetLabel(cell_id, 1);
//////		}
//////	}
//////	TopologicalRegularGrid::AllCellsIterator all_cell_iter(topo_grid);
//////	for (all_cell_iter.begin(); all_cell_iter.valid(); all_cell_iter.advance()) {
//////		INDEX_TYPE cell_id = all_cell_iter.value();
//////		if (TEST_PARTITION->GetLabel(cell_id) != 1) printf("ERROR!!!\n");
//////	}
//////	
////	timer.EndGlobal();
////	timer.PrintAll();
////	exit(1);
////
////
////
////	//printf("done with discrete function creation!\ncomputing steepest descent discrete gradient\n");
////	//DiscreteGradientLabeling* topo_grad = new DiscreteGradientLabeling(topo_grid);
////	//topo_grad->ClearAllGradient();
////	//TopologicalRegionGrowingSimpleGradientBuilder* topo_rg_builder =
////	//	new TopologicalRegionGrowingSimpleGradientBuilder(topo_func, topo_grid, topo_grad);
////	//topo_rg_builder->computeGradient();
////	//printf("done!\n");
////
////	//DenseLabeling<char>* testgrad = new DenseLabeling<char>(topo_grid->numCells());
////	//testgrad->SetAll(0);
////	//TopologicalRegularGrid::AllCellsIterator it(topo_grid);
////	//for (it.begin(); it.valid(); it.advance()) {
////	//	INDEX_TYPE cid = it.value();
////	//	if (topo_grad->getCritical(cid))
////	//		testgrad->SetLabel(cid, topo_grid->dimension(cid) + 1);
////	//}
////	//testgrad->OutputToFile("topograd.raw");
////
////
////
////
////	DiscreteGradientLabeling* topo_grad2 = new DiscreteGradientLabeling(topo_grid);
////	TopologicalConverventGradientBuilder* topo_cg_builder = new TopologicalConverventGradientBuilder(topo_func, topo_grid, topo_grad2);
////	topo_cg_builder->computeGradient();
////
////	//// OK FUNCTION SEEMS TO BE FINE!
////	//DenseLabeling<float>* testfunc = new DenseLabeling<float>(topo_grid->numCells());
////	//TopologicalRegularGrid::AllCellsIterator it(topo_grid);
////	//for (it.begin(); it.valid(); it.advance()) {
////	//	INDEX_TYPE cid = it.value();
////	//	testfunc->SetLabel(cid, topo_func->cellValue(cid));
////	//}
////	//testfunc->OutputToFile("topofunc.raw");
////
////	//testboundary->OutputToFile("boundaries.raw");
////	/// OK i think boundaries are correct!
////	//DenseLabeling<char>* testboundary = new DenseLabeling<char>(topo_grid->numCells());
////	//TopologicalRegularGrid::AllCellsIterator it(topo_grid);
////	//for (it.begin(); it.valid(); it.advance()) {
////	//	INDEX_TYPE cid = it.value();
////	//	testboundary->SetLabel(cid, topo_grid->boundaryValue(cid));
////	//}
////	//testboundary->OutputToFile("boundaries.raw");
////
////	char gradname2[2048];
////	sprintf(gradname2, "%s.grad", argv[4]);
////	topo_grad2->outputToFile(gradname2);
////
////	timer.EndGlobal();
////	timer.PrintAll();
//////
//////
//////	INT_TYPE PROBLEMSIZE = GlobalGrid->NumElements();
//////	GlobalGrid->load_image_from_file(argv[4]);
//////
//////
//////	float grad_threshold;
//////	float error_threshold;
//////	sscanf(argv[5], "%f", &grad_threshold);
//////	sscanf(argv[6], "%f", &error_threshold);
//////
//////	int numiter;
//////	sscanf(argv[7], "%d", &numiter);
//////
//////	int rkindex = 1;
//////	if (argc > 8) sscanf(argv[8], "%d", &rkindex);
//////
//////	printf("computing gradient using rk%d-indices...\n", rkindex);
//////	GlobalGrid->compute_grad_from_image(rkindex);
//////
//////	//int *sources = new int[PROBLEMSIZE];
//////	dests = new int[PROBLEMSIZE];
//////	certains = new int[PROBLEMSIZE];
//////	//#pragma acc data copyin(Gradient[0:PROBLEMSIZE]) copyout(sources[0:PROBLEMSIZE], dests[0:PROBLEMSIZE])
//////	//{
//////	// number of iterations
//////	//private(stepsize, xx, yy, SOURCE, itercount, tsource, res1, grad, gradL, res01, res02, grad01, error, orig_x, orig_y, ORIGIN, itercount2, dest_x, dest_y, DESTINATION)
//////	//#pragma acc parallel loop 
//////	std::vector<INT_TYPE> maxvec;
//////
//////	// assign steepest ascent and find maxima
//////	printf("finding maxima...\n");
//////#pragma omp parallel shared(dests, maxvec) 
//////	{
//////#pragma omp for schedule(dynamic)  nowait 
//////		for (int p = 0; p < PROBLEMSIZE; p++) {
//////			//initialize certains to no-assigment
//////			certains[p] = -1;
//////
//////			// now find steepest ascent partner
//////			Vec3i xxyyzz = GlobalGrid->XYZ3d(p);
//////
//////			// initially set all dests to be one higher
//////			INT_TYPE thighest = highest_neighbor(xxyyzz);
//////			dests[p] = thighest;
//////			if (thighest == p) {
//////#pragma omp critical
//////				{
//////					maxvec.push_back(p);
//////				}
//////			}
//////		}
//////	}
//////
//////
//////
//////	// expand steepest ascent unchangable volumes
////////	int nummaxs = maxvec.size();
////////	printf("found %d maxima!\nExpanding maxima certain regions\n", nummaxs);
////////
////////#pragma omp parallel shared(dests, maxvec) 
////////	{
////////#pragma omp for schedule(dynamic)  nowait 
////////		for (int m = 0; m < nummaxs; m++) {
////////			INT_TYPE maximum = maxvec[m];
////////			certains[maximum] = m;
////////			Expand_Lower_Neighborhood(maximum);
////////		}
////////	}
////////
////////	// START SECTION JUST FOR DEBUGGING
////////	printf("done expanding!\n");
////////	int mycounter = 0;
////////	for (int p = 0; p < PROBLEMSIZE; p++) {
////////		if (certains[p] >= 0) mycounter++;
////////	}
////////	printf("%f percent of cells identified as certain\n", ((float)mycounter / (float)PROBLEMSIZE) * 100.0);
////////	FILE* fout1 = fopen("dest_fist_expand.raw", "wb");
////////	//fwrite(sources, sizeof(int), PROBLEMSIZE, fout);
////////	fwrite(certains, sizeof(int), PROBLEMSIZE, fout1);
////////	fclose(fout1);
//////	// END SECTION JUST FOR DEBUGGIN
//////
////////#pragma omp parallel shared(dests, GlobalGrid) 
////////	{
////////#pragma omp for schedule(dynamic)  nowait 
//////		
//////	// STRUCTURE TO HOLD GRAD ASSIGNMENT
//////	unsigned char* cell_dgrad;
//////	// SET ALL TO UNASSIGNED = 0;
//////
//////	SET_ALL_CRITS_TO_ASSIGNED(maxvec, cell_dgrad); // NO CONCURRENCY CHECK NEEDED
//////
//////	// REPLACE THIS ORDERING WITH NEW ORDERING on ALL CELLS
//////	for (int p = 0; p < PROBLEMSIZE; p++) {
//////	// INSTEAD, FIRST EXTEND EXISTING 
//////
//////
//////		//for (int xx = 0; xx < X; xx +=skip) {
//////			//for (int yy = 0; yy < Y; yy+=skip) {
//////
//////
//////			// a large gradient step could really mess things up. // start small
//////			Vec3i startpoint_i = GlobalGrid->XYZ3d(p);
//////	
//////			// CHECK IF TH STARTPOINT HAS BEEN ASSIGNED;
//////			if (CHECK_ASSIGNED(p, cell_dgrad)) continue;
//////			
//////			//if (startpoint_i[0] == 0 && startpoint_i[1] == 0) { printf("doing "); startpoint_i.print_vi(); }
//////			FLOATTYPE tstartgradmag = GlobalGrid->sampleG(startpoint_i).Mag();
//////			FLOATTYPE stepsize = 0.5 / tstartgradmag; // don't start moving more than half a grid cell
//////			//int p = xx + yy * X;
//////
//////			Vec3d source_d = startpoint_i;
//////
//////			int itercount = 0;
//////			//Vec3d tsource_d = source_d;
//////
//////			bool hasearlytermination = false;
//////
//////			/// temps to accelerate later stuff, but for now make sure they are values that will not be reused
//////			Vec3d oldgrad(-1, -1, -1);
//////			FLOATTYPE oldgradL = -1;
//////			Vec3d oldsource(-1, -1, -1);
//////			Vec3i oldbase(-10, -10, -10);
//////			Vec3d surroundinggrad[8];
//////
//////			Vec3i oldnext1base(-1, -1, -1);
//////			Vec3d oldnext1cache[8];
//////			FLOATTYPE dist_now = GlobalGrid->DistToBoundary(source_d);
//////			while (itercount < numiter) {
//////
//////
//////
//////				itercount++;
//////
//////				// fill these in
//////				Vec3d grad;
//////				FLOATTYPE gradL;
//////				// do we need to resample grad?
//////				if (oldsource == source_d) {
//////					//no, we did not move point - probably changed stepsize
//////					grad = oldgrad;
//////					gradL = oldgradL;
//////				}
//////				else {
//////					// do early termination - only recheckcheck this when the sample point is different
//////					Vec3i closestcell;
//////					if (dist_now <= 2) {
//////						closestcell = GlobalGrid->Inbounds(source_d + 0.5); // do we need IntFloor? probably not
//////					}
//////					else {
//////						closestcell = source_d + 0.5;
//////					}
//////					INT_TYPE tid = GlobalGrid->Index3d(closestcell);
//////					if (certains[tid] >= 0) {
//////						dests[p] = dests[tid];
//////						hasearlytermination = true;
//////						break;
//////					}
//////
//////					// we moved the source point - do we have to re-sample the surrounding graident?
//////					oldsource = source_d;
//////					Vec3i newbase = source_d.IntFloor();
//////					if (!(newbase == oldbase)) {
//////
//////						// WE MOVE CELLS, RECORD WITH A GRADIENT ARROW
//////						if (CHECK_ASSIGNED(newbase, cell_dgrad)) continue;
//////						ADD_GRAD_ARROW(cell_dgrad, oldbase, newbase);
//////
//////
//////						oldbase = newbase;
//////						if (dist_now <= 2) dist_now = GlobalGrid->DistToBoundary(newbase); // can sometimes get bigger because we approximate with offset
//////						if (dist_now <= 2) {
//////							GlobalGrid->get_grad_surrounding(newbase, surroundinggrad);
//////						}
//////						else {
//////							GlobalGrid->get_grad_surrounding_no_boundary_check(newbase, surroundinggrad);
//////						}
//////					}
//////					grad = GlobalGrid->TriLinInterpGrad(source_d, surroundinggrad);
//////					gradL = grad.Mag();
//////
//////					if (gradL * stepsize > 0.5) stepsize = 0.5 / gradL; // never take too big steps!!
//////
//////					oldgrad = grad;
//////					oldgradL = gradL;
//////				}
//////
//////				// we are done integrating if we reach a really flat area
//////				if (gradL < grad_threshold) {
//////					//MAKE THINGS CRITICAL
//////					MAKE_CRITICAL(oldbase, cell_dgrad);
//////					break;
//////				}
//////
//////
//////				Vec3d next;
//////				Vec3d next0; 
//////				Vec3d next0grad;
//////				Vec3d next1;
//////				if (dist_now <= 2) {
//////					next = GlobalGrid->IStep(source_d, grad, stepsize);
//////
//////					// now try 1/2 stepsize to see if we are within error
//////					next0 = GlobalGrid->IStep(source_d, grad, stepsize*0.5);
//////
//////					// we need to take another 1/2 step - so we need another gradient computation
//////					//  but we first try to reuse the same surrounding samples, so we don't have to
//////					// start sampling again
//////					next0grad;
//////					if (next0.IntFloor() == oldbase) {
//////						next0grad = GlobalGrid->TriLinInterpGrad(next0, surroundinggrad);
//////					}
//////					else {
//////						next0grad = GlobalGrid->TriLinInterpGrad(next0);
//////					}
//////					next1 = GlobalGrid->IStep(next0, next0grad, stepsize*0.5);
//////				}
//////				else {
//////					next = GlobalGrid->IStep_no_boundary_check(source_d, grad, stepsize);
//////
//////					// now try 1/2 stepsize to see if we are within error
//////					next0 = GlobalGrid->IStep_no_boundary_check(source_d, grad, stepsize*0.5);
//////
//////					// we need to take another 1/2 step - so we need another gradient computation
//////					//  but we first try to reuse the same surrounding samples, so we don't have to
//////					// start sampling again
//////					FLOATTYPE nextstepsize = stepsize*0.5;
//////					next0grad;
//////					Vec3i next0int = next0.IntFloor();
//////					if (next0int == oldbase) {
//////						next0grad = GlobalGrid->TriLinInterpGrad(next0, surroundinggrad);
//////					}
//////					else {
//////
//////						if (! (oldnext1base == next0int)) {
//////							GlobalGrid->get_grad_surrounding_no_boundary_check(next0int, oldnext1cache);
//////						}
//////						next0grad = GlobalGrid->TriLinInterpGrad(next0, oldnext1cache);
//////					}
//////					FLOATTYPE next0gradmag = next0grad.Mag(); //avoid taking square root, so instead do comparison on squaring both sides
//////					if (nextstepsize  * next0gradmag  > 0.5) {
//////						// this will be a HUGE error anyway, so just goto next
//////						stepsize = stepsize * 0.5;
//////						continue;
//////					}
//////					next1 = GlobalGrid->IStep_no_boundary_check(next0, next0grad, nextstepsize);
//////				}
//////				
//////				FLOATTYPE error = (next1 - next).MagSq();
//////
//////				if (error < error_threshold) {
//////					
//////					FLOATTYPE actualstepingrid = gradL * stepsize;
//////					dist_now -= actualstepingrid;
//////					if (actualstepingrid >= 0.5) {
//////						// set stepsize so we dont move too fast
//////						stepsize = 0.5 / gradL;
//////					}
//////					else {
//////						stepsize = stepsize * 2 * 0.9;
//////					}
//////					source_d = next1;
//////				}
//////				else {
//////					stepsize = stepsize * 0.5;
//////					continue;
//////				}
//////			}
//////			if (!hasearlytermination) {
//////				/// AGAIN MAKE CRITICAL
//////				MAKE_CRITICAL(asf);
//////				
//////				Vec3i origin = GlobalGrid->Inbounds(source_d);
//////				int ORIGIN = GlobalGrid->Index3d(origin);
//////				dests[p] = ORIGIN;
//////			}
//////		}
//////
//////		//FILE* fout = fopen("out.raw", "wb");
//////
//////
//////		//for (int xx = 0; xx < X; xx +=skip) {
//////		//	for (int yy = 0; yy < Y; yy+=skip) {
//////		//		int p = xx + yy * X;
//////		//		int i = Sources[p].x;
//////		//		int j = Sources[p].y;
//////		//		if (i < 0) i = 0;
//////		//		if (i > X-1) i = X-1;
//////		//		if (j < 0) j = 0;
//////		//		if (j > Y -1) j = Y-1;
//////		//		int d = i + j * X;
//////		//		fwrite(&p, sizeof(int), 1, fout);
//////		//		fwrite(&d, sizeof(int), 1, fout);
//////		//	}
//////		//}
//////
//////
//////		//for (int xx = 0; xx < X; xx +=skip) {
//////		//	for (int yy = 0; yy < Y; yy+=skip) {
//////		//		int p = xx + yy * X;
//////		//		float fx = (float) xx;
//////		//		float fy = (float) yy;
//////
//////		//		fwrite(&fx, sizeof(float), 1, fout);
//////		//		fwrite(&fy, sizeof(float), 1, fout);
//////		//		float dx = (float) Sources[p].x;
//////		//		float dy = (float) Sources[p].y;
//////		//		fwrite(&dx, sizeof(float), 1, fout);
//////		//		fwrite(&dy, sizeof(float), 1, fout);
//////		//		vect2 pp; pp.x = fx; pp.y = fy;
//////		//		vect2 gg; m_grad(pp, gg);
//////		//		//printf("(%f, %f), (%f, %f) (%f, %f)\n", fx, fy, Sources[p].x, Sources[p].y, gg.x, gg.y); //, Sources[p].x, Sources[p].y);
//////		//	}
//////		//}
//////		//
//////
//////	//}
//////
//////	FILE* fout2 = fopen("dest_pre_clean.raw", "wb");
//////	//fwrite(sources, sizeof(int), PROBLEMSIZE, fout);
//////	fwrite(dests, sizeof(int), PROBLEMSIZE, fout2);
//////	fclose(fout2);
//////
//////	printf("doing cleanup\n");
//////	for (int i = 0; i < PROBLEMSIZE; i++) {
//////		int d = find(i, dests);
//////		dests[i] = d;
//////	}
//////
//////	FILE* fout = fopen("dest.raw", "wb");
//////	//fwrite(sources, sizeof(int), PROBLEMSIZE, fout);
//////	fwrite(dests, sizeof(int), PROBLEMSIZE, fout);
//////	fclose(fout);
//////
//////	time(&now);
//////	printf("Read file in %.f seconds\n", difftime(now, starttimer));
//////	//printf("success: %f, %f\n", p.x, p.y);
//////	return 1;
////#endif
////}
