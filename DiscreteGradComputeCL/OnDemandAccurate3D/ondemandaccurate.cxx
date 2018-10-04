#include "base/gi_ondemand_accurate_grad_builder.h"

int main(int argc, char** argv) {
	GInt::OndemandDiscreteGradientBuilder* builder = new GInt::OndemandDiscreteGradientBuilder();
	builder->BuildGradient(argc, argv);
	return 1;
}


//
//#include <vector>
//#include <set>
//#include <queue>
//#include <time.h>
//
//#include "base/gi_timing.h"
//#include "base/gi_topological_regular_grid_3d.h"
//#include "base/gi_isolated_region_remover.h"
//#include "base/gi_isolated_region_remover_masked.h"
//#include "base/gi_numeric_integrator_path_compressing.h"
//#include "base/gi_numeric_streamline_integrator_digitizing.h"
//#include "base/gi_timing.h"
//#include "base/gi_adaptive_euler_advector_2d.h"
//#include "base/gi_adaptive_euler_advector_3d.h"
//#include "base/gi_advection_checkers.h"
//#include "base/gi_advection_events.h"
//#include "base/gi_index_comparer.h"
//#include "base/gi_maxmin_vertex_labeling.h"
//#include "base/gi_conforming_discrete_gradient.h"
//#include "base/gi_robins_sliding_regular_grid.h"
//#include "base/gi_labeling_to_bounary_labeling.h"
//#include "base/gi_topological_gradient_using_algorithms.h"
//
//
////#include "base/gi_topological_explicit_mesh_function.h"
//
//
////#include "gi_topological_region_growing_simple_gradient_builder.h"
////#include "gi_topological_convergent_gradient_builder.h"
////#include "base/gi_robin_labeling.h"
////#include "base/gi_adaptive_in_quad_euler_advector.h"
////#include "gi_numeric_integrator_2d_restricted_expanding_region.h"
////#include "gi_topological_2d_restricted_expanding_regions.h"
//#include "base/gi_topological_gradient_using_algorithms.h"
////#include "base/gi_topological_regular_grid_restricted.h"
//#include "base/gi_isolated_region_remover.h"
////#include "base/gi_topological_utility_functions.h"
////#include "gi_morse_smale_complex_basic.h"
////#include "gi_morse_smale_complex_restricted.h"
////#include "gi_kdtree.h"
////#include "gi_msc_selectors.h"
////#include "base/gi_numeric_streamline_integrator.h"
////#include "gi_experimental.h"
////#include "gi_experimental2.h"
////#include "gi_experimental3.h"
////#include "gi_convergent_2dmanifold_gradient_builder.h"
////#include "gi_experimental5.h"
//#include "base/gi_bifiltration_pairing.h"
//#include "base/gi_topological_max_vertex_mesh_function.h"
//#include "base/gi_extrema_region_builder.h"
//#include "base/gi_numeric_integrator_path_compressing.h"
////#include "base/gi_fast_robins_noalloc.h"
//
//using namespace GInt;
//
//
////typedef IsolatedRegionRemover<ComparerASC> RegionRemoverTypeASC;
//
////#include "integrate3.hpp"
////#include "vector2.hpp"
////#include <stdio.h>
//
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
////	int nn = GlobalGrid->GatherExistingNeighborsSameBdry6(xyz, neighbors);
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
//using namespace GInt;
//
//
//typedef RegularGrid3D GridType;
//typedef RegularGridTrilinearFunction GridFuncType;
//
////typedef UncachedRegularGridTrilinearFunction GridFuncType;
//typedef TopologicalRegularGrid3D MeshType;
//typedef IndexCompareLessThan<GridFuncType> ComparerASC;
//typedef IndexCompareGreaterThan<GridFuncType> ComparerDSC;
////typedef MorseSmaleComplexBasic<FLOATTYPE, MeshType, MeshFuncType, GradType> MSCType;
////typedef NumericIntegratorExpandingRegionStopWithCutoff<AdaptiveEulerAdvector3D<-1>, ComparerASC> IntegratorTypeWC;
////typedef NumericIntegratorExpandingRegionStop<AdaptiveEulerAdvector3D<-1>, ComparerASC> IntegratorTypeASC;
//typedef NumericIntegratorPathCompressingToTerminal<AdaptiveEulerAdvector3D<GridFuncType, -1>, GridFuncType > IntegratorTypeASC;
//typedef NumericIntegratorPathCompressingToTerminal<AdaptiveEulerAdvector3D<GridFuncType, 1>, GridFuncType > IntegratorTypeDSC;
//typedef IsolatedRegionRemoverMasked<ComparerASC> RegionRemoverTypeASC;
//typedef IsolatedRegionRemoverMasked<ComparerDSC> RegionRemoverTypeDSC;
//typedef DigitizingNumericStreamlineIntegrator3dASC<MeshType, GridFuncType, AdaptiveEulerAdvector3D<GridFuncType, 1> > StreamlineIntegratorTypeASC;
//typedef DigitizingNumericStreamlineIntegrator3dDSC<MeshType, GridFuncType, AdaptiveEulerAdvector3D<GridFuncType, -1> > StreamlineIntegratorTypeDSC;
//typedef DiscreteGradientLabeling<MeshType> GradType;
////typedef UncachedMaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
////typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
//typedef RegularGridMaxMinVertexLabeling3D<MeshType, GridFuncType> MaxVLType;
//typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType, 4, 6> RobinsType;
//typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> MeshFuncType;
//typedef SlidingWindowRobinsNoalloc < RegularGrid3D, RegularGridTrilinearFunction, MeshType, MaxVLType, GradType> NewRobinsType;
//
//int X, Y, Z;
//int iteration_limit;
//int per_x, per_y, per_z;
//float error_threshold, gradient_threshold;
//std::string filename;
//int parallelism = -1;
//int outputdebug = 0;
//int hacknum = 1000;
//bool usecutoff = false;
//float g_pre_simp_threshold = 0.0f;
////// we need asc man 3 dsc man 3?
//int need_ASC_3 = false;
//int need_DSC_3 = false;
//int need_ASC_1 = false;
//int need_DSC_1 = false;
//int needsad = false;
//
//bool GetOptions(int argc, char** argv) {
//	if (argc < 11) { printf("Usage: X Y Z filename error_threshold grad_threshold maxnumiter needASC1 needDSC1 needASC2 needDSC2 PresimpThesh [parallelism=ompmaxnumthreads] [outputdebug=0] [integrationinteraltimer=0]\n"); return 0; }
//	sscanf(argv[1], "%d", &X);
//	sscanf(argv[2], "%d", &Y);
//	sscanf(argv[3], "%d", &Z);
//	filename = std::string(argv[4]);
//	sscanf(argv[5], "%f", &error_threshold);
//	sscanf(argv[6], "%f", &gradient_threshold);
//	sscanf(argv[7], "%d", &iteration_limit);
//	sscanf(argv[8], "%d", &need_ASC_1);
//	sscanf(argv[9], "%d", &need_DSC_1);
//	sscanf(argv[10], "%d", &need_ASC_3);
//	sscanf(argv[11], "%d", &need_DSC_3);
//	sscanf(argv[12], "%f", &g_pre_simp_threshold);
//	if (argc >= 14)
//		sscanf(argv[13], "%d", &parallelism);
//
//	// set remaining values
//	if (parallelism != -1) {
//		omp_set_num_threads(parallelism);
//	}
//
//	printf("dims=(%d,%d,%d)\nfile=%s\nintegration parameters: e=%f, gt=%f, il=%d\nondemandacc: a1=%d, d1=%d, a2=%d, d2=%d, ps=%f\npar=%d\n",
//		X, Y, Z, argv[4], error_threshold, gradient_threshold, iteration_limit, need_ASC_1, need_DSC_1, need_ASC_3, need_DSC_3, g_pre_simp_threshold, parallelism);
//
//}
//
//
//std::chrono::steady_clock::time_point g_start_time;
//std::chrono::steady_clock::time_point task_start_time;
//std::chrono::steady_clock::time_point prior_time;
//std::chrono::steady_clock::time_point now_time;
//FILE* gtiming;
//
//void StartTask() {
//	task_start_time = std::chrono::steady_clock::now();
//	printf("starting task\n");
//}
//
//std::vector<int> timings_sequence;
//
//void EndTaskAndRecord(const char* s) {
//	now_time = std::chrono::steady_clock::now();
//	// format: [global activity name] [task] [start] [end] [dration]
//	timings_sequence.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count());
//	fprintf(gtiming, "%s %d %d %d\n", s,
//		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
//		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
//		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
//	printf("TIMING: %s %d %d %d\n", s,
//		std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - g_start_time).count(),
//		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
//		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - task_start_time).count());
//}
//
//void bigEnd(int parallelism) {
//	for (auto i : timings_sequence)  {
//		printf("%d ", i);
//	}
//	printf("\n");
//}
//RegularGrid3D* g_grid;
//GridFuncType* g_rgt_func;
//MeshType *g_topo_grid;
//IntegratorTypeASC* g_num_integrator;
////IntegratorTypeWC* g_num_integrator_with_cutoff;
//RegionRemoverTypeASC* g_region_remover;
//#ifdef USE_REGION_CLEANER
//VertexLabelingToBoundaryLabeling<INDEX_TYPE>* g_edge_map;
//#else
//VertexLabelingToBoundaryLabeling<int, MaxVLType>* g_edge_map;
//#endif
//
//MeshFuncType* g_topo_func;
//GradType *base_grad;
////RobinsLabelingAlgorithm<MeshType, MeshFuncType> *g_robin_alg;
//TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>* g_topo_alg;
//StreamlineIntegratorTypeASC* g_digitizing_streamline_integrator_asc;
//StreamlineIntegratorTypeDSC* g_digitizing_streamline_integrator_dsc;
//MaxVLType* g_maxv_labeling;
//
////void CombineLabels() {
////	INDEX_TYPE num = g_topo_grid->numCells();
////
////#pragma omp parallel for
////	for (INDEX_TYPE i = 0; i < num; i++) {
////		char v1 = g_digitizing_streamline_integrator_asc->get_output()->GetLabel(i);
////		char v2 = g_edge_map->GetOutputLabels()->GetLabel(i);
////
////		g_edge_map->GetOutputLabels()->SetLabel(i, max(v1, v2));
////	}
//
//
//
//void ReIntegrateUpFrom2Saddles(GradType* base_grad) {
//
//	// FIRST gather all the critical 2-saddles from the discrete gradient
//	StartTask();
//	std::vector<INDEX_TYPE> topo_index_partition;
//	int num_threads;
//	std::vector<std::pair<float, INDEX_TYPE> > criticals;
//	printf("gothere 2\n");
//#pragma omp parallel
//	{
//#pragma omp single
//		{
//			num_threads = omp_get_num_threads();
//			ArrayIndexPartitioner::EvenChunkSplit(g_topo_grid->numCells(), num_threads, topo_index_partition);
//		}
//
//		int thread_num = omp_get_thread_num();
//		// in parallel go through and find all 2-saddles
//		std::vector<std::pair<float, INDEX_TYPE> > lcriticals;
//		MeshType::DCellsIterator face_iterator(g_topo_grid, 2, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
//		for (face_iterator.begin(); face_iterator.valid(); face_iterator.advance()) {
//			INDEX_TYPE cell_id = face_iterator.value();
//
//			if (base_grad->getCritical(cell_id)) {
//				std::pair<float, INDEX_TYPE> p(g_topo_func->cellValue(cell_id), cell_id);
//
//				lcriticals.push_back(p);
//
//			}
//		}
//#pragma omp critical
//		{
//			criticals.insert(criticals.end(), lcriticals.begin(), lcriticals.end());
//		}
//	}
//	printf("gothere 3\n");
//	EndTaskAndRecord("Topological GatherCritical2Saddles");
//
//	StartTask();
//	std::sort(criticals.begin(), criticals.end());
//	EndTaskAndRecord("Topological Sort2Saddles");
//	//for (auto p : criticals) printf("cp: %f %llu\n", p.first, p.second);
//
//	// so criticals is now sorted list, highest value last
//
//	INDEX_TYPE total_left = criticals.size() - 1;
//#pragma omp parallel 
//	{
//		while (true) {
//			INDEX_TYPE local_id;
//			//int tots;
//#pragma omp critical
//			{
//				local_id = total_left;
//				total_left -= 1;
//				//tots = rand() % 10000;
//				//printf("%d doing %llu %llu\n", omp_get_thread_num(), local_id, total_left);
//			}
//
//			if (local_id < 0) break;
//
//			//std::vector<int> fff;
//			//fff.clear();
//			//for (int i = 0; i < tots; i++) {
//			//	fff.push_back(i * i);
//			//}
//			//printf("fff size %d\n", fff.size());
//
//			INDEX_TYPE sad_id = criticals[local_id].second;
//			// 
//			//// now find each arc
//			std::vector<INDEX_TYPE> result;
//			std::queue<INDEX_TYPE> cell_queue;
//			cell_queue.push(sad_id);
//
//			// THIS IS A SUPER FAST WAY OF filling in gometry of arc... only 2 in each direction
//			// gather the 4 hexes along paths on either side of the critical saddle
//			std::set<INDEX_TYPE> cell_visited;
//			int counter = 4;
//			while (!cell_queue.empty() && counter >= 0) {
//				INDEX_TYPE current = cell_queue.front();
//				cell_queue.pop();
//
//				cell_visited.insert(current);
//				//result.push_back(current);
//
//				MeshType::CofacetsIterator cofacets(g_topo_grid);
//				for (cofacets.begin(current); cofacets.valid(); cofacets.advance()) {
//					INDEX_TYPE temp_id = cofacets.value();
//
//					if (base_grad->getCritical(temp_id) &&
//						cell_visited.count(temp_id) == 0) {
//						result.push_back(temp_id);
//						cell_visited.insert(temp_id);
//					}
//					else if (cell_visited.count(temp_id) == 0) {
//						INDEX_TYPE pair = base_grad->getPair(temp_id);
//						result.push_back(temp_id);
//						//result.push_back(pair);
//						cell_visited.insert(temp_id);
//						cell_visited.insert(pair);
//						cell_queue.push(pair);
//					}
//				}
//				counter--;
//			}
//
//			//printf("result size %d\n", result.size());
//			for (auto arc_hex_id : result) {
//				//printf("%llu ", arc_hex_id);
//				if (g_topo_grid->dimension(arc_hex_id) != 3) continue;
//				std::vector<Vec3d> points;
//				std::vector<INDEX_TYPE> dline;
//				Vec3d seed; g_topo_grid->centroid(arc_hex_id, seed);
//				seed = seed * 0.5; // back to grid coordinates
//				g_digitizing_streamline_integrator_asc->IntegrateStreamline(seed, points, dline);
//			}
//			//printf("\n");
//			g_digitizing_streamline_integrator_asc->set_label(sad_id);
//
//
//
//		}
//
//	}
//
//
//}
//
//void ReIntegrateDownFrom1Saddles(GradType* base_grad) {
//
//	// FIRST gather all the critical 2-saddles from the discrete gradient
//	StartTask();
//	std::vector<INDEX_TYPE> topo_index_partition;
//	int num_threads;
//	std::vector<std::pair<float, INDEX_TYPE> > criticals;
//	printf("gothere 2\n");
//#pragma omp parallel
//	{
//#pragma omp single
//		{
//			num_threads = omp_get_num_threads();
//			ArrayIndexPartitioner::EvenChunkSplit(g_topo_grid->numCells(), num_threads, topo_index_partition);
//		}
//
//		int thread_num = omp_get_thread_num();
//		// in parallel go through and find all 2-saddles
//		std::vector<std::pair<float, INDEX_TYPE> > lcriticals;
//		MeshType::DCellsIterator face_iterator(g_topo_grid, 1, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
//		for (face_iterator.begin(); face_iterator.valid(); face_iterator.advance()) {
//			INDEX_TYPE cell_id = face_iterator.value();
//
//			if (base_grad->getCritical(cell_id)) {
//				std::pair<float, INDEX_TYPE> p(-1 * g_topo_func->cellValue(cell_id), cell_id);
//
//				lcriticals.push_back(p);
//
//			}
//		}
//#pragma omp critical
//		{
//			criticals.insert(criticals.end(), lcriticals.begin(), lcriticals.end());
//		}
//	}
//	printf("gothere 3\n");
//	EndTaskAndRecord("Topological GatherCritical1Saddles");
//
//	StartTask();
//	std::sort(criticals.begin(), criticals.end());
//	EndTaskAndRecord("Topological Sort2Saddles");
//	//for (auto p : criticals) printf("cp: %f %llu\n", p.first, p.second);
//
//	// so criticals is now sorted list, highest value last
//
//	INDEX_TYPE total_left = criticals.size() - 1;
//#pragma omp parallel 
//	{
//		while (true) {
//			INDEX_TYPE local_id;
//			//int tots;
//#pragma omp critical
//			{
//				local_id = total_left;
//				total_left -= 1;
//				//tots = rand() % 10000;
//				//printf("%d doing %llu %llu\n", omp_get_thread_num(), local_id, total_left);
//			}
//
//			if (local_id < 0) break;
//
//			//std::vector<int> fff;
//			//fff.clear();
//			//for (int i = 0; i < tots; i++) {
//			//	fff.push_back(i * i);
//			//}
//			//printf("fff size %d\n", fff.size());
//
//			INDEX_TYPE sad_id = criticals[local_id].second;
//			// 
//			//// now find each arc
//			std::vector<INDEX_TYPE> result;
//			std::queue<INDEX_TYPE> cell_queue;
//			cell_queue.push(sad_id);
//
//			// THIS IS A SUPER FAST WAY OF filling in gometry of arc... only 2 in each direction
//			// gather the 4 vertices along paths on either side of the critical saddle
//			std::set<INDEX_TYPE> cell_visited;
//			int counter = 4;
//			while (!cell_queue.empty() && counter >= 0) {
//				INDEX_TYPE current = cell_queue.front();
//				cell_queue.pop();
//
//				cell_visited.insert(current);
//				//result.push_back(current);
//
//				MeshType::FacetsIterator facets(g_topo_grid);
//				for (facets.begin(current); facets.valid(); facets.advance()) {
//					INDEX_TYPE temp_id = facets.value();
//
//					if (base_grad->getCritical(temp_id) &&
//						cell_visited.count(temp_id) == 0) {
//						result.push_back(temp_id);
//						cell_visited.insert(temp_id);
//					}
//					else if (cell_visited.count(temp_id) == 0) {
//						INDEX_TYPE pair = base_grad->getPair(temp_id);
//						result.push_back(temp_id);
//						//result.push_back(pair);
//						cell_visited.insert(temp_id);
//						cell_visited.insert(pair);
//						cell_queue.push(pair);
//					}
//				}
//				counter--;
//			}
//
//			//printf("result size %d\n", result.size());
//			for (auto arc_vert_id : result) {
//				//printf("%llu ", arc_hex_id);
//				if (g_topo_grid->dimension(arc_vert_id) != 0) continue;
//				std::vector<Vec3d> points;
//				std::vector<INDEX_TYPE> dline;
//				Vec3d seed; g_topo_grid->centroid(arc_vert_id, seed);
//				seed = seed * 0.5; // back to grid coordinates
//				g_digitizing_streamline_integrator_dsc->IntegrateStreamline(seed, points, dline);
//			}
//			//printf("\n");
//			g_digitizing_streamline_integrator_dsc->set_label(sad_id);
//
//
//
//		}
//
//	}
//
//
//}
//
//
//
////-------------------------------------------------------------
////-------------------------------------------------------------
////-------------------------------------------------------------
//// DO LOCAL SIMPLIFICATION
////-------------------------------------------------------------
////-------------------------------------------------------------
////-------------------------------------------------------------
////-------------------------------------------------------------
//
//void BasicHardSimplifyGradient(GradType* grad, DenseLabeling<char>* restriction_labels) {
//
//
//	g_topo_alg = new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(g_topo_func, g_topo_grid, grad);
//	printf("before local cancellations:\n");
//	g_topo_alg->count_critical_points(4);
//	printf("checking for loops:\n");
//	g_topo_alg->CheckGradientForLoops();
//	printf("done\n");
//
//	int counts[3];
//	for (int i = 0; i < 3; i++) counts[i] = 0;
//	//for (int k = 0; k < 2; k++) {
//	MeshType::DCellsIterator eit(g_topo_grid, 0);
//	for (eit.begin(); eit.valid(); eit.advance()) {
//		INDEX_TYPE eid = eit.value();
//		//INDEX_TYPE v1 = g_maxv_labeling->Cell2HighestVertex(eid);
//		if (grad->getCritical(eid)) {
//			char lab1 = restriction_labels->GetLabel(eid);
//			std::vector<INDEX_TYPE> candidates;
//			MeshType::CofacetsIterator cfit(g_topo_grid);
//			for (cfit.begin(eid); cfit.valid(); cfit.advance()) {
//				INDEX_TYPE fid = cfit.value();
//
//				if (grad->getCritical(fid)) {
//					char lab2 = restriction_labels->GetLabel(fid);
//					if (lab1 == lab2) candidates.push_back(fid);
//				}
//			}
//			//if (candidates.size() == 1) {
//			for (auto fid : candidates) {
//				if (g_topo_grid->boundaryValue(eid) == g_topo_grid->boundaryValue(fid)
//					) {
//					counts[0]++;
//					grad->setPair(eid, fid);
//					grad->setPair(fid, eid);
//					break;
//					//printf("%d: %llu -> %llu : %d - %d\n",k, eid, fid, lab1, lab2);
//				}
//			}
//			//else if (lab1 != lab2) {
//			//	counts[g_topo_grid->dimension(eid)]++;
//			//	continue;
//			//}
//			//else {
//			//	if (g_topo_grid->boundaryValue(eid) == g_topo_grid->boundaryValue(fid)
//			//		//lab1 != lab2 &&
//			//		/*g_maxv_labeling->Before(g_topo_grid->VertexNumberFjkromCellID(g_robin_alg->lowest_vertex(fid)),
//			//		g_topo_grid->VertexNumberFromCellID(g_robin_alg->lowest_vertex(eid)))*/) {
//			//		base_grad->setPair(eid, fid);
//			//		base_grad->setPair(fid, eid);
//			//		break;
//			//		//printf("%d: %llu -> %llu : %d - %d\n",k, eid, fid, lab1, lab2);
//			//	}
//			//}
//		}
//		//}
//
//
//	}
//	MeshType::DCellsIterator eit2(g_topo_grid, 3);
//	for (eit2.begin(); eit2.valid(); eit2.advance()) {
//		INDEX_TYPE eid = eit2.value();
//		//INDEX_TYPE v1 = g_maxv_labeling->Cell2HighestVertex(eid);
//		if (grad->getCritical(eid)) {
//			char lab1 = restriction_labels->GetLabel(eid);
//			std::vector<INDEX_TYPE> candidates;
//			MeshType::FacetsIterator cfit(g_topo_grid);
//			for (cfit.begin(eid); cfit.valid(); cfit.advance()) {
//				INDEX_TYPE fid = cfit.value();
//
//				if (grad->getCritical(fid)) {
//					char lab2 = restriction_labels->GetLabel(fid);
//					if (lab1 == lab2) candidates.push_back(fid);
//				}
//			}
//			//if (candidates.size() == 1) {
//			for (auto fid : candidates) {
//				if (g_topo_grid->boundaryValue(eid) == g_topo_grid->boundaryValue(fid)
//					) {
//					counts[2]++;
//					grad->setPair(eid, fid);
//					grad->setPair(fid, eid);
//					break;
//					//printf("%d: %llu -> %llu : %d - %d\n",k, eid, fid, lab1, lab2);
//				}
//			}
//			//else if (lab1 != lab2) {
//			//	counts[g_topo_grid->dimension(eid)]++;
//			//	continue;
//			//}
//			//else {
//			//	if (g_topo_grid->boundaryValue(eid) == g_topo_grid->boundaryValue(fid)
//			//		//lab1 != lab2 &&
//			//		/*g_maxv_labeling->Before(g_topo_grid->VertexNumberFjkromCellID(g_robin_alg->lowest_vertex(fid)),
//			//		g_topo_grid->VertexNumberFromCellID(g_robin_alg->lowest_vertex(eid)))*/) {
//			//		base_grad->setPair(eid, fid);
//			//		base_grad->setPair(fid, eid);
//			//		break;
//			//		//printf("%d: %llu -> %llu : %d - %d\n",k, eid, fid, lab1, lab2);
//			//	}
//			//}
//		}
//		//}
//
//
//	}
//	MeshType::DCellsIterator eit3(g_topo_grid, 2);
//	for (eit3.begin(); eit3.valid(); eit3.advance()) {
//		INDEX_TYPE eid = eit3.value();
//		//INDEX_TYPE v1 = g_maxv_labeling->Cell2HighestVertex(eid);
//		if (grad->getCritical(eid)) {
//			char lab1 = restriction_labels->GetLabel(eid);
//			std::vector<INDEX_TYPE> candidates;
//			MeshType::FacetsIterator cfit(g_topo_grid);
//			for (cfit.begin(eid); cfit.valid(); cfit.advance()) {
//				INDEX_TYPE fid = cfit.value();
//
//				if (grad->getCritical(fid)) {
//					char lab2 = restriction_labels->GetLabel(fid);
//					if (lab1 == lab2) candidates.push_back(fid);
//				}
//			}
//			//if (candidates.size() == 1) {
//			for (auto fid : candidates) {
//				if (g_topo_grid->boundaryValue(eid) == g_topo_grid->boundaryValue(fid)
//					) {
//					counts[1]++;
//					grad->setPair(eid, fid);
//					grad->setPair(fid, eid);
//					break;
//					//printf("%d: %llu -> %llu : %d - %d\n",k, eid, fid, lab1, lab2);
//				}
//			}
//			//else if (lab1 != lab2) {
//			//	counts[g_topo_grid->dimension(eid)]++;
//			//	continue;
//			//}
//			//else {
//			//	if (g_topo_grid->boundaryValue(eid) == g_topo_grid->boundaryValue(fid)
//			//		//lab1 != lab2 &&
//			//		/*g_maxv_labeling->Before(g_topo_grid->VertexNumberFjkromCellID(g_robin_alg->lowest_vertex(fid)),
//			//		g_topo_grid->VertexNumberFromCellID(g_robin_alg->lowest_vertex(eid)))*/) {
//			//		base_grad->setPair(eid, fid);
//			//		base_grad->setPair(fid, eid);
//			//		break;
//			//		//printf("%d: %llu -> %llu : %d - %d\n",k, eid, fid, lab1, lab2);
//			//	}
//			//}
//		}
//		//}
//
//
//	}
//
//	//}
//
//	//}
//	for (int i = 0; i < 3; i++) printf("dim-%d crits from mismatch\n", counts[i]);
//
//	printf("after local cancellations:\n");
//	g_topo_alg->count_critical_points(4);
//
//	printf("chekcing for loops:\n");
//	g_topo_alg->CheckGradientForLoops();
//}
//
//
//
//
//void RecordGrad(GradType* grad, const char* gradname) {
//	printf("setting dim asc man\n");
//	g_topo_alg->setAscendingManifoldDimensions();
//	printf("outputting to file %s\n", gradname);
//	grad->outputToFile(gradname);
//	//return 1;
//}
//
//GradType* stageGrad[4];
//DenseLabeling<char>* stageLabel[4];
//
//int main(int argc, char** argv) {
//
//	// read command line options
//	GetOptions(argc, argv);
//
//	printf("ondemand accuracy: a1=%d d1=%d a2=%d d2=%d\n", need_ASC_1, need_DSC_1, need_ASC_3, need_DSC_3);
//
//	printf("SANITY\n");
//	// start timing overall algorithm
//	ThreadedTimer timer(1);
//	timer.StartGlobal();
//
//
//	char gradname[1024];
//	sprintf(gradname, "%s.grad", argv[4]);
//
//
//	// will write timing to this file
//	char timingname[2048];
//	sprintf(timingname, "%s.%03d.gtime.txt", argv[4], parallelism);
//	gtiming = fopen(timingname, "w");
//
//	// START IO ---------------------------
//	//StartTask();
//	g_grid = new RegularGrid3D(Vec3i(X, Y, Z), Vec3b(per_x, per_y, per_z));
//	g_rgt_func = new GridFuncType(g_grid);
//	g_rgt_func->LoadImageFromFile(filename.c_str());
//	//EndTaskAndRecord(Read Data");
//
//
//	g_start_time = std::chrono::steady_clock::now();
//	// CREATE MAX VERTEX LABELING
//	StartTask();
//	g_topo_grid = new MeshType(g_grid);
//
//	g_maxv_labeling = new MaxVLType(g_topo_grid, g_rgt_func);
//	g_maxv_labeling->ComputeOutput();
//	EndTaskAndRecord("Topological MaxVlabel");
//	// create a topology function
//	StartTask();
//	printf("topo function...\n");
//	g_topo_func = new MeshFuncType();
//	g_topo_func->setMeshAndFuncAndMaxVLabeling(g_topo_grid, g_rgt_func, g_maxv_labeling);
//	EndTaskAndRecord("Topological Function");
//
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	// DO FIRST DISCRETE GRADIENT COMPUTATION WITH NO RESTRICTION
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//
//	StartTask();
//	base_grad = new GradType(g_topo_grid);
//	base_grad->ClearAllGradient();
//	RobinsType* first_robins = new RobinsType(g_topo_grid, g_maxv_labeling, base_grad);
//	first_robins->ComputePairing();
//	EndTaskAndRecord("Topological Robins0");
//
//	//g_topo_alg = new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
//	//printf("after base first robins:\n");
//	//g_topo_alg->count_critical_points(4);
//
//	if (!(need_ASC_3 || need_DSC_3 || need_ASC_1 || need_DSC_1 || needsad)) {
//		g_topo_alg = new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
//		printf("no accuraccy needed, outputting\n");
//		RecordGrad(base_grad, gradname);
//		return 1;
//	}
//
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	// DO VOLUME ACCURATE GRADIENT COMPUTATION WITH NO RESTRICTION
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//
//	//-------------------------------------------------------------
//	// IF WE WANT NUMERIC ACCURACY, WE WILL NEED NUMEIC GRADIENT
//	//-------------------------------------------------------------
//
//	// Do gradient vectors computation from raw image data
//	printf("computing gradient\n");
//	StartTask();
//	g_rgt_func->ComputeGradFromImage(1);
//	//g_rgt_func->Negate();
//	EndTaskAndRecord("NumericalTracing GradCompute");
//
//	// CREATE RESTRICTION MAP - WILL NEED
//	DenseLabeling<char>* restriction_labels = new DenseLabeling<char>(g_topo_grid->numCells());
//	restriction_labels->SetAll(0);
//	// we will always need a constrained robins alg
//	RobinsType* constrained_robins = new RobinsType(g_topo_grid, g_maxv_labeling, restriction_labels, base_grad);
//
//
//
//	//-------------------------------------------------------------
//	// IF WE WANT ACCURATE 3-MANIFOLDS, COMPUTE SIMPLIFIED MAPS
//	//-------------------------------------------------------------
//
//	if (need_ASC_3 || need_DSC_3 || needsad) {
//		StartTask();
//		SimplifiedExtremumGraph<MeshType, MeshFuncType, GradType>* simplified_ext_graph =
//			new SimplifiedExtremumGraph<MeshType, MeshFuncType, GradType>(g_topo_grid, g_topo_func, base_grad);
//		if ((need_ASC_3 && need_DSC_3) || needsad) {
//			simplified_ext_graph->SetMode(SimplifiedExtremumGraph<MeshType, MeshFuncType, GradType>::EXTGRAPHMODE::BOTH);
//		}
//		else if (need_ASC_3) {
//			simplified_ext_graph->SetMode(SimplifiedExtremumGraph<MeshType, MeshFuncType, GradType>::EXTGRAPHMODE::MINS);
//		}
//		else {
//			simplified_ext_graph->SetMode(SimplifiedExtremumGraph<MeshType, MeshFuncType, GradType>::EXTGRAPHMODE::MAXS);
//		}
//		simplified_ext_graph->ComputeMinMapFromGradient(g_pre_simp_threshold);
//		EndTaskAndRecord("Topological SimplifiedGraph");
//		printf("done creating simplified extremum graph\n");
//
//		g_edge_map = new VertexLabelingToBoundaryLabeling<int, MaxVLType>(g_topo_grid, restriction_labels);
//		g_edge_map->InitializeFirst();
//
//		//UFMergeGraph<MeshFuncType>* asdf = new UFMergeGraph<MeshFuncType>();
//		// NOW BUILD A TERMINAL MAP FROM THE SIMPLIFIED MIN/MAX HIERARCHIES
//		// DO ASCENDING MANIFOLDS
//
//		if (need_ASC_3 || needsad) {
//			// a simplified map means that we only care about accurate boundaries for 
//			// extrema that persist above a given threshold. 
//			// per Julien's observation, we can compute the extremal simplification graphs
//			// without spending the cost of doing a full MS Complex (even though we have the gradient)
//
//
//			StartTask();
//         std::unordered_map<INDEX_TYPE, INT_TYPE> extmapASC;
//			for (auto p : simplified_ext_graph->mMinGraph->mCellIndexToListIdMap) {
//				extmapASC[p.first] = simplified_ext_graph->mMinGraph->Representative(p.second);
//			}
//			GridSimplifiedExtremalRegionBuilder<ComparerASC, GridFuncType, MeshType>* test_simp_reg_builder_asc =
//				new GridSimplifiedExtremalRegionBuilder<ComparerASC, GridFuncType, MeshType>(g_rgt_func, g_grid, g_topo_grid);
//			test_simp_reg_builder_asc->BeginIntegration(extmapASC);
//			EndTaskAndRecord("Topological SimExtRBASC");
//
//			//test_simp_reg_builder_asc->GetOutputLabels()->OutputToIntFile("extremal_asc.raw");
//
//			// now do actual numeric integration using simplified extremum regions map as target
//			StartTask();
//			IntegratorTypeASC* newintegrator_asc =
//				new IntegratorTypeASC(g_rgt_func, g_grid, error_threshold,
//				gradient_threshold, iteration_limit);
//			newintegrator_asc->BeginIntegration(test_simp_reg_builder_asc->GetIdMap(),
//				test_simp_reg_builder_asc->GetOutputLabels(), true);
//			//newintegrator_asc->GetOutputLabels()->OutputToIntFile("newintegrator_asc.raw");
//			EndTaskAndRecord("NumericalTracing NewIntegrationASC");
//			//test_simp_reg_builder_asc->GetOutputLabels()->OutputToIntFile("integrated_asc.raw");
//			// REMOVE DISCONNECTED COMPONENTS
//			StartTask();
//			IsolatedCCRegionRemoverNEW<ComparerASC, GridFuncType>* cleaner1_asc =
//				new IsolatedCCRegionRemoverNEW<ComparerASC, GridFuncType>(g_rgt_func, newintegrator_asc->GetOutputLabels());
//			printf("Removing Disconnected Regions\n");
//			cleaner1_asc->ComputeOutput();
//			printf("here2\n");
//			EndTaskAndRecord("Removed CleanerASC");
//			//newintegrator_asc->GetOutputLabels()->OutputToIntFile("newintegrator_cleaned_asc.raw");
//			
//			//// TEST FIX
//			//auto lab = test_simp_reg_builder_asc->GetOutputLabels();
//			//map<pair<int, int>, int> counter;
//			//for (INDEX_TYPE id = 0; id < g_grid->NumElements(); id++) {
//
//			//	Vec3l t_neighbors[6];
//			//	Vec3l t_coords = g_grid->XYZ3d(id);
//			//	int t_num_neighbors = g_grid->GatherExistingNeighborsSameBdry6(t_coords, t_neighbors);
//			//	int lab1 = lab->GetLabel(id);
//
//			//	for (int i = 0; i < t_num_neighbors; i++) {
//			//		INDEX_TYPE t_neighbor_vertex = g_grid->Index3d(t_neighbors[i]);
//			//		int lab2 = lab->GetLabel(t_neighbor_vertex);
//			//		if (lab1 != lab2) {
//			//			if (lab1 < lab2) {
//			//				pair<int, int> p(lab1, lab2);
//			//				if (counter.count(p) == 0) counter[p] = 1;
//			//				else (counter[p]++);
//			//			}
//			//			else {
//			//				pair<int, int> p(lab2, lab1);
//			//				if (counter.count(p) == 0) counter[p] = 1;
//			//				else (counter[p]++);
//			//			}
//			//		}
//			//	}
//
//
//			//}
//			//for (auto p : counter) {
//			//	if (p.first.first != -1 && p.second > 2000)
//			//		printf("<%d, %d>=%d\n", p.first.first, p.first.second, p.second);
//			//}
//			// add this guy's contribution to the restriction labeling
//			StartTask();
//			g_edge_map->ComputeMINBoundary(newintegrator_asc->GetOutputLabels());
//			EndTaskAndRecord("Topological EdgeMapASC");
//		}
//
//		if (need_DSC_3 || needsad) {
//
//			// DO DESCENDING MANIFOLDS
//			StartTask();
//         std::unordered_map<INDEX_TYPE, INT_TYPE> extmapDSC;
//			for (auto p : simplified_ext_graph->mMaxGraph->mCellIndexToListIdMap) {
//				extmapDSC[g_maxv_labeling->Cell2HighestVertex(p.first)] = simplified_ext_graph->mMaxGraph->Representative(p.second);
//			}
//			GridSimplifiedExtremalRegionBuilder<ComparerDSC, GridFuncType, MeshType>* test_simp_reg_builder_dsc =
//				new GridSimplifiedExtremalRegionBuilder<ComparerDSC, GridFuncType, MeshType>(g_rgt_func, g_grid, g_topo_grid);
//			test_simp_reg_builder_dsc->BeginIntegration(extmapDSC);
//			EndTaskAndRecord("Topological SimExtRBDSC");
//
//			//test_simp_reg_builder_dsc->GetOutputLabels()->OutputToIntFile("extremal_dsc.raw");
//
//
//			StartTask();
//
//			IntegratorTypeDSC* newintegrator_dsc =
//				new IntegratorTypeDSC(g_rgt_func, g_grid, error_threshold,
//				gradient_threshold, iteration_limit);
//			newintegrator_dsc->BeginIntegration(test_simp_reg_builder_dsc->GetIdMap(),
//				test_simp_reg_builder_dsc->GetOutputLabels(), true);
//
//			//newintegrator_dsc->GetOutputLabels()->OutputToIntFile("newintegrator_dsc.raw");
//			EndTaskAndRecord("NumericalTracing NewIntegrationDSC");
//			//test_simp_reg_builder_dsc->GetOutputLabels()->OutputToIntFile("integrated_dsc.raw");
//
//			// REMOVE DISCONNECTED COMPONENTS
//			StartTask();
//			IsolatedCCRegionRemoverNEW<ComparerDSC, GridFuncType>* cleaner1_dsc =
//				new IsolatedCCRegionRemoverNEW<ComparerDSC, GridFuncType>(g_rgt_func, newintegrator_dsc->GetOutputLabels());
//			printf("Removing Disconnected Regions\n");
//			cleaner1_dsc->ComputeOutput();
//			printf("here2\n");
//			EndTaskAndRecord("NumericalTracing CleaningDSC");
//			//newintegrator_dsc->GetOutputLabels()->OutputToIntFile("newintegrator_cleaned_dsc.raw");
//
//
//			//// TEST FIX
//			//auto lab = test_simp_reg_builder_dsc->GetOutputLabels();
//			//map<pair<int, int>, int> counter;
//			//for (INDEX_TYPE id = 0; id < g_grid->NumElements(); id++) {
//
//			//	Vec3l t_neighbors[6];
//			//	Vec3l t_coords = g_grid->XYZ3d(id);
//			//	int t_num_neighbors = g_grid->GatherExistingNeighborsSameBdry6(t_coords, t_neighbors);
//			//	int lab1 = lab->GetLabel(id);
//
//			//	for (int i = 0; i < t_num_neighbors; i++) {
//			//		INDEX_TYPE t_neighbor_vertex = g_grid->Index3d(t_neighbors[i]);
//			//		int lab2 = lab->GetLabel(t_neighbor_vertex);
//			//		if (lab1 != lab2) {
//			//			if (lab1 < lab2) {
//			//				pair<int, int> p(lab1, lab2);
//			//				if (counter.count(p) == 0) counter[p] = 1;
//			//				else (counter[p]++);
//			//			}
//			//			else {
//			//				pair<int, int> p(lab2, lab1);
//			//				if (counter.count(p) == 0) counter[p] = 1;
//			//				else (counter[p]++);
//			//			}
//			//		}
//			//	}
//
//
//			//}
//			//for (auto p : counter) {
//			//	if (p.first.first != -1 && p.second > 2000)
//			//		printf("<%d, %d>=%d\n", p.first.first, p.first.second, p.second);
//			//}
//
//			// add this guy's contribution to the restriction labeling
//			StartTask();
//			g_edge_map->ComputeMAXBoundary(newintegrator_dsc->GetOutputLabels(), g_maxv_labeling);
//			EndTaskAndRecord("Topological EdgeMapDSC");
//
//		}
//		//restriction_labels->OutputToFile("boundary_labels_after_3m.raw");
//		if (outputdebug) {
//			//restriction_labels->OutputToFile("boundary_labels_after_3m.raw");
//		}
//	}
//	
//
//	//-------------------------------------------------------------
//	// IF ALL WE NEED IS ACCURATE ASC/DSC 3-m then recompute grad and exit
//	//-------------------------------------------------------------
//	if (!(need_ASC_1 || need_DSC_1)) {
//		StartTask();
//		printf("redoing discrete gradient in changed boundarids\n");
//		std::vector<INDEX_TYPE> topo_index_partition;
//
//		int num_threads;
//#pragma omp parallel
//		{
//#pragma omp single
//		{
//			num_threads = omp_get_num_threads();
//			ArrayIndexPartitioner::EvenChunkSplit(g_topo_grid->numCells(), num_threads, topo_index_partition);
//		}
//
//			int thread_num = omp_get_thread_num();
//			INDEX_TYPE threadfixcont = 0;
//			// iterate over all vertices
//			MeshType::DCellsIterator verts(g_topo_grid, 0, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
//			for (verts.begin(); verts.valid(); verts.advance()) {
//				INDEX_TYPE vert_GI = verts.value();
//				bool hasdiff = false;
//				if (restriction_labels->GetLabel(vert_GI) > 0) hasdiff = true;
//				if (!hasdiff) {
//					MeshType::AdjacentCellsIterator edgeit(g_topo_grid);
//					//bool hasdiff = false;
//
//					for (edgeit.begin(vert_GI); edgeit.valid(); edgeit.advance()) {
//						INDEX_TYPE edge_GI = edgeit.value();
//						if (restriction_labels->GetLabel(edge_GI) > 0 &&
//							g_maxv_labeling->Cell2HighestVertex(edge_GI) == vert_GI) {
//							hasdiff = true;
//							break;
//						}
//					}
//				}
//				if (hasdiff) {
//					constrained_robins->ComputeLowerStar(vert_GI);
//					threadfixcont++;
//				}
//			}
//#pragma omp critical
//		{
//			printf("thread %d did %llu fixed vertices\n", thread_num, threadfixcont);
//		}
//		}
//		EndTaskAndRecord("Topological ConformingGrad");
//
//		//printf("asdf\n");
//		//BasicHardSimplifyGradient(base_grad, restriction_labels);
//
//		g_topo_alg = new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
//		//printf("after Second robins:\n");
//		//g_topo_alg->count_critical_points(4);
//		RecordGrad(base_grad, gradname);
//		return 1;
//	}
//
//
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	// DO LINES ACCURATE GRADIENT COMPUTATION 
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//	//-------------------------------------------------------------
//
//	if (need_ASC_1) {
//		StartTask();
//		int labeltarget = 4;
//		g_digitizing_streamline_integrator_asc = new StreamlineIntegratorTypeASC(g_grid, g_rgt_func, g_topo_grid, error_threshold, gradient_threshold, iteration_limit);
//		g_digitizing_streamline_integrator_asc->SetDigitizingTarget(restriction_labels, labeltarget);
//		ReIntegrateUpFrom2Saddles(base_grad);
//		EndTaskAndRecord("Numerical SaddleIntASC");
//	}
//	
//	if (need_DSC_1) {
//		StartTask();
//		int labeltarget = 5;
//		g_digitizing_streamline_integrator_dsc = new StreamlineIntegratorTypeDSC(g_grid, g_rgt_func, g_topo_grid, error_threshold, gradient_threshold, iteration_limit);
//		g_digitizing_streamline_integrator_dsc->SetDigitizingTarget(restriction_labels, labeltarget);
//		ReIntegrateDownFrom1Saddles(base_grad);
//		EndTaskAndRecord("Numerical SaddleIntDSC");
//	}
//
//	//restriction_labels->OutputToFile("boundary_labels_final.raw");
//
//
//	//-------------------------------------------------------------
//	// NOW REDO DISCRETE GRADIENT IN NEIGHBORHOOD
//	//-------------------------------------------------------------
//	StartTask();
//	std::vector<INDEX_TYPE> topo_index_partition;
//	int num_threads;
//#pragma omp parallel
//	{
//#pragma omp single
//		{
//			num_threads = omp_get_num_threads();
//			ArrayIndexPartitioner::EvenChunkSplit(g_topo_grid->numCells(), num_threads, topo_index_partition);
//		}
//
//		int thread_num = omp_get_thread_num();
//		INDEX_TYPE threadfixcont = 0;
//		// iterate over all vertices
//		MeshType::DCellsIterator verts(g_topo_grid, 0, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
//		for (verts.begin(); verts.valid(); verts.advance()) {
//			INDEX_TYPE vert_GI = verts.value();
//			bool hasdiff = false;
//			if (restriction_labels->GetLabel(vert_GI) > 0) hasdiff = true;
//			if (!hasdiff) {
//				MeshType::AdjacentCellsIterator cocells(g_topo_grid);
//				//bool hasdiff = false;
//
//				for (cocells.begin(vert_GI); cocells.valid(); cocells.advance()) {
//					INDEX_TYPE cocell_GI = cocells.value();
//					if (restriction_labels->GetLabel(cocell_GI) > 0 &&
//						g_maxv_labeling->Cell2HighestVertex(cocell_GI) == vert_GI) {
//						hasdiff = true;
//						break;
//					}
//				}
//			}
//			if (hasdiff) {
//				constrained_robins->ComputeLowerStar(vert_GI);
//				threadfixcont++;
//			}
//		}
//#pragma omp critical
//		{
//			printf("thread %d did %llu fixed vertices\n", thread_num, threadfixcont);
//		}
//	}
//	EndTaskAndRecord("Topological ConformingGrad");
//
//	printf("asdf\n");
//	//BasicHardSimplifyGradient(base_grad, restriction_labels);
//	//bigEnd(parallelism, );
//	//return 1;
//	g_topo_alg = new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
//	printf("after Second robins:\n");
//	g_topo_alg->count_critical_points(4);
//
//
//	RecordGrad(base_grad, gradname);
//	fclose(gtiming);
//	 
//
//	return 1;
//};
//
//
