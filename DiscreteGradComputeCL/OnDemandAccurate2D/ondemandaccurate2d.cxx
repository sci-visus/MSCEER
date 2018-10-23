/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/




int main(int argc, char** argv) {
	//GInt::OndemandDiscreteGradientBuilder* builder = new GInt::OndemandDiscreteGradientBuilder();
	//builder->BuildGradient(argc, argv);
	return 1;
}

//#include "base/gi_ondemand_accurate_grad_builder_2d.h"
//
//int main(int argc, char** argv) {
//	GInt::OndemandDiscreteGradientBuilder* builder = new GInt::OndemandDiscreteGradientBuilder();
//	builder->BuildGradient(argc, argv);
//	return 1;
//}
//
//
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
////#include "gi_topological_region_growing_simple_gradient_builder.h"
////#include "gi_topological_convergent_gradient_builder.h"
//#include "gi_robin_labeling.h"
//#include "gi_adaptive_in_quad_euler_advector.h"
////#include "gi_numeric_integrator_2d_restricted_expanding_region.h"
////#include "gi_topological_2d_restricted_expanding_regions.h"
//#include "gi_topological_gradient_using_algorithms.h"
//#include "gi_topological_regular_grid_restricted.h"
//#include "gi_isolated_region_remover.h"
//#include "gi_isolated_region_remover2.h"
//#include "gi_topological_utility_functions.h"
////#include "gi_morse_smale_complex_basic.h"
////#include "gi_morse_smale_complex_restricted.h"
////#include "gi_kdtree.h"
////#include "gi_msc_selectors.h"
//#include "gi_numeric_streamline_integrator.h"
//#include "gi_numeric_integrator_expanding_region_stop_filtered2.h"
//#include "gi_numeric_integrator_expanding_region_stop_filtered.h"
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
////#include "gi_topological_region_growing_simple_gradient_builder.h"
////#include "gi_topological_convergent_gradient_builder.h"
//#include "gi_robin_labeling.h"
//#include "gi_adaptive_in_quad_euler_advector.h"
////#include "gi_numeric_integrator_2d_restricted_expanding_region.h"
////#include "gi_topological_2d_restricted_expanding_regions.h"
//#include "gi_topological_gradient_using_algorithms.h"
//#include "gi_topological_regular_grid_restricted.h"
//#include "gi_isolated_region_remover.h"
//#include "gi_topological_utility_functions.h"
////#include "gi_morse_smale_complex_basic.h"
////#include "gi_morse_smale_complex_restricted.h"
////#include "gi_kdtree.h"
////#include "gi_msc_selectors.h"
//#include "gi_numeric_streamline_integrator.h"
////#include "gi_experimental.h"
////#include "gi_experimental2.h"
////#include "gi_experimental3.h"
////#include "gi_convergent_2dmanifold_gradient_builder.h"
////#include "gi_experimental5.h"
//#include "gi_bifiltration_pairing.h"
//#include "gi_topological_max_vertex_mesh_function.h"
//#include "gi_extrema_region_builder.h"
//#include "gi_numeric_integrator_path_compressing.h"
//#include "gi_fast_robins_noalloc.h"
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
//typedef RegularGrid2D GridType;
//typedef RegularGridBilinearFunction GridFuncType;
//
////typedef UncachedRegularGridTrilinearFunction GridFuncType;
//typedef TopologicalRegularGrid2D MeshType;
//typedef IndexCompareLessThan<GridFuncType> ComparerASC;
//typedef IndexCompareGreaterThan<GridFuncType> ComparerDSC;
////typedef MorseSmaleComplexBasic<FLOATTYPE, MeshType, MeshFuncType, GradType> MSCType;
////typedef NumericIntegratorExpandingRegionStopWithCutoff<AdaptiveEulerAdvector3D<-1>, ComparerASC> IntegratorTypeWC;
////typedef NumericIntegratorExpandingRegionStop<AdaptiveEulerAdvector3D<-1>, ComparerASC> IntegratorTypeASC;
//typedef DigitizingNumericStreamlineIntegrator2dASC<MeshType, GridFuncType, AdaptiveEulerAdvector2D<GridFuncType, 1> > StreamlineIntegratorTypeASC;
//typedef DigitizingNumericStreamlineIntegrator2dDSC<MeshType, GridFuncType, AdaptiveEulerAdvector2D<GridFuncType, -1> > StreamlineIntegratorTypeDSC;
//typedef DiscreteGradientLabeling<MeshType> GradType;
////typedef UncachedMaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
////typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
//typedef RegularGridMaxMinVertexLabeling2D<MeshType, GridFuncType> MaxVLType;
//typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType, 4, 6> RobinsType;
//typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> MeshFuncType;
////typedef SlidingWindowRobinsNoalloc < RegularGrid2D, RegularGridTrilinearFunction, MeshType, MaxVLType, GradType> NewRobinsType;
//
//int X, Y;
//int iteration_limit;
//int per_x, per_y;
//float error_threshold, gradient_threshold;
//std::string filename;
//int parallelism = -1;
//int outputdebug = 0;
//int hacknum = 1000;
//bool usecutoff = false;
//float g_pre_simp_threshold = 0.0f;
////// we need asc man 3 dsc man 3?
//
//int need_ASC_1 = false;
//int need_DSC_1 = false;
//int needsad = false;
//
//bool GetOptions(int argc, char** argv) {
//	if (argc < 11) { printf("Usage: X Y filename error_threshold grad_threshold maxnumiter needASC1 needDSC1 PresimpThesh [parallelism=ompmaxnumthreads] [outputdebug=0] [integrationinteraltimer=0]\n"); return 0; }
//	sscanf(argv[1], "%d", &X);
//	sscanf(argv[2], "%d", &Y);
//	filename = std::string(argv[3]);
//	sscanf(argv[4], "%f", &error_threshold);
//	sscanf(argv[5], "%f", &gradient_threshold);
//	sscanf(argv[6], "%d", &iteration_limit);
//	sscanf(argv[7], "%d", &need_ASC_1);
//	sscanf(argv[8], "%d", &need_DSC_1);
//
//	sscanf(argv[9], "%f", &g_pre_simp_threshold);
//	if (argc >= 10)
//		sscanf(argv[10], "%d", &parallelism);
//
//	// set remaining values
//	if (parallelism != -1) {
//		omp_set_num_threads(parallelism);
//	}
//
//	printf("dims=(%d,%d)\nfile=%s\nintegration parameters: e=%f, gt=%f, il=%d\nondemandacc: a1=%d, d1=%d, ps=%f\npar=%d\n",
//		X, Y, argv[3], error_threshold, gradient_threshold, iteration_limit, need_ASC_1, need_DSC_1,  g_pre_simp_threshold, parallelism);
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
//GridType* g_grid;
//GridFuncType* g_rgt_func;
//MeshType *g_topo_grid;
//
//#ifdef USE_REGION_CLEANER
//VertexLabelingToBoundaryLabeling<INDEX_TYPE>* g_edge_map;
//#else
//VertexLabelingToBoundaryLabeling<int, MaxVLType>* g_edge_map;
//#endif
//
//GradType *base_grad;
////RobinsLabelingAlgorithm<MeshType, MeshFuncType> *g_robin_alg;
////TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>* g_topo_alg;
//StreamlineIntegratorTypeASC* g_digitizing_streamline_integrator_asc;
//StreamlineIntegratorTypeDSC* g_digitizing_streamline_integrator_dsc;
//MaxVLType* g_maxv_labeling;
//MeshFuncType* g_topo_func;
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
//#define OUTPUT_DEBUG_2LINES
//
//void ReIntegrateUpFromSaddles(GradType* base_grad) {
//
//
//#ifdef OUTPUT_DEBUG_2LINES
//	FILE* fout = fopen("test_uplines.bin", "wb");
//	int countwrite = 0;
//	int writelength = 0;
//#endif
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
//			int counter = 10;
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
//				if (g_topo_grid->dimension(arc_hex_id) != 2) continue;
//				std::vector<Vec2d> points;
//				std::vector<INDEX_TYPE> dline;
//				Vec2d seed; g_topo_grid->centroid(arc_hex_id, seed);
//				seed = seed * 0.5; // back to grid coordinates
//				g_digitizing_streamline_integrator_asc->IntegrateStreamline(seed, points, dline);
//#ifdef OUTPUT_DEBUG_2LINES
//#pragma omp critical
//				{
//					int num = points.size();
//					fwrite(&num, sizeof(int), 1, fout);
//					fwrite(points.data(), sizeof(Vec2d), points.size(), fout);
//					countwrite++; writelength += points.size();
//				}
//#endif
//
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
//#ifdef OUTPUT_DEBUG_2LINES
//	printf("countwrite = %d, length = %d\n", countwrite, writelength);
//	fclose(fout);
//#endif
//
//}
//
//void ReIntegrateDownFromSaddles(GradType* base_grad) {
//#ifdef OUTPUT_DEBUG_2LINES
//	FILE* fout = fopen("test_dnlines.bin", "wb");
//#endif
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
//			int counter = 10;
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
//				std::vector<Vec2d> points;
//				std::vector<INDEX_TYPE> dline;
//				Vec2d seed; g_topo_grid->centroid(arc_vert_id, seed);
//				seed = seed * 0.5; // back to grid coordinates
//				g_digitizing_streamline_integrator_dsc->IntegrateStreamline(seed, points, dline);
//#ifdef OUTPUT_DEBUG_2LINES
//#pragma omp critical
//				{
//					int num = points.size();
//					fwrite(&num, sizeof(int), 1, fout);
//					fwrite(points.data(), sizeof(Vec2d), points.size(), fout);
//				}
//#endif
//			}
//			//printf("\n");
//			g_digitizing_streamline_integrator_dsc->set_label(sad_id);
//
//
//
//		}
//#ifdef OUTPUT_DEBUG_2LINES
//
//
//
//		for (auto p : criticals) {
//			INDEX_TYPE id = p.second;
//			g_digitizing_streamline_integrator_dsc->get_output()->SetLabel(id, 4);
//		}
//
//		fclose(fout);
//#endif
//	}
//
//
//}
//
//
//
//
//
//void RecordGrad(GradType* grad, const char* gradname) {
//	printf("setting dim asc man\n");
//	//g_topo_alg->setAscendingManifoldDimensions();
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
//	printf("ondemand accuracy: a1=%d d1=%d\n", need_ASC_1, need_DSC_1);
//
//	// start timing overall algorithm
//	ThreadedTimer timer(1);
//	timer.StartGlobal();
//
//
//	char gradname[1024];
//	sprintf(gradname, "%s.grad", argv[3]);
//
//
//	// will write timing to this file
//	char timingname[2048];
//	sprintf(timingname, "%s.%03d.gtime.txt", argv[3], parallelism);
//	gtiming = fopen(timingname, "w");
//
//	// START IO ---------------------------
//	//StartTask();
//	g_grid = new GridType(Vec2i(X, Y), Vec2b(per_x, per_y));
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
//	if (!(need_ASC_1 || need_DSC_1)) {
//		//g_topo_alg = new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
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
//		int labeltarget = 1;
//		g_digitizing_streamline_integrator_asc = new StreamlineIntegratorTypeASC(g_grid, g_rgt_func, g_topo_grid, error_threshold, gradient_threshold, iteration_limit);
//		g_digitizing_streamline_integrator_asc->SetDigitizingTarget(restriction_labels, labeltarget);
//		ReIntegrateUpFromSaddles(base_grad);
//		EndTaskAndRecord("Numerical SaddleIntASC");
//	}
//	
//	if (need_DSC_1) {
//		StartTask();
//		int labeltarget = 2;
//		g_digitizing_streamline_integrator_dsc = new StreamlineIntegratorTypeDSC(g_grid, g_rgt_func, g_topo_grid, error_threshold, gradient_threshold, iteration_limit);
//		g_digitizing_streamline_integrator_dsc->SetDigitizingTarget(restriction_labels, labeltarget);
//		ReIntegrateDownFromSaddles(base_grad);
//		EndTaskAndRecord("Numerical SaddleIntDSC");
//	}
//
//	restriction_labels->OutputToFile("labs.bin");
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
//	//g_topo_alg = new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
//	//printf("after Second robins:\n");
//	//g_topo_alg->count_critical_points(4);
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
