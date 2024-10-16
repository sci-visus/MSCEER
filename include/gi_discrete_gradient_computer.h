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
#include "gi_modified_robins.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_extrema_region_builder.h"
#include "gi_numeric_integrator_path_compressing.h"
#include "gi_fast_robins_noalloc.h"

namespace GInt {

	namespace Accurate2D {
		
		// the classes we use for 2d ondemand accurate gradient construction
		typedef RegularGrid2D GridType;
		typedef RegularGridBilinearFunction GridFuncType;
		typedef TopologicalRegularGrid2D MeshType;
		typedef IndexCompareLessThan<GridFuncType> ComparerASC;
		typedef IndexCompareGreaterThan<GridFuncType> ComparerDSC;
		typedef DigitizingNumericStreamlineIntegrator2dASC<MeshType, GridFuncType, AdaptiveEulerAdvector2D<GridFuncType, 1> > StreamlineIntegratorTypeASC;
		typedef DigitizingNumericStreamlineIntegrator2dDSC<MeshType, GridFuncType, AdaptiveEulerAdvector2D<GridFuncType, -1> > StreamlineIntegratorTypeDSC;
		typedef DiscreteGradientLabeling<MeshType> GradType;
		typedef RegularGridMaxMinVertexLabeling2D<MeshType, GridFuncType> MaxVLType;
		typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType, 4, 6> RobinsType;
		typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> MeshFuncType;
		
		class DiscreteGradientBuilder {
		protected:
			// parameters for building gradient
			int m_X, m_Y;
			int m_iteration_limit;
			int m_per_x, m_per_y;
			float m_error_threshold, m_gradient_threshold;
			std::string m_raw_filename;
			int m_parallelism; // = -1;
			int m_need_ASC_1;// = false;
			int m_need_DSC_1;// = false;		
			bool m_verbose;

			// for timing
			std::chrono::steady_clock::time_point m_start_time;

			struct TimedTask {
				std::chrono::steady_clock::time_point start_time;
				std::chrono::steady_clock::time_point task_start_time;
				std::chrono::steady_clock::time_point end_time;
				bool m_verbose;
				std::string taskname;
				TimedTask(std::chrono::steady_clock::time_point start, bool verbose = false) {
					start_time = start;
					m_verbose = verbose;
				}

				std::chrono::steady_clock::time_point StartTask(std::string name = "") {
					task_start_time = std::chrono::steady_clock::now();
					taskname = name;
					printf("starting task %s\n", name.c_str());
					return task_start_time;
				}
				void EndTask(std::string name, bool verbose_override = false) {
					end_time = std::chrono::steady_clock::now();
					taskname = name;
					if (verbose_override || m_verbose) {
						// format: [global activity name] [task] [start] [end] [dration]
						printf("TIMING: %s started: %d cumulative_elapsed: %d task_time: %d\n", taskname.c_str(),
							std::chrono::duration_cast<std::chrono::milliseconds>(task_start_time - start_time).count(),
							std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count(),
							std::chrono::duration_cast<std::chrono::milliseconds>(end_time - task_start_time).count());
					}
				}
			};

			std::vector<TimedTask> m_timing_seqence;

			// internal data structures
			GridType* g_grid;
			GridFuncType* g_rgt_func;
			MeshType* g_topo_grid;
			VertexLabelingToBoundaryLabeling<int, MaxVLType>* g_edge_map;
			GradType* base_grad;
			TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>* g_topo_alg;
			StreamlineIntegratorTypeASC* g_digitizing_streamline_integrator_asc;
			StreamlineIntegratorTypeDSC* g_digitizing_streamline_integrator_dsc;
			MaxVLType* g_maxv_labeling;
			MeshFuncType* g_topo_func;

			// meat of accuracy in 2d is tracing lines
			void ReIntegrateUpFromSaddles(GradType* base_grad) {


#ifdef OUTPUT_DEBUG_2LINES
				FILE* fout = fopen("test_uplines.bin", "wb");
				int countwrite = 0;
				int writelength = 0;
#endif

				// FIRST gather all the critical 2-saddles from the discrete gradient
				auto task0 = TimedTask(m_start_time, true); task0.StartTask("Topological GatherCritical2Saddles");
				std::vector<INDEX_TYPE> topo_index_partition;
				int num_threads;
				std::vector<std::pair<float, INDEX_TYPE> > criticals;
	
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
				task0.EndTask("Topological GatherCritical2Saddles");

				auto task1 = TimedTask(m_start_time, true); task1.StartTask("Topological Sort2Saddles");
				std::sort(criticals.begin(), criticals.end());
				task1.EndTask("Topological Sort2Saddles");
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
				auto task2 = TimedTask(m_start_time, true); task2.StartTask("Topological GatherCritical1Saddles");
				std::vector<INDEX_TYPE> topo_index_partition;
				int num_threads;
				std::vector<std::pair<float, INDEX_TYPE> > criticals;
				//printf("gothere 2\n");
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
				//printf("gothere 3\n");
				task2.EndTask("Topological GatherCritical1Saddles");

				auto task3 = TimedTask(m_start_time, true); task3.StartTask("Topological Sort2Saddles");
				std::sort(criticals.begin(), criticals.end());
				task3.EndTask("Topological Sort2Saddles");
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


			void RecordGrad(GradType* grad, const char* gradname) {
				printf("outputting to file %s\n", gradname);
				grad->outputToFile(gradname);
				//return 1;
			}

			bool m_needs_grid;
			bool m_needs_gfunc;
			bool m_needs_mesh;

		public:
			DiscreteGradientBuilder() {
				// set default options
				m_X = 0; m_Y = 0;
				m_iteration_limit = 500;
				m_per_x = 0; m_per_y = 0;
				m_error_threshold = 0.0001; m_gradient_threshold = 0.0;
				m_parallelism = -1; // = -1;
				m_raw_filename = "None";
				m_need_ASC_1= false;
				m_need_DSC_1= false;
				m_verbose = false;

				m_needs_grid = true;
				m_needs_gfunc = true;
				m_needs_mesh = true;
			}

			void SetGridAndFunc(GridType* grid, GridFuncType* func) {
				m_needs_grid = false;
				m_needs_gfunc = false;
				g_grid = grid;
				g_rgt_func = func;
				m_X = g_grid->XY()[0];
				m_Y = g_grid->XY()[1];
			}

			void SetTopoMesh(MeshType* mesh) {
				m_needs_mesh = false;
				g_topo_grid = mesh;
			}

			void SetFloadArrayAndDims(int x, int y, float* values) {
				m_needs_grid = false;
				m_needs_gfunc = false;
				g_grid = new GridType(Vec2i(x, y), Vec2b(m_per_x, m_per_y));
				g_rgt_func = new GridFuncType(g_grid, values);
				m_X = g_grid->XY()[0];
				m_Y = g_grid->XY()[1];
			}

			void SetRawFileNameAndDims(int x, int y, std::string rawname) {
				m_raw_filename = rawname;
				m_X = x;
				m_Y = y;
			}
			
			void SetNeededAccuracy(bool ascending, bool descending) {
				m_need_ASC_1 = ascending;
				m_need_DSC_1 = descending;
			}

			void SetParallelism(int threads) {
				m_parallelism = threads;
			}

			MaxVLType* GetMaxVertLabeleing() {
				return g_maxv_labeling;
			}

			MeshFuncType* GetMeshFunc() {
				return g_topo_func;
			}
			
			MeshType* GetTopoMesh() {
				return g_topo_grid;
			}

			GridType* GetGrid() {
				return g_grid;
			}

			GridFuncType* GetGridFunc() {
				return g_rgt_func;
			}

			GradType* GetGrad() {
				return base_grad;
			}

			void ComputeDiscreteGradient() {
				
				m_start_time = std::chrono::steady_clock::now();

				if (m_needs_grid)
					g_grid = new GridType(Vec2i(m_X, m_Y), Vec2b(m_per_x, m_per_y));
				if (m_needs_gfunc) {
					g_rgt_func = new GridFuncType(g_grid);
					if (m_raw_filename == "None") {
						printf("DiscreteGradientBuilder::ComputeDiscreteGradient() needs either filename or existing grid and function!\n");
						return;
					}
					g_rgt_func->LoadImageFromFile(m_raw_filename.c_str());
				}
				


				// CREATE MAX VERTEX LABELING
				auto task4 = TimedTask(m_start_time, true); task4.StartTask("Topological MaxVlabel");
				g_topo_grid = new MeshType(g_grid);

				g_maxv_labeling = new MaxVLType(g_topo_grid, g_rgt_func);
				g_maxv_labeling->ComputeOutput();
				task4.EndTask("Topological MaxVlabel");
				// create a topology function
				auto task5 = TimedTask(m_start_time, true); task5.StartTask("Topological Function");
				//printf("topo function...\n");
				g_topo_func = new MeshFuncType();
				g_topo_func->setMeshAndFuncAndMaxVLabeling(g_topo_grid, g_rgt_func, g_maxv_labeling);
				task5.EndTask("Topological Function");
				auto task6 = TimedTask(m_start_time, true); task6.StartTask("Topological Robins0");
				base_grad = new GradType(g_topo_grid);
				base_grad->ClearAllGradient();
				RobinsType* first_robins = new RobinsType(g_topo_grid, g_maxv_labeling, base_grad);
				first_robins->ComputePairing();
				task6.EndTask("Topological Robins0");

				g_topo_alg = new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
				//printf("after base first robins:\n");
				//g_topo_alg->count_critical_points(3);

				if (!(m_need_ASC_1 || m_need_DSC_1)) {
					//g_topo_alg = new TopologicalGradientUsingAlgorithms<MeshType, MeshFuncType, GradType>(g_topo_func, g_topo_grid, base_grad);
					printf("no accuraccy needed, outputting\n");

				}
				else {
					// Do gradient vectors computation from raw image data
					//printf("computing gradient\n");
					auto task7 = TimedTask(m_start_time, true); task7.StartTask("NumericalTracing GradCompute");
					g_rgt_func->ComputeGradFromImage(1);
					//g_rgt_func->Negate();
					task7.EndTask("NumericalTracing GradCompute");

					// CREATE RESTRICTION MAP - WILL NEED
					DenseLabeling<char>* restriction_labels = new DenseLabeling<char>(g_topo_grid->numCells());
					restriction_labels->SetAll(0);
					// we will always need a constrained robins alg
					RobinsType* constrained_robins = new RobinsType(g_topo_grid, g_maxv_labeling, restriction_labels, base_grad);

					//-------------------------------------------------------------
					//-------------------------------------------------------------
					//-------------------------------------------------------------
					// DO LINES ACCURATE GRADIENT COMPUTATION 
					//-------------------------------------------------------------
					//-------------------------------------------------------------
					//-------------------------------------------------------------
					//-------------------------------------------------------------

					if (m_need_ASC_1) {
						auto task8 = TimedTask(m_start_time, true); task8.StartTask("Numerical SaddleIntASC");
						int labeltarget = 1;
						g_digitizing_streamline_integrator_asc = new StreamlineIntegratorTypeASC(g_grid, g_rgt_func, g_topo_grid, m_error_threshold, m_gradient_threshold, m_iteration_limit);
						g_digitizing_streamline_integrator_asc->SetDigitizingTarget(restriction_labels, labeltarget);
						ReIntegrateUpFromSaddles(base_grad);
						task8.EndTask("Numerical SaddleIntASC");
						delete g_digitizing_streamline_integrator_asc;
					}

					if (m_need_DSC_1) {
						auto task9 = TimedTask(m_start_time, true); task9.StartTask("Numerical SaddleIntDSC");
						int labeltarget = 2;
						g_digitizing_streamline_integrator_dsc = new StreamlineIntegratorTypeDSC(g_grid, g_rgt_func, g_topo_grid, m_error_threshold, m_gradient_threshold, m_iteration_limit);
						g_digitizing_streamline_integrator_dsc->SetDigitizingTarget(restriction_labels, labeltarget);
						ReIntegrateDownFromSaddles(base_grad);
						task9.EndTask("Numerical SaddleIntDSC");
						delete g_digitizing_streamline_integrator_dsc;
					}
#ifdef OUTPUT_DEBUG_2LINES
					restriction_labels->OutputToFile("labs.bin");
#endif

					//-------------------------------------------------------------
					// NOW REDO DISCRETE GRADIENT IN NEIGHBORHOOD
					//-------------------------------------------------------------
					auto task10 = TimedTask(m_start_time, true); task10.StartTask("Topological ConformingGrad");
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
					task10.EndTask("Topological ConformingGrad");
					delete restriction_labels;
				}
				auto task11 = TimedTask(m_start_time, true); task11.StartTask("Set Dim of Asc Manifolds");
				g_topo_alg->setAscendingManifoldDimensions();
				task11.EndTask("Set Dim of Asc Manifolds");

			}

			void WriteGrad(const char* filename) {
				RecordGrad(base_grad, filename);
			}

		};



	};






};