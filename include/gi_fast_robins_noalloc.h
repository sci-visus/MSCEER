#ifndef BAST_ROBINS_NOALLOC_H
#define BAST_ROBINS_NOALLOC_H

#include "gi_labeling.h"
#include "gi_discrete_gradient_labeling.h"
#include "gi_topological_simplicial_complex.h"
#include "gi_bifiltration_pairing.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_max_vertex_labeling.h"


//#define SANITY_CHECKS

namespace GInt {

	//#define DEBUGPARALLEL
	//#define DEBUGPARALLEL

	template<class GridType, class GridFuncType, class MeshType, class MaxVLType, class GradType>
	class SlidingWindowRobinsNoalloc
	{
	protected:

		GridType* mGrid;
		GridFuncType* mFunc;
		MeshType* mMesh;
		MaxVLType* mMaxVL;
		DenseLabeling<char>* mResLabel;
		GradType* mGrad;
		MyRobinsNoalloc<MeshType, MaxVLType, GradType, 5, 4>* mStandardRobins;
		INDEX_TYPE m_data_27_offsets[27];

#ifdef DEBUGPARALLEL
		int* db_counter;
#endif

		struct MESH_CONTEXT {
			Explicit3x3x3SmallRegularGrid* small_mesh;
			RegularGrid3D* small_grid;
			RegularGridTrilinearFunction* small_grid_func;
			RegularGridMaxMinVertexLabeling3D<Explicit3x3x3SmallRegularGrid, RegularGridTrilinearFunction>* small_mesh_maxmin_labeling;
			INDEX_TYPE small_mesh_id_to_big_mesh_id[125];
			FLOATTYPE small_grid_values[27];
		};

		struct small_INDEX_vector {
			INDEX_TYPE vec[27];
			int size;

			void push_back(INDEX_TYPE val) {
				vec[size] = val;
				++size;
			}

			small_INDEX_vector() : size(0) {}
			const INDEX_TYPE& operator[](int i) const { return vec[i]; }
			INDEX_TYPE& operator[](int i) { return vec[i]; }

		};

		struct cell_pairing {
			int num_missing;
			INDEX_TYPE pair;
			bool paired;
			cell_pairing() {
				//printf("SHOULD NEVER CALL\n");
				pair = -1;
				paired = false;
				num_missing = 0;
			}

		};

		struct myStaticMap {
			INDEX_TYPE sm_ids_in_lstar[27];
			int in_lstar[125];
			cell_pairing cell_pairings[125];
			int size;

			void push_back(INDEX_TYPE small_mesh_id) {
				sm_ids_in_lstar[size] = small_mesh_id;
				in_lstar[small_mesh_id] = 1;
				++size;
			}

			// initialize all in_star to 0
			myStaticMap() : size(0), in_lstar{} {}

			int is_in_lstar(INDEX_TYPE id) const {
				return in_lstar[id] != 0;
			}

		};

		void printStaticMapState(myStaticMap& m) {
			for (int i = 0; i < m.size; i++) {
				INDEX_TYPE id = m.sm_ids_in_lstar[i];
				int inlstar = m.in_lstar[id];
				cell_pairing& cp = m.cell_pairings[id];
				printf("%lld -> %lld, nm=%d, inlst=%d, pinlst=%d, p=%d, ref=%d\n",
					id, cp.pair, cp.num_missing, inlstar, m.in_lstar[cp.pair], cp.paired,
					m.cell_pairings[cp.pair].pair == id);
			}
		}


		//std::queue<INDEX_TYPE> readytogo;

		void decrementCofacets(const MESH_CONTEXT& mc, INDEX_TYPE id, myStaticMap& small_mesh_cell_pairings) const {
			typename Explicit3x3x3SmallRegularGrid::CofacetsIterator cfit(mc.small_mesh);
			for (cfit.begin(id); cfit.valid(); cfit.advance()) {
				INDEX_TYPE cid = cfit.value();
				if (small_mesh_cell_pairings.is_in_lstar(cid))
					small_mesh_cell_pairings.cell_pairings[cid].num_missing--;
			}
		}

		bool is_steeper(const MESH_CONTEXT& mc, INDEX_TYPE vn_1, INDEX_TYPE vn_2) const {
			return mc.small_mesh_maxmin_labeling->Before(vn_1, vn_2);
		}

	public:
		INDEX_TYPE lowest_vertex(const MESH_CONTEXT& mc, INDEX_TYPE cid) const {
			return mc.small_mesh->VertexNumberFromCellID(mc.small_mesh_maxmin_labeling->Cell2LowestVertex(cid));
		}
	protected:


		INDEX_TYPE PickLowestCandidate(const MESH_CONTEXT& mc, small_INDEX_vector& cands, myStaticMap& small_mesh_cell_pairings) const {
			if (cands.size == 1) return cands[0];

			INDEX_TYPE curr_lowest_id = cands[0];
			INDEX_TYPE lv_vid = lowest_vertex(mc, curr_lowest_id);
			for (int i = 1; i < cands.size; i++) {
				INDEX_TYPE olv_vid = lowest_vertex(mc, cands[i]);
				if (mc.small_mesh_maxmin_labeling->Before(olv_vid, lv_vid)) {
					lv_vid = olv_vid;
					curr_lowest_id = cands[i];
				}
			}
			return curr_lowest_id;
		}



		void HomotopyExpand(const MESH_CONTEXT& mc, small_INDEX_vector& lstar_cell_sm_ids, myStaticMap& small_mesh_cell_pairings) const {

			// first push all small mesh ids into the cell pairings map
			small_INDEX_vector list_of_d_cells_sm_ids[4];
			for (int i = 0; i < lstar_cell_sm_ids.size; i++) {
				INDEX_TYPE small_mesh_id = lstar_cell_sm_ids[i];
				small_mesh_cell_pairings.push_back(small_mesh_id);
			}
			// count number of facets that are in the lower star
			for (int i = 0; i < small_mesh_cell_pairings.size; i++) {
				INDEX_TYPE small_mesh_cell_id = small_mesh_cell_pairings.sm_ids_in_lstar[i];

				list_of_d_cells_sm_ids[mc.small_mesh->dimension(small_mesh_cell_id)].push_back(small_mesh_cell_id);
				typename Explicit3x3x3SmallRegularGrid::FacetsIterator small_mesh_facets_iterator(mc.small_mesh);
				for (small_mesh_facets_iterator.begin(small_mesh_cell_id); small_mesh_facets_iterator.valid(); small_mesh_facets_iterator.advance()) {
					INDEX_TYPE small_mesh_facet_id = small_mesh_facets_iterator.value();
					if (small_mesh_cell_pairings.is_in_lstar(small_mesh_facet_id)) small_mesh_cell_pairings.cell_pairings[small_mesh_cell_id].num_missing++;
				}
			}


			// if there is a vertex we should pick steepest descent
			if (list_of_d_cells_sm_ids[0].size > 0) {
				// we should only have one
#ifdef SANITY_CHECKS
				if (list_of_d_cells_sm_ids[0].size > 1) {
					printf("ERROR: too many vertices %d\n", list_of_d_cells_sm_ids[0].size);
					printf("\n");
				}
#endif
				// the id of the vertex is simply the first element of the list
				INDEX_TYPE sm_vertex_id = list_of_d_cells_sm_ids[0][0];

				// if there are no edges, then make the vertex critical, else pair with an edge
				if (list_of_d_cells_sm_ids[1].size == 0) {
					// make vertex critical
					small_mesh_cell_pairings.cell_pairings[sm_vertex_id].pair = sm_vertex_id;
					small_mesh_cell_pairings.cell_pairings[sm_vertex_id].paired = true;
					decrementCofacets(mc, sm_vertex_id, small_mesh_cell_pairings);
				}
				else {
					// to pair with the steepest down edge, we want to look through the list
					INDEX_TYPE sm_lowest_edge_id;
					// if there is only one edge, it's easy, pick that!
					if (list_of_d_cells_sm_ids[1].size == 1) {
						// just pair with only option 
						sm_lowest_edge_id = list_of_d_cells_sm_ids[1][0];
					}
					else {
						// find minimal edge
						sm_lowest_edge_id = list_of_d_cells_sm_ids[1][0]; // set to first
						INDEX_TYPE temp_lowest_vertex_vn =
							mc.small_mesh->VertexNumberFromCellID(mc.small_mesh_maxmin_labeling->Cell2LowestVertex(sm_lowest_edge_id));

#ifdef SANITY_CHECKS
						if (temp_lowest_vertex_vn == sm_vertex_id) {
							printf("ERROR: how the heck can the lowest vertex of an edge be its lstar thingy\n");
						}
#endif
						for (int i = 1; i < list_of_d_cells_sm_ids[1].size; i++) {
							INDEX_TYPE other_edge_id = list_of_d_cells_sm_ids[1][i];
							INDEX_TYPE temp_other_vertex_vn =
								mc.small_mesh->VertexNumberFromCellID(mc.small_mesh_maxmin_labeling->Cell2LowestVertex(other_edge_id));

							if (is_steeper(mc, temp_other_vertex_vn, temp_lowest_vertex_vn)) {
								sm_lowest_edge_id = other_edge_id;
								temp_lowest_vertex_vn = temp_other_vertex_vn;
							}
						}


					}

					// pair in direction of steepest descent
					small_mesh_cell_pairings.cell_pairings[sm_vertex_id].pair = sm_lowest_edge_id;
					small_mesh_cell_pairings.cell_pairings[sm_vertex_id].paired = true;
					small_mesh_cell_pairings.cell_pairings[sm_lowest_edge_id].pair = sm_vertex_id;
					small_mesh_cell_pairings.cell_pairings[sm_lowest_edge_id].paired = true;
					decrementCofacets(mc, sm_vertex_id, small_mesh_cell_pairings);
					decrementCofacets(mc, sm_lowest_edge_id, small_mesh_cell_pairings);

				}
			}



			for (int i = 0; i < 4; i++) {
				//while (!sorted.empty()) {



				// logic is we need to process every cell of dimension i
				// until all have been processed, first try to pair
				// if no pairing was successful, make one critical and repeat
				int num_processed = 0;
				int total_to_process = 0;
				for (int j = 0; j < list_of_d_cells_sm_ids[i].size; j++) {
					INDEX_TYPE i_cell_id = list_of_d_cells_sm_ids[i][j];
					if (!small_mesh_cell_pairings.cell_pairings[i_cell_id].paired) total_to_process++;
				}
				while (num_processed < total_to_process) {

					int start_num_proc = num_processed;
					// try to pair as many as possible
					for (int j = 0; j < list_of_d_cells_sm_ids[i].size; j++) {
						INDEX_TYPE i_cell_id = list_of_d_cells_sm_ids[i][j];
						if (small_mesh_cell_pairings.cell_pairings[i_cell_id].paired) continue; // already paired
#ifdef DEBUG_PARALLEL
						if (small_mesh_cell_pairings.cell_pairings[i_cell_id].num_missing > 0) {
							printf("ERROR: should never get here1\n");
						}
#endif

						small_INDEX_vector candidates;
						typename Explicit3x3x3SmallRegularGrid::CofacetsIterator cfit(mc.small_mesh);
						for (cfit.begin(i_cell_id); cfit.valid(); cfit.advance()) {
							INDEX_TYPE cfid = cfit.value();

							if (!small_mesh_cell_pairings.is_in_lstar(cfid)) continue; // not in our lower star
#ifdef DEBUG_PARALLEL
							if (small_mesh_cell_pairings.cell_pairings[cfid].paired) {
								printf("ERROR: should never get here2\n");
							}
#endif
							if (small_mesh_cell_pairings.cell_pairings[cfid].num_missing == 1) {
								// pair lstar_cell_sm_ids
								candidates.push_back(cfid);

							}
						}

						//if (candidates.size() > 1) printf("got here candidates: %d\n", candidates.size());
						if (candidates.size > 0) {
							INDEX_TYPE cfid = PickLowestCandidate(mc, candidates, small_mesh_cell_pairings);
							small_mesh_cell_pairings.cell_pairings[i_cell_id].pair = cfid;
							small_mesh_cell_pairings.cell_pairings[i_cell_id].paired = true;
							small_mesh_cell_pairings.cell_pairings[cfid].pair = i_cell_id;
							small_mesh_cell_pairings.cell_pairings[cfid].paired = true;
							decrementCofacets(mc, i_cell_id, small_mesh_cell_pairings);
							decrementCofacets(mc, cfid, small_mesh_cell_pairings);

							num_processed++;
							break;
						}
					}

					if (start_num_proc == num_processed) {
						// then no more pairs were possible 
						small_INDEX_vector candidates;
						for (int j = 0; j < list_of_d_cells_sm_ids[i].size; j++) {
							INDEX_TYPE i_cell_id = list_of_d_cells_sm_ids[i][j];
							if (small_mesh_cell_pairings.cell_pairings[i_cell_id].paired) continue; // already paired
							// make one critical and break
							candidates.push_back(i_cell_id);
						}

						INDEX_TYPE id = PickLowestCandidate(mc, candidates, small_mesh_cell_pairings);
						//asdf want to make lowest critical!?
						small_mesh_cell_pairings.cell_pairings[id].pair = id;
						small_mesh_cell_pairings.cell_pairings[id].paired = true;
						decrementCofacets(mc, id, small_mesh_cell_pairings);
						num_processed++;

					}
				}

			}
			//printf("out:");
			//for (auto c : small_mesh_cell_pairings) {
			//	if (c.second.pair == c.second.id) printf(" (%d:%llu)", mc.small_mesh->dimension(c.second.id), c.second.pair);
			//	if (mc.small_mesh->dimension(c.second.id) < mc.small_mesh->dimension(c.second.pair)) 
			//		printf(" (%d:%llu->%d:%llu)", mc.small_mesh->dimension(c.second.id), c.second.id, mc.small_mesh->dimension(c.second.pair), c.second.pair);
			//}
			//printf("\n");
		}

		void init() {
			mStandardRobins = new MyRobinsNoalloc<MeshType, MaxVLType, GradType, 5, 4>(mMesh, mMaxVL, NULL, mGrad);
			// get data offsets for a 27 neighborhood
			INDEX_TYPE t_did111 = mGrid->Index3d(Vec3l(1, 1, 1));
			int t_pos = 0;
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						m_data_27_offsets[t_pos++] = mGrid->Index3d(Vec3l(i + 1, j + 1, k + 1)) - t_did111;
					}
				}
			}
		}
	public:


		SlidingWindowRobinsNoalloc(
			GridType* grid,
			GridFuncType* grid_func,
			MeshType* mesh,
			MaxVLType* label1,
			GradType* grad) : mGrid(grid), mFunc(grid_func), mMesh(mesh), mMaxVL(label1), mResLabel(NULL), mGrad(grad) {
			init();
		}
		~SlidingWindowRobinsNoalloc() {
			delete mStandardRobins;
		}



		void ComputeLowerStar(const MESH_CONTEXT& mc, INDEX_TYPE small_mesh_vertex_id) {
			small_INDEX_vector lower_star_list[1];
			// now add all lower star to restriciton sets 
			typename Explicit3x3x3SmallRegularGrid::AdjacentCellsIterator star(mc.small_mesh);
			for (star.begin(small_mesh_vertex_id); star.valid(); star.advance()) {
				INDEX_TYPE small_mesh_vertex_neighbor = star.value();

				// discard a cell if its highest vertex is NOT the vertex, hence not part of lower star
				INDEX_TYPE highest_small_mesh_vertex_id = mc.small_mesh_maxmin_labeling->Cell2HighestVertex(small_mesh_vertex_neighbor);
				if (highest_small_mesh_vertex_id != small_mesh_vertex_id) continue; // not in lower star of f1
				lower_star_list[0].push_back(small_mesh_vertex_neighbor);
#ifdef DEBUGPARALLEL
#pragma omp critical
				{
					if (omp_get_thread_num() > 0) printf("here\n");
					db_counter[mc.small_mesh_id_to_big_mesh_id[small_mesh_vertex_neighbor]]++;
				}
#endif

			}

			// now do homotopy expand on each subset!
			// do homotopy expand
			myStaticMap small_mesh_cell_pairings;
			HomotopyExpand(mc, lower_star_list[0], small_mesh_cell_pairings);
#ifdef SANITY_CHECKS
			printStaticMapState(small_mesh_cell_pairings);
			printf("\n");
#endif
			for (int j = 0; j < small_mesh_cell_pairings.size; j++) {
				INDEX_TYPE small_mesh_id = small_mesh_cell_pairings.sm_ids_in_lstar[j];
				INDEX_TYPE small_mesh_id_pair = small_mesh_cell_pairings.cell_pairings[small_mesh_id].pair;

				INDEX_TYPE big_mesh_id = mc.small_mesh_id_to_big_mesh_id[small_mesh_id];
				INDEX_TYPE big_mesh_id_pair = mc.small_mesh_id_to_big_mesh_id[small_mesh_id_pair];
				// NOW GO BACK TO GLOBAL?
				if (big_mesh_id == big_mesh_id_pair) {
					mGrad->setCritical(big_mesh_id, true);
					mGrad->setAssigned(big_mesh_id, 1);
#ifdef DEBUGPARALLEL2
#pragma omp critical
					{
						db_counter[big_mesh_id]++;
					}
#endif
				}
				else {

					// SANITY CHECKS
#ifdef SANITY_CHECKS
					// check small mesh sanity
					if (mc.small_mesh_maxmin_labeling->Cell2HighestVertex(small_mesh_id) !=
						mc.small_mesh_maxmin_labeling->Cell2HighestVertex(small_mesh_id_pair)) {
						printf("whoathere\n");
					}
					if (mMaxVL->Cell2HighestVertex(big_mesh_id) !=
						mMaxVL->Cell2HighestVertex(big_mesh_id_pair)) {
						printf("whoathere\n");
					}
#endif
					// END SANITY CHECKS

					mGrad->setPair(big_mesh_id, big_mesh_id_pair);
					mGrad->setPair(big_mesh_id_pair, big_mesh_id);
					mGrad->setAssigned(big_mesh_id, 1);
					mGrad->setAssigned(big_mesh_id_pair, 1);
#ifdef DEBUGPARALLEL2
#pragma omp critical
					{
						db_counter[big_mesh_id]++;
						db_counter[big_mesh_id_pair]++;
					}
#endif
				}
			}
		}



		void ComputePairing() {

			std::chrono::steady_clock::time_point now_time = std::chrono::steady_clock::now();
			std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

			// dimensions of the mesh
			Vec3l big_mesh_xyz = mMesh->XYZ();
			int lstars_count = 0;

#ifdef DEBUGPARALLEL
			db_counter = new int[mMesh->numCells()];
			memset(db_counter, 0, sizeof(int) * mMesh->numCells());
#endif


			// START PARALLEL WORK
#pragma omp parallel
			{
				std::vector<INDEX_TYPE> topo_index_partition;
				int num_threads;
				num_threads = omp_get_num_threads();
				ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, topo_index_partition);
				int thread_num = omp_get_thread_num();

				// these coordinates are INCLUSIVE - which means do start and end
				INDEX_TYPE thread_start_id = topo_index_partition[thread_num];
				INDEX_TYPE thread_end_id = topo_index_partition[thread_num + 1] - 1;
				Vec3l start_coord, end_coord;
				mMesh->cellid2Coords(thread_start_id, start_coord);
				mMesh->cellid2Coords(thread_end_id, end_coord); // get inclusive coord
				//#pragma omp critical 
				//				{
				//					printf("thread %d doing:\n\t", thread_num);
				//					start_coord.PrintInt(); printf("\t");
				//					end_coord.PrintInt();
				//				}
				// iterate over all vertices
				MESH_CONTEXT mc; // gather the pointers rather than have to pass a million items
				// only need maxvl labeling and function values
				// and maybe reslabel
				mc.small_mesh = new Explicit3x3x3SmallRegularGrid();
				mc.small_grid = new RegularGrid3D(Vec3l(3, 3, 3), Vec3b(0, 0, 0));


				mc.small_grid_func = new RegularGridTrilinearFunction(mc.small_grid, mc.small_grid_values); // wrapper for our values
				// place to store our local copy of the max/min vertices for each cell
				mc.small_mesh_maxmin_labeling = new RegularGridMaxMinVertexLabeling3D<Explicit3x3x3SmallRegularGrid, RegularGridTrilinearFunction>(mc.small_mesh, mc.small_grid_func);
				mc.small_mesh_maxmin_labeling->HACK_init();

				const INDEX_TYPE kernel_baseid = mc.small_mesh->coords2Cellid(Vec3l(2, 2, 2));
				const INDEX_TYPE kernel_data_baseid = mc.small_grid->Index3d(Vec3l(1, 1, 1));


				int kstart = start_coord[2];
				if (kstart % 2 == 1) kstart--; // kstart cannot start on an odd number, if it is odd, start on prior?
				if (kstart == 0) kstart = 2;
				int kend = end_coord[2];
				if (kend == big_mesh_xyz[2] - 1) kend = big_mesh_xyz[2] - 2;
				// NOW DO ALL INTERIOR VERTICES
#ifdef DEBUGPARALLEL
#pragma omp critical 
				{
					printf("thread %d doing actual k: [%d:%d]\n", thread_num, kstart, kend);
				}
#endif

				for (int k = kstart; k <= kend; k += 2) { // do parallel division of work
					const int d_k = k >> 1; // data k

					int jstart = 2;
					int jend = big_mesh_xyz[1] - 1;
					//if (k == kstart) {
					//	jstart = start_coord[1];
					//} 
					//if (k == kend) {
					//	jend = end_coord[1] - 2;
					//}

					for (int j = jstart; j < jend; j += 2) {
						const int d_j = j >> 1; // data j
						const INDEX_TYPE baseid_nox = mMesh->coords2Cellid(Vec3l(0, j, k));
						const INDEX_TYPE data_baseid_nox = mGrid->Index3d(Vec3l(0, d_j, d_k));

						int istart = 2;
						int iend = big_mesh_xyz[0] - 1;
						//if (k == kstart && j == start_coord[1]) {
						//	istart = start_coord[0];
						//}
						//if (k == kend && j == end_coord[1]) {
						//	iend = end_coord[0] - 2;
						//}

						for (int i = istart; i < iend; i += 2) {
							const int d_i = i >> 1; // data i
							const INDEX_TYPE baseid = baseid_nox + i;
							if (baseid < thread_start_id || baseid > thread_end_id) continue;
							const INDEX_TYPE data_baseid = data_baseid_nox + d_i;

#ifdef DEBUGPARALLEL
#pragma omp critical 
							{
								printf("thread %d doing actual %d,%d,%d\n", thread_num, i,j,k);
							}
#endif

#ifdef SANITY_CHECKS

							this->mStandardRobins->ComputeLowerStar(baseid);
							INDEX_TYPE pre_pair = mGrad->getPair(baseid);
							INDEX_TYPE pre_ppair = mGrad->getPair(pre_pair);

							BYTE_TYPE GRADS[27];

#endif
							// so for each vertex FIRST copy in the values
							// we can optimize this later to do less global lookups
							for (int pos = 0; pos < 27; pos++) {
								//int sd_nid = mc.small_mesh->get27NeighborOffset(pos) + kernel_baseid;
								INDEX_TYPE big_mesh_vertex_id = mMesh->get27NeighborOffset(pos) + baseid;
								INDEX_TYPE big_grid_vertex_data_id = m_data_27_offsets[pos] + data_baseid;
								INDEX_TYPE kernel_vertex_id = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos);
								//INDEX_TYPE kernel_data_nid = kernel_data_baseid + m_data_27_offsets[pos]; // this should just = pos??
#ifdef SANITY_CHECKS						
								if (this->mMaxVL->Cell2HighestVertex(big_mesh_vertex_id) == baseid) {
									GRADS[pos] = mGrad->getAsChar(big_mesh_vertex_id);
									mGrad->clearGrad(big_mesh_vertex_id);
								}
#endif
								mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id] = big_mesh_vertex_id;

								mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id, this->mMaxVL->GetUncompressedMaxVal(big_mesh_vertex_id));
								mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id, this->mMaxVL->GetUncompressedMinVal(big_mesh_vertex_id));

								mc.small_grid_values[pos] = this->mFunc->SampleImage(big_grid_vertex_data_id);
							}

							ComputeLowerStar(mc, kernel_baseid);

#ifdef SANITY_CHECKS						
							for (int pos = 0; pos < 27; pos++) {
								//int sd_nid = mc.small_mesh->get27NeighborOffset(pos) + kernel_baseid;
								INDEX_TYPE big_mesh_vertex_id = mMesh->get27NeighborOffset(pos) + baseid;
								//INDEX_TYPE kernel_data_nid = kernel_data_baseid + m_data_27_offsets[pos]; // this should just = pos??
								if (this->mMaxVL->Cell2HighestVertex(big_mesh_vertex_id) == baseid) {
									BYTE_TYPE comp = mGrad->getAsChar(big_mesh_vertex_id);
									if (comp != GRADS[pos]) {
										printf("Error %d != %d\n", comp, GRADS[pos]);
									}
								}
							}
#endif

#ifdef SANITY_CHECKS
							INDEX_TYPE post_pair = mGrad->getPair(baseid);
							lstars_count++;

#endif

						}
					}
				}
			}
			//printf("INTERIOR: new robins1 in  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			now_time = std::chrono::steady_clock::now();
			//lstars_count = 0;
			// DO Z Plane Boundaries
#pragma omp parallel for
			for (int j = 0; j < big_mesh_xyz[1]; j += 2) {
				for (auto k : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[2] - 1 })) { // do parallel division of work
					for (int i = 0; i < big_mesh_xyz[0]; i += 2) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d Z boundaries\n", lstars_count);
			//int tmp = lstars_count;
			//lstars_count = 0;

			// DO Y Plane Boundaries
#pragma omp parallel for
			for (int k = 2; k < big_mesh_xyz[2] - 2; k += 2) { // smaller range since we did k = 0 and k = xyz[2]-1
				for (auto j : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[1] - 1 })) {
					for (int i = 0; i < big_mesh_xyz[0]; i += 2) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d Y boundaries\n", lstars_count);
			//int tmp2 = lstars_count;
			//lstars_count = 0;
			// DO X Plane Boundaries
#pragma omp parallel for
			for (int k = 2; k < big_mesh_xyz[2] - 2; k += 2) { // smaller range since we did k = 0 and k = xyz[2]-1
				for (int j = 2; j < big_mesh_xyz[1] - 2; j += 2) { // again smaller range
					for (auto i : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[0] - 1 })) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d X boundaries\n", lstars_count);
			//printf("BOUNDARY: new robins1 in  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			now_time = std::chrono::steady_clock::now();

			//printf("new robins1 %d lower stars in  %dms\n", lstars_count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			//now_time = std::chrono::steady_clock::now();
			//lstars_count = 0;
			//for (int k = 2; k < big_mesh_xyz[2] - 3; k += 2) { // do parallel division of work
			//	for (int j = 2; j < big_mesh_xyz[1] - 3; j += 2) {
			//		for (int i = 2; i < big_mesh_xyz[0] - 3; i += 2) {
			//			const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));

			//			INDEX_TYPE pre_pair = mGrad->getPair(baseid);

			//			this->mStandardRobins->ComputeLowerStar(baseid);
			//			lstars_count++;
			//			INDEX_TYPE post_pair = mGrad->getPair(baseid);

			//			if (pre_pair != post_pair) {
			//				printf("asdasdf\n");
			//			}

			//		}
			//	}
			//}

			//printf("old robins1 %d lower stars in  %dms\n", lstars_count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());

#ifdef DEBUGPARALLEL
			FILE* fout = fopen("test_out.raw", "wb");
			fwrite(db_counter, sizeof(int), mMesh->numCells(), fout);
			fclose(fout);
			//printf("doing seen checks!\n");
			//for (INDEX_TYPE i = 0; i < mMesh->numCells(); i++) {
			//	
			//	if (mMesh->boundaryValue(i) == 0 && db_counter[i] != 1) {
			//		printf("index %lld seen %d times: ", i, db_counter[i]);
			//		Vec3l c;
			//		mMesh->cellid2Coords(i, c);
			//		c.PrintInt();
			//	}

			//}
			printf("done seen checks\n");


#endif


		}
		//
		//			std::vector<INDEX_TYPE> topo_index_partition;
		//			int num_threads;
		//#pragma omp parallel
		//			{
		//#pragma omp single
		//				{
		//					num_threads = omp_get_num_threads();
		//					ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, topo_index_partition);
		//				}
		//
		//				int thread_num = omp_get_thread_num();
		//
		//				// iterate over all vertices
		//				typename MeshType::DCellsIterator verts(mMesh, 0, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
		//				for (verts.begin(); verts.valid(); verts.advance()){
		//					INDEX_TYPE small_mesh_vertex_id = verts.value();
		//
		//					ComputeLowerStar(small_mesh_vertex_id);
		//
		//
		//				}
		//			}
		//		}
		//		//DenseLabeling<INDEX_TYPE>* GetLabeling() { return mPairs; }


		void ComputePairing_slidingSerial() {
			// dimensions of the mesh
			Vec3l big_mesh_xyz = mMesh->XYZ();
			//int lstars_count = 0;

			{
				std::vector<INDEX_TYPE> topo_index_partition;
				int num_threads = 1;
				ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, topo_index_partition);
				int thread_num = 0;

				// these coordinates are INCLUSIVE - which means do start and end
				INDEX_TYPE thread_start_id = topo_index_partition[thread_num];
				INDEX_TYPE thread_end_id = topo_index_partition[thread_num + 1] - 1;
				Vec3l start_coord, end_coord;
				mMesh->cellid2Coords(thread_start_id, start_coord);
				mMesh->cellid2Coords(thread_end_id, end_coord); // get inclusive coord
				// iterate over all vertices
				MESH_CONTEXT mc; // gather the pointers rather than have to pass a million items
				// only need maxvl labeling and function values
				// and maybe reslabel
				mc.small_mesh = new Explicit3x3x3SmallRegularGrid();
				mc.small_grid = new RegularGrid3D(Vec3l(3, 3, 3), Vec3b(0, 0, 0));


				mc.small_grid_func = new RegularGridTrilinearFunction;
				mc.small_grid_func->RegularGridTrilinearFunctionSerial(mc.small_grid, mc.small_grid_values); // wrapper for our values
				// place to store our local copy of the max/min vertices for each cell
				mc.small_mesh_maxmin_labeling = new RegularGridMaxMinVertexLabeling3D<Explicit3x3x3SmallRegularGrid, RegularGridTrilinearFunction>(mc.small_mesh, mc.small_grid_func);
				mc.small_mesh_maxmin_labeling->HACK_init();

				const INDEX_TYPE kernel_baseid = mc.small_mesh->coords2Cellid(Vec3l(2, 2, 2));
				const INDEX_TYPE kernel_data_baseid = mc.small_grid->Index3d(Vec3l(1, 1, 1));


				int kstart = start_coord[2];
				if (kstart % 2 == 1) kstart--; // kstart cannot start on an odd number, if it is odd, start on prior?
				if (kstart == 0) kstart = 2;
				int kend = end_coord[2];
				if (kend == big_mesh_xyz[2] - 1) kend = big_mesh_xyz[2] - 2;
				// NOW DO ALL INTERIOR VERTICES

				for (int k = kstart; k <= kend; k += 2) {
					const int d_k = k >> 1; // data k

					int jstart = 2;
					int jend = big_mesh_xyz[1] - 1;
					//if (k == kstart) {
					//	jstart = start_coord[1];
					//} 
					//if (k == kend) {
					//	jend = end_coord[1] - 2;
					//}

					for (int j = jstart; j < jend; j += 2) {
						const int d_j = j >> 1; // data j
						const INDEX_TYPE baseid_nox = mMesh->coords2Cellid(Vec3l(0, j, k));
						const INDEX_TYPE data_baseid_nox = mGrid->Index3d(Vec3l(0, d_j, d_k));

						int istart = 2;
						int iend = big_mesh_xyz[0] - 1;
						//if (k == kstart && j == start_coord[1]) {
						//	istart = start_coord[0];
						//}
						//if (k == kend && j == end_coord[1]) {
						//	iend = end_coord[0] - 2;
						//}
						if (istart >= iend) continue;
						
						// DO FIRST WINDOW - COPY ALL ELEMENTS

						int i = istart;
						const int d_i_0 = i >> 1; // data i
						const INDEX_TYPE baseid_0 = baseid_nox + i;
						if (!(baseid_0 < thread_start_id || baseid_0 > thread_end_id)) {
							const INDEX_TYPE data_baseid_0 = data_baseid_nox + d_i_0;

							// so for each vertex FIRST copy in the values
							// we can optimize this later to do less global lookups
							for (int pos = 0; pos < 27; pos++) {
								//int sd_nid = mc.small_mesh->get27NeighborOffset(pos) + kernel_baseid;
								INDEX_TYPE big_mesh_vertex_id = mMesh->get27NeighborOffset(pos) + baseid_0;
								INDEX_TYPE big_grid_vertex_data_id = m_data_27_offsets[pos] + data_baseid_0;
								INDEX_TYPE kernel_vertex_id = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos);
								//INDEX_TYPE kernel_data_nid = kernel_data_baseid + m_data_27_offsets[pos]; // this should just = pos??

								mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id] = big_mesh_vertex_id;
								mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id, this->mMaxVL->GetUncompressedMaxVal(big_mesh_vertex_id));
								mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id, this->mMaxVL->GetUncompressedMinVal(big_mesh_vertex_id));
								mc.small_grid_values[pos] = this->mFunc->SampleImage(big_grid_vertex_data_id);
							}

							ComputeLowerStar(mc, kernel_baseid);
						}

						
						istart += 2;

						for (i = istart; i < iend; i += 2) {
							const int d_i = i >> 1; // data i
							const INDEX_TYPE baseid = baseid_nox + i;
							if (baseid < thread_start_id || baseid > thread_end_id) continue;
							const INDEX_TYPE data_baseid = data_baseid_nox + d_i;

							for (int pos = 0; pos < 27; pos += 3) {
								//int sd_nid = mc.small_mesh->get27NeighborOffset(pos) + kernel_baseid;
								//INDEX_TYPE kernel_data_nid = kernel_data_baseid + m_data_27_offsets[pos]; // this should just = pos??
								INDEX_TYPE kernel_vertex_id_0 = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos);
								INDEX_TYPE kernel_vertex_id_next = kernel_vertex_id_0 + 2;
								mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_0] = mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_next];
								mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id_0, mc.small_mesh_maxmin_labeling->GetUncompressedMaxVal(kernel_vertex_id_next));
								mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id_0, mc.small_mesh_maxmin_labeling->GetUncompressedMinVal(kernel_vertex_id_next));


								INDEX_TYPE kernel_vertex_id_1 = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos + 1);
								INDEX_TYPE big_mesh_vertex_id_1 = mMesh->get27NeighborOffset(pos + 1) + baseid;
								mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_1] = big_mesh_vertex_id_1;
								mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id_1, this->mMaxVL->GetUncompressedMaxVal(big_mesh_vertex_id_1));
								mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id_1, this->mMaxVL->GetUncompressedMinVal(big_mesh_vertex_id_1));

								INDEX_TYPE kernel_vertex_id_2 = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos + 2);
								INDEX_TYPE big_mesh_vertex_id_2 = mMesh->get27NeighborOffset(pos + 2) + baseid;
								mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_2] = big_mesh_vertex_id_2;
								mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id_2, this->mMaxVL->GetUncompressedMaxVal(big_mesh_vertex_id_2));
								mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id_2, this->mMaxVL->GetUncompressedMinVal(big_mesh_vertex_id_2));

							}




							mc.small_grid_values[0] = mc.small_grid_values[1];
							mc.small_grid_values[1] = mc.small_grid_values[2];
							mc.small_grid_values[2] = this->mFunc->SampleImage(m_data_27_offsets[2] + data_baseid);

							mc.small_grid_values[3] = mc.small_grid_values[4];
							mc.small_grid_values[4] = mc.small_grid_values[5];
							mc.small_grid_values[5] = this->mFunc->SampleImage(m_data_27_offsets[5] + data_baseid);

							mc.small_grid_values[6] = mc.small_grid_values[7];
							mc.small_grid_values[7] = mc.small_grid_values[8];
							mc.small_grid_values[8] = this->mFunc->SampleImage(m_data_27_offsets[8] + data_baseid);

							mc.small_grid_values[9] = mc.small_grid_values[10];
							mc.small_grid_values[10] = mc.small_grid_values[11];
							mc.small_grid_values[11] = this->mFunc->SampleImage(m_data_27_offsets[11] + data_baseid);

							mc.small_grid_values[12] = mc.small_grid_values[13];
							mc.small_grid_values[13] = mc.small_grid_values[14];
							mc.small_grid_values[14] = this->mFunc->SampleImage(m_data_27_offsets[14] + data_baseid);

							mc.small_grid_values[15] = mc.small_grid_values[16];
							mc.small_grid_values[16] = mc.small_grid_values[17];
							mc.small_grid_values[17] = this->mFunc->SampleImage(m_data_27_offsets[17] + data_baseid);

							mc.small_grid_values[18] = mc.small_grid_values[19];
							mc.small_grid_values[19] = mc.small_grid_values[20];
							mc.small_grid_values[20] = this->mFunc->SampleImage(m_data_27_offsets[20] + data_baseid);

							mc.small_grid_values[21] = mc.small_grid_values[22];
							mc.small_grid_values[22] = mc.small_grid_values[23];
							mc.small_grid_values[23] = this->mFunc->SampleImage(m_data_27_offsets[23] + data_baseid);

							mc.small_grid_values[24] = mc.small_grid_values[25];
							mc.small_grid_values[25] = mc.small_grid_values[26];
							mc.small_grid_values[26] = this->mFunc->SampleImage(m_data_27_offsets[26] + data_baseid);

							ComputeLowerStar(mc, kernel_baseid);

						}

						

					}
				}

				delete mc.small_mesh;
				delete mc.small_grid;
				delete mc.small_grid_func;
				delete mc.small_mesh_maxmin_labeling;
			}
			//lstars_count = 0;
			// DO Z Plane Boundaries
			for (int j = 0; j < big_mesh_xyz[1]; j += 2) {
				for (auto k : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[2] - 1 })) { // do parallel division of work
					for (int i = 0; i < big_mesh_xyz[0]; i += 2) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d Z boundaries\n", lstars_count);
			//int tmp = lstars_count;
			//lstars_count = 0;

			// DO Y Plane Boundaries
			for (int k = 2; k < big_mesh_xyz[2] - 2; k += 2) { // smaller range since we did k = 0 and k = xyz[2]-1
				for (auto j : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[1] - 1 })) {
					for (int i = 0; i < big_mesh_xyz[0]; i += 2) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d Y boundaries\n", lstars_count);
			//int tmp2 = lstars_count;
			//lstars_count = 0;
			// DO X Plane Boundaries
			for (int k = 2; k < big_mesh_xyz[2] - 2; k += 2) { // smaller range since we did k = 0 and k = xyz[2]-1
				for (int j = 2; j < big_mesh_xyz[1] - 2; j += 2) { // again smaller range
					for (auto i : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[0] - 1 })) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d X boundaries\n", lstars_count);

			//printf("new robins1 %d lower stars in  %dms\n", lstars_count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			//now_time = std::chrono::steady_clock::now();
			//lstars_count = 0;
			//for (int k = 2; k < big_mesh_xyz[2] - 3; k += 2) { // do parallel division of work
			//	for (int j = 2; j < big_mesh_xyz[1] - 3; j += 2) {
			//		for (int i = 2; i < big_mesh_xyz[0] - 3; i += 2) {
			//			const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));

			//			INDEX_TYPE pre_pair = mGrad->getPair(baseid);

			//			this->mStandardRobins->ComputeLowerStar(baseid);
			//			lstars_count++;
			//			INDEX_TYPE post_pair = mGrad->getPair(baseid);

			//			if (pre_pair != post_pair) {
			//				printf("asdasdf\n");
			//			}

			//		}
			//	}
			//}

			//printf("old robins1 %d lower stars in  %dms\n", lstars_count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
		}

		void ComputePairing_slidingSlices() {

			std::chrono::steady_clock::time_point now_time = std::chrono::steady_clock::now();
			std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

			// dimensions of the mesh
			Vec3l big_mesh_xyz = mMesh->XYZ();
			//int lstars_count = 0;

#ifdef DEBUGPARALLEL
			db_counter = new int[mMesh->numCells()];
			memset(db_counter, 0, sizeof(int) * mMesh->numCells());
#endif


			// START PARALLEL WORK
#pragma omp parallel
			{
				std::vector<INDEX_TYPE> topo_index_partition;
				int num_threads;
				num_threads = omp_get_num_threads();
				ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, topo_index_partition);
				int thread_num = omp_get_thread_num();

				// these coordinates are INCLUSIVE - which means do start and end
				INDEX_TYPE thread_start_id = topo_index_partition[thread_num];
				INDEX_TYPE thread_end_id = topo_index_partition[thread_num + 1] - 1;
				Vec3l start_coord, end_coord;
				mMesh->cellid2Coords(thread_start_id, start_coord);
				mMesh->cellid2Coords(thread_end_id, end_coord); // get inclusive coord
				//#pragma omp critical 
				//				{
				//					printf("thread %d doing:\n\t", thread_num);
				//					start_coord.PrintInt(); printf("\t");
				//					end_coord.PrintInt();
				//				}
				// iterate over all vertices
				MESH_CONTEXT mc; // gather the pointers rather than have to pass a million items
				// only need maxvl labeling and function values
				// and maybe reslabel
				mc.small_mesh = new Explicit3x3x3SmallRegularGrid();
				mc.small_grid = new RegularGrid3D(Vec3l(3, 3, 3), Vec3b(0, 0, 0));


				mc.small_grid_func = new RegularGridTrilinearFunction(mc.small_grid, mc.small_grid_values); // wrapper for our values
				// place to store our local copy of the max/min vertices for each cell
				mc.small_mesh_maxmin_labeling = new RegularGridMaxMinVertexLabeling3D<Explicit3x3x3SmallRegularGrid, RegularGridTrilinearFunction>(mc.small_mesh, mc.small_grid_func);
				mc.small_mesh_maxmin_labeling->HACK_init();

				const INDEX_TYPE kernel_baseid = mc.small_mesh->coords2Cellid(Vec3l(2, 2, 2));
				const INDEX_TYPE kernel_data_baseid = mc.small_grid->Index3d(Vec3l(1, 1, 1));


				int kstart = start_coord[2];
				if (kstart % 2 == 1) kstart--; // kstart cannot start on an odd number, if it is odd, start on prior?
				if (kstart == 0) kstart = 2;
				int kend = end_coord[2];
				if (kend == big_mesh_xyz[2] - 1) kend = big_mesh_xyz[2] - 2;
				// NOW DO ALL INTERIOR VERTICES
#ifdef DEBUGPARALLEL
#pragma omp critical 
				{
					printf("thread %d doing actual k: [%d:%d]\n", thread_num, kstart, kend);
				}
#endif

				for (int k = kstart; k <= kend; k += 2) { // do parallel division of work
					const int d_k = k >> 1; // data k

					int jstart = 2;
					int jend = big_mesh_xyz[1] - 1;
					//if (k == kstart) {
					//	jstart = start_coord[1];
					//} 
					//if (k == kend) {
					//	jend = end_coord[1] - 2;
					//}

					for (int j = jstart; j < jend; j += 2) {
						const int d_j = j >> 1; // data j
						const INDEX_TYPE baseid_nox = mMesh->coords2Cellid(Vec3l(0, j, k));
						const INDEX_TYPE data_baseid_nox = mGrid->Index3d(Vec3l(0, d_j, d_k));

						int istart = 2;
						int iend = big_mesh_xyz[0] - 1;
						//if (k == kstart && j == start_coord[1]) {
						//	istart = start_coord[0];
						//}
						//if (k == kend && j == end_coord[1]) {
						//	iend = end_coord[0] - 2;
						//}
						if (istart >= iend) continue;
						
						// DO FIRST WINDOW - COPY ALL ELEMENTS

						int i = istart;
						const int d_i_0 = i >> 1; // data i
						const INDEX_TYPE baseid_0 = baseid_nox + i;
						if (!(baseid_0 < thread_start_id || baseid_0 > thread_end_id)) {
							const INDEX_TYPE data_baseid_0 = data_baseid_nox + d_i_0;

							// so for each vertex FIRST copy in the values
							// we can optimize this later to do less global lookups
							for (int pos = 0; pos < 27; pos++) {
								//int sd_nid = mc.small_mesh->get27NeighborOffset(pos) + kernel_baseid;
								INDEX_TYPE big_mesh_vertex_id = mMesh->get27NeighborOffset(pos) + baseid_0;
								INDEX_TYPE big_grid_vertex_data_id = m_data_27_offsets[pos] + data_baseid_0;
								INDEX_TYPE kernel_vertex_id = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos);
								//INDEX_TYPE kernel_data_nid = kernel_data_baseid + m_data_27_offsets[pos]; // this should just = pos??

								mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id] = big_mesh_vertex_id;
								mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id, this->mMaxVL->GetUncompressedMaxVal(big_mesh_vertex_id));
								mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id, this->mMaxVL->GetUncompressedMinVal(big_mesh_vertex_id));
								mc.small_grid_values[pos] = this->mFunc->SampleImage(big_grid_vertex_data_id);
							}

							ComputeLowerStar(mc, kernel_baseid);
						}

						
						istart += 2;

						for (i = istart; i < iend; i += 2) {
							const int d_i = i >> 1; // data i
							const INDEX_TYPE baseid = baseid_nox + i;
							if (baseid < thread_start_id || baseid > thread_end_id) continue;
							const INDEX_TYPE data_baseid = data_baseid_nox + d_i;

							for (int pos = 0; pos < 27; pos += 3) {
								//int sd_nid = mc.small_mesh->get27NeighborOffset(pos) + kernel_baseid;
								//INDEX_TYPE kernel_data_nid = kernel_data_baseid + m_data_27_offsets[pos]; // this should just = pos??
								INDEX_TYPE kernel_vertex_id_0 = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos);
								INDEX_TYPE kernel_vertex_id_next = kernel_vertex_id_0 + 2;
								mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_0] = mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_next];
								mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id_0, mc.small_mesh_maxmin_labeling->GetUncompressedMaxVal(kernel_vertex_id_next));
								mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id_0, mc.small_mesh_maxmin_labeling->GetUncompressedMinVal(kernel_vertex_id_next));


								INDEX_TYPE kernel_vertex_id_1 = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos + 1);
								INDEX_TYPE big_mesh_vertex_id_1 = mMesh->get27NeighborOffset(pos + 1) + baseid;
								mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_1] = big_mesh_vertex_id_1;
								mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id_1, this->mMaxVL->GetUncompressedMaxVal(big_mesh_vertex_id_1));
								mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id_1, this->mMaxVL->GetUncompressedMinVal(big_mesh_vertex_id_1));

								INDEX_TYPE kernel_vertex_id_2 = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos + 2);
								INDEX_TYPE big_mesh_vertex_id_2 = mMesh->get27NeighborOffset(pos + 2) + baseid;
								mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_2] = big_mesh_vertex_id_2;
								mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id_2, this->mMaxVL->GetUncompressedMaxVal(big_mesh_vertex_id_2));
								mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id_2, this->mMaxVL->GetUncompressedMinVal(big_mesh_vertex_id_2));

							}




							mc.small_grid_values[0] = mc.small_grid_values[1];
							mc.small_grid_values[1] = mc.small_grid_values[2];
							mc.small_grid_values[2] = this->mFunc->SampleImage(m_data_27_offsets[2] + data_baseid);

							mc.small_grid_values[3] = mc.small_grid_values[4];
							mc.small_grid_values[4] = mc.small_grid_values[5];
							mc.small_grid_values[5] = this->mFunc->SampleImage(m_data_27_offsets[5] + data_baseid);

							mc.small_grid_values[6] = mc.small_grid_values[7];
							mc.small_grid_values[7] = mc.small_grid_values[8];
							mc.small_grid_values[8] = this->mFunc->SampleImage(m_data_27_offsets[8] + data_baseid);

							mc.small_grid_values[9] = mc.small_grid_values[10];
							mc.small_grid_values[10] = mc.small_grid_values[11];
							mc.small_grid_values[11] = this->mFunc->SampleImage(m_data_27_offsets[11] + data_baseid);

							mc.small_grid_values[12] = mc.small_grid_values[13];
							mc.small_grid_values[13] = mc.small_grid_values[14];
							mc.small_grid_values[14] = this->mFunc->SampleImage(m_data_27_offsets[14] + data_baseid);

							mc.small_grid_values[15] = mc.small_grid_values[16];
							mc.small_grid_values[16] = mc.small_grid_values[17];
							mc.small_grid_values[17] = this->mFunc->SampleImage(m_data_27_offsets[17] + data_baseid);

							mc.small_grid_values[18] = mc.small_grid_values[19];
							mc.small_grid_values[19] = mc.small_grid_values[20];
							mc.small_grid_values[20] = this->mFunc->SampleImage(m_data_27_offsets[20] + data_baseid);

							mc.small_grid_values[21] = mc.small_grid_values[22];
							mc.small_grid_values[22] = mc.small_grid_values[23];
							mc.small_grid_values[23] = this->mFunc->SampleImage(m_data_27_offsets[23] + data_baseid);

							mc.small_grid_values[24] = mc.small_grid_values[25];
							mc.small_grid_values[25] = mc.small_grid_values[26];
							mc.small_grid_values[26] = this->mFunc->SampleImage(m_data_27_offsets[26] + data_baseid);

							ComputeLowerStar(mc, kernel_baseid);

						}

						

					}
				}
			}
			//printf("INTERIOR: new robins1 in  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			now_time = std::chrono::steady_clock::now();
			//lstars_count = 0;
			// DO Z Plane Boundaries
#pragma omp parallel for
			for (int j = 0; j < big_mesh_xyz[1]; j += 2) {
				for (auto k : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[2] - 1 })) { // do parallel division of work
					for (int i = 0; i < big_mesh_xyz[0]; i += 2) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d Z boundaries\n", lstars_count);
			//int tmp = lstars_count;
			//lstars_count = 0;

			// DO Y Plane Boundaries
#pragma omp parallel for
			for (int k = 2; k < big_mesh_xyz[2] - 2; k += 2) { // smaller range since we did k = 0 and k = xyz[2]-1
				for (auto j : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[1] - 1 })) {
					for (int i = 0; i < big_mesh_xyz[0]; i += 2) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d Y boundaries\n", lstars_count);
			//int tmp2 = lstars_count;
			//lstars_count = 0;
			// DO X Plane Boundaries
#pragma omp parallel for
			for (int k = 2; k < big_mesh_xyz[2] - 2; k += 2) { // smaller range since we did k = 0 and k = xyz[2]-1
				for (int j = 2; j < big_mesh_xyz[1] - 2; j += 2) { // again smaller range
					for (auto i : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[0] - 1 })) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d X boundaries\n", lstars_count);
			//printf("BOUNDARY: new robins1 in  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			now_time = std::chrono::steady_clock::now();

			//printf("new robins1 %d lower stars in  %dms\n", lstars_count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			//now_time = std::chrono::steady_clock::now();
			//lstars_count = 0;
			//for (int k = 2; k < big_mesh_xyz[2] - 3; k += 2) { // do parallel division of work
			//	for (int j = 2; j < big_mesh_xyz[1] - 3; j += 2) {
			//		for (int i = 2; i < big_mesh_xyz[0] - 3; i += 2) {
			//			const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));

			//			INDEX_TYPE pre_pair = mGrad->getPair(baseid);

			//			this->mStandardRobins->ComputeLowerStar(baseid);
			//			lstars_count++;
			//			INDEX_TYPE post_pair = mGrad->getPair(baseid);

			//			if (pre_pair != post_pair) {
			//				printf("asdasdf\n");
			//			}

			//		}
			//	}
			//}

			//printf("old robins1 %d lower stars in  %dms\n", lstars_count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());

#ifdef DEBUGPARALLEL
			FILE* fout = fopen("test_out.raw", "wb");
			fwrite(db_counter, sizeof(int), mMesh->numCells(), fout);
			fclose(fout);
			//printf("doing seen checks!\n");
			//for (INDEX_TYPE i = 0; i < mMesh->numCells(); i++) {
			//	
			//	if (mMesh->boundaryValue(i) == 0 && db_counter[i] != 1) {
			//		printf("index %lld seen %d times: ", i, db_counter[i]);
			//		Vec3l c;
			//		mMesh->cellid2Coords(i, c);
			//		c.PrintInt();
			//	}

			//}
			printf("done seen checks\n");


#endif


		}

		void ComputePairing_sliding() {

			std::chrono::steady_clock::time_point now_time = std::chrono::steady_clock::now();
			std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

			// dimensions of the mesh
			Vec3l big_mesh_xyz = mMesh->XYZ();
			//int lstars_count = 0;

#ifdef DEBUGPARALLEL
			db_counter = new int[mMesh->numCells()];
			memset(db_counter, 0, sizeof(int) * mMesh->numCells());
#endif


			// START PARALLEL WORK
#pragma omp parallel
			{
				std::vector<INDEX_TYPE> topo_index_partition;
				int num_threads;
				num_threads = omp_get_num_threads();
				ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, topo_index_partition);
				int thread_num = omp_get_thread_num();

				// these coordinates are INCLUSIVE - which means do start and end
				INDEX_TYPE thread_start_id = topo_index_partition[thread_num];
				INDEX_TYPE thread_end_id = topo_index_partition[thread_num + 1] - 1;
				Vec3l start_coord, end_coord;
				mMesh->cellid2Coords(thread_start_id, start_coord);
				mMesh->cellid2Coords(thread_end_id, end_coord); // get inclusive coord
				//#pragma omp critical 
				//				{
				//					printf("thread %d doing:\n\t", thread_num);
				//					start_coord.PrintInt(); printf("\t");
				//					end_coord.PrintInt();
				//				}
				// iterate over all vertices
				MESH_CONTEXT mc; // gather the pointers rather than have to pass a million items
				// only need maxvl labeling and function values
				// and maybe reslabel
				mc.small_mesh = new Explicit3x3x3SmallRegularGrid();
				mc.small_grid = new RegularGrid3D(Vec3l(3, 3, 3), Vec3b(0, 0, 0));


				mc.small_grid_func = new RegularGridTrilinearFunction(mc.small_grid, mc.small_grid_values); // wrapper for our values
				// place to store our local copy of the max/min vertices for each cell
				mc.small_mesh_maxmin_labeling = new RegularGridMaxMinVertexLabeling3D<Explicit3x3x3SmallRegularGrid, RegularGridTrilinearFunction>(mc.small_mesh, mc.small_grid_func);
				mc.small_mesh_maxmin_labeling->HACK_init();

				const INDEX_TYPE kernel_baseid = mc.small_mesh->coords2Cellid(Vec3l(2, 2, 2));
				const INDEX_TYPE kernel_data_baseid = mc.small_grid->Index3d(Vec3l(1, 1, 1));


				int kstart = start_coord[2];
				if (kstart % 2 == 1) kstart--; // kstart cannot start on an odd number, if it is odd, start on prior?
				if (kstart == 0) kstart = 2;
				int kend = end_coord[2];
				if (kend == big_mesh_xyz[2] - 1) kend = big_mesh_xyz[2] - 2;
				// NOW DO ALL INTERIOR VERTICES
#ifdef DEBUGPARALLEL
#pragma omp critical 
				{
					printf("thread %d doing actual k: [%d:%d]\n", thread_num, kstart, kend);
				}
#endif

				for (int k = kstart; k <= kend; k += 2) { // do parallel division of work
					const int d_k = k >> 1; // data k

					int jstart = 2;
					int jend = big_mesh_xyz[1] - 1;
					//if (k == kstart) {
					//	jstart = start_coord[1];
					//} 
					//if (k == kend) {
					//	jend = end_coord[1] - 2;
					//}

					for (int j = jstart; j < jend; j += 2) {
						const int d_j = j >> 1; // data j
						const INDEX_TYPE baseid_nox = mMesh->coords2Cellid(Vec3l(0, j, k));
						const INDEX_TYPE data_baseid_nox = mGrid->Index3d(Vec3l(0, d_j, d_k));

						int istart = 2;
						int iend = big_mesh_xyz[0] - 1;
						//if (k == kstart && j == start_coord[1]) {
						//	istart = start_coord[0];
						//}
						//if (k == kend && j == end_coord[1]) {
						//	iend = end_coord[0] - 2;
						//}
						if (istart >= iend) continue;
						
						bool first_window = true;
						for (int i = istart; i < iend; i += 2) {
							const int d_i = i >> 1; // data i
							const INDEX_TYPE baseid = baseid_nox + i;
							if (baseid < thread_start_id || baseid > thread_end_id) continue;
							const INDEX_TYPE data_baseid = data_baseid_nox + d_i;

							// DO FIRST WINDOW - COPY ALL ELEMENTS
							if (first_window) {
								// so for each vertex FIRST copy in the values
								// we can optimize this later to do less global lookups
								for (int pos = 0; pos < 27; pos++) {
									//int sd_nid = mc.small_mesh->get27NeighborOffset(pos) + kernel_baseid;
									INDEX_TYPE big_mesh_vertex_id = mMesh->get27NeighborOffset(pos) + baseid;
									INDEX_TYPE big_grid_vertex_data_id = m_data_27_offsets[pos] + data_baseid;
									INDEX_TYPE kernel_vertex_id = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos);
									//INDEX_TYPE kernel_data_nid = kernel_data_baseid + m_data_27_offsets[pos]; // this should just = pos??

									mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id] = big_mesh_vertex_id;
									mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id, this->mMaxVL->GetUncompressedMaxVal(big_mesh_vertex_id));
									mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id, this->mMaxVL->GetUncompressedMinVal(big_mesh_vertex_id));
									mc.small_grid_values[pos] = this->mFunc->SampleImage(big_grid_vertex_data_id);
								}

								first_window = false;
							} else {
								for (int pos = 0; pos < 27; pos += 3) {
									//int sd_nid = mc.small_mesh->get27NeighborOffset(pos) + kernel_baseid;
									//INDEX_TYPE kernel_data_nid = kernel_data_baseid + m_data_27_offsets[pos]; // this should just = pos??
									INDEX_TYPE kernel_vertex_id_0 = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos);
									INDEX_TYPE kernel_vertex_id_next = kernel_vertex_id_0 + 2;
									mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_0] = mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_next];
									mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id_0, mc.small_mesh_maxmin_labeling->GetUncompressedMaxVal(kernel_vertex_id_next));
									mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id_0, mc.small_mesh_maxmin_labeling->GetUncompressedMinVal(kernel_vertex_id_next));

									INDEX_TYPE kernel_vertex_id_1 = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos + 1);
									INDEX_TYPE big_mesh_vertex_id_1 = mMesh->get27NeighborOffset(pos + 1) + baseid;
									mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_1] = big_mesh_vertex_id_1;
									mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id_1, this->mMaxVL->GetUncompressedMaxVal(big_mesh_vertex_id_1));
									mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id_1, this->mMaxVL->GetUncompressedMinVal(big_mesh_vertex_id_1));

									INDEX_TYPE kernel_vertex_id_2 = kernel_baseid + mc.small_mesh->get27NeighborOffset(pos + 2);
									INDEX_TYPE big_mesh_vertex_id_2 = mMesh->get27NeighborOffset(pos + 2) + baseid;
									mc.small_mesh_id_to_big_mesh_id[kernel_vertex_id_2] = big_mesh_vertex_id_2;
									mc.small_mesh_maxmin_labeling->SetUncompressedMaxVal(kernel_vertex_id_2, this->mMaxVL->GetUncompressedMaxVal(big_mesh_vertex_id_2));
									mc.small_mesh_maxmin_labeling->SetUncompressedMinVal(kernel_vertex_id_2, this->mMaxVL->GetUncompressedMinVal(big_mesh_vertex_id_2));

								}

								mc.small_grid_values[0] = mc.small_grid_values[1];
								mc.small_grid_values[1] = mc.small_grid_values[2];
								mc.small_grid_values[2] = this->mFunc->SampleImage(m_data_27_offsets[2] + data_baseid);

								mc.small_grid_values[3] = mc.small_grid_values[4];
								mc.small_grid_values[4] = mc.small_grid_values[5];
								mc.small_grid_values[5] = this->mFunc->SampleImage(m_data_27_offsets[5] + data_baseid);

								mc.small_grid_values[6] = mc.small_grid_values[7];
								mc.small_grid_values[7] = mc.small_grid_values[8];
								mc.small_grid_values[8] = this->mFunc->SampleImage(m_data_27_offsets[8] + data_baseid);

								mc.small_grid_values[9] = mc.small_grid_values[10];
								mc.small_grid_values[10] = mc.small_grid_values[11];
								mc.small_grid_values[11] = this->mFunc->SampleImage(m_data_27_offsets[11] + data_baseid);

								mc.small_grid_values[12] = mc.small_grid_values[13];
								mc.small_grid_values[13] = mc.small_grid_values[14];
								mc.small_grid_values[14] = this->mFunc->SampleImage(m_data_27_offsets[14] + data_baseid);

								mc.small_grid_values[15] = mc.small_grid_values[16];
								mc.small_grid_values[16] = mc.small_grid_values[17];
								mc.small_grid_values[17] = this->mFunc->SampleImage(m_data_27_offsets[17] + data_baseid);

								mc.small_grid_values[18] = mc.small_grid_values[19];
								mc.small_grid_values[19] = mc.small_grid_values[20];
								mc.small_grid_values[20] = this->mFunc->SampleImage(m_data_27_offsets[20] + data_baseid);

								mc.small_grid_values[21] = mc.small_grid_values[22];
								mc.small_grid_values[22] = mc.small_grid_values[23];
								mc.small_grid_values[23] = this->mFunc->SampleImage(m_data_27_offsets[23] + data_baseid);

								mc.small_grid_values[24] = mc.small_grid_values[25];
								mc.small_grid_values[25] = mc.small_grid_values[26];
								mc.small_grid_values[26] = this->mFunc->SampleImage(m_data_27_offsets[26] + data_baseid);
							}

							ComputeLowerStar(mc, kernel_baseid);
						}
					}
				}
				
				delete mc.small_mesh;
				delete mc.small_grid;
				delete mc.small_grid_func;
				delete mc.small_mesh_maxmin_labeling;
			}
			//printf("INTERIOR: new robins1 in  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			now_time = std::chrono::steady_clock::now();
			//lstars_count = 0;
			// DO Z Plane Boundaries
#pragma omp parallel for
			for (int j = 0; j < big_mesh_xyz[1]; j += 2) {
				for (auto k : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[2] - 1 })) { // do parallel division of work
					for (int i = 0; i < big_mesh_xyz[0]; i += 2) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d Z boundaries\n", lstars_count);
			//int tmp = lstars_count;
			//lstars_count = 0;

			// DO Y Plane Boundaries
#pragma omp parallel for
			for (int k = 2; k < big_mesh_xyz[2] - 2; k += 2) { // smaller range since we did k = 0 and k = xyz[2]-1
				for (auto j : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[1] - 1 })) {
					for (int i = 0; i < big_mesh_xyz[0]; i += 2) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d Y boundaries\n", lstars_count);
			//int tmp2 = lstars_count;
			//lstars_count = 0;
			// DO X Plane Boundaries
#pragma omp parallel for
			for (int k = 2; k < big_mesh_xyz[2] - 2; k += 2) { // smaller range since we did k = 0 and k = xyz[2]-1
				for (int j = 2; j < big_mesh_xyz[1] - 2; j += 2) { // again smaller range
					for (auto i : std::vector<INDEX_TYPE>({ 0, big_mesh_xyz[0] - 1 })) {
						const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));
						this->mStandardRobins->ComputeLowerStar(baseid);
						//lstars_count++;
					}
				}
			}
			//printf("did %d X boundaries\n", lstars_count);
			//printf("BOUNDARY: new robins1 in  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			now_time = std::chrono::steady_clock::now();

			//printf("new robins1 %d lower stars in  %dms\n", lstars_count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			//now_time = std::chrono::steady_clock::now();
			//lstars_count = 0;
			//for (int k = 2; k < big_mesh_xyz[2] - 3; k += 2) { // do parallel division of work
			//	for (int j = 2; j < big_mesh_xyz[1] - 3; j += 2) {
			//		for (int i = 2; i < big_mesh_xyz[0] - 3; i += 2) {
			//			const INDEX_TYPE baseid = mMesh->coords2Cellid(Vec3l(i, j, k));

			//			INDEX_TYPE pre_pair = mGrad->getPair(baseid);

			//			this->mStandardRobins->ComputeLowerStar(baseid);
			//			lstars_count++;
			//			INDEX_TYPE post_pair = mGrad->getPair(baseid);

			//			if (pre_pair != post_pair) {
			//				printf("asdasdf\n");
			//			}

			//		}
			//	}
			//}

			//printf("old robins1 %d lower stars in  %dms\n", lstars_count, std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());

#ifdef DEBUGPARALLEL
			FILE* fout = fopen("test_out.raw", "wb");
			fwrite(db_counter, sizeof(int), mMesh->numCells(), fout);
			fclose(fout);
			//printf("doing seen checks!\n");
			//for (INDEX_TYPE i = 0; i < mMesh->numCells(); i++) {
			//	
			//	if (mMesh->boundaryValue(i) == 0 && db_counter[i] != 1) {
			//		printf("index %lld seen %d times: ", i, db_counter[i]);
			//		Vec3l c;
			//		mMesh->cellid2Coords(i, c);
			//		c.PrintInt();
			//	}

			//}
			printf("done seen checks\n");


#endif


		}



	};





}
#endif
