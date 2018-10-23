/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef EXTREMAL_REGION_BUILDER
#define EXTREMAL_REGION_BUILDER

#include <set>
#include <queue>
#include <vector>
#include <stack>
#include <unordered_map>
#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_labeling.h"
//#include "gi_regular_grid.h"
//#include "gi_regular_grid_trilinear_function.h"
#include "gi_timing.h"
#include "gi_union_find_labeling.h"
//#include "gi_topological_regular_grid.h"
#include "gi_array_index_partition.h"
#include "omp.h"
//#include "gi_numeric_integrator_expanding_region_stop_filtered2.h"
//#define OUTPUTINTERMEDIATE

namespace GInt {





	template<class FUNC_TYPE, class MESH_TYPE>
	class UFMergeGraph {
	protected:
		FUNC_TYPE* mFunc;
		MESH_TYPE* mMesh;
		
		struct ExtremumSet {
			INDEX_TYPE extremum_cell_index;
			INT_TYPE representative_list_id;
			typename FUNC_TYPE::DType extremum_fval;
			ExtremumSet(INDEX_TYPE gid, INT_TYPE lid, typename FUNC_TYPE::DType val) : extremum_cell_index(gid), representative_list_id(lid), extremum_fval(val) {}
			void printme() const {
				printf("extremum_cell_index=%llu, rep=%d, val=%f\n", extremum_cell_index, representative_list_id, extremum_fval);
			}
		};

		struct MergeArc {
			INDEX_TYPE saddle_cell_index;
			INT_TYPE extremum1_list_id;
			INT_TYPE extremum2_list_id;
			typename FUNC_TYPE::DType merge_persistence;
			typename FUNC_TYPE::DType merge_saddle_value;
			void printme() const {
				printf("saddle_cell_index=%llu, e1=%d, e2=%d, p=%f, s=%f\n", saddle_cell_index, extremum1_list_id, extremum2_list_id, merge_persistence, merge_saddle_value);
			}

			bool operator()(const MergeArc& _Left, const MergeArc& _Right) {
				if (_Left.merge_persistence < _Right.merge_persistence) return false;
				if (_Left.merge_persistence > _Right.merge_persistence) return true;
				return _Left.saddle_cell_index > _Right.saddle_cell_index;
			}

		};
		int countcancels;
	public:
      std::unordered_map<INDEX_TYPE, INT_TYPE> mCellIndexToListIdMap;
      std::vector<ExtremumSet> mExtrema;
      std::priority_queue<MergeArc, std::vector< MergeArc>, MergeArc > mMergesToCancel;

		void MakeSet(INDEX_TYPE gid) {
			//printf("%llu gid\n", gid);
			INT_TYPE lid = mExtrema.size();
			//printf("%d lid\n", lid);
			mExtrema.push_back(ExtremumSet{ gid, lid, mFunc->cellValue(gid) });
			//printf("extrema size %d\n", mExtrema.size());
         //printf("value = %f\n", mExtrema[mExtrema.size()-1].extremum_fval);
			mCellIndexToListIdMap[gid] = lid;
			//printf("idmap size %d\n", mCellIndexToListIdMap.size());
		}

		INT_TYPE FindRepresentative(INT_TYPE lid) {
			ExtremumSet& e = mExtrema[lid];
			//printf("%d->%d\n", lid, e.representative_list_id);
			if (e.representative_list_id == lid) return lid;
			INT_TYPE i = FindRepresentative(e.representative_list_id);
			e.representative_list_id = i;
			return i;
		}

		bool mDoMinHierarchy;


		void MergeByVal(INT_TYPE lid1, INT_TYPE lid2) {
			ExtremumSet& e1 = mExtrema[lid1];
			ExtremumSet& e2 = mExtrema[lid2];
			// if reverse order is false, set the lower one to point
			// to the higher one - so reverse order true used for basins
			if (mFunc->lessThan(e1.extremum_cell_index, e2.extremum_cell_index) ^ mDoMinHierarchy) {
				e1.representative_list_id = e2.representative_list_id;
				//printf("%f -> %f\n", e1.extremum_fval, e2.extremum_fval);
			}
			else {
				e2.representative_list_id = e1.representative_list_id;
				//printf("%f -> %f\n", e2.extremum_fval, e1.extremum_fval);
			}

		}

		void PerformCancel(MergeArc& m) {
			if (m.merge_persistence > gThreshold) {
#ifdef 	DEBUGGSERB		
				printf("threshold too great %f %f\n", m.merge_persistence, gThreshold);
#endif
				return;
			}

			//printf("e1=%d, e2=%d\n", m.extremum1_list_id, m.extremum2_list_id);
			INT_TYPE current_extrep1_lid = FindRepresentative(m.extremum1_list_id);
			INT_TYPE current_extrep2_lid = FindRepresentative(m.extremum2_list_id);
			//printf("found %d %d\n", current_extrep1_lid, current_extrep2_lid);
			if (current_extrep1_lid == current_extrep2_lid) return; // do nothing for loops
			ExtremumSet& current_extrep1 = mExtrema[current_extrep2_lid];
			ExtremumSet& current_extrep2 = mExtrema[current_extrep1_lid];

			// if this arc is no longer valid then insert the new arc
			// if the current_extrep2 is not the same then we need to modify
			if (m.extremum1_list_id != current_extrep1_lid ||
				m.extremum2_list_id != current_extrep2_lid) {
				// reinsert?
				typename FUNC_TYPE::DType d1 = abs(m.merge_saddle_value - current_extrep1.extremum_fval);
				typename FUNC_TYPE::DType d2 = abs(m.merge_saddle_value - current_extrep2.extremum_fval);
				if (d1 < d2) {
					if (d1 > gThreshold) {
#ifdef 	DEBUGGSERB		
						printf("thresha too big %f, %f, %f\n", d1, d2, gThreshold);
#endif
						return;
					}
					m.merge_persistence = d1;
					m.extremum1_list_id = current_extrep2_lid;
					m.extremum2_list_id = current_extrep1_lid;
					mMergesToCancel.push(m);
					//struct MergeArc {
					//	INDEX_TYPE saddle_cell_index;
					//	INT_TYPE extremum1_list_id;
					//	INT_TYPE headExtId;
					//	typename FUNC_TYPE::DType merge_persistence;
					//	typename FUNC_TYPE::DType merge_saddle_value;
				}
				else {
					if (d2 > gThreshold) {
#ifdef 	DEBUGGSERB		
						printf("threshb too big %f, %f, %f\n", d1, d2, gThreshold);
#endif
						return;
					}
					m.merge_persistence = d2;
					m.extremum1_list_id = current_extrep1_lid;
					m.extremum2_list_id = current_extrep2_lid;
					mMergesToCancel.push(m);

				}
				//printf("merged pers %f\n", m.merge_persistence);
				return;
			}
			//if (m.merge_saddle_value == -1) {
			//	printf("pers = %f, sad = %f\n", m.merge_persistence, m.merge_saddle_value);
			//	m.printme();
			//	current_extrep1.printme();
			//	current_extrep2.printme();
			//	printf("frommesh: dim=%d, val=%f, boundary=%d\n", mMesh->dimension(m.saddle_cell_index), mFunc->cellValue(m.saddle_cell_index), mMesh->boundaryValue(m.saddle_cell_index));
			//	MESH_TYPE::FacetsIterator fit(mMesh);
			//	for (fit.begin(m.saddle_cell_index); fit.valid(); fit.advance()) {
			//		printf(" --> %f\n", mFunc->cellValue(fit.value()));
			//	}
			//}
			//return;
			// otherwise actually do the merge
			MergeByVal(current_extrep1_lid, current_extrep2_lid);
			countcancels++;
		}


		typename FUNC_TYPE::DType gThreshold;

	public:
		enum Direction { MAXIMAL, MINIMAL };
		UFMergeGraph(FUNC_TYPE* func, MESH_TYPE* mesh, Direction d ) : mFunc(func), mMesh(mesh) {
			if (d == MINIMAL) {
				mDoMinHierarchy = true;
			}
			else {
				mDoMinHierarchy = false;
			}
			countcancels = 0;
		}
		INT_TYPE NumExtrema() const {
			return this->mExtrema.size();
		}

		INT_TYPE Representative(INT_TYPE nid) const {
			return this->mExtrema[nid].representative_list_id;
		}

		INDEX_TYPE TopoIndex(INT_TYPE nid) const {
			return this->mExtrema[nid].extremum_cell_index;

		}

		void AddNode(INDEX_TYPE key) {
			MakeSet(key);
		}

		void AddArc(INDEX_TYPE a, INDEX_TYPE b, INDEX_TYPE s) {
			if (a == b) return; // ignore these
			
			// ignore connectors that span different boundary types
			DIM_TYPE bva = mMesh->boundaryValue(a);
			DIM_TYPE bvb = mMesh->boundaryValue(b);
			//if (bvb != bva) return;
			DIM_TYPE bvs = mMesh->boundaryValue(s);
			//if (bva != bvs) return;

			// so this could be a valid cancellation so add the mergearc
			MergeArc m;
			m.saddle_cell_index = s;
			m.extremum1_list_id = mCellIndexToListIdMap[a];
			m.extremum2_list_id = mCellIndexToListIdMap[b];
			m.merge_saddle_value = mFunc->cellValue(s);
			//printf("sadval = %f\n", m.merge_saddle_value);
			// reinsert?
			typename FUNC_TYPE::DType d1 =  mExtrema[m.extremum1_list_id].extremum_fval - m.merge_saddle_value ;
			typename FUNC_TYPE::DType d2 =  mExtrema[m.extremum2_list_id].extremum_fval - m.merge_saddle_value ;
         if (d1 < 0) d1 *= -1;
         if (d2 < 0) d2 *= -1;
			if (d1 < d2) {
				m.merge_persistence = d1;
			}
			else {
				m.merge_persistence = d2;
			}
		//	printf("found %f pers d1=%f d2=%f s=%f e1=%f e2=%f\n", m.merge_persistence, d1, d2, m.merge_saddle_value, mExtrema[m.extremum1_list_id].extremum_fval, mExtrema[m.extremum2_list_id].extremum_fval);
			mMergesToCancel.push(m);

		}


		int SimplifyToThreshold(typename FUNC_TYPE::DType threshold) {
			gThreshold = threshold;
			
			///printf("simplifying to %f\n", threshold);

			while (!mMergesToCancel.empty()) {
				MergeArc edge = mMergesToCancel.top();
				mMergesToCancel.pop();
				//printf("got %d left\n", mMergesToCancel.size());
				PerformCancel(edge);
			}
			//printf("gothere - siplified now flattening\n");
			// flatten UF so all subsequent queries are CONST
			for (INT_TYPE esid = 0; esid < mExtrema.size(); esid++) {
				FindRepresentative(esid);
			}
			//printf("done!\n");
			//printf("performed %d merges\n", countcancels);
			return countcancels;
		}

		

	};

	// we will use FUNC_TYPE to access index comparator
	template<class MESH_TYPE, class FUNC_TYPE, class GRAD_TYPE>
	class SimplifiedExtremumGraph {


	public:
		// global context
		GRAD_TYPE* mGrad;
		MESH_TYPE* mMesh;
		FUNC_TYPE* mFunc;
		UFMergeGraph<FUNC_TYPE, MESH_TYPE>* mMinGraph;
		UFMergeGraph<FUNC_TYPE, MESH_TYPE>* mMaxGraph;
		bool do_mins;
		bool do_maxs;


		INDEX_TYPE rec_td( INDEX_TYPE cellid) const {
			//printf("a %llu d=%d\n", cellid, mMesh->dimension(cellid));
			INDEX_TYPE current = cellid;

			if (mGrad->getCritical(current)) return current;
			//printf("b\n");
			INDEX_TYPE pair = mGrad->getPair(current);

			//printf("c pair = %llu, d=%d\n", pair, mMesh->dimension(pair));
			if (mGrad->getCritical(pair)) return pair; // should NEVER happen
			//printf("c.1\n");
			typename MESH_TYPE::FacetsIterator facets(mMesh);
			facets.begin(pair);
			INDEX_TYPE next = facets.value();
			//printf("next = %llu\n", next);
			if (next == current) {
				facets.advance();
				next = facets.value();
			}
			//printf("next = %llu\n", next);

			//printf("d\n");
			return rec_td(next);
		}

		void trace_down_1saddle(INDEX_TYPE start, INDEX_TYPE& min1, INDEX_TYPE& min2) const {
			
			typename MESH_TYPE::FacetsIterator facets(mMesh);
			facets.begin(start);

			min1 = rec_td(facets.value());
			facets.advance();
			if (!facets.valid()) {
				// error
				printf("should never get here in tarce down 1saddle\n");
				min2 = min1; return;
			}

			min2 = rec_td(facets.value());
			return;
		}



		INDEX_TYPE rec_tu(INDEX_TYPE cellid) const {
			INDEX_TYPE current = cellid;

			if (mGrad->getCritical(current)) return current;

			INDEX_TYPE pair = mGrad->getPair(current);

			if (mGrad->getCritical(pair)) return pair; // should NEVER happen

			typename MESH_TYPE::CofacetsIterator cofacets(mMesh);
			cofacets.begin(pair);
			INDEX_TYPE next = cofacets.value();
			if (next == current) {
				cofacets.advance();
				if (!cofacets.valid()) {
					printf("WHOATHERE\n");
				}
				next = cofacets.value();
			}


			return rec_tu(next);
		}

		void trace_up_2saddle(const INDEX_TYPE& start, INDEX_TYPE& max1, INDEX_TYPE& max2) const {
			if (mMesh->boundaryValue(start) != 0) {
				max1 = max2 = -1;
				return;
			}
			typename MESH_TYPE::CofacetsIterator cofacets(mMesh);
			cofacets.begin(start);

			max1 = rec_tu(cofacets.value());
			cofacets.advance();
			if (!cofacets.valid()) {
				// error
				printf("should never get here in tarce up 2saddle\n");
				max2 = max1; return;
			}

			max2 = rec_tu(cofacets.value());
			return;
		}

	public:


		SimplifiedExtremumGraph(MESH_TYPE* mesh, FUNC_TYPE* func, GRAD_TYPE* grad) :
			mGrad(grad), mMesh(mesh), mFunc(func){
			mMinGraph = new UFMergeGraph < FUNC_TYPE, MESH_TYPE >(func, mesh,
				UFMergeGraph < FUNC_TYPE, MESH_TYPE >::MINIMAL);
			mMaxGraph = new UFMergeGraph < FUNC_TYPE, MESH_TYPE >(func, mesh,
				UFMergeGraph < FUNC_TYPE, MESH_TYPE >::MAXIMAL);
		}
		enum EXTGRAPHMODE { NONE, BOTH, MINS, MAXS };
		void SetMode(EXTGRAPHMODE m) {
			switch (m) {
			case NONE:
				this->do_maxs = false;
				this->do_mins = false;
				return;
			case BOTH:
				this->do_maxs = true;
				this->do_mins = true;
				return;
			case MINS:
				this->do_maxs = false;
				this->do_mins = true;
				return;
			case MAXS:
				this->do_maxs = true;
				this->do_mins = false;
				return;
			}
		}

		void ComputeMinMapFromGradient(typename FUNC_TYPE::DType THRESHOLD) {
         std::vector<INDEX_TYPE> saddles1;
         std::vector<INDEX_TYPE> saddles2;
			std::vector<INDEX_TYPE> topo_index_partition;
			int num_threads;

#pragma omp parallel
			{
				// START PARALLEL CONSTRUCTION

				// divide index space
#pragma omp single
		{
			num_threads = omp_get_num_threads();
			ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, topo_index_partition);
		}
#pragma omp barrier


				// SCAN ALL CELLS COLLECT MINIMA, MAXIMA, 1SADDLES, 2SADDLES
				int thread_num = omp_get_thread_num();
				// in parallel go through and find all 2-saddles
				std::vector<INDEX_TYPE> lminima;
				std::vector<INDEX_TYPE> l1saddles;
				std::vector<INDEX_TYPE> lmaxima;
				std::vector<INDEX_TYPE> l2saddles;
//#pragma omp critical
//				{
//					printf("thread %d doing index %d-%d\n", thread_num, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
//				}
				typename MESH_TYPE::AllCellsIterator cellit(mMesh, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (cellit.begin(); cellit.valid(); cellit.advance()) {
					INDEX_TYPE cell_id = cellit.value();
					if (mGrad->getCritical(cell_id)) {
						DIM_TYPE d = mMesh->dimension(cell_id);
						if (do_mins) {
							if (d == 0) {

								lminima.push_back(cell_id);
							}
							else if (d == 1) {
								l1saddles.push_back(cell_id);
							}
						}
						if (do_maxs) {
							if (d == 2) {

								l2saddles.push_back(cell_id);
							}
							else if (d == 3) {
								lmaxima.push_back(cell_id);
							}
						}
					}
				}
#pragma omp barrier

				// COMBINE CRITICAL POINT LISTS
#pragma omp critical
				{
					//printf("gothere mins %d, saddles %d\n", lminima.size(), l1saddles.size());

					for (auto id : lminima) {
						mMinGraph->AddNode(id);
					}
					//printf("hothere\n");
					saddles1.insert(saddles1.end(), l1saddles.begin(), l1saddles.end() );
					//printf("asdfasdfasdfasdf\n");
				}
#pragma omp critical
				{
					//printf("gothere maxs\n");
					for (auto id : lmaxima) {
						mMaxGraph->AddNode(id);
					}
					saddles2.insert(saddles2.end(), l2saddles.begin(), l2saddles.end() );
				}
#pragma omp barrier

				int num1s = saddles1.size();
#pragma omp single
				{
					//printf("gothere2 \n");
					ArrayIndexPartitioner::EvenChunkSplit(num1s, num_threads, topo_index_partition);
				}
#pragma omp barrier
				thread_num = omp_get_thread_num();
				for (int index = topo_index_partition[thread_num]; index < topo_index_partition[thread_num + 1]; index++) {

					INDEX_TYPE sad1 = saddles1[index];
					INDEX_TYPE min1, min2;
					
					trace_down_1saddle(sad1, min1, min2);
					if (min1 != min2) {

#pragma omp critical
						{		
							//printf("adding %llu, %llu, %llu\n", sad1, min1, min2);
							mMinGraph->AddArc(min1, min2, sad1);
							//printf("added!\n");
						}
					}
				}
#pragma omp barrier

				int num2s = saddles2.size();
#pragma omp single
				{
					//printf("gothere3 \n");
					ArrayIndexPartitioner::EvenChunkSplit(num2s, num_threads, topo_index_partition);
				}
#pragma omp barrier
				thread_num = omp_get_thread_num();
				for (int index = topo_index_partition[thread_num]; index < topo_index_partition[thread_num + 1]; index++) {

					INDEX_TYPE sad2 = saddles2[index];
					INDEX_TYPE max1, max2;
					trace_up_2saddle(sad2, max1, max2);
					if (max1 != max2) {
#pragma omp critical
						{
							mMaxGraph->AddArc(max1, max2, sad2);
						}
					}
				}
#pragma omp barrier
#pragma omp single
				{
					printf("built minmaxgraph:\n mingraph = %d nodes, %d arcs\n maxgraph = %d nodes, %d arcs\n",
						this->mMinGraph->NumExtrema(), this->mMinGraph->mMergesToCancel.size(),
						this->mMaxGraph->NumExtrema(), this->mMaxGraph->mMergesToCancel.size());
				}
#pragma omp single
				{
					int numc = mMinGraph->SimplifyToThreshold(THRESHOLD);
					printf("MinGraph did %d cancellations\n", numc);
				}
#pragma omp single
				{
					int numc = mMaxGraph->SimplifyToThreshold(THRESHOLD);
					printf("MaxGraph did %d cancellations\n", numc);
				}

			} // END OMP PARALLEL

		}
		// phase 0: (build nodes)
		// -- gather extrema
		// -- for each extremum, set itself as the representative_list_id
		// Phase 1: (build directed graph)
		// -- for each saddle s, take less extreme node nl, if |v(s) - v(nl)| < t
		// where t is max simp thresh, if s is "closer" to nl than nl's current saddle
		// replace nl's current saddle with s
		// Phase 2: (union-find representatives)

		//need map<INDEX_TYPE> -> int to find extrema when tracing paths
		//need struct arc { saddle_value; more_extreme_node_id; }


	};





















	//typedef IndexCompareLessThan Comparer;
	template< class Comparer>
	class GridExtremalRegionBuilder {

	protected:
		DenseLabeling<int>* m_dest_label;
		RegularGrid3D* m_grid;
		RegularGridTrilinearFunction* m_func;

		Vec3i m_xyz;
		Vec3b m_periodic;
		inline bool AComesBeforeB(INDEX_TYPE a, INDEX_TYPE b) const {
			return mCompare->Compare(a, b);
		}

		bool IsExtremeVertexIn6Neighborhood(INDEX_TYPE id) const {
			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->XYZ3d(id);
			int t_num_neighbors = m_grid->GatherExistingNeighborsSameBdry6(t_coords, t_neighbors);

			INDEX_TYPE t_current_lowest = id;
			for (int i = 0; i < t_num_neighbors; i++) {
				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
				if (AComesBeforeB(t_neighbor_vertex, t_current_lowest)) {
					return false;
				}
			}
			return true;
		}
		void Enqueue_Later_Neighbors(Vec3l xyz, std::priority_queue<INDEX_TYPE, std::vector<INDEX_TYPE>, Comparer > &expansion, std::set<INDEX_TYPE>&enqueued_set) {
			INDEX_TYPE currentVID = m_grid->Index3d(xyz);

			Vec3l neighbors[6];
			int num_neighbors = m_grid->GatherExistingNeighborsSameBdry6(xyz, neighbors);
			for (int i = 0; i < num_neighbors; i++) {
				INDEX_TYPE neighborVID = m_grid->Index3d(neighbors[i]);
				if (m_dest_label->GetLabel(neighborVID) == -1 && AComesBeforeB(currentVID, neighborVID) && enqueued_set.count(neighborVID) == 0) {
					enqueued_set.insert(neighborVID);
					expansion.push(neighborVID);
				}

			}
		}

		// look at neighborhood of currentVID, for each neighbor, if it's "earlier"
		// check if it has been assigned- if not, then this point is unassigned
		// if it is assigned, check that each point has same assignment,
		// if there are no neighbors, return its original label
		int InspectPriorRegions(INDEX_TYPE currentVID) {
			INDEX_TYPE neighborVID;
			int extremal_certain = m_dest_label->GetLabel(currentVID);
			bool has_extremal = false;

			Vec3l neighbors[6];
			int num_neighbors = m_grid->GatherExistingNeighborsSameBdry6(m_grid->XYZ3d(currentVID), neighbors);
			for (int i = 0; i < num_neighbors; i++) {
				INDEX_TYPE neighborVID = m_grid->Index3d(neighbors[i]);
				if (AComesBeforeB(neighborVID, currentVID)) {
					if (m_dest_label->GetLabel(neighborVID) == -1) return -1; // if a extremal one is uncertain, we are uncertain
					if (!has_extremal) {
						extremal_certain = m_dest_label->GetLabel(neighborVID);
						has_extremal = true;
					}
					else {
						if (extremal_certain != m_dest_label->GetLabel(neighborVID)) return -1;
					}
				}
			}

			//if (!has_extremal) {
			//	printf("ERROR should never get here2\n");
			//	Vec3l coords = m_grid->XYZ3d(currentVID); coords.PrintInt();

			//	return -1;
			//}
			return extremal_certain;

		}

		void Expand_Lower_Neighborhood(INDEX_TYPE startid, int start_label) {
			Vec3l xyz = m_grid->XYZ3d(startid);
			std::set<INDEX_TYPE> enqueued_set;

			INDEX_TYPE currentVID = startid;
			// the natural ordering using the < operator on pairs will give us the highest
			// element first, simulating region growing from high to low
			std::priority_queue<INDEX_TYPE, std::vector<INDEX_TYPE>, Comparer > growing_front(*mCompare);
			enqueued_set.insert(startid);
			m_dest_label->SetLabel(startid, start_label);

			Enqueue_Later_Neighbors(xyz, growing_front, enqueued_set);

			while (!growing_front.empty()) {

				INDEX_TYPE currid = growing_front.top();
				growing_front.pop();

				int cellvale = InspectPriorRegions(currid);
				// find extremals

				// cellvalue >=0 indicates that there is certainty here, so lets expand
				if (cellvale >= 0) {
					m_dest_label->SetLabel(currid, cellvale);
					Enqueue_Later_Neighbors(m_grid->XYZ3d(currid), growing_front, enqueued_set);
				}

			}
		}


		Comparer* mCompare;
		std::vector<INDEX_TYPE> mExtrema;

	public:

		GridExtremalRegionBuilder(RegularGridTrilinearFunction* func, RegularGrid3D* grid) :
			m_xyz(func->GetGrid()->XYZ()), m_periodic(func->GetGrid()->Periodic()), m_func(func), m_grid(grid) {
			mCompare = new Comparer(func);
		}

		DenseLabeling<int>* GetOutputLabels() { return m_dest_label; }
		RegularGrid3D* GetGrid() { return m_grid; }
		RegularGridTrilinearFunction* GetFunction() { return m_func; }
		const std::vector<INDEX_TYPE>& GetExtrema() const { return mExtrema; }


		void BeginIntegration(bool verbose = false) {

			ThreadedTimer gtimer(1);
			gtimer.StartGlobal();

			m_dest_label = new DenseLabeling<int>(m_grid->NumElements());

			const INDEX_TYPE t_num_vertices = m_grid->NumElements();

#pragma omp parallel for
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				m_dest_label->SetLabel(i, -1);

				if (IsExtremeVertexIn6Neighborhood(i)) {
#pragma omp critical
					{
						mExtrema.push_back(i);
					}
				}
			}

			int num_extrema = mExtrema.size();
			printf("Found and expanding %d regions\n", num_extrema);
#pragma omp parallel 
			{
#pragma omp for schedule(dynamic)  nowait
				for (int m = 0; m < num_extrema; m++) {
					INDEX_TYPE maximum = mExtrema[m];
					Expand_Lower_Neighborhood(maximum, m);
				}
			}
		}


	};


//#define DEBUGGSERB

	//typedef IndexCompareLessThan Comparer;
	template< class Comparer, class GridFuncType, class TopoGridType>
	class GridSimplifiedExtremalRegionBuilder {

	protected:
		DenseLabeling<int>* m_dest_label;
		RegularGrid3D* m_grid;
		GridFuncType* m_func;

		inline bool AComesBeforeB(INDEX_TYPE a, INDEX_TYPE b) const {
			return mCompare->Compare(a, b);
		}

		bool IsExtremeVertexIn6Neighborhood(INDEX_TYPE id) const {
			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->XYZ3d(id);
			int t_num_neighbors = m_grid->GatherExistingNeighborsSameBdry6(t_coords, t_neighbors);

			INDEX_TYPE t_current_lowest = id;
			for (int i = 0; i < t_num_neighbors; i++) {
				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
				if (AComesBeforeB(t_neighbor_vertex, t_current_lowest)) {
					return false;
				}
			}
			return true;
		}
		void Enqueue_Later_Neighbors(Vec3l xyz, std::priority_queue<INDEX_TYPE, std::vector<INDEX_TYPE>, Comparer > &expansion, std::set<INDEX_TYPE>&enqueued_set) const {
			INDEX_TYPE currentVID = m_grid->Index3d(xyz);

			Vec3l neighbors[6];
			int num_neighbors = m_grid->GatherExistingNeighborsAll6(xyz, neighbors); // we want to grow towards the middle?
			for (int i = 0; i < num_neighbors; i++) {
				INDEX_TYPE neighborVID = m_grid->Index3d(neighbors[i]);
				if (m_dest_label->GetLabel(neighborVID) == -1 && 
					AComesBeforeB(currentVID, neighborVID) && 
					enqueued_set.count(neighborVID) == 0) {
#ifdef 	DEBUGGSERB		
					testout->SetLabel(neighborVID, 1);
#endif
					enqueued_set.insert(neighborVID);
					expansion.push(neighborVID);
				}

			}
		}

		// look at neighborhood of currentVID, for each neighbor, if it's "earlier"
		// check if it has been assigned- if not, then this point is unassigned
		// if it is assigned, check that each point has same assignment,
		// if there are no neighbors, return its original label
#ifdef 	DEBUGGSERB		
		int InspectPriorRegions(INDEX_TYPE currentVID, std::set<INDEX_TYPE>& enqueued_set)
#else
		int InspectPriorRegions(INDEX_TYPE currentVID)
#endif
		{
			INDEX_TYPE neighborVID;
			int extremal_certain = m_dest_label->GetLabel(currentVID);
			bool has_extremal = false;

			Vec3l neighbors[6];
			int num_neighbors = m_grid->GatherExistingNeighborsAll6(m_grid->XYZ3d(currentVID), neighbors);
			for (int i = 0; i < num_neighbors; i++) {
				INDEX_TYPE neighborVID = m_grid->Index3d(neighbors[i]);
				if (AComesBeforeB(neighborVID, currentVID)) {
					if (m_dest_label->GetLabel(neighborVID) == -1) {
#ifdef 	DEBUGGSERB		
						//printf("%f 's before neighbor %f has -1 label : c%d, n%d\n", 
						//	m_func->SampleImage(currentVID), 
						//	m_func->SampleImage(neighborVID), 
						//	enqueued_set.count(currentVID), 
						//	enqueued_set.count(neighborVID));
						//auto v1 = m_grid->XYZ3d(currentVID);
						//auto v2 = m_grid->XYZ3d(neighborVID);
						//v1.PrintInt();
						//v2.PrintInt();
						
#endif
						return -1; // if a extremal one is uncertain, we are uncertain
					}
					if (!has_extremal) {
						extremal_certain = m_dest_label->GetLabel(neighborVID);
						has_extremal = true;
					}
					else {
						if (extremal_certain != m_dest_label->GetLabel(neighborVID)) {
						

							return -1;

						}
					}
				}
			}

			//if (!has_extremal) {
			//	printf("ERROR should never get here2\n");
			//	Vec3l coords = m_grid->XYZ3d(currentVID); coords.PrintInt();

			//	return -1;
			//}

			return extremal_certain;

		}

#ifdef DEBUGGSERB
		DenseLabeling<int>* testout;
#endif

		void Expand_Lower_Neighborhood(const std::vector<INDEX_TYPE>& startids, int start_label) {
			std::set<INDEX_TYPE> enqueued_set;




			std::priority_queue<INDEX_TYPE, std::vector<INDEX_TYPE>, Comparer > growing_front(*mCompare);
			for (auto startid : startids) {
				//if (! this->IsExtremeVertexIn6Neighborhood(startid)) {
				//	printf("NONMAX START ID\n");
				//}
				// the natural ordering using the < operator on pairs will give us the highest
				// element first, simulating region growing from high to low
				enqueued_set.insert(startid);
				m_dest_label->SetLabel(startid, start_label);
#ifdef DEBUGGSERB
				testout->SetLabel(startid, 2);
#endif

			}

#ifdef 	DEBUGGSERB		
			//int numvalid = 0;
			//for (INDEX_TYPE i = 0; i < m_grid->NumElements(); i++) {
			//	if (this->IsExtremeVertexIn6Neighborhood(i)) {
			//		if (m_dest_label->GetLabel(i) == -1) {
			//			printf("UNFOUND EXTREMA!! %llu\n", i);
			//		}
			//		numvalid++;// printf("extremum validated: %llu\n", i);
			//	}
			//	else {
			//		//growing_front.push(i);
			//		//enqueued_set.insert(i);
			//	}

			//}
			//printf("validated %d extrema\n", numvalid);

			//while (!growing_front.empty()) {
			//	INDEX_TYPE currid = growing_front.top();

			//	growing_front.pop();
			//	int cellvale = InspectPriorRegions(currid, enqueued_set);
			//	if (cellvale >= 0) {
			//		m_dest_label->SetLabel(currid, cellvale);
			//		testout->SetLabel(currid, 2);
			//	}
			//}
			//return;

#endif

			for (auto startid : startids) {

#ifdef 	DEBUGGSERB		
				INDEX_TYPE TESTID = InspectPriorRegions(startid, enqueued_set);
				if (TESTID == -1) printf("WHOATHERE %ll %llu\n", TESTID, startid);
#endif

				Vec3l xyz = m_grid->XYZ3d(startid);
				Enqueue_Later_Neighbors(xyz, growing_front, enqueued_set);
			}
#ifdef 	DEBUGGSERB		
			printf("growing front size %d\n", growing_front.size());
#endif

			while (!growing_front.empty()) {

				INDEX_TYPE currid = growing_front.top();

				growing_front.pop();
				//printf("asdf %f\n", m_func->SampleImage(currid));

#ifdef DEBUGGSERB
				int cellvale = InspectPriorRegions(currid, enqueued_set);
#else
				int cellvale = InspectPriorRegions(currid);
#endif
				// find extremals
				//continue;
				// cellvalue >=0 indicates that there is certainty here, so lets expand
				if (cellvale >= 0) {
					m_dest_label->SetLabel(currid, cellvale);
#ifdef DEBUGGSERB
					testout->SetLabel(currid, 2);
#endif
					Enqueue_Later_Neighbors(m_grid->XYZ3d(currid), growing_front, enqueued_set);
				}

			}
		}

      std::unordered_map<INT_TYPE, std::vector<INDEX_TYPE> > mIdToVIndexMap;

		Comparer* mCompare;
		//std::vector<INDEX_TYPE> mExtrema;
		TopoGridType* mMesh;
	public:

	
		GridSimplifiedExtremalRegionBuilder(GridFuncType* func, RegularGrid3D* grid, TopoGridType* tgrid) :
			m_func(func), m_grid(grid), mMesh(tgrid) {
			mCompare = new Comparer(func);
#ifdef DEBUGGSERB
			testout = new DenseLabeling<int>(grid->NumElements());
#endif
		}

		DenseLabeling<int>* GetOutputLabels() { return m_dest_label; }
		RegularGrid3D* GetGrid() { return m_grid; }
		GridFuncType* GetFunction() { return m_func; }
		const std::vector<INDEX_TYPE>& GetExtrema() const { return this->mExtrema; }
		const std::unordered_map<INT_TYPE, std::vector<INDEX_TYPE> >& GetIdMap() {
			return mIdToVIndexMap;
		}

		// the extremum map takes in the GID of a VERTEX (actually any cell but it gets the vertex number by just
		// rounding down. to get the right vertex in a cell, pass in a map which uses max labeling to replace each
		// cell with its highest vertex.
		void BeginIntegration(const std::unordered_map<INDEX_TYPE, INT_TYPE>& extremaGIDtoLabelMap, bool verbose = false) {


#ifdef 	DEBUGGSERB		
			printf("gothere1 %d\n", extremaGIDtoLabelMap.size());
			mCompare->PrintRule();
			testout->SetAll(0);
#endif
			ThreadedTimer gtimer(1);
			gtimer.StartGlobal();

			m_dest_label = new DenseLabeling<int>(m_grid->NumElements());
			mIdToVIndexMap.clear();
			const INDEX_TYPE t_num_vertices = m_grid->NumElements();
			m_dest_label->SetAll(-1);
			printf("starting expansion\n");
//#pragma omp parallel for
//			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
//				m_dest_label->SetLabel(i, -1);
//
//				if (IsExtremeVertexIn6Neighborhood(i)) {
//#pragma omp critical
//					{
//						mExtrema.push_back(i);
//					}
//				}
//			}
#ifdef 	DEBUGGSERB		
			printf("gothere2\n");
#endif
         std::vector<INT_TYPE> uniqueids;
			for (auto p : extremaGIDtoLabelMap) {
				if (mIdToVIndexMap.count(p.second) == 0) uniqueids.push_back(p.second);
				mIdToVIndexMap[p.second].push_back(mMesh->VertexNumberFromCellID(p.first));
			}
#ifdef 	DEBUGGSERB		
			for (auto vv : mIdToVIndexMap) {
				printf("%d: %d ids\n", vv.first, vv.second.size());
			}
#endif
			//for (auto m : mExtrema) {
			//	INDEX_TYPE tgridid = mMesh->CellIDFromVertexNumber(m);
			//	if (extremaGIDtoLabelMap.count(tgridid) == 0) {
			//		printf("WHOA THERE EXTMAP DOES NOT CONTAIN CP\n");
			//	}
			//}

			int num_unique = uniqueids.size();
			printf("Found and expanding %d regions\n", num_unique);
#pragma omp parallel 
			{
#pragma omp for schedule(dynamic)  nowait
				for (int m = 0; m < num_unique; m++) {
					Expand_Lower_Neighborhood(mIdToVIndexMap[uniqueids[m]], m);
				}
			}
#ifdef 	DEBUGGSERB		
			printf("gothere3\n");
			testout->OutputToIntFile("testoutasdfasdfasdf.raw");
#endif
		}


	};

}
#endif
