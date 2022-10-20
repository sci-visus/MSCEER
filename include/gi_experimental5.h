#ifndef GI_EXPERIMENTAL5_H
#define GI_EXPERIMENTAL5_H

#include <vector>
#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_labeling.h"
#include "gi_regular_grid.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_topological_regular_grid.h"
#include "gi_discrete_gradient_labeling.h"
#include "gi_topological_utility_functions.h"
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
#include "gi_numeric_integrator_expanding_region_stop.h"


//#include "concurrentqueue.h"
#include <thread>
#include <mutex>

#include "gi_experimental.h"
#include "gi_experimental3b.h"
#include "gi_convergent_2dmanifold_gradient_builder.h"

//#define OUTPUTINTERMEDIATE

// This file is a safe space to put experimental classes and functionality.

namespace GInt {



	template< class Advector, class Comparer>
	class NumericIntegratorExpandingRegionStopWithCutoff {

	protected:

		//struct bridge {
		//	INDEX_TYPE labela;
		//	INDEX_TYPE labelb;
		//	float value;
		//};
		//

		DenseLabeling<INDEX_TYPE>* m_destinations;
		DenseLabeling<int>* m_certains;
		RegularGrid3D* m_grid;
		RegularGridTrilinearFunction* m_func;
		int m_num_iterations_left;

		Vec3i m_xyz;
		Vec3b m_periodic;
		FLOATTYPE m_error_threshold;
		FLOATTYPE m_gradient_threshold;

		char* m_fname;

		//struct cmp{
		//	bool operator()(INDEX_TYPE a, INDEX_TYPE b) {
		//		return mCompare->Compare(a, b);
		//	}
		//};

		bool IsExtremeVertexIn6Neighborhood(INDEX_TYPE id) const {
			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->Coords(id);
			int t_num_neighbors = m_grid->GatherExistingNeighborsSameBdry6(t_coords, t_neighbors);

			INDEX_TYPE t_current_lowest = id;
			for (int i = 0; i < t_num_neighbors; i++) {
				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
				if (mCompare->Compare(t_neighbor_vertex, t_current_lowest)) {
					return false;
				}
			}
			return true;
		}
		void Enqueue_Later_Neighbors(Vec3l xyz, std::priority_queue<INDEX_TYPE, std::vector<INDEX_TYPE>, Comparer > &expansion, std::set<INDEX_TYPE>&seen) {
			INDEX_TYPE tid = m_grid->Index3d(xyz);

			Vec3l neighbors[6];
			int nn = m_grid->GatherExistingNeighborsSameBdry6(xyz, neighbors);
			for (int i = 0; i < nn; i++) {
				INDEX_TYPE tneg = m_grid->Index3d(neighbors[i]);
				if (m_certains->GetLabel(tneg) == -1 && mCompare->Compare(tid, tneg) && seen.count(tneg) == 0) {
					seen.insert(tneg);
					expansion.push(tneg);
				}

			}
		}

		int Inspect_Higher_Certains(INDEX_TYPE tid) {
			INDEX_TYPE tneg;
			int extremal_certain = m_certains->GetLabel(tid);
			bool has_extremal = false;

			Vec3l neighbors[6];
			int nn = m_grid->GatherExistingNeighborsSameBdry6(m_grid->Coords(tid), neighbors);
			for (int i = 0; i < nn; i++) {
				INDEX_TYPE tneg = m_grid->Index3d(neighbors[i]);
				if (mCompare->Compare(tneg, tid)) {
					if (m_certains->GetLabel(tneg) < 0) return -1; // if a extremal one is uncertain, we are uncertain
					if (!has_extremal) {
						extremal_certain = m_certains->GetLabel(tneg);
						has_extremal = true;
					}
					else {
						if (extremal_certain != m_certains->GetLabel(tneg)) return -1;
					}
				}
			}

			if (!has_extremal) {
				printf("ERROR should never get here\n");
				return -1;
			}
			return extremal_certain;

		}

		void Expand_Lower_Neighborhood(INDEX_TYPE startid) {
			Vec3l xyz = m_grid->Coords(startid);
			std::set<INDEX_TYPE> seen;

			INDEX_TYPE tid = startid;
			// the natural ordering using the < operator on pairs will give us the highest
			// element first, simulating region growing from high to low
			std::priority_queue<INDEX_TYPE, std::vector<INDEX_TYPE>, Comparer > growing_front(*mCompare);
			seen.insert(startid);
			Enqueue_Later_Neighbors(xyz, growing_front, seen);

			while (!growing_front.empty()) {

				INDEX_TYPE currid = growing_front.top();
				growing_front.pop();

				int cellvale = Inspect_Higher_Certains(currid);
				// find extremals

				// cellvalue >=0 indicates that there is certainty here, so lets expand 
				if (cellvale >= 0) {
					m_certains->SetLabel(currid, cellvale);
					m_destinations->SetLabel(currid, startid);
					Enqueue_Later_Neighbors(m_grid->Coords(currid), growing_front, seen);
				}

			}
		}


		Comparer* mCompare;
		std::unordered_set<INDEX_TYPE> m_extrema;
	public:

		FLOATTYPE mCutoff;

		NumericIntegratorExpandingRegionStopWithCutoff(RegularGridTrilinearFunction* func, RegularGrid3D* grid, FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int interation_limit, FLOATTYPE cutoff) :
			mCutoff(cutoff), m_num_iterations_left(interation_limit), m_xyz(func->GetGrid()->XYZ()), m_periodic(func->GetGrid()->Periodic()), m_func(func), m_grid(grid),
			m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold) {
			mCompare = new Comparer(func);
		}

		DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }
		RegularGrid3D* GetGrid() { return m_grid; }
		RegularGridTrilinearFunction* GetFunction() { return m_func; }
		std::unordered_set<INDEX_TYPE>& GetExtrema() { return m_extrema; }
		void BeginIntegration() {

			printf(" -- Performing numeric integration for volume assignment...\n");

			//m_func->ComputeGradFromImage(m_rkindex);

			m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
			m_certains = new DenseLabeling<int>(m_grid->NumElements());
			const INDEX_TYPE t_num_vertices = m_grid->NumElements();
			AdvectionChecker* inside_voxel_critical_advection_checker = new TerminateNearCertain(m_certains, m_grid);
			AdvectionChecker* no_check = new NoTermination();//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);

			printf("   -- finding extrema to terminate integral lines...");
			fflush(stdout);

			std::vector<INDEX_TYPE> extrema;

			// set all potential extrema, so we terminate near them
#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				
				if (mCompare->Compare(m_func->SampleImage(i), mCutoff)) {

					m_certains->SetLabel(i, t_num_vertices + 1);
					m_destinations->SetLabel(i, i);
				}
				else {
					m_certains->SetLabel(i, -1);
					if (IsExtremeVertexIn6Neighborhood(i)) {
						m_destinations->SetLabel(i, i);
#pragma omp critical
						{
							extrema.push_back(i);
							m_extrema.insert(i);
						}
					}
					else {
						m_destinations->SetLabel(i, -1);
					}
				}

			}
			printf(" found %d extrema!\n", extrema.size());
			printf("   -- expanding extrema certain regions...");
			fflush(stdout);

			int num_extrema = extrema.size();
#pragma omp parallel shared(extrema) 
			{
#pragma omp for schedule(dynamic)  nowait 
				for (int m = 0; m < num_extrema; m++) {
					INDEX_TYPE maximum = extrema[m];
					m_certains->SetLabel(maximum, m);
					Expand_Lower_Neighborhood(maximum);
					m_func->SetGradExplicit(maximum, Vec3d(0, 0, 0));
				}
			}

			printf(" expansion done!\n", extrema.size());
			printf("   -- doing numerical integration first pass with path compression...");
			fflush(stdout);

			//m_destinations->OutputToIntFile("dest.raw");
			int t1, t2; t1 = t2 = 0;
#pragma omp parallel 
			{
				Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_critical_advection_checker);
				std::vector<INDEX_TYPE> t_path;
				t_path.reserve(100);
#pragma omp for schedule(guided)  nowait 
				for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
					// early skip if this is already a maximum
					if (m_certains->GetLabel(i) != -1 || m_destinations->GetLabel(i) != -1) {
						continue;
					}

					t_path.clear();
					t_path.push_back(i);
					Vec3l t_coords = m_grid->Coords(i); // get the coordinates of the poitn
					//if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
					Vec3d t_current_point = t_coords;
					int t_num_iterations_left = m_num_iterations_left;

					bool t_continue = true;

					while (t_continue) {
						Vec3d t_next_point;
						ADVECTION_EVENT t_return_code;
						if (m_grid->DistToBoundary(t_coords) <= 1) {
							t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
							t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
							t1++;
						}
						else {
							t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
							t_coords = (t_current_point + 0.5);
							t2++;
						}
						INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
						t_path.push_back(t_next_id);
						// if we terminated or hit a critical point, then we are done
						if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
							t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
							t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
							t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
							INDEX_TYPE t_dest_label = m_destinations->GetLabel(t_next_id);
							int t_certain_label = m_certains->GetLabel(t_next_id);
							//#pragma omp critical
							{
								for (int j = 0; j < t_path.size(); j++) {
									m_destinations->SetLabel(t_path[j], t_dest_label);
									//m_certains->SetLabel(t_path[j], t_certain_label);
								}
							}
							t_continue = false;
						}
					}
				}
#pragma omp single
			{
				printf(" done!\n");
				printf("   -- checking unambiguous voxels...");
				fflush(stdout);
			}
#pragma omp for schedule(static)  
				for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
					// early skip if this is already a maximum
					if (m_certains->GetLabel(i) != -1) {
						continue;
					}

					Vec3l t_coords = m_grid->Coords(i); // get the coordinates of the poitn
					//if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
					Vec3l negs[6];

					int nn = m_grid->GatherExistingNeighborsSameBdry6(t_coords, negs);

					bool needsupdate = false;
					INDEX_TYPE t_dest_label = m_destinations->GetLabel(i);
					for (int j = 0; j < nn; j++)  {
						if (m_destinations->GetLabel(m_grid->Index3d(negs[j])) != t_dest_label) {
							needsupdate = true;
							break;
						}
					}

					if (needsupdate) continue;

					m_certains->SetLabel(i, 1);
				}
#pragma omp single
			{
				printf(" done!\n");
				printf("   -- re-integrating boundaries that may have ended too soon...");
				fflush(stdout);
			}
#pragma omp for schedule(guided)  nowait 
				for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
					// early skip if this is already a maximum
					if (m_certains->GetLabel(i) != -1) {
						continue;
					}

					Vec3l t_coords = m_grid->Coords(i); // get the coordinates of the poitn
					//if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);



					Vec3d t_current_point = t_coords;
					int t_num_iterations_left = 500;
					bool t_continue = true;

					while (t_continue) {
						Vec3d t_next_point;
						ADVECTION_EVENT t_return_code;
						if (m_grid->DistToBoundary(t_coords) <= 1) {
							t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
							t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
							t1++;
						}
						else {
							t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
							t_coords = (t_current_point + 0.5);
							t2++;
						}
						INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
						t_path.push_back(t_next_id);
						// if we terminated or hit a critical point, then we are done
						if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
							t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
							t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
							t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
							m_destinations->SetLabel(i, m_destinations->GetLabel(t_next_id));
							t_continue = false;
						}
					}
				}

				//if (m_grid->DistToBoundary()) {}
				//if (IsLowestVertexIn6Neighborhood(i)) {
				//	m_destinations->SetLabel(i, i);
				//}
				//else {
				//	m_destinations->SetLabel(i, -1);
				//}



				//printf("%d boundary, %d noboundary\n", t1, t2);
			}
			printf(" done!\n");
			printf(" -- Done numeric integration of 3d volume!\n");

		}



	};


	template <class PayloadType, class DependencyGraphType>
	class myCustomMap2 {
	public:
		PayloadType** mItems;
		std::unordered_map<INDEX_TYPE, INDEX_TYPE> mMap;
		PayloadType defitem;
		myCustomMap2(DependencyGraphType* graph){
			mItems = new PayloadType*[graph->NumIds()];
			INDEX_TYPE item = 0;
			typename DependencyGraphType::vertex_iterator it(graph);
			for (it.begin(); it.valid(); it.advance()) {
				INDEX_TYPE id = it.value();
				mMap[id] = item++;
			}
			memset(mItems, 0, sizeof(PayloadType*)* graph->NumIds());

		}
		~myCustomMap2(){
			delete[] mItems;
		}

		INDEX_TYPE count(INDEX_TYPE id) {
			return mItems[mMap[id]] != NULL;
		}

		PayloadType& operator[](INDEX_TYPE id) {
			if (count(id) == 0) {
				printf("whoa\n");
				return defitem;
			}
			return  *(mItems[mMap[id]]);
		}

		void insert(INDEX_TYPE id, PayloadType& p) {
			mItems[mMap[id]] = new PayloadType(p);
		}
		void erase(INDEX_TYPE id) {
			if (count(id) != 0) {
				delete mItems[mMap[id]];
				mItems[mMap[id]] = NULL;
			}
			//printf("re-erasing\n");
		}
	};

//	template<class PayloadType>
//	class myCustomMap {
//
//	public:
//		omp_lock_t* mLocks;
//		std::unordered_map<INDEX_TYPE, PayloadType>* mMaps;
//		myCustomMap(INDEX_TYPE nume){
//			INDEX_TYPE mysize = (nume >> 10) + 1;
//			mLocks = new omp_lock_t[mysize];
//			mMaps = new std::unordered_map<INDEX_TYPE, PayloadType>[mysize];
//#pragma omp parallel for
//			for (int i = 0; i < mysize; i++){
//				omp_init_lock(&(mLocks[i]));
//			}
//		}
//		~myCustomMap(){
//			delete[] mLocks;
//			delete[] mMaps;
//		}
//
//		INDEX_TYPE count(INDEX_TYPE id) {
//			return mMaps[id >> 10].count(id);
//		}
//
//		PayloadType& operator[](INDEX_TYPE id) {
//			return mMaps[id >> 10][id];
//		}
//
//		void insert(INDEX_TYPE id, PayloadType& p) {
//			INDEX_TYPE lid = id >> 10;
//			omp_set_lock(&(mLocks[lid]));
//			mMaps[lid][id] = p;
//			omp_unset_lock(&(mLocks[lid]));
//		}
//		void erase(INDEX_TYPE id) {
//			INDEX_TYPE lid = id >> 10;
//			omp_set_lock(&(mLocks[lid]));
//			mMaps[lid].erase(id);
//			omp_unset_lock(&(mLocks[lid]));
//
//		}
//	};
	template <class PayloadType, class DependencyGraphType>
	class TopologicalGraphExpandingRegionsBTiming {
	public:
		//typedef BoundaryEdgeGraph DependencyGraphType;
		//typedef int PayloadType;

		// VARIABLES AND TYPE DECLARATIONS
		typedef std::unordered_set<INDEX_TYPE> index_set_type;
		typedef std::unordered_map<INDEX_TYPE, PayloadType> payload_map_type;
		DependencyGraphType* mD;
		DenseLabeling<char>* mGloballyProcessedLabel; // debug
		DenseLabeling<char>* mIndepTempProcLabel; // for now use int because it's atomic - dont know if byte is...


		class IdComparer {
		public:
			DependencyGraphType* dep;
			IdComparer(DependencyGraphType* g) : dep(g) {}

			bool operator()(INDEX_TYPE a, INDEX_TYPE b){
				return dep->Before(b, a);
			}
		};
		IdComparer m_id_comparer;
		typedef std::priority_queue<INDEX_TYPE, vector<INDEX_TYPE>, IdComparer > SORTED_ID_QUEUE_TYPE;
		struct RegionState {
			RegionState(IdComparer comp) : unprocessed_band(comp) {}
			INDEX_TYPE region_id;
			INDEX_TYPE bottleneck_id;
			//SORTEDPAIRSQUEUE unprocessed_band;
			SORTED_ID_QUEUE_TYPE unprocessed_band;
			index_set_type seen_set;
			index_set_type processed_set;
			payload_map_type prior_data;

		};
		std::map<INDEX_TYPE, RegionState*> m_bottlenecks_map;
		myCustomMap2< PayloadType, DependencyGraphType> m_prior_data;

		enum WorkType { MERGE_AND_EXPAND, EXPAND_ONLY, PREP_FOR_GLOBAL, NONE };
		struct WorkRequest {
			WorkType worktype;
			INDEX_TYPE start_cell;
			RegionState* region1;
			RegionState* region2;
		};


		bool IsProcessed(INDEX_TYPE id) const { return mGloballyProcessedLabel->GetLabel(id) != 0; }
		void SetProcessed(INDEX_TYPE id, char value) { mGloballyProcessedLabel->SetLabel(id, value); }
		inline bool AIsGreaterThanB(INDEX_TYPE a, INDEX_TYPE b) const {
			return mD->Before(b, a);
		}

		// ONLY CALL THIS FROM THREAD SAFE LOCATION, i.e. while other threads are NOT assigning ids
		bool IsUNASSIGNEDMinInNeighborhood(INDEX_TYPE id) const {
			// return false if the edge is not part of our permissable set
			if (IsProcessed(id)) {
				return false;
			}
			typename DependencyGraphType::neighbor_iterator nit(mD);
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (IsProcessed(neighbor_id)) continue;
				// else compare the values;
				if (AIsGreaterThanB(id, neighbor_id)) return false;
			}
			return true;
		}

		
		void Enqueue_Higher_Neighbors(INDEX_TYPE id,
			/*FLOATTYPE t_original_edge_value,*/
			/*SORTEDPAIRSQUEUE &expansion,*/ SORTED_ID_QUEUE_TYPE& expansion,
			index_set_type&seen_set) const {

			// to get surrounding edges, look at surrounding quads, and all their edges
			typename DependencyGraphType::neighbor_iterator nit(mD);
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				// make sure it is permissable and not already in our traversal
				if (!IsProcessed(neighbor_id) && // inserted --
					seen_set.count(neighbor_id) == 0) {
					// i think we can ignore the value test? - no we can't don't want to "spill over"
					// into another basin - only insert strictly lower stuff!
					//FLOATTYPE t_neg_value = m_func->cellValue(neighbor_id);
					if (AIsGreaterThanB(neighbor_id, id)) {
						seen_set.insert(neighbor_id);
						//expansion.push(FIPAIR(t_neg_value, neighbor_id));
						expansion.push(neighbor_id);
					}
				}

			}
		}

		// returns true if all my PRIOR neighbors (in this case Lower because of low to high traversal)
		// are either in the GLOBALLY processed set or in my LOCALLY processed set
		// return false if there exists a neighbor that has lower value but is not globally or locally processed!!
		// this means that it will be an obstruction point
		//bool AllLowerNeighborsAlreadyProcessedLocallyOrGlobally(INDEX_TYPE primal_edge_topo_id, FLOATTYPE original_edge_value, index_set_type& processed_set) const {
		bool AllLowerNeighborsAlreadyProcessedLocallyOrGlobally(INDEX_TYPE primal_edge_topo_id, index_set_type& processed_set) const {
			typename DependencyGraphType::neighbor_iterator nit(mD);
			for (nit.begin(primal_edge_topo_id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (IsProcessed(neighbor_id))
					continue;
				// if our neighbor is lower, check it is processed for sure
				if (AIsGreaterThanB(primal_edge_topo_id, neighbor_id)) {
					if (processed_set.count(neighbor_id) == 0) {
						// then its not certain so we can't be certain
						return false;
					}
				}

			}
			return true;
		}

		bool HasNoLowerNeighbors(INDEX_TYPE id) {
			typename DependencyGraphType::neighbor_iterator nit(mD);
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (AIsGreaterThanB(id, neighbor_id)) return false;
			}
			return true;
		}

		bool AllNeighborsAlreadyProcessedGlobally(INDEX_TYPE primal_edge_topo_id) const {
			typename DependencyGraphType::neighbor_iterator nit(mD);
			for (nit.begin(primal_edge_topo_id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (mIndepTempProcLabel->GetLabel(neighbor_id) == 0)
				return false;
			}
			return true;
		}
		bool AllHigherNeighborsAlreadyProcessedLocallyOrGlobally(INDEX_TYPE primal_edge_topo_id, index_set_type& processed_set) const {
			typename DependencyGraphType::neighbor_iterator nit(mD);
			for (nit.begin(primal_edge_topo_id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (IsProcessed(neighbor_id))
					continue;
				if (AIsGreaterThanB(neighbor_id, primal_edge_topo_id)) {
					if (processed_set.count(neighbor_id) == 0) {
						return false;
					}
				}

			}
			return true;
		}
		RegionState* MergeRegions(RegionState* a, RegionState* b) const {
			RegionState* aa;
			RegionState* bb;
			if (a->seen_set.size() > b->seen_set.size()) {
				aa = a;
				bb = b;
			}
			else  {
				aa = b;
				bb = a;
			}
			// insert unprocessed band first, so we can check duplicates
			while (!bb->unprocessed_band.empty()) {
				if (aa->seen_set.count(bb->unprocessed_band.top()) == 0) aa->unprocessed_band.push(bb->unprocessed_band.top());
				bb->unprocessed_band.pop();

			}
			aa->processed_set.insert(bb->processed_set.begin(), bb->processed_set.end());
			aa->prior_data.insert(bb->prior_data.begin(), bb->prior_data.end()); // here could be optimized to be sparse
			aa->seen_set.insert(bb->seen_set.begin(), bb->seen_set.end());
			//printf("tops = aa %d bb %d\n", aa->unprocessed_band.top(), bb->unprocessed_band.top());


			aa->region_id = (aa->region_id > bb->region_id ? aa->region_id + 1 : bb->region_id + 1);
			delete bb;
			//printf("regionid = %d\n", a->region_id);
			return aa;
		}

		//void Sparsify_Region(RegionState* region, unordered_set<INDEX_TYPE>& lp) const {
		//	for (auto it = lp.begin(); it != lp.end(); it++) {
		//		INDEX_TYPE id = *it;
		//		if (region->prior_data.count(id) == 0) continue;
		//		bool all_local = true;
		//		typename DependencyGraphType::neighbor_iterator nit(mD);
		//		for (nit.begin(id); nit.valid(); nit.advance()) {
		//			INDEX_TYPE neighbor_id = nit.value();
		//			if (region->processed_set.count(neighbor_id) == 0 && !IsProcessed(neighbor_id)){
		//				all_local = false;
		//				break;
		//			}
		//		}
		//		if (all_local) {
		//			region->prior_data.erase(id);
		//		}

		//	}
		//}

		RegionState* Expand_MergeRegions(RegionState* region)  {

			RegionState* local_region = region;

			//INDEX_TYPE primal_edge_topo_id = region->unprocessed_band.top().second;
			INDEX_TYPE primal_edge_topo_id = region->unprocessed_band.top();
			// the natural ordering using the < operator on pairs will give us the highest
			// element first, simulating region growing from high to low
			//unordered_set<INDEX_TYPE> t_locally_processed;
			//SORTEDPAIRSQUEUE growing_front(region->unprocessed_band);
			//region->unprocessed_band = SORTEDPAIRSQUEUE(); // clear the old delayed front
			SORTED_ID_QUEUE_TYPE growing_front(region->unprocessed_band);
			region->unprocessed_band = SORTED_ID_QUEUE_TYPE(m_id_comparer); // clear the old delayed front
			//std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> >  new_delayed_front (growing_front);

			while (!growing_front.empty()) {

				//INDEX_TYPE currid = growing_front.top().second;
				INDEX_TYPE currid = growing_front.top();
				//FLOATTYPE currid_value = growing_front.top().first;
				growing_front.pop();

				// cellvalue >=0 indicates that there is certainty here, so lets expand 
				//if (AllLowerNeighborsAlreadyProcessedLocallyOrGlobally(currid, currid_value, local_region->processed_set)) {
				if (AllLowerNeighborsAlreadyProcessedLocallyOrGlobally(currid, local_region->processed_set)) {
					if (local_region->processed_set.count(currid) != 0) printf("reprocessing2 %d\n", currid);

					local_region->processed_set.insert(currid);
					//t_locally_processed.insert(currid);

					mIndepTempProcLabel->SetLabel(currid, omp_get_thread_num() + 1);
					m_thread_work_count[omp_get_thread_num()]++;
					// do work
					local_region->prior_data[currid] = DoWork(currid, local_region, 1);


					//Enqueue_Higher_Neighbors(currid, currid_value, growing_front, local_region->seen_set);
					Enqueue_Higher_Neighbors(currid, growing_front, local_region->seen_set);
				}
				else {
					//region->unprocessed_band.push(FIPAIR(currid_value, currid));
					region->unprocessed_band.push(currid);
				}

			}
			if (!local_region->unprocessed_band.empty()) {
				//local_region->bottleneck_id = region->unprocessed_band.top().second;
				local_region->bottleneck_id = region->unprocessed_band.top();
			}

			//Sparsify_Region(local_region, t_locally_processed);
			return local_region;
		}

		// does a flood fill upwards starting from base_id
		// -- uses the global processed as context
		// first creates a new local region to represent this local flood fill
		// 
		RegionState* ExpandRegionUpwardsFlood(INDEX_TYPE base_id)  {

			RegionState* local_region = new RegionState(m_id_comparer);

			local_region->processed_set.clear();
			INDEX_TYPE t_id = base_id;
			local_region->seen_set.insert(base_id);
			local_region->processed_set.insert(base_id);

			mIndepTempProcLabel->SetLabel(base_id, omp_get_thread_num() + 1); // debug?
			m_thread_work_count[omp_get_thread_num()]++;
			local_region->prior_data[base_id] = DoInitialWork(base_id, local_region); // this represents the meat of the flood fill approach - do the base case of the recurrence relation

			//FLOATTYPE original_edge_value = m_func->cellValue(base_id);
			//SORTEDPAIRSQUEUE growing_front;
			SORTED_ID_QUEUE_TYPE growing_front(m_id_comparer);
			//Enqueue_Higher_Neighbors(base_id, original_edge_value, growing_front, local_region->seen_set);
			Enqueue_Higher_Neighbors(base_id, growing_front, local_region->seen_set);

			while (!growing_front.empty()) {

				//INDEX_TYPE currid = growing_front.top().second;
				INDEX_TYPE currid = growing_front.top();
				//FLOATTYPE currid_value = growing_front.top().first;
				growing_front.pop();

				// cellvalue >=0 indicates that there is certainty here, so lets expand 
				//if (AllLowerNeighborsAlreadyProcessedLocallyOrGlobally(currid, currid_value, local_region->processed_set)) {
				if (AllLowerNeighborsAlreadyProcessedLocallyOrGlobally(currid, local_region->processed_set)) {
					if (local_region->processed_set.count(currid) != 0) printf("reprocessing %d\n", currid);
					local_region->processed_set.insert(currid);
					mIndepTempProcLabel->SetLabel(currid, omp_get_thread_num() + 1); // debug?
					m_thread_work_count[omp_get_thread_num()]++;
					/// DO COMBO WORK
					local_region->prior_data[currid] = DoWork(currid, local_region, 0); // this represents the meat of the flood fill approach - do the inductive case of the recurrence relation

					//Enqueue_Higher_Neighbors(currid, currid_value, growing_front, local_region->seen_set);
					Enqueue_Higher_Neighbors(currid, growing_front, local_region->seen_set);
				}
				else {
					//local_region->unprocessed_band.push(FIPAIR(currid_value, currid));
					local_region->unprocessed_band.push(currid);
				}

			}
			if (!local_region->unprocessed_band.empty()) {
				//local_region->bottleneck_id = local_region->unprocessed_band.top().second;
				local_region->bottleneck_id = local_region->unprocessed_band.top();
			}
			local_region->region_id = 1; // debug?
			//Sparsify_Region(local_region);
			return local_region;
		}

		TopologicalGraphExpandingRegionsBTiming(
			//TopologicalExplicitDenseMeshFunction* func,
			DependencyGraphType* graph) :
			//m_func(func),
			mD(graph),
			m_id_comparer(graph),
			m_prior_data(graph)
		{
		}
		~TopologicalGraphExpandingRegionsBTiming() {
			delete mIndepTempProcLabel;
			delete mGloballyProcessedLabel;
		}
		//TopologicalExplicitDenseMeshFunction* GetFunction() { return m_func; }



		virtual PayloadType DoWork(INDEX_TYPE id, RegionState* region, int asdf) = 0;
		//	//std::vector<INDEX_TYPE> negs;
		//	//int mval = 0;
		//	typename DependencyGraphType::neighbor_iterator nit(mD);
		//	for (nit.begin(id); nit.valid(); nit.advance()) {
		//		INDEX_TYPE neighbor_id = nit.value();
		//		if (AIsGreaterThanB(id, neighbor_id)) {
		//			//negs.push_back(neighbor_id);
		//			// check
		//			if (mIndepTempProcLabel->GetLabel(neighbor_id) == 0) printf("WHOathere\n");

		//			if (region->prior_data.count(neighbor_id) != 0) {
		//				int otherval = region->prior_data[neighbor_id];
		//			}
		//			else if (m_prior_data.count(neighbor_id) != 0) {
		//				int otherval = m_prior_data[neighbor_id];
		//			}
		//			else {
		//				printf("whoa, prior data not in local or global %d \n", asdf);
		//			}
		//			//if (AllHigherNeighborsAlreadyProcessedLocallyOrGlobally(neighbor_id, region->processed_set)) {
		//			//	region->prior_data.erase(neighbor_id);
		//			//}

		//		}
		//	}
		//	region->prior_data[id] = 1;
		//}

		virtual PayloadType DoInitialWork(INDEX_TYPE id, RegionState* region) = 0;
		//	region->prior_data[id] = 0;
		//}

		// deletes prior data that cannot be accessed by anyone who is not processed
		// -- super conservative for now, only removes data that is completely surrounded
		//    IN THIS REGION - ignores globally processed things
		void Sparsify_Region(RegionState* region) const {
			set<INDEX_TYPE> keys;
			int count = 0;
			for (auto it = region->prior_data.begin(); it != region->prior_data.end(); it++) keys.insert((*it).first);
			for (auto it = keys.begin(); it != keys.end(); it++) {
				//for (auto it = region->processed_set.begin(); it != region->processed_set.end(); it++) {

				INDEX_TYPE id = *it;
				if (region->prior_data.count(id) == 0) continue;
				bool all_local = true;
				typename DependencyGraphType::neighbor_iterator nit(mD);
				for (nit.begin(id); nit.valid(); nit.advance()) {
					INDEX_TYPE neighbor_id = nit.value();
					if (region->processed_set.count(neighbor_id) == 0 && !IsProcessed(neighbor_id)){
						all_local = false;
						break;
					}
				}
				if (all_local) {
					region->prior_data.erase(id);
					count++;
				}

			}
			//printf("remoed %d of %d \n", count, region->prior_data.size());
		}



#ifdef DEBUGDELAY
		std::vector<INDEX_TYPE> delayed_cells; // debug
		std::vector<INDEX_TYPE> critical_primal_edges;// debug
#endif

		//std::stack<WorkRequest> m_work_to_do;

		struct pq_comp {
			int priority;
			WorkRequest w;
			bool operator()(pq_comp& p1, pq_comp& p2) {
				return p1.priority > p2.priority;
			}
		};
		std::priority_queue< pq_comp, std::vector<pq_comp>, pq_comp > pq_work_request;
		std::stack<WorkRequest> m_leftover_work;
		std::stack<RegionState*> m_leftover_regions;
		//std::map<INDEX_TYPE, INDEX_TYPE> bottleneckcount;
		std::map<INDEX_TYPE, std::vector<INDEX_TYPE> > bottleneckregions;
		int* m_thread_work_count;
		std::map<INDEX_TYPE, RegionState*> m_bottleneckpoint_to_region_map;


		void DoSingleton() {

			std::vector<INDEX_TYPE> mins;
			typename DependencyGraphType::vertex_iterator primal_edge_iterator(mD);
			for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
				INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

				if (IsUNASSIGNEDMinInNeighborhood(primal_edge_topo_id)) {
					mins.push_back(primal_edge_topo_id);
				}
			}

			RegionState* local_region = new RegionState(m_id_comparer);
			SORTED_ID_QUEUE_TYPE growing_front(m_id_comparer);

			local_region->processed_set.clear();

			for (int i = 0; i < mins.size(); i++) {
				INDEX_TYPE base_id = mins[i];
				local_region->seen_set.insert(base_id);
				local_region->processed_set.insert(base_id);

				//mIndepTempProcLabel->SetLabel(base_id, omp_get_thread_num() + 1); // debug?
				//m_thread_work_count[omp_get_thread_num()]++;
				INDEX_TYPE t_id = base_id;
				local_region->prior_data[base_id] = DoInitialWork(base_id, local_region); // this represents the meat of the flood fill approach - do the base case of the recurrence relation

				Enqueue_Higher_Neighbors(base_id, growing_front, local_region->seen_set);
			}
			//FLOATTYPE original_edge_value = m_func->cellValue(base_id);
			//SORTEDPAIRSQUEUE growing_front;
			//Enqueue_Higher_Neighbors(base_id, original_edge_value, growing_front, local_region->seen_set);
			while (!growing_front.empty()) {

				//INDEX_TYPE currid = growing_front.top().second;
				INDEX_TYPE currid = growing_front.top();
				//FLOATTYPE currid_value = growing_front.top().first;
				growing_front.pop();

				// cellvalue >=0 indicates that there is certainty here, so lets expand 
				//if (AllLowerNeighborsAlreadyProcessedLocallyOrGlobally(currid, currid_value, local_region->processed_set)) {
				if (AllLowerNeighborsAlreadyProcessedLocallyOrGlobally(currid, local_region->processed_set)) {
					if (local_region->processed_set.count(currid) != 0) printf("reprocessing %d\n", currid);
					local_region->processed_set.insert(currid);
					mIndepTempProcLabel->SetLabel(currid, omp_get_thread_num() + 1); // debug?
					m_thread_work_count[omp_get_thread_num()]++;
					/// DO COMBO WORK
					local_region->prior_data[currid] = DoWork(currid, local_region, 0); // this represents the meat of the flood fill approach - do the inductive case of the recurrence relation

					//Enqueue_Higher_Neighbors(currid, currid_value, growing_front, local_region->seen_set);
					Enqueue_Higher_Neighbors(currid, growing_front, local_region->seen_set);
				}
				else {
					//local_region->unprocessed_band.push(FIPAIR(currid_value, currid));
					local_region->unprocessed_band.push(currid);
				}

			}
			if (!local_region->unprocessed_band.empty()) {
				//local_region->bottleneck_id = local_region->unprocessed_band.top().second;
				local_region->bottleneck_id = local_region->unprocessed_band.top();
			}
			local_region->region_id = 1; // debug?
			//Sparsify_Region(local_region);
			//return local_region;


		}
		int hacknum;
		int m_num_threads;
		bool AfterLocalWork(RegionState* local_region) {

			// if this region ended, just continue
			if (local_region->unprocessed_band.empty()) {
#pragma omp critical
								{
									m_leftover_regions.push(local_region);
								}
				return true;
			}

			// otherwise should we add another merge?
			// should we even try to keep merging this?
			if (local_region->seen_set.size() > hacknum) {
#pragma omp critical
					{
						//printf("delaying seen set %d\n", local_region->seen_set.size());
						pq_comp wr;
						wr.priority = local_region->seen_set.size();
						wr.w = WorkRequest{ PREP_FOR_GLOBAL, 0, local_region, NULL };
						pq_work_request.push(wr);
					}
				return true;
			}

			// begin critical update to merge list and work to do
#pragma omp critical
			{
				if (m_bottleneckpoint_to_region_map.count(local_region->bottleneck_id)) {
					RegionState* other_region = m_bottleneckpoint_to_region_map[local_region->bottleneck_id];
					m_bottleneckpoint_to_region_map.erase(local_region->bottleneck_id);
					if (pq_work_request.size() < m_num_threads) {
						//keep_looping = false;
						pq_comp wr;
						wr.priority = local_region->seen_set.size();
						wr.w = WorkRequest{ PREP_FOR_GLOBAL, 0, local_region, NULL };
						pq_work_request.push(wr);
						pq_comp wr2;
						wr2.priority = other_region->seen_set.size();
						wr2.w = WorkRequest{ PREP_FOR_GLOBAL, 0, other_region, NULL };
					}
					else {
						pq_comp wr;
						wr.priority = local_region->seen_set.size() + other_region->seen_set.size();
						wr.w = WorkRequest{ MERGE_AND_EXPAND, 0, local_region, other_region };
						pq_work_request.push(wr);

					}
				}
				else {
					m_bottleneckpoint_to_region_map[local_region->bottleneck_id] = local_region;
				}
			} // end critical update
			return true;

		}


		enum TimingType {
			TI_GLOBALIO,
			TI_IEXPAND,
			TI_INTEGRATION,
			TI_WAIT,
			TI_CLEANUP,
			TI_EXPAND_REGION,
			TI_MERGE_REGION,
			TI_SYNCHRONIZE,
			TI_ROBINS
		};


		void BeginIntegration(ThreadedTimer* timer) {
#pragma omp parallel
			{
				int thread_num = omp_get_thread_num();
				timer->RecordThreadActivity(thread_num, TI_EXPAND_REGION);

			}
			printf("staring overall\nCreating labelings\n");
			mIndepTempProcLabel = new DenseLabeling<char>(mD->IndexSpaceSize());
			mIndepTempProcLabel->SetAll(0);
			mGloballyProcessedLabel = new DenseLabeling<char>(mD->IndexSpaceSize());
			mGloballyProcessedLabel->SetAll(0);
#pragma omp parallel
			{
				int thread_num = omp_get_thread_num();
				timer->RecordThreadActivity(thread_num, TI_WAIT);

			}
			//const INDEX_TYPE t_num_vertices = m_grid->NumElements();

			//AdvectionChecker* inside_voxel_critical_advection_checker = new TerminateNearCertain(m_certains, m_grid);
			//AdvectionChecker* no_check = new NoTermination();//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);

			m_thread_work_count = new int[256];
			memset(m_thread_work_count, 0, 256 * sizeof(int));
			//int num_threads;

//
//#pragma omp parallel
//			{
//#pragma omp single
//				{
//					m_num_threads = omp_get_num_threads();
//				}
//			}
			//if (m_num_threads == 1) {
			//	DoSingleton();
			//	return;
			//}
			//mTimings = new vector<TimingReport>[m_num_threads];


			// set all potential critical_primal_edges, so we terminate near them
			std::vector<INDEX_TYPE> topo_index_partition;
			printf("starting rounds\n"); 
			for (int rounds = 0; rounds < 100; rounds++) {

				// ADD IN ALL THE CURRENT MINIMA ----------------------------------
#pragma omp parallel
			{
#pragma omp single
				{
					m_num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mD->NumIds(), m_num_threads, topo_index_partition);
					printf("starting CP detect\n");
				}

				int thread_num = omp_get_thread_num();
				timer->RecordThreadActivity(thread_num, TI_EXPAND_REGION);
					typename DependencyGraphType::vertex_iterator primal_edge_iterator(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
					INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();
					bool pushback = IsUNASSIGNEDMinInNeighborhood(primal_edge_topo_id);
					if (pushback) {
						//timer->RecordThreadActivity(thread_num, TI_WAIT);
#pragma omp critical
						{
							pq_comp wr;
							wr.priority = 1;
							wr.w = WorkRequest{ EXPAND_ONLY, primal_edge_topo_id, NULL, NULL };
							pq_work_request.push(wr);

							//mIndepTempProcLabel->SetLabel(primal_edge_topo_id, thread_num + 1);
						}
					}
				}
				timer->RecordThreadActivity(thread_num, TI_WAIT);

			}
				// END ADD IN ALL THE CURRENT MINIMA ----------------------------------
				printf("round %d: found %d mins\n", pq_work_request.size());
				printf("doing expand\n");
				// MAIN WHILE LOOP: WHILE WORK IN THIS ROUND =================================	
#pragma omp parallel
				{
					int thread_num = omp_get_thread_num();
					bool keep_looping = true;
					WorkRequest work;
					while (keep_looping) { // begin while
#pragma omp critical
						{ // begin critical
							if (pq_work_request.empty()) {
								keep_looping = false;
							}
							else {
								work = pq_work_request.top().w;
								//printf("doing p %d\n", pq_work_request.top().priority);
								//work = m_work_to_do.front();
								pq_work_request.pop();
							}
						} // end critical

						if (keep_looping) {
							if (work.worktype == EXPAND_ONLY) {
								timer->RecordThreadActivity(thread_num, TI_EXPAND_REGION);
								RegionState* local_region = ExpandRegionUpwardsFlood(work.start_cell);
								//Sparsify_Region(local_region);
								timer->RecordThreadActivity(thread_num, TI_WAIT);
								keep_looping = AfterLocalWork(local_region);
								if (!keep_looping) delete local_region;

							}
							else if (work.worktype == MERGE_AND_EXPAND) {
								timer->RecordThreadActivity(thread_num, TI_MERGE_REGION);
								RegionState* merged = MergeRegions(work.region1, work.region2);
								timer->RecordThreadActivity(thread_num, TI_EXPAND_REGION);
								RegionState* modified = Expand_MergeRegions(merged);
								timer->RecordThreadActivity(thread_num, TI_WAIT);
								//Sparsify_Region(modified);
								keep_looping = AfterLocalWork(modified);
								if (!keep_looping) delete modified;
							}
							else if (work.worktype == PREP_FOR_GLOBAL) {
								timer->RecordThreadActivity(thread_num, TI_EXPAND_REGION);
								Sparsify_Region(work.region1);
								timer->RecordThreadActivity(thread_num, TI_WAIT);

#pragma omp critical
								{
									m_leftover_regions.push(work.region1);
								}
							}

						} // end if(keep_looping)
					} // end while (keep_looping) 		
				} // end parallel region

				// END MAIN WHILE LOOP  ======================================================	

				// SYNCHRONIZATION STEP
				for (int i = 0; i < m_num_threads; i++) {
					printf("thread %d did %d work\n", i, m_thread_work_count[i]);
				}
				int countassigned = 0;
				int countunassigned = 0;
				vector<INDEX_TYPE> to_erase;
				printf("collecting toerase\n");

#pragma omp parallel
				{
#pragma omp single
				{
					m_num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mD->NumIds(), m_num_threads, topo_index_partition);
				}
					int tempcountassigned = 0;
					int tempcountunassigned = 0;
					int thread_num = omp_get_thread_num();
					vector<INDEX_TYPE> myvectoerase;
					timer->RecordThreadActivity(thread_num, TI_SYNCHRONIZE);
					typename DependencyGraphType::vertex_iterator primal_edge_iterator(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
					for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
						INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

						if (mIndepTempProcLabel->GetLabel(primal_edge_topo_id) == 0) {
							tempcountunassigned++;
							//printf("ERROR: edge that is part of boundar is not certain!\n");
						}
						else {
							if (m_prior_data.count(primal_edge_topo_id) > 0 && AllNeighborsAlreadyProcessedGlobally(primal_edge_topo_id)) {
								myvectoerase.push_back(primal_edge_topo_id);
							}
							mGloballyProcessedLabel->SetLabel(primal_edge_topo_id, 1);
							tempcountassigned++;
						}
					}
#pragma omp critical
					{
						to_erase.insert(to_erase.end(), myvectoerase.begin(), myvectoerase.end());
						countassigned += tempcountassigned;
						countunassigned += tempcountunassigned;
					}
					timer->RecordThreadActivity(thread_num, TI_WAIT);

				}

				timer->RecordThreadActivity(omp_get_thread_num(), TI_SYNCHRONIZE);
				printf("round%d: %d assigned, %d unassigned\n", rounds, countassigned, countunassigned);
				if (countunassigned == 0) {
					rounds = 100;
					return;

				}
				printf("doing serial gather regions\n");

				vector<RegionState*> gather_regions;

				while (!pq_work_request.empty()) {
					WorkRequest w = pq_work_request.top().w;
					//WorkRequest w = m_work_to_do.front();
					pq_work_request.pop();
					if (w.worktype == MERGE_AND_EXPAND) {
						gather_regions.push_back(w.region1);
						gather_regions.push_back(w.region2);						
					}
				}
				while (!m_leftover_work.empty()) {
					WorkRequest w = m_leftover_work.top();
					m_leftover_work.pop();
					if (w.worktype == MERGE_AND_EXPAND) {
						gather_regions.push_back(w.region1);
						gather_regions.push_back(w.region2);
					}
				}
				while (!m_leftover_regions.empty()) {
					RegionState* r = m_leftover_regions.top();
					m_leftover_regions.pop();
					gather_regions.push_back(r);
				}
				for (auto it = m_bottleneckpoint_to_region_map.begin(); it != m_bottleneckpoint_to_region_map.end(); it++) {
					RegionState* r = (*it).second;
					gather_regions.push_back(r);
				}
				
				int numregs = gather_regions.size();
#pragma omp parallel
				{
					int thread_num = omp_get_thread_num();
					timer->RecordThreadActivity(thread_num, TI_SYNCHRONIZE);

				}
				printf("doing parallel sparsify  regions\n");

#pragma omp parallel for schedule(dynamic)
				for (int i = 0; i < numregs; i++) {
					Sparsify_Region(gather_regions[i]);

				}
#pragma omp parallel
				{
					int thread_num = omp_get_thread_num();
					timer->RecordThreadActivity(thread_num, TI_SYNCHRONIZE);
				} // end parallel counting

				printf("doing serial record and cleanup\n");

				//timer->RecordThreadActivity(omp_get_thread_num(), TI_SYNCHRONIZE);
#pragma omp parallel for
				for (int i = 0; i < numregs; i++) {
					for (auto it = gather_regions[i]->prior_data.begin(); it != gather_regions[i]->prior_data.end(); it++) {
						m_prior_data.insert((*it).first, (*it).second);
					}
					delete gather_regions[i];
				}
#pragma omp parallel for
				for (int i = 0; i < to_erase.size(); i++) {
					m_prior_data.erase(to_erase[i]);
				}
				////char outname[2048];
				////sprintf(outname, "tmp_out_%d.raw", rounds);
				////mGloballyProcessedLabel->OutputToFile(outname);

				//printf("worktodo.size = %d, leftoverwork.size = %d, leftoverregions.size = %d, unpariedregions.size = %d\n", pq_work_request.size(), m_leftover_work.size(), m_leftover_regions.size(), m_bottleneckpoint_to_region_map.size());

				//while (!pq_work_request.empty()) {
				//	WorkRequest w = pq_work_request.top().w;
				//	//WorkRequest w = m_work_to_do.front();
				//	pq_work_request.pop();
				//	if (w.worktype == MERGE_AND_EXPAND) {
				//		m_prior_data.insert(w.region1->prior_data.begin(), w.region1->prior_data.end());
				//		delete w.region1;
				//		m_prior_data.insert(w.region2->prior_data.begin(), w.region2->prior_data.end());
				//		delete w.region2;
				//	}
				//}
				//while (!m_leftover_work.empty()) {
				//	WorkRequest w = m_leftover_work.top();
				//	m_leftover_work.pop();
				//	if (w.worktype == MERGE_AND_EXPAND) {
				//		m_prior_data.insert(w.region1->prior_data.begin(), w.region1->prior_data.end());
				//		delete w.region1;
				//		m_prior_data.insert(w.region2->prior_data.begin(), w.region2->prior_data.end());
				//		delete w.region2;
				//	}
				//}
				//while (!m_leftover_regions.empty()) {
				//	RegionState* r = m_leftover_regions.top();
				//	m_leftover_regions.pop();
				//	m_prior_data.insert(r->prior_data.begin(), r->prior_data.end());
				//	delete r;
				//}
				//for (auto it = m_bottleneckpoint_to_region_map.begin(); it != m_bottleneckpoint_to_region_map.end(); it++) {
				//	RegionState* r = (*it).second;
				//	m_prior_data.insert(r->prior_data.begin(), r->prior_data.end());
				//	delete r;
				//}
				//printf("global priors.size() = %d\n", m_prior_data.size());

				m_bottleneckpoint_to_region_map.clear();
				timer->RecordThreadActivity(omp_get_thread_num(), TI_WAIT);

			} // end for rounds
		} // end begin integration function


		void BeginIntegration() {

			mIndepTempProcLabel = new DenseLabeling<char>(mD->IndexSpaceSize());
			mIndepTempProcLabel->SetAll(0);
			mGloballyProcessedLabel = new DenseLabeling<char>(mD->IndexSpaceSize());
			mGloballyProcessedLabel->SetAll(0);

			m_thread_work_count = new int[256];
			memset(m_thread_work_count, 0, 256 * sizeof(int));

			std::vector<INDEX_TYPE> topo_index_partition;

			for (int rounds = 0; rounds < 100; rounds++) {

#ifdef OUTPUTINTERMEDIATE
				vector< INDEX_TYPE> crits;
#endif				
				// ADD IN ALL THE CURRENT MINIMA ----------------------------------
#pragma omp parallel
			{
#pragma omp single
				{
					m_num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mD->NumIds(), m_num_threads, topo_index_partition);
				}



				int thread_num = omp_get_thread_num();
				//timer->RecordThreadActivity(thread_num, TI_EXPAND_REGION);
				typename DependencyGraphType::vertex_iterator primal_edge_iterator(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
					INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();
					bool pushback = IsUNASSIGNEDMinInNeighborhood(primal_edge_topo_id);
					if (pushback) {
						//timer->RecordThreadActivity(thread_num, TI_WAIT);
#pragma omp critical
						{
							pq_comp wr;
							wr.priority = 1;
							wr.w = WorkRequest{ EXPAND_ONLY, primal_edge_topo_id, NULL, NULL };
							pq_work_request.push(wr);

							//mIndepTempProcLabel->SetLabel(primal_edge_topo_id, thread_num + 1);
#ifdef OUTPUTINTERMEDIATE
							crits.push_back(primal_edge_topo_id);
#endif
						}
					}
				}
				//timer->RecordThreadActivity(thread_num, TI_WAIT);

			}
				printf("round %d: found %d mins\n",rounds,  pq_work_request.size());

				// END ADD IN ALL THE CURRENT MINIMA ----------------------------------
#ifdef OUTPUTINTERMEDIATE
				char nnames2[2048];
				sprintf(nnames2, "crit_ids_%d.txt", rounds);
		
				FILE* pfile2 = fopen(nnames2, "w");
				
				for (auto it = crits.begin(); it != crits.end(); it++) {
					INDEX_TYPE id = *it;
					
						fprintf(pfile2, "%d\n", id);
						
				}
				fclose(pfile2);

#endif

				// MAIN WHILE LOOP: WHILE WORK IN THIS ROUND =================================	
#pragma omp parallel
				{
					int thread_num = omp_get_thread_num();
					bool keep_looping = true;
					WorkRequest work;
					while (keep_looping) { // begin while
#pragma omp critical
						{ // begin critical
							if (pq_work_request.empty()) {
								keep_looping = false;
							}
							else {
								work = pq_work_request.top().w;
								//printf("doing p %d\n", pq_work_request.top().priority);
								//work = m_work_to_do.front();
								pq_work_request.pop();
							}
						} // end critical

						if (keep_looping) {
							if (work.worktype == EXPAND_ONLY) {
								//timer->RecordThreadActivity(thread_num, TI_EXPAND_REGION);
								RegionState* local_region = ExpandRegionUpwardsFlood(work.start_cell);
								//Sparsify_Region(local_region);
								//timer->RecordThreadActivity(thread_num, TI_WAIT);
								keep_looping = AfterLocalWork(local_region);
								if (!keep_looping) delete local_region;
							}
							else if (work.worktype == MERGE_AND_EXPAND) {
								//timer->RecordThreadActivity(thread_num, TI_MERGE_REGION);
								RegionState* merged = MergeRegions(work.region1, work.region2);
								//timer->RecordThreadActivity(thread_num, TI_EXPAND_REGION);
								RegionState* modified = Expand_MergeRegions(merged);
								//timer->RecordThreadActivity(thread_num, TI_WAIT);
								//Sparsify_Region(modified);
								keep_looping = AfterLocalWork(modified);
								if (!keep_looping) delete modified;
							}
							else if (work.worktype == PREP_FOR_GLOBAL) {
								//timer->RecordThreadActivity(thread_num, TI_EXPAND_REGION);
								Sparsify_Region(work.region1);
								//timer->RecordThreadActivity(thread_num, TI_WAIT);

#pragma omp critical
								{
									m_leftover_regions.push(work.region1);
								}
							}

						} // end if(keep_looping)
					} // end while (keep_looping) 		
				} // end parallel region

				// END MAIN WHILE LOOP  ======================================================	

				// SYNCHRONIZATION STEP
				for (int i = 0; i < m_num_threads; i++) {
					printf("thread %d did %d work\n", i, m_thread_work_count[i]);
				}
#ifdef OUTPUTINTERMEDIATE
				char nname[2048];
				sprintf(nname, "exapand_certains_%d.raw", rounds);
				mGloballyProcessedLabel->OutputToFile(nname);
				sprintf(nname, "exapand_indep_%d.raw", rounds);
				mIndepTempProcLabel->OutputToFile(nname);
				sprintf(nname, "gpriordata_a_%d.txt", rounds);
				FILE* pfile = fopen(nname, "w");
				typename DependencyGraphType::vertex_iterator it(mD);
				for (it.begin(); it.valid(); it.advance()) {
					INDEX_TYPE id = it.value();
					if (m_prior_data.count(id) != 0) {
						fprintf(pfile, "%d\n", id);
					}
				}
				fclose(pfile);

#endif

				int countassigned = 0;
				int countunassigned = 0;
				vector<INDEX_TYPE> to_erase;

#pragma omp parallel
				{
#pragma omp single
				{
					m_num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mD->NumIds(), m_num_threads, topo_index_partition);
				}
					int tempcountassigned = 0;
					int tempcountunassigned = 0;
					int thread_num = omp_get_thread_num();
					vector<INDEX_TYPE> myvectoerase;
					//timer->RecordThreadActivity(thread_num, TI_SYNCHRONIZE);
					typename DependencyGraphType::vertex_iterator primal_edge_iterator(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
					for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
						INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

						if (mIndepTempProcLabel->GetLabel(primal_edge_topo_id) == 0) {
							tempcountunassigned++;
							//printf("ERROR: edge that is part of boundar is not certain!\n");
						}
						else {
							if (m_prior_data.count(primal_edge_topo_id) > 0 && AllNeighborsAlreadyProcessedGlobally(primal_edge_topo_id)) {
								myvectoerase.push_back(primal_edge_topo_id);
							}
							mGloballyProcessedLabel->SetLabel(primal_edge_topo_id, 1);
							tempcountassigned++;
						}
					}
#pragma omp critical
					{
						to_erase.insert(to_erase.end(), myvectoerase.begin(), myvectoerase.end());
						countassigned += tempcountassigned;
						countunassigned += tempcountunassigned;
					}
					//timer->RecordThreadActivity(thread_num, TI_WAIT);

				} // end parallel counting

				//timer->RecordThreadActivity(omp_get_thread_num(), TI_SYNCHRONIZE);
				printf("round%d: %d assigned, %d unassigned\n", rounds, countassigned, countunassigned);
				if (countunassigned == 0) {
					rounds = 100;
					return;

				}

				vector<RegionState*> gather_regions;

				while (!pq_work_request.empty()) {
					WorkRequest w = pq_work_request.top().w;
					//WorkRequest w = m_work_to_do.front();
					pq_work_request.pop();
					if (w.worktype == MERGE_AND_EXPAND) {
						gather_regions.push_back(w.region1);
						gather_regions.push_back(w.region2);
					}
				}
				while (!m_leftover_work.empty()) {
					WorkRequest w = m_leftover_work.top();
					m_leftover_work.pop();
					if (w.worktype == MERGE_AND_EXPAND) {
						gather_regions.push_back(w.region1);
						gather_regions.push_back(w.region2);
					}
				}
				while (!m_leftover_regions.empty()) {
					RegionState* r = m_leftover_regions.top();
					m_leftover_regions.pop();
					gather_regions.push_back(r);
				}
				for (auto it = m_bottleneckpoint_to_region_map.begin(); it != m_bottleneckpoint_to_region_map.end(); it++) {
					RegionState* r = (*it).second;
					gather_regions.push_back(r);
				}

				int numregs = gather_regions.size();

#pragma omp parallel for schedule(dynamic)
				for (int i = 0; i < numregs; i++) {
					Sparsify_Region(gather_regions[i]);

				}
				//printf("doing serial record and cleanup\n");

				//timer->RecordThreadActivity(omp_get_thread_num(), TI_SYNCHRONIZE);
#pragma omp parallel for
				for (int i = 0; i < numregs; i++) {
					for (auto it = gather_regions[i]->prior_data.begin(); it != gather_regions[i]->prior_data.end(); it++) {
						m_prior_data.insert((*it).first, (*it).second);
					}
					delete gather_regions[i];
				}
#ifdef OUTPUTINTERMEDIATE

				sprintf(nname, "gpriordata_b_%d.txt", rounds);
				pfile = fopen(nname, "w");
				for (it.begin(); it.valid(); it.advance()) {
					INDEX_TYPE id = it.value();
					if (m_prior_data.count(id) != 0) {
						fprintf(pfile, "%d\n", id);
					}
				}
				fclose(pfile);

#endif

#pragma omp parallel for
				for (int i = 0; i < to_erase.size(); i++) {
					m_prior_data.erase(to_erase[i]);
				}

				//printf("global priors.size() = %d\n", m_prior_data.size());

				m_bottleneckpoint_to_region_map.clear();
				//timer->RecordThreadActivity(omp_get_thread_num(), TI_WAIT);

			} // end for rounds
		} // end begin integration function


	};








}
#endif
