#ifndef GI_EXPERIMENTAL3_H
#define GI_EXPERIMENTAL3_H

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

//#include "concurrentqueue.h"
#include <thread>
#include <mutex>

#include "gi_experimental.h"
// This file is a safe space to put experimental classes and functionality.

namespace GInt {



	template <class PayloadType, class DependencyGraphType>
	class TopologicalGraphExpandingRegions {
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
		std::unordered_map<INDEX_TYPE, PayloadType> m_prior_data;
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
				if (AllLowerNeighborsAlreadyProcessedLocallyOrGlobally(currid,  local_region->processed_set)) {
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

		TopologicalGraphExpandingRegions(
			//TopologicalExplicitDenseMeshFunction* func,
			DependencyGraphType* graph) :
			//m_func(func),
			mD(graph),
			m_id_comparer(graph)
		{
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
					if (region->processed_set.count(neighbor_id) == 0  && ! IsProcessed(neighbor_id)){
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

		std::stack<WorkRequest> m_work_to_do;
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


		void BeginIntegration() {
			mIndepTempProcLabel = new DenseLabeling<char>(mD->IndexSpaceSize());
			mIndepTempProcLabel->SetAll(0);
			mGloballyProcessedLabel = new DenseLabeling<char>(mD->IndexSpaceSize());
			mGloballyProcessedLabel->SetAll(0);

			//const INDEX_TYPE t_num_vertices = m_grid->NumElements();

			//AdvectionChecker* inside_voxel_critical_advection_checker = new TerminateNearCertain(m_certains, m_grid);
			//AdvectionChecker* no_check = new NoTermination();//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);

			m_thread_work_count = new int[256];
			memset(m_thread_work_count, 0, 256 * sizeof(int));
			int num_threads;


#pragma omp parallel
			{
#pragma omp single
				{
					num_threads = omp_get_num_threads();
				}
			}
			if (num_threads == 1) {
				DoSingleton();
				return;


			}


			// set all potential critical_primal_edges, so we terminate near them
			std::vector<INDEX_TYPE> topo_index_partition;
//#define PRINTPROGRESS
#ifdef PRINTPROGRESS
			int progresscount = 0;
			printf(" -- progress %d\n", progresscount++);
#endif
#pragma omp parallel
			{
#pragma omp single
				{
					num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mD->IndexSpaceSize(), num_threads, topo_index_partition);
				}

				int thread_num = omp_get_thread_num();
				typename DependencyGraphType::vertex_iterator primal_edge_iterator(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
					INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

					if (IsUNASSIGNEDMinInNeighborhood(primal_edge_topo_id)) {
#pragma omp critical
						{
							m_work_to_do.push(WorkRequest{ EXPAND_ONLY, primal_edge_topo_id, NULL, NULL });
#ifdef DEBUGDELAY
							critical_primal_edges.push_back(primal_edge_topo_id);
#endif
							//mIndepTempProcLabel->SetLabel(primal_edge_topo_id, thread_num + 1);
						}
					}
				}

			}
#ifdef PRINTPROGRESS
			printf(" -- progress %d\n", progresscount++);
#endif
#ifdef DEBUGDELAY
			printf("found %d critical_primal_edges!\nExpanding critical_primal_edges certain regions\n", critical_primal_edges.size());
			FILE* fmerge = fopen("merge.txt", "w");
#endif

#ifdef DOTIMINGS
			std::vector< std::pair<std::chrono::steady_clock::time_point, int> > timings;
			timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
#endif
#pragma omp parallel
			{
				bool keep_looping = true;
				WorkRequest work;
				while (keep_looping) { // begin while
#pragma omp critical
						{ // begin critical
							if (m_work_to_do.empty()) {
								keep_looping = false;
							}
							else {
								work = m_work_to_do.top();
								//work = m_work_to_do.front();
								m_work_to_do.pop();
#ifdef DOTIMINGS
								timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
#endif
							}
#ifdef PRINTPROGRESS
							printf(" -- progress %d %d\n", progresscount++, work.worktype);
#endif
						} // end critical

					if (keep_looping) {
						if (work.worktype == EXPAND_ONLY) {

																	// debug
							int thread_num = omp_get_thread_num();														// debug

							RegionState* local_region = ExpandRegionUpwardsFlood(work.start_cell);
							//Sparsify_Region(local_region);

							if (!local_region->unprocessed_band.empty()) {
#pragma omp critical
								{ // begin critical update to merge list and work to do
									if (m_bottleneckpoint_to_region_map.count(local_region->bottleneck_id) > 0) {
										RegionState* other_region = m_bottleneckpoint_to_region_map[local_region->bottleneck_id];
										m_bottleneckpoint_to_region_map.erase(local_region->bottleneck_id);
										m_work_to_do.push(WorkRequest{ MERGE_AND_EXPAND, 0, local_region, other_region });
#ifdef DOTIMINGS
										timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
#endif

									}
									else {
										m_bottleneckpoint_to_region_map[local_region->bottleneck_id] = local_region;
									}

#ifdef DEBUGDELAY
									// BEGIN DEBUG REGION
									SORTEDPAIRSQUEUE tmp_front(local_region->unprocessed_band);
									while (!tmp_front.empty()) {
										delayed_cells.push_back(tmp_front.top().second);
										tmp_front.pop();
									}
#endif
									// END DEBUG REGION

									//while (!local_region->unprocessed_band.empty()) {											// debug
									//	delayed_cells.push_back(local_region->unprocessed_band.top().second);					// debug
									//	local_region->unprocessed_band.pop();													// debug
									//}																						// debug
								} // end critical update
							}
							else {
#pragma omp critical
								{
									m_leftover_regions.push(local_region);
								}
							}

						}
						else if (work.worktype == MERGE_AND_EXPAND) {
							RegionState* merged = MergeRegions(work.region1, work.region2);
							RegionState* modified = Expand_MergeRegions(merged);
							//Sparsify_Region(modified);
							// only add more work if there IS work to be added - could end at a max!
							if (!modified->unprocessed_band.empty()/* && modified->seen_set.size() < 1000*/){
#pragma omp critical
								{ // begin critical update to merge list and work to do
									if (m_bottleneckpoint_to_region_map.count(modified->bottleneck_id)) {
										RegionState* other_region = m_bottleneckpoint_to_region_map[modified->bottleneck_id];
										m_bottleneckpoint_to_region_map.erase(modified->bottleneck_id);
										if (m_work_to_do.size() < num_threads) {
											keep_looping = false;
											m_leftover_work.push(WorkRequest{ MERGE_AND_EXPAND, 0, modified, other_region });
										}
										else if (m_work_to_do.size() < num_threads) {
											keep_looping = false;
											m_leftover_work.push(WorkRequest{ MERGE_AND_EXPAND, 0, modified, other_region });
										}
										else {
											m_work_to_do.push(WorkRequest{ MERGE_AND_EXPAND, 0, modified, other_region });
#ifdef DOTIMINGS
											timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
#endif
#ifdef DEBUGDELAY
											fprintf(fmerge, "pushing merge point %d %d %d %d %d\n", (int)modified->bottleneck_id, modified->seen_set.size(), other_region->seen_set.size(), m_work_to_do.size(), omp_get_thread_num());
#endif
										}
									}
									else {
										m_bottleneckpoint_to_region_map[modified->bottleneck_id] = modified;
									}

#ifdef DEBUGDELAY
									// debug
									SORTEDPAIRSQUEUE tmp_front(modified->unprocessed_band);
									while (!tmp_front.empty()) {
										delayed_cells.push_back(tmp_front.top().second);
										tmp_front.pop();
									}
#endif
								} // end critical update
							}
							else {
								// end if delayed front not empty
#pragma omp critical
								{
									m_leftover_regions.push(modified);
								}
							}
						}
					} // end if(keep_looping)
				} // end while (keep_looping) 		
#ifdef DOTIMINGS
				timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
#endif
			} // end parallel region
#ifdef PRINTPROGRESS
			printf(" -- progress b %d\n", progresscount++);
#endif
			for (int i = 0; i < num_threads; i++) {
				printf("thread %d did %d work\n", i, m_thread_work_count[i]);
			}

			printf("worktodo.size = %d, leftoverwork.size = %d, leftoverregions.size = %d, unpariedregions.size = %d\n", m_work_to_do.size(), m_leftover_work.size(), m_leftover_regions.size(), m_bottleneckpoint_to_region_map.size());

			int countassigned = 0;
			int countunassigned = 0;
#pragma omp parallel
			{
#pragma omp single
				{
					num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mD->IndexSpaceSize(), num_threads, topo_index_partition);
				}
				int tempcountassigned = 0;
				int tempcountunassigned = 0;
				int thread_num = omp_get_thread_num();
				typename DependencyGraphType::vertex_iterator primal_edge_iterator(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
					INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();
					if (mIndepTempProcLabel->GetLabel(primal_edge_topo_id) == 0) {
						tempcountunassigned++;
						//printf("ERROR: edge that is part of boundar is not certain!\n");
					}
					else {
						tempcountassigned++;
					}
				}
#pragma omp critical
					{
						countassigned += tempcountassigned;
						countunassigned += tempcountunassigned;
					}

			} // end parallel counting

			
			printf("%d assigned, %d unassigned\n", countassigned, countunassigned);
			while (!m_work_to_do.empty()) {
				WorkRequest w = m_work_to_do.top();
				//WorkRequest w = m_work_to_do.front();
				m_work_to_do.pop();
				if (w.worktype == MERGE_AND_EXPAND) {
					m_prior_data.insert(w.region1->prior_data.begin(), w.region1->prior_data.end());
					delete w.region1;
					m_prior_data.insert(w.region2->prior_data.begin(), w.region2->prior_data.end());
					delete w.region2;
				}
			}
			while (!m_leftover_work.empty()) {
				WorkRequest w = m_leftover_work.top();
				m_leftover_work.pop();
				if (w.worktype == MERGE_AND_EXPAND) {
					m_prior_data.insert(w.region1->prior_data.begin(), w.region1->prior_data.end());
					delete w.region1;
					m_prior_data.insert(w.region2->prior_data.begin(), w.region2->prior_data.end());
					delete w.region2;
				}
			}
			while (!m_leftover_regions.empty()) {
				RegionState* r = m_leftover_regions.top();
				m_leftover_regions.pop();
				m_prior_data.insert(r->prior_data.begin(), r->prior_data.end());
				delete r;
			}
			for (auto it = m_bottleneckpoint_to_region_map.begin(); it != m_bottleneckpoint_to_region_map.end(); it++) {
				RegionState* r = (*it).second;
				m_prior_data.insert(r->prior_data.begin(), r->prior_data.end());
				delete r;
			}
			printf("global priors.size() = %d\n", m_prior_data.size());

			m_bottleneckpoint_to_region_map.clear();
#ifdef PRINTPROGRESS
			printf(" -- progress b %d\n", progresscount++);
#endif
			for (int rounds = 0; rounds < 100; rounds++) {

				// this simulates a global "merge" - we merge the processed labels
#pragma omp parallel
			{

				int thread_num = omp_get_thread_num();
				typename DependencyGraphType::vertex_iterator primal_edge_iterator(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
					INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

					if (mIndepTempProcLabel->GetLabel(primal_edge_topo_id) > 0) {
						mGloballyProcessedLabel->SetLabel(primal_edge_topo_id, 1);
					}
				}

			}

#ifdef DEBUGDELAY
				critical_primal_edges.clear();
#endif

#pragma omp parallel
				{

					int thread_num = omp_get_thread_num();
					typename DependencyGraphType::vertex_iterator primal_edge_iterator(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
					for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
						INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();

						if (IsUNASSIGNEDMinInNeighborhood(primal_edge_topo_id)) {
#pragma omp critical
						{
							m_work_to_do.push(WorkRequest{ EXPAND_ONLY, primal_edge_topo_id, NULL, NULL });
							//printf("found critical2\n");
#ifdef DEBUGDELAY
							critical_primal_edges.push_back(primal_edge_topo_id);
#endif
							//mIndepTempProcLabel->SetLabel(primal_edge_topo_id, thread_num + 1);
						}
						}
					}

				} // end parallel min finding
#ifdef DEBUGDELAY
				printf("found %d critical_primal_edges!\nExpanding critical_primal_edges certain regions\n", critical_primal_edges.size());
#endif

#pragma omp parallel
				{
					bool keep_looping = true;
					WorkRequest work;
					while (keep_looping) { // begin while
#pragma omp critical
						{ // begin critical
							if (m_work_to_do.empty()) {
								keep_looping = false;
							}
							else {
								work = m_work_to_do.top();
								//work = m_work_to_do.front();
								m_work_to_do.pop();
#ifdef DOTIMINGS
								timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
#endif
							}
						} // end critical

						if (keep_looping) {
							if (work.worktype == EXPAND_ONLY) {

								std::set<INDEX_TYPE> local_certains;														// debug
								int thread_num = omp_get_thread_num();														// debug

								RegionState* local_region = ExpandRegionUpwardsFlood(work.start_cell);
								//Sparsify_Region(local_region);

								if (!local_region->unprocessed_band.empty()) {
#pragma omp critical
								{ // begin critical update to merge list and work to do
									if (m_bottleneckpoint_to_region_map.count(local_region->bottleneck_id) > 0) {
										RegionState* other_region = m_bottleneckpoint_to_region_map[local_region->bottleneck_id];
										m_bottleneckpoint_to_region_map.erase(local_region->bottleneck_id);
										m_work_to_do.push(WorkRequest{ MERGE_AND_EXPAND, 0, local_region, other_region });
#ifdef DOTIMINGS
										timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
#endif
									}
									else {
										m_bottleneckpoint_to_region_map[local_region->bottleneck_id] = local_region;
									}

#ifdef DEBUGDELAY
									// BEGIN DEBUG REGION
									SORTEDPAIRSQUEUE tmp_front(local_region->unprocessed_band);
									while (!tmp_front.empty()) {
										delayed_cells.push_back(tmp_front.top().second);
										tmp_front.pop();
									}
									// END DEBUG REGION
#endif

									//while (!local_region->unprocessed_band.empty()) {											// debug
									//	delayed_cells.push_back(local_region->unprocessed_band.top().second);					// debug
									//	local_region->unprocessed_band.pop();													// debug
									//}																						// debug
								} // end critical update
								}
								else {
#pragma omp critical
								{
									m_leftover_regions.push(local_region);
								}
								}

							}
							else if (work.worktype == MERGE_AND_EXPAND) {
								RegionState* merged = MergeRegions(work.region1, work.region2);
								RegionState* modified = Expand_MergeRegions(merged);
								//Sparsify_Region(modified);
								// only add more work if there IS work to be added - could end at a max!
								if (!modified->unprocessed_band.empty() /*&& modified->seen_set.size() < 1000*/){
#pragma omp critical
								{ // begin critical update to merge list and work to do
									if (m_bottleneckpoint_to_region_map.count(modified->bottleneck_id)) {
										RegionState* other_region = m_bottleneckpoint_to_region_map[modified->bottleneck_id];
										m_bottleneckpoint_to_region_map.erase(modified->bottleneck_id);
										if (m_work_to_do.size() < num_threads) {
											keep_looping = false;
											m_leftover_work.push(WorkRequest{ MERGE_AND_EXPAND, 0, modified, other_region });

										}
										else {
											m_work_to_do.push(WorkRequest{ MERGE_AND_EXPAND, 0, modified, other_region });
#ifdef DOTIMINGS
											timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
#endif
#ifdef DEBUGDELAY
											fprintf(fmerge, "pushing merge point %d %d %d %d %d\n", (int)modified->bottleneck_id, modified->seen_set.size(), other_region->seen_set.size(), m_work_to_do.size(), omp_get_thread_num());
#endif
										}
									}
									else {
										m_bottleneckpoint_to_region_map[modified->bottleneck_id] = modified;
									}

#ifdef DEBUGDELAY
									// debug
									SORTEDPAIRSQUEUE tmp_front(modified->unprocessed_band);
									while (!tmp_front.empty()) {
										delayed_cells.push_back(tmp_front.top().second);
										tmp_front.pop();
									}
#endif
								} // end critical update
								}
								else{
#pragma omp critical
								{
									m_leftover_regions.push(modified);
								}
								}
							}
						} // end if(keep_looping)
					} // end while (keep_looping) 		
#ifdef DOTIMINGS
					timings.push_back(std::pair<std::chrono::steady_clock::time_point, int>(std::chrono::steady_clock::now(), m_work_to_do.size()));
#endif
				} // end parallel region

				for (int i = 0; i < num_threads; i++) {
					printf("thread %d did %d work\n", i, m_thread_work_count[i]);
				}
				countassigned = 0;
				countunassigned = 0;
#pragma omp parallel
				{
#pragma omp single
				{
					num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mD->IndexSpaceSize(), num_threads, topo_index_partition);
				}
					int tempcountassigned = 0;
					int tempcountunassigned = 0;
					int thread_num = omp_get_thread_num();
					typename DependencyGraphType::vertex_iterator primal_edge_iterator(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
					for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
						INDEX_TYPE primal_edge_topo_id = primal_edge_iterator.value();
						if (mIndepTempProcLabel->GetLabel(primal_edge_topo_id) == 0) {
							tempcountunassigned++;
							//printf("ERROR: edge that is part of boundar is not certain!\n");
						}
						else {
							tempcountassigned++;
						}
					}
#pragma omp critical
					{
						countassigned += tempcountassigned;
						countunassigned += tempcountunassigned;
					}

				} // end parallel counting
				printf("round%d: %d assigned, %d unassigned\n", rounds, countassigned, countunassigned);

				if (countunassigned == 0) {
					rounds = 100;
					continue;
				}
				printf("worktodo.size = %d, leftoverwork.size = %d, leftoverregions.size = %d, unpariedregions.size = %d\n", m_work_to_do.size(), m_leftover_work.size(), m_leftover_regions.size(), m_bottleneckpoint_to_region_map.size());

				while (!m_work_to_do.empty()) {
					WorkRequest w = m_work_to_do.top();
					//WorkRequest w = m_work_to_do.front();
					m_work_to_do.pop();
					if (w.worktype == MERGE_AND_EXPAND) {
						m_prior_data.insert(w.region1->prior_data.begin(), w.region1->prior_data.end());
						delete w.region1;
						m_prior_data.insert(w.region2->prior_data.begin(), w.region2->prior_data.end());
						delete w.region2;
					}
				}
				while (!m_leftover_work.empty()) {
					WorkRequest w = m_leftover_work.top();
					m_leftover_work.pop();
					if (w.worktype == MERGE_AND_EXPAND) {
						m_prior_data.insert(w.region1->prior_data.begin(), w.region1->prior_data.end());
						delete w.region1;
						m_prior_data.insert(w.region2->prior_data.begin(), w.region2->prior_data.end());
						delete w.region2;
					}
				}
				while (!m_leftover_regions.empty()) {
					RegionState* r = m_leftover_regions.top();
					m_leftover_regions.pop();
					m_prior_data.insert(r->prior_data.begin(), r->prior_data.end());
					delete r;
				}
				for (auto it = m_bottleneckpoint_to_region_map.begin(); it != m_bottleneckpoint_to_region_map.end(); it++) {
					RegionState* r = (*it).second;
					m_prior_data.insert(r->prior_data.begin(), r->prior_data.end());
					delete r;
				}
				printf("global priors.size() = %d\n", m_prior_data.size());

				m_bottleneckpoint_to_region_map.clear();
			}
#ifdef DEBUGDELAY
			fclose(fmerge);
#endif
			////			// THIS WILL ALL CHANGE -> one more indirection needed
			////			std::map<INDEX_TYPE, INDEX_TYPE>* pmaps =
			////				new std::map<INDEX_TYPE, INDEX_TYPE>[topo_index_partition.size() - 1];
			////
			////			INDEX_TYPE num_critical_primal_edges = critical_primal_edges.size();
			////#pragma omp parallel for schedule(dynamic) 
			////			for (INDEX_TYPE i = 0; i < num_critical_primal_edges; i++) {
			////				std::set<INDEX_TYPE> local_certains;
			////				int thread_num = omp_get_thread_num();
			////				RegionState* local_region = Expand_Lower_Neighborhood2(critical_primal_edges[i]);
			////
			////				for (auto it = local_region->processed_set.begin(); it != local_region->processed_set.end(); it++) {
			////
			////					pmaps[thread_num][*it] = critical_primal_edges[i];
			////					mIndepTempProcLabel->SetLabel(*it, mIndepTempProcLabel->GetLabel(*it) + 1);
			////				}
			////#pragma omp critical
			////				{
			////					while (!local_region->unprocessed_band.empty()){
			////						delayed_cells.push_back(local_region->unprocessed_band.top().second);
			////						local_region->unprocessed_band.pop();
			////					}
			////				}
			////			}



		} // end begin integration function

	};

}
#endif
