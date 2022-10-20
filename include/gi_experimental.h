#ifndef GI_EXPERIMENTAL_H
#define GI_EXPERIMENTAL_H

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
#include "gi_array_index_partition.h"
#include "gi_graphs.h"


//#include "concurrentqueue.h"
#include <omp.h>
#include <thread>
#include <mutex>
#include <map>

// This file is a safe space to put experimental classes and functionality.

namespace GInt {




	// these are the methods and iterators used for static polymorphism - in the templated implementation of a graph traversal algorithm
	// any graph MUST implement ALL these methods with the same function prototype,
	// as well as have a vertex_iterator (iterate over all nodes of the graph) and neighbor_iterator (iterate over all nodes adjacent to a 
	// given node). furthermore, all function prototypes must match as well


	class TopoVertexGraph {
	protected:
		TopologicalRegularGrid3D* mGrid;
		TopologicalExplicitDenseMeshFunction<TopologicalRegularGrid3D, float>* mFunc;
	public:
		TopoVertexGraph(TopologicalRegularGrid3D* grid, TopologicalExplicitDenseMeshFunction<TopologicalRegularGrid3D, float>* func) : mFunc(func), mGrid(grid) {}
		//INDEX_TYPE NumElements() { return mGrid->numCells(0); }

		class vertex_iterator {
		protected:
			TopologicalRegularGrid3D::DCellsIterator mit;
		public:
			vertex_iterator(TopoVertexGraph* g) : mit(g->mGrid, 0) {}
			void begin() { mit.begin(); }
			void advance() { mit.advance(); }
			bool valid() { return mit.valid(); }
			INDEX_TYPE value() const { return mit.value(); }

		};

		class neighbor_iterator {
		protected:
			TopologicalRegularGrid3D::CofacetsIterator mit;
			INDEX_TYPE mid;
		public:
			neighbor_iterator(TopoVertexGraph* g) : mit(g->mGrid) {}
			void begin(INDEX_TYPE id) { mit.begin(id); mid = id; }
			void advance() { mit.advance(); }
			bool valid() { return mit.valid(); }
			INDEX_TYPE value() const { return 2 * mit.value() - mid; } // double offset // 0 1 2 3 4 5, mid = 2, 3*2-2=4, 1*2-2=0, // mid = 3, 2*2-3=1, 2*4-3=5

		};

		bool Before(INDEX_TYPE a, INDEX_TYPE b) {
			return mFunc->lessThan(a, b);
		}
	};

	// skips anything where restricted = 0
	class BoundaryEdgeGraph {
	protected:
		TopologicalRegularGridRestricted* mGrid;
		TopologicalExplicitDenseMeshFunction<TopologicalRegularGridRestricted, float>* mFunc;
	public:
		BoundaryEdgeGraph(TopologicalRegularGridRestricted* grid, TopologicalExplicitDenseMeshFunction<TopologicalRegularGridRestricted, float>* func) : mFunc(func), mGrid(grid) {}
		//INDEX_TYPE NumElements() { return mGrid->numCells(1); }
		INDEX_TYPE IndexSpaceSize() { return mGrid->numCells(); }
		class vertex_iterator {
		protected:
			TopologicalRegularGridRestricted::DCellsIterator mit;
			TopologicalRegularGridRestricted* mGrid;
		public:
			vertex_iterator(BoundaryEdgeGraph* g) : mit(g->mGrid, 1), mGrid(g->mGrid) {}
			vertex_iterator(BoundaryEdgeGraph* g, INDEX_TYPE start, INDEX_TYPE end) : mit(g->mGrid, 1, start, end), mGrid(g->mGrid) {}
			void begin() {
				mit.begin();
				if (mit.valid() && mGrid->restrictionLabel(mit.value()) == 0) advance(); // this hardcoded 1 depends on the output of the labelign - an abstraction weakness
			}
			void advance() {
				mit.advance();
				while (mit.valid() && mGrid->restrictionLabel(mit.value()) == 0) mit.advance();
			}
			bool valid() {
				return mit.valid();
			}
			INDEX_TYPE value() const {
				return mit.value();
			}

		};

		/// this is the tough one
		class neighbor_iterator {
		protected:
			TopologicalRegularGridRestricted::CofacetsIterator cit;
			TopologicalRegularGridRestricted::FacetsIterator fit;
			TopologicalRegularGridRestricted* mGrid;

			bool try_advance_facets() {
				if (!fit.valid()) return false;
				fit.advance();
				return fit.valid();
			}

			bool try_advance_cofacets() {
				if (!cit.valid()) return false;
				cit.advance();
				return cit.valid();
			}

			void invalidate() {
				mValid = false;
			}
			void validate() {
				mValid = true;
			}

			void next_edge() {
				while (true) {
					if (try_advance_facets()) return validate();
					if (try_advance_cofacets()) {
						fit.begin(cit.value());
						if (fit.valid()) return validate();
					}
					else {
						return invalidate();
					}
				}
			}

			bool test() {
				if (value() == mOrigin || mGrid->restrictionLabel(value()) == 0) return false;
				return true;
			}

			bool mValid;
			INDEX_TYPE mOrigin;
			//INDEX_TYPE mid;
		public:
			neighbor_iterator(BoundaryEdgeGraph* g) : cit(g->mGrid), fit(g->mGrid), mGrid(g->mGrid) {}
			void begin(INDEX_TYPE id) {
				mOrigin = id;
				mValid = true;
				cit.begin(mOrigin);
				if (!cit.valid()) {
					return invalidate();
				}
				fit.begin(cit.value());
				if (!fit.valid() || !test()) advance();
			}

			void advance() {

				next_edge();
				while (valid() && !test()) next_edge();
			}
			bool valid() { return mValid; }
			INDEX_TYPE value() const { return fit.value(); } // double offset // 0 1 2 3 4 5, mid = 2, 3*2-2=4, 1*2-2=0, // mid = 3, 2*2-3=1, 2*4-3=5

		};

		bool Before(INDEX_TYPE a, INDEX_TYPE b) {
			return mFunc->lessThan(a, b);
		}
	};

	//class testconnectivity {
	//public:
	//	BoundaryEdgeGraph* mB;
	//	vector<pair<INDEX_TYPE, INDEX_TYPE> > mv;
	//	testconnectivity(TopologicalRegularGridRestricted* grid, TopologicalExplicitDenseMeshFunction<TopologicalRegularGrid3D, float>* func) {
	//		mB = new BoundaryEdgeGraph(grid, func);

	//		printf("running connecticity tests...\n");
	//		BoundaryEdgeGraph::vertex_iterator vit(mB);
	//		for (vit.begin(); vit.valid(); vit.advance()) {

	//			set<INDEX_TYPE> negs;
	//			BoundaryEdgeGraph::neighbor_iterator nit(mB);
	//			for (nit.begin(vit.value()); nit.valid(); nit.advance()) {
	//				INDEX_TYPE neighbor_id = nit.value();
	//				if (neighbor_id == vit.value()) printf("ERROR: CONNECTIVITY: neighbor is base id\n");
	//				if (negs.count(neighbor_id) > 0) printf("ERROR: CONNECTIVITY: neighbor muliply stuffed\n");
	//				negs.insert(neighbor_id);

	//				bool has_reflexive = false;
	//				BoundaryEdgeGraph::neighbor_iterator nit2(mB);
	//				for (nit2.begin(neighbor_id); nit2.valid(); nit2.advance()) {
	//					INDEX_TYPE neighbor_id2 = nit2.value();

	//					if (neighbor_id2 == vit.value()) has_reflexive = true;


	//				}

	//				if (!has_reflexive)  {
	//					printf("ERROR: CONNECTIVITY: %d has neighbor %d, but not reflexively\n", vit.value(), neighbor_id);
	//				}
	//				mv.push_back(pair<INDEX_TYPE, INDEX_TYPE>(vit.value(), nit.value()));
	//			}



	//		}


	//	}


	//};



	// skips anything where restricted = 0
	template <class MeshType, class TopoFuncType>
	class BoundaryEdgeGraphSparse {
	protected:
		MeshType* mGrid;
		TopoFuncType* mFunc;
		vector<INDEX_TYPE> mIds;
	public:
		BoundaryEdgeGraphSparse(MeshType* grid, TopoFuncType* func) : mFunc(func), mGrid(grid) {
		}
		
		void Init() {
			std::vector<INDEX_TYPE> topo_index_partition;
			int num_threads;
#pragma omp parallel
				{
#pragma omp single
				{
					num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mGrid->numCells(), num_threads, topo_index_partition);
				}
					int thread_num = omp_get_thread_num();
					vector<INDEX_TYPE> my_ids;
					//timer->RecordThreadActivity(thread_num, TI_SYNCHRONIZE);
					typename MeshType::DCellsIterator primal_edge_iterator(mGrid, 1, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
					for (primal_edge_iterator.begin(); primal_edge_iterator.valid(); primal_edge_iterator.advance()) {
						INDEX_TYPE id = primal_edge_iterator.value();
						if (mGrid->restrictionLabel(id) != 0) my_ids.push_back(id);
					}
#pragma omp critical
					{
						mIds.insert(mIds.end(), my_ids.begin(), my_ids.end());
					}

				} // end parallel counting
		}
		INDEX_TYPE NumIds() { return mIds.size(); }
		//INDEX_TYPE NumElements() { return mGrid->numCells(1); }
		INDEX_TYPE IndexSpaceSize() { return mGrid->numCells(); }
		class vertex_iterator {
		protected:
			BoundaryEdgeGraphSparse<MeshType, TopoFuncType>* g;
			INDEX_TYPE startid;
			INDEX_TYPE currid;
			INDEX_TYPE endid;
		public:
			vertex_iterator(BoundaryEdgeGraphSparse<MeshType,TopoFuncType>* g) : g(g) { startid = 0; endid = g->NumIds(); }
			vertex_iterator(BoundaryEdgeGraphSparse<MeshType, TopoFuncType>* g, INDEX_TYPE start, INDEX_TYPE end) : g(g), startid(start), endid(end) {}
			void begin() {
				currid = startid;
			}
			void advance() {
				currid++;
			}
			bool valid() {
				return currid < endid;
			}
			INDEX_TYPE value() const {
				return g->mIds[currid];
			}

		};

		/// this is the tough one
		class neighbor_iterator {
		protected:
			typename MeshType::CofacetsIterator cit;
			typename MeshType::FacetsIterator fit;
			MeshType* mGrid;

			bool try_advance_facets() {
				if (!fit.valid()) return false;
				fit.advance();
				return fit.valid();
			}

			bool try_advance_cofacets() {
				if (!cit.valid()) return false;
				cit.advance();
				return cit.valid();
			}

			void invalidate() {
				mValid = false;
			}
			void validate() {
				mValid = true;
			}

			void next_edge() {
				while (true) {
					if (try_advance_facets()) return validate();
					if (try_advance_cofacets()) {
						fit.begin(cit.value());
						if (fit.valid()) return validate();
					}
					else {
						return invalidate();
					}
				}
			}

			bool test() {
				if (value() == mOrigin || mGrid->restrictionLabel(value()) == 0) return false;
				return true;
			}

			bool mValid;
			INDEX_TYPE mOrigin;
			//INDEX_TYPE mid;
		public:
			neighbor_iterator(BoundaryEdgeGraphSparse<MeshType, TopoFuncType>* g) : cit(g->mGrid), fit(g->mGrid), mGrid(g->mGrid) {}
			void begin(INDEX_TYPE id) {
				mOrigin = id;
				mValid = true;
				cit.begin(mOrigin);
				if (!cit.valid()) {
					return invalidate();
				}
				fit.begin(cit.value());
				if (!fit.valid() || !test()) advance();
			}

			void advance() {

				next_edge();
				while (valid() && !test()) next_edge();
			}
			bool valid() { return mValid; }
			INDEX_TYPE value() const { return fit.value(); } // double offset // 0 1 2 3 4 5, mid = 2, 3*2-2=4, 1*2-2=0, // mid = 3, 2*2-3=1, 2*4-3=5

		};

		bool Before(INDEX_TYPE a, INDEX_TYPE b) {
			return mFunc->lessThan(a, b);
		}
	};

	//class testconnectivity {
	//public:
	//	BoundaryEdgeGraph* mB;
	//	vector<pair<INDEX_TYPE, INDEX_TYPE> > mv;
	//	testconnectivity(TopologicalRegularGridRestricted* grid, TopologicalExplicitDenseMeshFunction* func) {
	//		mB = new BoundaryEdgeGraph(grid, func);

	//		printf("running connecticity tests...\n");
	//		BoundaryEdgeGraph::vertex_iterator vit(mB);
	//		for (vit.begin(); vit.valid(); vit.advance()) {

	//			set<INDEX_TYPE> negs;
	//			BoundaryEdgeGraph::neighbor_iterator nit(mB);
	//			for (nit.begin(vit.value()); nit.valid(); nit.advance()) {
	//				INDEX_TYPE neighbor_id = nit.value();
	//				if (neighbor_id == vit.value()) printf("ERROR: CONNECTIVITY: neighbor is base id\n");
	//				if (negs.count(neighbor_id) > 0) printf("ERROR: CONNECTIVITY: neighbor muliply stuffed\n");
	//				negs.insert(neighbor_id);

	//				bool has_reflexive = false;
	//				BoundaryEdgeGraph::neighbor_iterator nit2(mB);
	//				for (nit2.begin(neighbor_id); nit2.valid(); nit2.advance()) {
	//					INDEX_TYPE neighbor_id2 = nit2.value();

	//					if (neighbor_id2 == vit.value()) has_reflexive = true;


	//				}

	//				if (!has_reflexive)  {
	//					printf("ERROR: CONNECTIVITY: %d has neighbor %d, but not reflexively\n", vit.value(), neighbor_id);
	//				}
	//				mv.push_back(pair<INDEX_TYPE, INDEX_TYPE>(vit.value(), nit.value()));
	//			}



	//		}


	//	}


	//};



	// general class to compare elements
	class Ordering {
	public:
		virtual bool  Before(INDEX_TYPE a, INDEX_TYPE b) { return false; }
	};
	// a StrongOrdering means we can compare ANY two elements in the index space and get a well defined result
	class StrongOrdering : public Ordering {};

	// a RelativeOrdering means only comparisons of adjacent elements return a valid result - however we can still
	// compare non-adjacent things but the result is not guaranteed to be consistent
	class RelativeOrdering : public Ordering {};

	//
	class DepenencyGraph : public GraphInterface, public RelativeOrdering {	};

	class Work {
	public:
		void DoWork(INDEX_TYPE a);

	};


	// the WaveFront class begins a propagation accross a dependency graph
	// it is a purely local structure, an independent piece of propagation. it can tell the difference between a local 
	// dependency and global one, as a local one will exist in its priors
	// Again, this is simply a wavefront interface, which will be used by an actual propagaion algorithm

	//// IMPLEMENTATION OF WAVEFRONT
	// have one unordered map to store PRIORS, thier payloads, as well as a count to how many WAITING things reference them
	// have one unordered map to store WAITINIG, and count how many unresolved PRIORS they have. 
	// have one list of ready-to-go items (stack or queue)	
	// general class to compare elements

	class WaveFrontInterface {
	protected:
		DepenencyGraph* mGlobalGraph;

	public:
		WaveFrontInterface(DepenencyGraph* global_graph) : mGlobalGraph(global_graph) {}

		void setSeed(INDEX_TYPE a) {} // initialize at a seed point

		bool getLowestAdjacent(INDEX_TYPE& ret) {} // gets A local minimum (i.e. before all its neighbors) that is adjacent to but not before the priors

		bool isLocal(INDEX_TYPE id) {} // returns true if all neighbors that are "before" are part of priors

		//	bool getNextToProcess(INDEX_TYPE& ret) {} // get an index of an element that is ready to go - i.e. ALL its predecessors exist in this wavefront

		void setElementAsProcessed(INDEX_TYPE a) {} // removes any PRIORS that have no more depentant ACTIVES

		// does a "deep" split of this wavefront - into "independent" wavefronts - does memory copy to new items
		// -- finds "seeds" among all WAITING and creates a new wavefront
		// -- for each "seed" finds all reachable that are in WAITING, and adds to seeds wavefront
		// -- for each new wavefront, find all PRIORS and add to the wavefront 
		// -- finally, update counts for each! - maybe can set PRIOR count when adding
		void deepSplitWaveFront(vector<WaveFrontInterface*> res) {}

		// does a "deep" merge and deletes the other wavefront - inserts items and updates counts
		// -- add all other's PRIORS, adding up counts of PRIORS that already existed
		// -- add all other's WAITING, subtracting from neighboring PRIORS if the WAITING id already existed
		void mergeAndConsumeWaveFront(WaveFrontInterface* other) {}

	};

	class randoqueue {
	protected:


	public:
		stack<INDEX_TYPE> ms;
		queue<INDEX_TYPE> mq;
		int which;
		void push(INDEX_TYPE id) {
			if (rand() % 2 == 0) {
				ms.push(id);
			}
			else {
				mq.push(id);
			}
		}
		INDEX_TYPE top() {
			if (mq.empty()) {
				which = 0;
			}
			else if (ms.empty()) {
				which = 1;
			}
			else {
				which = rand() % 2;
			}
			if (which == 0) {
				return ms.top();
			}
			else {
				return mq.front();
			}
		}

		void pop() {
			if (which == 0) {
				return ms.pop();
			}
			else {
				return mq.pop();
			}

		}

		INDEX_TYPE size() {
			return ms.size() + mq.size();
		}


	};

	//template <class PAYLOAD_TYPE>
	typedef int PAYLOAD_TYPE;
	typedef BoundaryEdgeGraph DEPGRAPH_TYPE;
	class HashingWaveFront {
	public:

		static int count;

		int my_wave_id;

		DEPGRAPH_TYPE* mGlobalGraph;

		struct PriorElement {
			INDEX_TYPE id;
			BYTE_TYPE count_unprocessed_afters;
			PAYLOAD_TYPE payload;

			PriorElement() {
				count_unprocessed_afters = 0;
			}
			PriorElement(INDEX_TYPE id, PAYLOAD_TYPE p) : id(id), payload(p) {
				count_unprocessed_afters = 0;
			}
			//~PriorElement() {
			//	delete payload; // for now this also deletes payload
			//}
		};

		map<INDEX_TYPE, PriorElement > mPriors;
		map<INDEX_TYPE, BYTE_TYPE> mWaiting;
#define QUEQUE
#ifdef RANDQ
		randoqueue mReady2Go;
#endif
#ifdef QUEQUE
		//queue<INDEX_TYPE> mReady2Go;
		set<INDEX_TYPE> mReady2Go;
#endif
#ifdef STACKQ
		stack<INDEX_TYPE> mReady2Go;
#endif

		bool ABeforeB(INDEX_TYPE a, INDEX_TYPE b) {
			return mGlobalGraph->Before(a, b);
		}


		//=====================================

		BYTE_TYPE CountAftersNotInWaiting(INDEX_TYPE my_id) {
			BYTE_TYPE counthigher = 0;
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(my_id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (ABeforeB(my_id, neighbor_id) && (InWaiting(neighbor_id) || mReady2Go.count(neighbor_id) != 0)) {
					if (InWaiting(neighbor_id)) printf("unserviced in waiting\n");
					if (mReady2Go.count(neighbor_id) != 0) printf("unserviced in ready2go\n");

					counthigher++;
				}
			}
			return counthigher;
		}

		void UpdatePriorCountOrDelete(INDEX_TYPE id) {
			if (!InPriors(id)) {
				printf("ERROR, why do we get here\n");
				return;
			}
			PriorElement& e = mPriors[id];
			if (e.count_unprocessed_afters == 1) {
				//printf("erasing...\n");
				if (CountAftersNotInWaiting(id) != 0) {
					printf("ERROR: CountAftersNotInWaiting is not zero!!!\n"); // debug
				}
				mPriors.erase(id);
				return;
			}
			// means that a "prior" is waiting for one less! - when this reaches zero, we delete the prior
			//printf("size = %d\n", e.count_unprocessed_afters);
			e.count_unprocessed_afters--;
			//if (CountAftersNotInWaiting(id) != e.count_unprocessed_afters) printf("asdf;laksjdf;laskjdlkj\n");
		}
		void displayNeighborhood(INDEX_TYPE id) {
			//vector<pair<char, INDEX_TYPE> > befores;
			printf("w%d - ", my_wave_id);
			if (InPriors(id)) {
				printf("id %d = p %d:\n", id, mPriors[id].count_unprocessed_afters);
			}
			else if (InWaiting(id)) {
				printf("id %d = w %d:\n", id, mWaiting[id]);
			}
			else if (mReady2Go.count(id) != 0) {
				printf("id %d = r 0:\n", id);
			}
			else if (InProcessed(id)) {
				printf("id %d = processed and removed %d\n", id);
			}
			else {
				printf("id %d = unclassified %d:\n", id);
			}
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (ABeforeB(neighbor_id, id)) {
					if (InPriors(neighbor_id)) {
						printf("\t<-- p %d[%d]\n", neighbor_id, mPriors[neighbor_id].count_unprocessed_afters);
					}
					else if (InWaiting(neighbor_id)) {
						printf("\t<-- w %d[%d]\n", neighbor_id, mWaiting[neighbor_id]);
					}
					else if (mReady2Go.count(neighbor_id) != 0) {
						printf("\t<-- r %d[0]\n", neighbor_id);
					}
					else if (InProcessed(neighbor_id)) {
						printf("\t<-- processed and removed %d\n", neighbor_id);
					}
					else {
						printf("\t<-- missing %d\n", neighbor_id);
					}
				}
			}

			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (!ABeforeB(neighbor_id, id)) {
					if (InPriors(neighbor_id)) {
						printf("\t--> p %d[%d]\n", neighbor_id, mPriors[neighbor_id].count_unprocessed_afters);
					}
					else if (InWaiting(neighbor_id)) {
						printf("\t--> w %d[%d]\n", neighbor_id, mWaiting[neighbor_id]);
					}
					else if (mReady2Go.count(neighbor_id) != 0) {
						printf("\t--> r %d[0]\n", neighbor_id);
					}
					else if (InProcessed(neighbor_id)) {
						printf("\t--> processed and removed %d\n", neighbor_id);
					}
					else {
						printf("\t--> missing %d\n", neighbor_id);
					}
				}
			}


		}


		BYTE_TYPE CountBeforesNotInPrior(INDEX_TYPE my_id) {
			BYTE_TYPE countlower = 0;
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(my_id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (ABeforeB(neighbor_id, my_id)) {
					if (!InPriors(neighbor_id)) {
						countlower++;
					}
					//else {
					//	mPriors[neighbor_id].count_unprocessed_afters++;
					//}
				}
			}
			return countlower;
		}

		void AddEnqueuedIdToPriorCount(INDEX_TYPE my_id, INDEX_TYPE fromid) {
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(my_id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (neighbor_id != fromid && ABeforeB(neighbor_id, my_id)) {
					if (InPriors(neighbor_id)) {
						mPriors[neighbor_id].count_unprocessed_afters++;
					}
				}
			}
		}
		//void InsertIntoWaiting(INDEX_TYPE id, BYTE_TYPE count) {
		//	DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
		//	for (nit.begin(my_id); nit.valid(); nit.advance()) {
		//		INDEX_TYPE neighbor_id = nit.value();
		//		if (ABeforeB(neighbor_id, my_id) && !InPriors(neighbor_id)) {
		//			countlower++;
		//		}
		//	}
		//	return countlower;
		//}




		void EnqueueID(INDEX_TYPE my_id, INDEX_TYPE fromid) {
			// count the number of vertices of my_id that are lower and not already processed
			BYTE_TYPE countlower = CountBeforesNotInPrior(my_id);
			// if there are none, then it means the vertex is ready to go!
			if (countlower == 0) {
				//mReady2Go.push(my_id);
				if (InUnprocessedFront(my_id)) printf("SDF:KLSDJFSDL:KFJSD:KLFJSL:DKFJS:DKLJF\n");
				mReady2Go.insert(my_id);
			}
			else {
				// otherwise, we have to make this vertex wait
				if (InWaiting(my_id)) { // debug
					printf("ERROROR: %d already in waiting and is added again!\n"); // debug
				} // debug
				if (InReady(my_id) != 0) {
					printf("ERROROR: %d already in ready and is added again!\n"); // debug
				} // debug

				AddEnqueuedIdToPriorCount(my_id, fromid);
				mWaiting[my_id] = countlower;
			}
		}

		void UpdateNeighborWaitingCountOrInsert(INDEX_TYPE id, INDEX_TYPE fromid) {
			// id is already in our waiting set, so check if its ready or decrement its counter
			if (InReady(id)) printf("showa neighbor of just processed already ready\n");
			if (InWaiting(id)) {
				// if it only was waiting on current node, then we can add to ready-to-go
				if (mWaiting[id] == 1) {
					mWaiting.erase(id);
					printf("\t... erasing %d\n", id);
					//mReady2Go.push(id);
					if (mReady2Go.count(id) != 0) printf("SDF:KLSDJFSDL:KFJSD:KLFJSL:DKFJS:DKLJF\n");
					mReady2Go.insert(id);

					// debug
					if (CountBeforesNotInPrior(id) != 0) printf("ERRORORORORORRO: asdfasdfasf\n"); // debug
				}
				else {
					mWaiting[id]--;
				}
			}
			else {
				// this never existed so we have to enqueue it
				EnqueueID(id, fromid);
			}
		}

		set<INDEX_TYPE> mHackProcessed;

		bool InWaiting(INDEX_TYPE id) { return mWaiting.count(id) != 0; }
		bool InUnprocessedFront(INDEX_TYPE id) { return InWaiting(id) || InReady(id); }
		bool InProcessed(INDEX_TYPE id) { return mHackProcessed.count(id) != 0; }
		bool InPriors(INDEX_TYPE id) { return mPriors.count(id) != 0; }
		bool InReady(INDEX_TYPE id) { return mReady2Go.count(id) != 0; }
		void MarkAsReady(INDEX_TYPE id) { mReady2Go.insert(id); }
		void RemoveFromReady(INDEX_TYPE id) { mReady2Go.erase(id); }


	public:

		//// start with a legal state
		////-- each prior has count of afters depending on it
		////-- each waiting has a count of unprocessed befores, this is always > 0, else moved to readytogo
		////-- there is a "ready to go" list where each node has ALL its befores in priors
		//// to add a seed
		//// first check mpriors - if it's there, then simply replace the contents -- although this should never happen
		//// next check waiting - if its there, then "process" the element 
		////
		//// difference between processing elements and seeds is that a seed already has its payload - not true
		//// they both have the payload!
		//// actually the payload is the RESULT of a processing, so waiting and ready2go don't even have them.
		//// so priors need to know how many afters are depending on them
		//// waiting only need to know how many befores they depend on
		//// ready to go is just an index
		HashingWaveFront(DEPGRAPH_TYPE* global_graph) : mGlobalGraph(global_graph) { my_wave_id = count++; }


		void SetSeed(INDEX_TYPE id, PAYLOAD_TYPE& p) {
			mReady2Go.insert(id); // this is a faster hack than EnqueueID(id); MarkAsReady(id);
			SetAsProcessed(id, p);
		}

		void SetAsProcessed(INDEX_TYPE id, PAYLOAD_TYPE& p) {
			// this should not be in either the waiting or the priors
			if (InWaiting(id) || InPriors(id) || !InReady(id)) printf("WHOATHEREANOTHERNELLY\n");// debug
			RemoveFromReady(id);
			if (CountBeforesNotInPrior(id) != 0) printf("processed missing a prior\n");
			// simply add it to the priors list
			PriorElement ne(id, p);
			mPriors[id] = ne;
			mHackProcessed.insert(id);
			// now update all neighbors
			// first decrement all priors that no longer are needed

			// increment first
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (!ABeforeB(neighbor_id, id)) {
					// something that comes after THIS id should not already be in the priors list
					if (InPriors(neighbor_id)) printf("ANOTHERWHOAWHAT2\n"); // debug
					UpdateNeighborWaitingCountOrInsert(neighbor_id, id);
					// for each neighbor that is after id, we will have to keep this prior around until they are all processed
					mPriors[id].count_unprocessed_afters++;
				}
			}

			//DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				// pick out the befores
				if (ABeforeB(neighbor_id, id)) {
					// this should never happen
					if (InUnprocessedFront(neighbor_id) != 0) printf("ANOTHERWHOAWHAT\n");//debug
					// if it exists decrement or delete i
					UpdatePriorCountOrDelete(neighbor_id);
				}
			}
		}

		bool FindNextToProcess(INDEX_TYPE& ret) {
			// if there are no readytogo, then return
			if (mReady2Go.size() == 0) {
				return false;
			}
#ifdef QUEQUE
			//ret = mReady2Go.front();
			ret = *mReady2Go.begin();
#else			
			ret = mReady2Go.top();
#endif
			//mReady2Go.pop();
			//mReady2Go.erase(mReady2Go.begin());

			// debug
			int actualcount = 0;
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(ret); nit.valid(); nit.advance()) {
				INDEX_TYPE nid = nit.value();
				if (ABeforeB(nid, ret)) {
					if (!InPriors(nid)) {
						printf("ERROR: returning item from readytogo, missing prior: %d -> %d \n", ret, nid);
					}
				}
			}
			return true;
		}

		PAYLOAD_TYPE& GetPayload(INDEX_TYPE id) {
			return mPriors[id].payload;
		}


		// create a new wavefront that has all the waiting elements reachable from base_id
		// as well as all the priors needed
		// incorrect if readytogo is not empty
		HashingWaveFront* createFromIncompleteSeed(INDEX_TYPE base_id) {
			if (mReady2Go.size() != 0) printf("WHOA could be problem\n"); // debug
			HashingWaveFront* newwave = new HashingWaveFront(mGlobalGraph);

			newwave->mHackProcessed.insert(mHackProcessed.begin(), mHackProcessed.end());
			// use a stack to do dfs in the dependency graph starting from base_id
			// we will traverse on all items from waiting that are reachable, hence depend on base_id
			stack<INDEX_TYPE> front;
			// start growing by adding to the front
			front.push(base_id);
			// keep a set of seen things so we only traverse each item once
			unordered_set<INDEX_TYPE> seen; // set of seen things
			while (!front.empty()) {
				INDEX_TYPE front_id = front.top();
				front.pop();
				// if it has already been processed, skip
				if (seen.count(front_id) > 0) continue;
				// mark as processed
				seen.insert(front_id);
				// debug - should never be already in the waiting list
				if (newwave->InWaiting(front_id)) printf("WHOATHERENELLY\n");
				// now insert it, so far with no dependants
				newwave->mWaiting[front_id] = 0;
				// now we will add dependent neighbors to the traversal,
				// and add all priors to this wavefront
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				for (nit.begin(front_id); nit.valid(); nit.advance()) {
					INDEX_TYPE neighbor_id = nit.value();
					// add higher neighbors to the traversal only if they are already
					// in the base wavefront - not to be confused with NEWWAVE's waiting
					if (ABeforeB(front_id, neighbor_id)) {
						if (InWaiting(neighbor_id)) {
							front.push(neighbor_id);
						}
					}
					else {
						// otherwise the neighbor is lower
						if (!InPriors(neighbor_id)) {
							// if there was no prior element to it in this wave, it will be waiting for one in the new
							// wave as well, so increment id of new wave waiting id.

							if (InProcessed(neighbor_id)) { printf("ERROR: there is a processed lower neighbor to a waiting item that is not in priors\n"); }

							newwave->mWaiting[front_id]++;
						}
						else {
							// it is WAS found in the priors of the original wavefront, then add it to 
							// the newwave's priors, or if its already there, then increment the counter
							if (!newwave->InPriors(neighbor_id)) {
								newwave->mPriors[neighbor_id] = mPriors[neighbor_id];
								newwave->mPriors[neighbor_id].count_unprocessed_afters = 1;
							}
							else {
								newwave->mPriors[neighbor_id].count_unprocessed_afters++;
							}
						}
					}
				}
			}
			//debug
			if (!newwave->ConsistencyCheck()) printf("newwave inconsistency\n"); // debug

			return newwave;

		}


		// does a "deep" split of this wavefront - into "independent" wavefronts - does memory copy to new items
		// -- finds "seeds" among all WAITING and creates a new wavefront
		// -- for each "seed" finds all reachable that are in WAITING, and adds to seeds wavefront
		// -- for each new wavefront, find all PRIORS and add to the wavefront 
		// -- finally, update counts for each! - maybe can set PRIOR count when adding
		INDEX_TYPE FindLowestSeed() {

			bool haslowest = false;
			INDEX_TYPE lowest;
			// check all waitings - any that have ONLY external dependency becomes a bottleneck point
			for (auto it = mWaiting.begin(); it != mWaiting.end(); it++) {
				INDEX_TYPE waitid = (*it).first;
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				// if it has an internal waiting then it is NOT a seed
				bool has_internal = false;
				for (nit.begin(waitid); nit.valid(); nit.advance()) {
					INDEX_TYPE neighbor_id = nit.value();
					//if (mGlobalGraph->Before(nfid, waitid) &&
					//	mPriors.count(nfid) == 0 &&
					//	mWaiting.count(nfid) == 0 &&
					//	tSeeds.count(nfid) == 0) {
					//	//printf("found obstruction\n");
					//	//obstructions.push_back(createFromSeed(nfid));
					//	obstructions.push_back(pair<INDEX_TYPE, HashingWaveFront*>(waitid, createFromSeed(waitid)));

					//}				

					// look for things that ONLY depend on an external thing
					if (ABeforeB(neighbor_id, waitid) &&
						InWaiting(neighbor_id)) {
						has_internal = true;
						//printf("found obstruction\n");
						//obstructions.push_back(createFromSeed(nfid));
					}
				}
				if (!has_internal) {
					if (haslowest) {
						if (mGlobalGraph->Before(waitid, lowest)) {
							lowest = waitid;
						}
					}
					else {
						lowest = waitid;
						haslowest = true;
					}
				}

			}
			return lowest;
		}


		// does a "deep" split of this wavefront - into "independent" wavefronts - does memory copy to new items
		// -- finds "seeds" among all WAITING and creates a new wavefront
		// -- for each "seed" finds all reachable that are in WAITING, and adds to seeds wavefront
		// -- for each new wavefront, find all PRIORS and add to the wavefront 
		// -- finally, update counts for each! - maybe can set PRIOR count when adding
		void deepSplitWaveFront(vector<HashingWaveFront*>& res, vector < pair<INDEX_TYPE, HashingWaveFront*> > & obstructions) {
			//unordered_set<INDEX_TYPE> tSeeds;
			//			while (mReady2Go.size() > 0) {
			//#ifdef QUEQUE
			//				INDEX_TYPE seed = mReady2Go.front();
			//#else
			//				INDEX_TYPE seed = mReady2Go.top();
			//#endif
			//				mReady2Go.pop();
			//				tSeeds.insert(seed);
			//				res.push_back(createFromSeed(seed));
			//
			//			}

			// check all waitings - any that have ONLY external dependency becomes a bottleneck point
			for (auto it = mWaiting.begin(); it != mWaiting.end(); it++) {
				INDEX_TYPE waitid = (*it).first;
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				// if it has an internal waiting then it is NOT a seed
				bool has_internal = false;
				for (nit.begin(waitid); nit.valid(); nit.advance()) {
					INDEX_TYPE neighbor_id = nit.value();
					//if (mGlobalGraph->Before(nfid, waitid) &&
					//	mPriors.count(nfid) == 0 &&
					//	mWaiting.count(nfid) == 0 &&
					//	tSeeds.count(nfid) == 0) {
					//	//printf("found obstruction\n");
					//	//obstructions.push_back(createFromSeed(nfid));
					//	obstructions.push_back(pair<INDEX_TYPE, HashingWaveFront*>(waitid, createFromSeed(waitid)));

					//}				

					// look for things that ONLY depend on an external thing
					if (ABeforeB(neighbor_id, waitid) &&
						InWaiting(neighbor_id)) {
						has_internal = true;
						//printf("found obstruction\n");
						//obstructions.push_back(createFromSeed(nfid));
					}
				}
				if (!has_internal) {
					obstructions.push_back(pair<INDEX_TYPE, HashingWaveFront*>(waitid, createFromIncompleteSeed(waitid)));
				}

			}
		}

		// does a "deep" merge and deletes the other wavefront - inserts items and updates counts
		// -- add all other's PRIORS, adding up counts of PRIORS that already existed
		// -- add all other's WAITING, subtracting from neighboring PRIORS if the WAITING id already existed
		void mergeAndConsumeWaveFront(HashingWaveFront* other) {
			//debug
			if (other->mReady2Go.size() != 0) printf("other ready to go not empty...\n"); // debug
			if (mReady2Go.size() != 0) printf("ready to go not empty...\n"); // debug
			for (auto it = other->mPriors.begin(); it != other->mPriors.end(); it++) { // debug
				if (mPriors.count((*it).first) != 0) printf("asdf already here...\n"); // debug
			} // debug
			// add priors data structures
			mHackProcessed.insert(other->mHackProcessed.begin(), other->mHackProcessed.end());

			mPriors.insert(other->mPriors.begin(), other->mPriors.end());


			mWaiting.insert(other->mWaiting.begin(), other->mWaiting.end());
			//		for (auto it = other->mPriors.begin(); it != other->mPriors.end(); it++) {
			//			INDEX_TYPE id = (*it).first;
			//			PriorElement& p = (*it).second;
			//			if (mPriors.count(id) == 0) {
			//				mPriors[id] = p;
			//			}
			////			else {
			//	//			mPriors[id].count_unprocessed_afters += p.count_unprocessed_afters;
			////			}
			//		}
			//		// set all priors to zero
			for (auto it = mPriors.begin(); it != mPriors.end(); it++) {
				(*it).second.count_unprocessed_afters = 0;
			}
			// now add all the waitings
			vector<INDEX_TYPE> newreadytogo;
			for (auto waiting_iterator = mWaiting.begin(); waiting_iterator != mWaiting.end(); waiting_iterator++) {
				INDEX_TYPE waiting_id = waiting_iterator->first;
				int tmpcount = 0; //debug
				waiting_iterator->second = 0; // how many this item is waiting on 
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				for (nit.begin(waiting_id); nit.valid(); nit.advance()) {
					INDEX_TYPE nid = nit.value();
					// for everything that is not a before with prior payload computed, count it up
					if (ABeforeB(nid, waiting_id)) {
						if (!InPriors(nid)) {
							waiting_iterator->second++;
							tmpcount++;
						}
					}
				}
				if (waiting_iterator->second == 0) {
					//mReady2Go.push(waiting_id);
					if (mReady2Go.count(waiting_id) != 0) printf("SDF:KLSDJFSDL:KFJSD:KLFJSL:DKFJS:DKLJF\n");
					mReady2Go.insert(waiting_id);

					newreadytogo.push_back(waiting_id);
					// debug
					if (CountBeforesNotInPrior(waiting_id) != 0) printf("ERORORORORO: has %d missing priors...\n", CountBeforesNotInPrior(waiting_id));
				}
				// now fill in ready-to-go
			}

			/// count number of unprocessed later neighbors
			for (auto prior_it = mPriors.begin(); prior_it != mPriors.end(); prior_it++) {
				//if (CountAftersNotInWaiting((*prior_it).first) != 0) printf("AHA! -------------------------\n");
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				for (nit.begin(prior_it->first); nit.valid(); nit.advance()) {
					INDEX_TYPE neighbor_id = nit.value();
					if (ABeforeB(prior_it->first, neighbor_id)) {
						if (InUnprocessedFront(neighbor_id)){
							prior_it->second.count_unprocessed_afters++;
						}
						//else if (! InProcessed(neighbor_id)) {
						//	printf("NOT_IN_UNPROCESSED FRONT:");
						//	printf("D1: "); displayNeighborhood(prior_it->first);
						//	printf("D2: "); displayNeighborhood(neighbor_id);
						//}
					}
				}
			}

			for (auto it = newreadytogo.begin(); it != newreadytogo.end(); it++) {
				mWaiting.erase(*it);
				// debug
				if (InPriors(*it)) printf("ERROR: %d already in priors...\n", *it);
			}
			if (!ConsistencyCheck()) printf("after merge not consistent\n");

		}

		int FrontSize() { return mWaiting.size() + mReady2Go.size(); }

		bool ConsistencyCheck() {
			bool valid = true;
			//printf("ready2gosize[%d] = %d\n", this, mReady2Go.size());
			//queue<INDEX_TYPE> tmpq = mReady2Go;
			//set<INDEX_TYPE> tmpq;
			// check ready2go
			set<INDEX_TYPE> r2gids;
			//while (!tmpq.empty()) {
			//	//printf("gothere\n");
			//	INDEX_TYPE id = tmpq.front();
			//	tmpq.pop();
			for (auto it = mReady2Go.begin(); it != mReady2Go.end(); it++){
				INDEX_TYPE id = *it;
				r2gids.insert(id);

				// ready to go must not be in waiting or priors
				if (mPriors.count(id) != 0) {
					printf("ERROR: ready-to-go id %d exists in priors\n", id); valid = false;
				}
				if (mWaiting.count(id) != 0) {
					printf("ERROR: ready-to-go id %d exists in waiting\n", id); valid = false;
				}

				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				for (nit.begin(id); nit.valid(); nit.advance()) {
					INDEX_TYPE nid = nit.value();
					if (mGlobalGraph->Before(nid, id)) {
						// must exist in priors
						if (mPriors.count(nid) == 0) {
							printf("ERROR: ready-to-go id %d missing %d in priors\n", id, nid); valid = false;
						}
						if (mWaiting.count(nid) != 0) {
							printf("ERROR: ready-to-go id %d has waiting %d in befores\n", id, nid); valid = false;
						}
					}
					else {
						if (mPriors.count(nid) != 0) {
							printf("ERROR: ready-to-go id %d exists with AFTER %d in priors\n", id, nid); valid = false;
						}
					}
				}

			}
			//printf("ready2gosize[%d] = %d\n", this, mReady2Go.size());

			// check that each waiting has right number of priors
			for (auto it = mPriors.begin(); it != mPriors.end(); it++) {

				INDEX_TYPE id = (*it).first;

				if (r2gids.count(id) != 0) {
					printf("ERROR: prior id %d exists in r2gids\n", id); valid = false;
				}
				if (mWaiting.count(id) != 0) {
					printf("ERROR: prior id %d exists in waiting\n", id); valid = false;
				}

				int tcountwaiting = 0;
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				for (nit.begin(id); nit.valid(); nit.advance()) {
					INDEX_TYPE nid = nit.value();
					if (mGlobalGraph->Before(nid, id)) {
						if (mWaiting.count(nid) != 0) {
							printf("ERROR: prior id %d exists with BEFORE %d in waiting\n", id, nid); valid = false;
						}
					}
					else {
						//// disable this check because priors a, b , a->b, a->c, b->d, d->c
						//// if a b are in priors, and d is processed, b can be removed, while "a" still has dependent c
						//if (mWaiting.count(nid) == 0 &&
						//	mPriors.count(nid) == 0 &&
						//	r2gids.count(nid) == 0) {
						//	printf("ERROR: prior id %d exists but its AFTER %d does not exist\n", id, nid); return false;
						//}
						// all ids in priors have been processed, so only count unprocessed stuff;
						if (mWaiting.count(nid) != 0 ||
							r2gids.count(nid) != 0){
							tcountwaiting++;
						}
					}
				}

				if ((*it).second.count_unprocessed_afters != tcountwaiting) {
					printf("ERROR: prior id %d has %d in count, %d in actual\n", id, (*it).second.count_unprocessed_afters, tcountwaiting); valid = false;
					displayNeighborhood(id);
				}
			}

			// check that each waiting has right number of priors
			for (auto it = mWaiting.begin(); it != mWaiting.end(); it++) {
				INDEX_TYPE id = (*it).first;

				if (r2gids.count(id) != 0) {
					printf("ERROR: waiting id %d exists in r2gids\n", id); valid = false;
				}
				if (mPriors.count(id) != 0) {
					printf("ERROR: waiting id %d exists in priors\n", id); valid = false;
				}

				int tcountwaiting = 0;
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				for (nit.begin(id); nit.valid(); nit.advance()) {
					INDEX_TYPE nid = nit.value();
					if (mGlobalGraph->Before(nid, id)) {

						//// this has a saddle point! this test does not work - since we keep ALL prior processed
						//// at an annulus this will break
						//if (InProcessed(nid) && !(InPriors(nid) || InUnprocessedFront(nid))) {
						//	printf("ERROR: waiting id %d has BEFORE %d that does not exist infront\n", id, nid);
						//	printf("DETAIL 1: ");
						//	displayNeighborhood(id);
						//	printf("DETAIL 2: ");
						//	displayNeighborhood(nid);
						//	valid = false;
						//}

						if (mPriors.count(nid) == 0) {
							tcountwaiting++;
						}



					}
					else {
						if (mPriors.count(nid) != 0) {
							printf("ERROR: waiting id %d exists with AFTER %d in priors\n", id, nid); valid = false;
						}
						if (r2gids.count(nid) != 0) {
							printf("ERROR: waiting id %d exists with AFTER %d in r2go\n", id, nid); valid = false;
						}
					}
				}

				if ((*it).second != tcountwaiting) {
					printf("ERROR: waiting id %d has %d in count, %d in actual\n", id, (*it).second, tcountwaiting); valid = false;
				}
			}
			return valid;
		}


	};






	class TestingDependencyGraphParallelTraversal {
	protected:
		DEPGRAPH_TYPE* mD;
		bool has_explicit_thred_num;
		int mNumThreads;
	public:
		TestingDependencyGraphParallelTraversal(DEPGRAPH_TYPE* d) : mD(d), has_explicit_thred_num(false) {


		}
		void SetNumThreads(int i) { mNumThreads = i; }
		void doStuff() {

			// create a new wavefront for every "minimum"		
			unordered_map<INDEX_TYPE, HashingWaveFront*> waves;


			std::vector<INDEX_TYPE> topo_index_partition;
			INDEX_TYPE tTotalNumToDo = 0;
			int num_threads = has_explicit_thred_num;
#pragma omp parallel
			{
#pragma omp single
				{
					if (!has_explicit_thred_num) num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mD->IndexSpaceSize(), num_threads, topo_index_partition);
				}

				int thread_num = omp_get_thread_num();
				INDEX_TYPE tLocalNumToDo = 0;
				DEPGRAPH_TYPE::vertex_iterator vit(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (vit.begin(); vit.valid(); vit.advance()) {
					INDEX_TYPE id = vit.value();
					tLocalNumToDo++;
					bool ismin = true;
					DEPGRAPH_TYPE::neighbor_iterator nit(mD);
					for (nit.begin(id); nit.valid(); nit.advance()) {
						if (mD->Before(nit.value(), id)) {
							ismin = false;
							break;
						}
					}
					if (ismin) {
						//(*ip) = 0;
#pragma omp critical
						{
							int zero = 0;
							HashingWaveFront* w = new HashingWaveFront(mD);
							printf("created w%d\n", w->my_wave_id);
							w->SetSeed(id, zero);

							waves[id] = w;
							printf("thread %d doing consistency check %d\n", thread_num, id);
							w->ConsistencyCheck();
						}
					}
				}
#pragma omp critical
				{
					printf("thread %d did %d ids\n", thread_num, tLocalNumToDo);
					tTotalNumToDo += tLocalNumToDo;
				}
			}
			printf("total number of edges: %d\n", tTotalNumToDo);

			char fname[2048];

			set<INDEX_TYPE> global_seen;
			INDEX_TYPE tTotalProcCount = 0;
			int round = 0;
			bool bcont = true;
			while (bcont) {

				if (round > 30) return;
				printf("doint round %d, parallelism = %d\n", round, waves.size());

				sprintf(fname, "test_s_%d.txt", round);
				printf("writing temp file %s\n", fname);
				FILE* ftest = fopen(fname, "w");

				printf("ROUND %d consitency check -- a\n", round);
				for (auto hit = waves.begin(); hit != waves.end(); hit++) {

					HashingWaveFront* h = (*hit).second;


					if (!h->ConsistencyCheck()) return;

					//printf("doing %d\n", (*hit).first);
					for (auto it = h->mPriors.begin(); it != h->mPriors.end(); it++) {
						fprintf(ftest, "p %d %d\n", (*hit).first, (*it).first);
					}
					for (auto it = h->mWaiting.begin(); it != h->mWaiting.end(); it++) {
						fprintf(ftest, "w %d %d\n", (*hit).first, (*it).first);
					}
					fprintf(ftest, "r %d %d\n", (*hit).first, (*hit).first);
				}
				round++;
				fclose(ftest);
				printf("done\n");


				//// do first expansion
				int tLocalProcCount = 0;
				int randocount = 0;
				for (auto it = waves.begin(); it != waves.end(); it++) {
					HashingWaveFront* h = (*it).second;
					printf("expanding w%d\n", h->my_wave_id);
					if (!h->ConsistencyCheck()) {
						printf("returning point 1...\n");
						return;
					}
					int counter = 0;
					bool cont = true;
					while (cont) {
						INDEX_TYPE id;
						cont = h->FindNextToProcess(id);
						if (!cont) break;
						if (global_seen.count(id) != 0) printf("WHOA %d already processed!\n", id);
						if (h->InPriors(id) || h->InWaiting(id)) printf("WHOA exists in priorowaiasdfk\n");
						// this is where work goes for payload
						tLocalProcCount++;

						if (!h->ConsistencyCheck()) {
							printf("returning point 1a...\n");
							return;
						}

						h->displayNeighborhood(id);
						h->SetAsProcessed(id, randocount);
						h->displayNeighborhood(id);
						global_seen.insert(id);
						counter++;
						randocount++;
						if (!h->ConsistencyCheck()) {
							printf("returning point 2 after %d processed... %d\n", counter, id);
							return;
						}

					}
				}

				sprintf(fname, "test_s_%d.txt", round);
				printf("writing temp file %s\n", fname);
				ftest = fopen(fname, "w");

				printf("ROUND %d consitency check -- b\n", round);
				for (auto hit = waves.begin(); hit != waves.end(); hit++) {

					HashingWaveFront* h = (*hit).second;
					if (!h->ConsistencyCheck()) return;
					//printf("doing %d\n", i);
					for (auto it = h->mPriors.begin(); it != h->mPriors.end(); it++) {
						fprintf(ftest, "p %d %d\n", (*hit).first, (*it).first);
					}
					for (auto it = h->mWaiting.begin(); it != h->mWaiting.end(); it++) {
						fprintf(ftest, "w %d %d\n", (*hit).first, (*it).first);
					}
					fprintf(ftest, "r %d %d\n", (*hit).first, (*hit).first);
				}
				round++;
				fclose(ftest);
				printf("done\n");


				tTotalProcCount += tLocalProcCount;
				printf("did %d of %d...\n", tTotalProcCount, tTotalNumToDo);

				if (tTotalProcCount >= tTotalNumToDo) {
					printf("did %d of %d, returning\n", tTotalProcCount, tTotalNumToDo);
					bcont = false;
					return;
				}

				vector < pair<INDEX_TYPE, HashingWaveFront*> > ohG;
				//for (auto it = waves.begin(); it != waves.end(); it++) {
				//	HashingWaveFront* h = (*it).second;
				//	vector<HashingWaveFront*> th;
				//	vector < pair<INDEX_TYPE, HashingWaveFront*> > oh;
				//	printf("splitting w%d size=%d\n", h->my_wave_id, h->mWaiting.size() + h->mReady2Go.size());
				//	h->deepSplitWaveFront(th, oh);
				//	delete h;
				//	if (th.size() > 0)  {
				//		printf("whoa should not get here\n");
				//	}

				//	for (auto sit = oh.begin(); sit != oh.end(); sit++) {
				//		printf("  into w%d at %d\n", (*sit).second->my_wave_id, (*sit).first);
				//		if (!(*sit).second->ConsistencyCheck()) return;
				//	}

				//	ohG.insert(ohG.end(), oh.begin(), oh.end());
				//}

				for (auto it = waves.begin(); it != waves.end(); it++) {
					HashingWaveFront* h = (*it).second;
					if (h->FrontSize() == 0) continue;
					INDEX_TYPE obid = h->FindLowestSeed();

					printf(" lowest = %d w%d size=%d\n", obid, h->my_wave_id, h->mWaiting.size() + h->mReady2Go.size());
					ohG.push_back(pair<INDEX_TYPE, HashingWaveFront*>(obid, h));
				}
				sprintf(fname, "test_s_%d.txt", round);
				printf("writing temp file %s\n", fname);
				ftest = fopen(fname, "w");

				printf("ROUND %d consitency check -- c\n", round);
				for (auto hit = ohG.begin(); hit != ohG.end(); hit++) {

					HashingWaveFront* h = (*hit).second;
					if (!h->ConsistencyCheck()) return;
					h->ConsistencyCheck();
					//printf("doing %d\n", i);
					for (auto it = h->mPriors.begin(); it != h->mPriors.end(); it++) {
						fprintf(ftest, "p %d %d\n", (*hit).first, (*it).first);
					}
					for (auto it = h->mWaiting.begin(); it != h->mWaiting.end(); it++) {
						fprintf(ftest, "w %d %d\n", (*hit).first, (*it).first);
					}
					fprintf(ftest, "r %d %d\n", (*hit).first, (*hit).first);
				}
				round++;
				fclose(ftest);
				printf("done\n");


				//vector<HashingWaveFront*> nexttogo;
				unordered_map<INDEX_TYPE, HashingWaveFront*> merges;
				set<HashingWaveFront*> seenh;
				printf("about to merge..\n");
				for (auto mit = ohG.begin(); mit != ohG.end(); mit++) {

					INDEX_TYPE mid = (*mit).first;
					HashingWaveFront* h = (*mit).second;
					if (seenh.count(h) != 0) printf("ERRORORORORRORO: h seen\n");
					if (merges.count(mid) != 0) {
						printf("merging w%d and w%d at point %d\n", merges[mid]->my_wave_id, h->my_wave_id, mid);
						merges[mid]->mergeAndConsumeWaveFront(h);
						delete h;
					}
					else {
						merges[mid] = h;
					}

				}

				////vector<HashingWaveFront*> nexttogo;
				//unordered_map<INDEX_TYPE, HashingWaveFront*> merges;
				//printf("about to merge..\n");
				//for (auto it = waves.begin(); it != waves.end(); it++) {
				//	HashingWaveFront* h = (*it).second;
				//	vector<HashingWaveFront*> th;
				//	vector < pair<INDEX_TYPE, HashingWaveFront*> > oh;
				//	h->deepSplitWaveFront(th, oh);
				//	if (th.size() > 0)  {
				//		printf("whoa should not get here\n");
				//	}

				//	for (auto mit = oh.begin(); mit != oh.end(); mit++) {
				//		INDEX_TYPE mid = (*mit).first;
				//		HashingWaveFront* h = (*mit).second;

				//		if (merges.count(mid) != 0) {
				//			merges[mid]->mergeAndConsumeWaveFront(h);
				//		}
				//		else {
				//			merges[mid] = h;
				//		}
				//	}
				//}

				waves.clear();
				waves = merges;
				printf("done  merge..\n");
				//round++;

				if (waves.size() == 0) return;


			}


			////// do first expansion
			//int randocount = 0;
			//for (auto it = waves.begin(); it != waves.end(); it++) {
			//	HashingWaveFront* h = (*it).second;
			//	bool cont = true;
			//	while (cont) {
			//		INDEX_TYPE id;
			//		cont = h->GetNextToProcess(id);
			//		if (!cont) break;
			//		// this is where work goes for payload
			//		h->SetAsProcessed(id, randocount);
			//		randocount++;

			//	}
			//}


			//vector<HashingWaveFront*> h2;
			//vector < pair<INDEX_TYPE, HashingWaveFront*> > oh2;
			//for (auto it = waves.begin(); it != waves.end(); it++) {
			//	HashingWaveFront* h = (*it).second;
			//	vector<HashingWaveFront*> th;
			//	vector < pair<INDEX_TYPE, HashingWaveFront*> > oh;
			//	h->deepSplitWaveFront(th, oh);
			//	if (th.size() > 0)  {
			//		printf("whoa should not get here\n");
			//	}

			//	for (auto mit = oh.begin(); mit != oh.end(); mit++) {
			//		if (merges.count((*mit).first) != 0) {
			//			(*mit).second->mergeAndConsumeWaveFront(merges[(*mit).first]);
			//			merges.erase((*mit).first);
			//			nexttogo.push_back((*mit).second);
			//		}
			//	}


			//	oh2.insert(oh2.end(), oh.begin(), oh.end());
			//	
			//}

			//FILE* ftest = fopen("testids_s.txt", "w");

			//for (auto hit = oh2.begin(); hit != oh2.end(); hit++) {

			//	HashingWaveFront* h = (*hit).second;
			//	//printf("doing %d\n", i);
			//	for (auto it = h->mPriors.begin(); it != h->mPriors.end(); it++) {
			//		fprintf(ftest, "p %d %d\n", (*hit).first, (*it).first);
			//	}
			//	for (auto it = h->mWaiting.begin(); it != h->mWaiting.end(); it++) {
			//		fprintf(ftest, "w %d %d\n", (*hit).first, (*it).first);
			//	}
			//	fprintf(ftest, "r %d %d\n", (*hit).first, (*hit).first);
			//}
			//fclose(ftest);
		}









	};




	class ExplicitWaveFrontTry2 {
	public:

		ExplicitWaveFrontTry2(DEPGRAPH_TYPE* global_graph) : mGlobalGraph(global_graph) {
#pragma omp critical 
			{
				my_wave_id = count++;
			}
		}

//		ExplicitWaveFrontTry2(DEPGRAPH_TYPE* global_graph, vector<ExplicitWaveFrontTry2*>& inputs) : mGlobalGraph(global_graph) { 
//#pragma omp critical 
//			{
//				my_wave_id = count++;
//			}
//
//			for (auto it = inputs.begin(); it != inputs.end(); it++) {
//				ExplicitWaveFrontTry2* other = (*it);
//				mPriors.insert(other->mPriors.begin(), other->mPriors.end());
//				mWaiting.insert(other->mWaiting.begin(), other->mWaiting.end());
//				mReady2Go.insert(other->mReady2Go.begin(), other->mReady2Go.end());
//				delete other;
//
//
//			}
//			vector<INDEX_TYPE> newreadytogo;
//			for (auto waiting_iterator = mWaiting.begin(); waiting_iterator != mWaiting.end(); waiting_iterator++) {
//				INDEX_TYPE waiting_id = *waiting_iterator;
//
//				if (CountBeforesNotInPrior(waiting_id) == 0) {
//					mReady2Go.insert(waiting_id);
//					newreadytogo.push_back(waiting_id);
//				}
//			}
//
//			for (auto it = newreadytogo.begin(); it != newreadytogo.end(); it++) {
//				mWaiting.erase(*it);
//			}
//
//
//		}
		
		
		
		
		
		
		struct PriorElement {
			PAYLOAD_TYPE payload;

			PriorElement() {
			}
			PriorElement(PAYLOAD_TYPE p) : payload(p) {
			}
		};


		// MEMBER VARIABLES

		DEPGRAPH_TYPE* mGlobalGraph;
		static int count;
		int my_wave_id;
		map<INDEX_TYPE, PriorElement > mPriors;
		set<INDEX_TYPE> mWaiting;
		set<INDEX_TYPE> mReady2Go;

		bool ABeforeB(INDEX_TYPE a, INDEX_TYPE b) { return mGlobalGraph->Before(a, b); }
		bool InWaiting(INDEX_TYPE id) { return mWaiting.count(id) != 0; }
		bool InUnprocessedFront(INDEX_TYPE id) { return InWaiting(id) || InReady(id); }
		bool InPriors(INDEX_TYPE id) { return mPriors.count(id) != 0; }
		bool InReady(INDEX_TYPE id) { return mReady2Go.count(id) != 0; }
		void MarkAsReady(INDEX_TYPE id) { mReady2Go.insert(id); }
		void RemoveFromReady(INDEX_TYPE id) { mReady2Go.erase(id); }
		static DenseLabeling<char>* mProcessed;

		BYTE_TYPE CountBeforesNotInPrior(INDEX_TYPE my_id) {
			BYTE_TYPE countlower = 0;
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(my_id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (ABeforeB(neighbor_id, my_id)) {
					if (!InPriors(neighbor_id)) {
						countlower++;
					}
				}
			}
			return countlower;
		}


		BYTE_TYPE CountUnprocessedAfters(INDEX_TYPE my_id) {
			BYTE_TYPE countlower = 0;
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(my_id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (ABeforeB(my_id, neighbor_id)) {
					if (InUnprocessedFront(neighbor_id)) {
						//if (! IsProcessed(neighbor_id)) {

						countlower++;
					}
				}
			}
			return countlower;
		}

		void EnqueueID(INDEX_TYPE my_id, INDEX_TYPE fromid) {
			BYTE_TYPE countlower = CountBeforesNotInPrior(my_id);
			if (countlower == 0) {
				mReady2Go.insert(my_id);
			}
			else {
				mWaiting.insert(my_id);
			}
		}

		void UpdateNeighborWaitingCountOrInsert(INDEX_TYPE id, INDEX_TYPE fromid) {
			if (InWaiting(id)) {
				if (CountBeforesNotInPrior(id) == 0) {
					mWaiting.erase(id);
					mReady2Go.insert(id);
				}
			}
			else if (!InReady(id)) {
				// this never existed so we have to enqueue it
				EnqueueID(id, fromid);
			}
		}

		void UpdatePriorCountOrDelete(INDEX_TYPE id) {
			if (CountUnprocessedAfters(id) == 0) {
				mPriors.erase(id);
				return;
			}
		}

		void SetAsProcessed(INDEX_TYPE id, PAYLOAD_TYPE& p) {
			RemoveFromReady(id);
			PriorElement ne(p);
			mPriors[id] = ne;
			mProcessed->SetLabel(id, 1);

			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (ABeforeB(id, neighbor_id)) {
					UpdateNeighborWaitingCountOrInsert(neighbor_id, id);
				}
			}
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				// pick out the befores
				if (ABeforeB(neighbor_id, id)) {
					UpdatePriorCountOrDelete(neighbor_id);
				}
			}
		}
		void SetSeed(INDEX_TYPE id, PAYLOAD_TYPE& p) {
			mReady2Go.insert(id); // this is a faster hack than EnqueueID(id); MarkAsReady(id);
			SetAsProcessed(id, p);
		}

		bool FindNextToProcess(INDEX_TYPE& ret) {
			if (mReady2Go.size() == 0) {
				return false;
			}
			ret = *mReady2Go.begin();

			// debug
			//DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			//for (nit.begin(ret); nit.valid(); nit.advance()) {
			//	INDEX_TYPE neighbor_id = nit.value();
			//	if (ABeforeB(neighbor_id, ret)) {
			//		if (!InPriors(neighbor_id)) {
			//			printf("WHOA processing element before it has all priors!\n");
			//		}
			//	}
			//}
			return true;
		}

		ExplicitWaveFrontTry2* createFromIncompleteSeed(INDEX_TYPE base_id, set<INDEX_TYPE>& marked) {
			ExplicitWaveFrontTry2* newwave;
#pragma omp critical 
			{
				newwave = new ExplicitWaveFrontTry2(mGlobalGraph);
			}
			stack<INDEX_TYPE> front;
			front.push(base_id);
			while (!front.empty()) {
				INDEX_TYPE front_id = front.top();
				front.pop();
				if (newwave->InUnprocessedFront(front_id)) continue;
				if (InWaiting(front_id)) {
					newwave->mWaiting.insert(front_id);
				}
				else if (InReady(front_id)) {
					newwave->mReady2Go.insert(front_id);
				}
				else {
					printf(" SJKLJJFJFJFJFJFJ dont know how i got here\n");
				}
				marked.insert(front_id);
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				for (nit.begin(front_id); nit.valid(); nit.advance()) {
					INDEX_TYPE neighbor_id = nit.value();
					if (ABeforeB(front_id, neighbor_id)) {
						if (InUnprocessedFront(neighbor_id)) {
							front.push(neighbor_id);
						}
					}
					else {
						if (InPriors(neighbor_id) && !newwave->InPriors(neighbor_id)) {
							newwave->mPriors[neighbor_id] = mPriors[neighbor_id];
						}
					}
				}
			}
			return newwave;
		}

		bool IsObstructed(INDEX_TYPE waitid) {
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			bool has_unassigned_lower = false;
			for (nit.begin(waitid); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (ABeforeB(neighbor_id, waitid)) {
					if (InUnprocessedFront(neighbor_id)) {
						return false;
					}
					//if (!InPriors(neighbor_id)) {
					//	has_unassigned_lower = true;
					//}
				}
			}
			return true;// has_unassigned_lower;
		}

		bool IsProcessed(INDEX_TYPE id) { return mProcessed->GetLabel(id) == 1; }
		bool IsRelativeMin(INDEX_TYPE id) {
			if (IsProcessed(id)) { printf("whoa min test procesed\n");  return false; }
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			bool has_unassigned_lower = false;
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (!IsProcessed(neighbor_id) && ABeforeB(neighbor_id, id)) {
					return false;
				}
			}
			return true;
		}

		void deepSplitWaveFront(vector < pair<INDEX_TYPE, ExplicitWaveFrontTry2*> > & obstructions) {
			set<INDEX_TYPE> marked;
			vector<INDEX_TYPE> res;
			FindObstructions(res);
			for (auto it = res.begin(); it != res.end(); it++) {
				INDEX_TYPE waitid = (*it);
				//if (IsObstructed(waitid) != IsRelativeMin(waitid)) printf("ERROR: relative mismatch %d %d \n", IsObstructed(waitid), IsRelativeMin(waitid));
				//if (IsObstructed(waitid)) {
				obstructions.push_back(pair<INDEX_TYPE, ExplicitWaveFrontTry2*>(waitid, createFromIncompleteSeed(waitid, marked)));
				//}

			}
			for (auto it = mWaiting.begin(); it != mWaiting.end(); it++) {

				if (marked.count(*it) == 0) printf("WHOA ERRORO ROROROROR not whole watiting is covered!!\n");
			}

		}


		void mergeAndConsumeWaveFront(ExplicitWaveFrontTry2* other) {
			mPriors.insert(other->mPriors.begin(), other->mPriors.end());
			mWaiting.insert(other->mWaiting.begin(), other->mWaiting.end());
			mReady2Go.insert(other->mReady2Go.begin(), other->mReady2Go.end());
			delete other;

			vector<INDEX_TYPE> newreadytogo;
			for (auto waiting_iterator = mWaiting.begin(); waiting_iterator != mWaiting.end(); waiting_iterator++) {
				INDEX_TYPE waiting_id = *waiting_iterator;

				if (CountBeforesNotInPrior(waiting_id) == 0) {
					mReady2Go.insert(waiting_id);
					newreadytogo.push_back(waiting_id);
				}
			}

			for (auto it = newreadytogo.begin(); it != newreadytogo.end(); it++) {
				mWaiting.erase(*it);
			}

		}

		int FrontSize() { return mWaiting.size() + mReady2Go.size(); }

		INDEX_TYPE FindLowestSeed() {

			vector<INDEX_TYPE> res;
			FindObstructions(res);
			bool haslowest = false;
			INDEX_TYPE lowest;
			for (auto it = res.begin(); it != res.end(); it++) {
				INDEX_TYPE waitid = (*it);

				//if (IsObstructed(waitid)) {
				if (haslowest) {
					if (ABeforeB(waitid, lowest)) {
						lowest = waitid;
					}
				}
				else {
					lowest = waitid;
					haslowest = true;
				}
				//}

			}
			return lowest;
		}
		int FindObstructions(vector<INDEX_TYPE>& res) {
			for (auto it = mWaiting.begin(); it != mWaiting.end(); it++) {
				INDEX_TYPE waitid = (*it);
				if (IsObstructed(waitid)) {
					res.push_back(waitid);
				}
			}
			res.insert(res.end(), mReady2Go.begin(), mReady2Go.end());
			return res.size();
		}


	};



	class TestingDependencyGraphParallelTraversalTry2Explicit {
	protected:
		DEPGRAPH_TYPE* mD;
		bool has_explicit_thred_num;
		int mNumThreads;
	public:
		TestingDependencyGraphParallelTraversalTry2Explicit(DEPGRAPH_TYPE* d) : mD(d), has_explicit_thred_num(false) {


		}
		void SetNumThreads(int i) { mNumThreads = i; }
		void OutputState(int round, vector<pair<INDEX_TYPE, ExplicitWaveFrontTry2*> >& waves) {
			char fname[2048];
			printf("doint round %d, parallelism = %d\n", round, waves.size());
			sprintf(fname, "test_s_%d.txt", round);
			printf("writing temp file %s\n", fname);
			FILE* ftest = fopen(fname, "w");
			for (auto hit = waves.begin(); hit != waves.end(); hit++) {

				ExplicitWaveFrontTry2* h = (*hit).second;


				for (auto it = h->mPriors.begin(); it != h->mPriors.end(); it++) {
					fprintf(ftest, "p %d %d\n", (*hit).first, (*it).first);
				}
				for (auto it = h->mWaiting.begin(); it != h->mWaiting.end(); it++) {
					if (h->IsObstructed((*it))) {
						fprintf(ftest, "o %d %d\n", (*hit).first, (*it));
					}
					else {
						fprintf(ftest, "w %d %d\n", (*hit).first, (*it));
					}

					DEPGRAPH_TYPE::neighbor_iterator nit(h->mGlobalGraph);
					for (nit.begin((*it)); nit.valid(); nit.advance()) {
						INDEX_TYPE neighbor_id = nit.value();
						if (h->ABeforeB((*it), neighbor_id)) {
							fprintf(ftest, "e %d %d\n", (*it), neighbor_id);
						}
						else {
							fprintf(ftest, "e %d %d\n", neighbor_id, (*it));
						}
					}


				}
				for (auto it = h->mReady2Go.begin(); it != h->mReady2Go.end(); it++) {
					fprintf(ftest, "r %d %d\n", (*it), *it);
					DEPGRAPH_TYPE::neighbor_iterator nit(h->mGlobalGraph);
					for (nit.begin((*it)); nit.valid(); nit.advance()) {
						INDEX_TYPE neighbor_id = nit.value();
						if (h->ABeforeB((*it), neighbor_id)) {
							fprintf(ftest, "e %d %d\n", (*it), neighbor_id);
						}
						else {
							fprintf(ftest, "e %d %d\n", neighbor_id, (*it));
						}
					}

				}

			}
			for (INDEX_TYPE i = 0; i < mD->IndexSpaceSize(); i++) {
				if (ExplicitWaveFrontTry2::mProcessed->GetLabel(i) == 1) {
					fprintf(ftest, "d %d %d\n", i, i);
				}
			}
			round++;
			fclose(ftest);
			printf("done\n");

		}

		void OutputState(int round, unordered_map<INDEX_TYPE, ExplicitWaveFrontTry2*>& waves) {
			char fname[2048];
			printf("doint round %d, parallelism = %d\n", round, waves.size());
			sprintf(fname, "test_s_%d.txt", round);
			printf("writing temp file %s\n", fname);
			FILE* ftest = fopen(fname, "w");
			for (auto hit = waves.begin(); hit != waves.end(); hit++) {

				ExplicitWaveFrontTry2* h = (*hit).second;


				for (auto it = h->mPriors.begin(); it != h->mPriors.end(); it++) {
					fprintf(ftest, "p %d %d\n", (*hit).first, (*it).first);
				}
				for (auto it = h->mWaiting.begin(); it != h->mWaiting.end(); it++) {
					if (h->IsObstructed((*it))) {
						fprintf(ftest, "o %d %d\n", (*hit).first, (*it));
					}
					else {
						fprintf(ftest, "w %d %d\n", (*hit).first, (*it));
					}

					DEPGRAPH_TYPE::neighbor_iterator nit(h->mGlobalGraph);
					for (nit.begin((*it)); nit.valid(); nit.advance()) {
						INDEX_TYPE neighbor_id = nit.value();
						if (h->ABeforeB((*it), neighbor_id)) {
							fprintf(ftest, "e %d %d\n", (*it), neighbor_id);
						}
						else {
							fprintf(ftest, "e %d %d\n", neighbor_id, (*it));
						}
					}


				}
				for (auto it = h->mReady2Go.begin(); it != h->mReady2Go.end(); it++) {
					fprintf(ftest, "r %d %d\n", (*it), *it);
					DEPGRAPH_TYPE::neighbor_iterator nit(h->mGlobalGraph);
					for (nit.begin((*it)); nit.valid(); nit.advance()) {
						INDEX_TYPE neighbor_id = nit.value();
						if (h->ABeforeB((*it), neighbor_id)) {
							fprintf(ftest, "e %d %d\n", (*it), neighbor_id);
						}
						else {
							fprintf(ftest, "e %d %d\n", neighbor_id, (*it));
						}
					}

				}

			}
			for (INDEX_TYPE i = 0; i < mD->IndexSpaceSize(); i++) {
				if (ExplicitWaveFrontTry2::mProcessed->GetLabel(i) == 1) {
					fprintf(ftest, "d %d %d\n", i, i);
				}
			}
			round++;
			fclose(ftest);
			printf("done\n");

		}

		void doStuff() {

			ExplicitWaveFrontTry2::mProcessed = new DenseLabeling<char>(mD->IndexSpaceSize());
			ExplicitWaveFrontTry2::mProcessed->SetAll(0);

			// create a new wavefront for every "minimum"		
			unordered_map<INDEX_TYPE, ExplicitWaveFrontTry2*> waves;
			unordered_map<INDEX_TYPE, ExplicitWaveFrontTry2*> global_waves;

			std::vector<INDEX_TYPE> topo_index_partition;
			INDEX_TYPE tTotalNumToDo = 0;
			int num_threads = has_explicit_thred_num;
#pragma omp parallel
			{
#pragma omp single
				{
					if (!has_explicit_thred_num) num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mD->IndexSpaceSize(), num_threads, topo_index_partition);
				}

				int thread_num = omp_get_thread_num();
				INDEX_TYPE tLocalNumToDo = 0;
				DEPGRAPH_TYPE::vertex_iterator vit(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (vit.begin(); vit.valid(); vit.advance()) {
					INDEX_TYPE id = vit.value();
					tLocalNumToDo++;
					bool ismin = true;
					DEPGRAPH_TYPE::neighbor_iterator nit(mD);
					for (nit.begin(id); nit.valid(); nit.advance()) {
						if (mD->Before(nit.value(), id)) {
							ismin = false;
							break;
						}
					}
					if (ismin) {
						//(*ip) = 0;
#pragma omp critical
						{
							int zero = 0;
							ExplicitWaveFrontTry2* w = new ExplicitWaveFrontTry2(mD);
							printf("created w%d\n", w->my_wave_id);
							w->SetSeed(id, zero);
							waves[id] = w;
						}
					}
				}
#pragma omp critical
				{
					printf("thread %d did %d ids\n", thread_num, tLocalNumToDo);
					tTotalNumToDo += tLocalNumToDo;
				}
			}
			printf("total number of edges: %d\n", tTotalNumToDo);


			set<INDEX_TYPE> global_seen;
			INDEX_TYPE tTotalProcCount = 0;
			int round = 0;
			bool bcont = true;
			while (bcont) {

				if (round > 30) return;

				OutputState(round++, waves);

				//// do first expansion
				int tLocalProcCount = 0;
				int randocount = 0;
				for (auto it = waves.begin(); it != waves.end(); it++) {
					ExplicitWaveFrontTry2* h = (*it).second;
					printf("-- expanding w%d\n", h->my_wave_id);
					bool cont = true;
					while (cont) {
						INDEX_TYPE id;
						cont = h->FindNextToProcess(id);
						if (!cont) break;
						tLocalProcCount++;

						if (waves.count(id) != 0 && waves[id] != h) {
							printf("about to process w%d obstruction %d in w%d\n", waves[id]->my_wave_id, id, h->my_wave_id);
							break;
						}

						h->SetAsProcessed(id, randocount);
					}
				}

				//OutputState(round++, waves);


				tTotalProcCount += tLocalProcCount;
				printf("did %d of %d...\n", tTotalProcCount, tTotalNumToDo);

				if (tTotalProcCount >= tTotalNumToDo) {
					printf("did %d of %d, returning\n", tTotalProcCount, tTotalNumToDo);
					bcont = false;
					return;
				}

				printf("waves.size() = %d\n", waves.size());
				vector<pair<INDEX_TYPE, ExplicitWaveFrontTry2*> > splits;
				for (auto it = waves.begin(); it != waves.end(); it++) {

					ExplicitWaveFrontTry2* h = (*it).second;
					printf("-- splitting w%d\n", h->my_wave_id);
					vector<INDEX_TYPE> obs;
					h->FindObstructions(obs);
					int ii = splits.size();
					if (obs.size() > 1) {
						if (waves.size() > 100) {
							splits.push_back(pair<INDEX_TYPE, ExplicitWaveFrontTry2*>(h->FindLowestSeed(), h));
						}
						else {
							h->deepSplitWaveFront(splits);
							for (; ii < splits.size(); ii++) {
								printf("-- split w%d -> w%d %d\n", h->my_wave_id, splits[ii].second->my_wave_id, splits[ii].first);
							}
							delete h;
						}
					}
					else if (obs.size() == 1) {
						printf("-- splitpassthrough w%d %d\n", h->my_wave_id, obs[0]);

						splits.push_back(pair<INDEX_TYPE, ExplicitWaveFrontTry2*>(obs[0], h));
					}
					else {
						printf("-- deleting w%d\n", h->my_wave_id);
						delete h;
					}
				}
				printf("done split -> %d\n", splits.size());
				//OutputState(round++, splits);
				//waves.clear();
				//if (splits.size() == 0) return;
				//printf("merging...\n");
				//for (int i = 1; i < splits.size(); i++) {
				//	splits[0].second->mergeAndConsumeWaveFront(splits[i].second);
				//}
				//waves.clear();
				//waves[splits[0].first] = splits[0].second;
				unordered_map<INDEX_TYPE, ExplicitWaveFrontTry2*> mergemap;
				for (auto it = splits.begin(); it != splits.end(); it++) {
					INDEX_TYPE id = (*it).first;
					ExplicitWaveFrontTry2* h = (*it).second;
					if (h->FrontSize() == 0) continue;

					if (mergemap.count(id) > 0) {
						printf("-- merge %d w%d <== w%d\n", id, mergemap[id]->my_wave_id, h->my_wave_id);

						mergemap[id]->mergeAndConsumeWaveFront(h);
					}
					else {
						printf("-- adding to merge %d w%d\n", id, h->my_wave_id);
						mergemap[id] = h;
					}				//if (h->FindLowestSeed() != id) printf("lowest = %d, seed = %d\n", h->FindLowestSeed(), id);
					//waves[id] = h;
				}
				waves = mergemap;



				printf("done  merge..\n");
				//round++;

				if (waves.size() == 0) return;


			}
		}
	};


	class TestingDependencyGraphParallelTraversalTry3Explicit {
	protected:
		DEPGRAPH_TYPE* mD;
		bool has_explicit_thred_num;
		int mNumThreads;
	public:
		TestingDependencyGraphParallelTraversalTry3Explicit(DEPGRAPH_TYPE* d) : mD(d), has_explicit_thred_num(false) {


		}
		void SetNumThreads(int i) { mNumThreads = i; has_explicit_thred_num = true; }
		void OutputState(int round, vector<pair<INDEX_TYPE, ExplicitWaveFrontTry2*> >& waves) {
			char fname[2048];
			printf("doint round %d, parallelism = %d\n", round, waves.size());
			sprintf(fname, "test_s_%d.txt", round);
			printf("writing temp file %s\n", fname);
			FILE* ftest = fopen(fname, "w");
			for (auto hit = waves.begin(); hit != waves.end(); hit++) {

				ExplicitWaveFrontTry2* h = (*hit).second;


				for (auto it = h->mPriors.begin(); it != h->mPriors.end(); it++) {
					fprintf(ftest, "p %d %d\n", (*hit).first, (*it).first);
				}
				for (auto it = h->mWaiting.begin(); it != h->mWaiting.end(); it++) {
					if (h->IsObstructed((*it))) {
						fprintf(ftest, "o %d %d\n", (*hit).first, (*it));
					}
					else {
						fprintf(ftest, "w %d %d\n", (*hit).first, (*it));
					}

					DEPGRAPH_TYPE::neighbor_iterator nit(h->mGlobalGraph);
					for (nit.begin((*it)); nit.valid(); nit.advance()) {
						INDEX_TYPE neighbor_id = nit.value();
						if (h->ABeforeB((*it), neighbor_id)) {
							fprintf(ftest, "e %d %d\n", (*it), neighbor_id);
						}
						else {
							fprintf(ftest, "e %d %d\n", neighbor_id, (*it));
						}
					}


				}
				for (auto it = h->mReady2Go.begin(); it != h->mReady2Go.end(); it++) {
					fprintf(ftest, "r %d %d\n", (*it), *it);
					DEPGRAPH_TYPE::neighbor_iterator nit(h->mGlobalGraph);
					for (nit.begin((*it)); nit.valid(); nit.advance()) {
						INDEX_TYPE neighbor_id = nit.value();
						if (h->ABeforeB((*it), neighbor_id)) {
							fprintf(ftest, "e %d %d\n", (*it), neighbor_id);
						}
						else {
							fprintf(ftest, "e %d %d\n", neighbor_id, (*it));
						}
					}

				}

			}
			for (INDEX_TYPE i = 0; i < mD->IndexSpaceSize(); i++) {
				if (ExplicitWaveFrontTry2::mProcessed->GetLabel(i) == 1) {
					fprintf(ftest, "d %d %d\n", i, i);
				}
			}
			round++;
			fclose(ftest);
			printf("done\n");

		}

		void OutputState(int round, unordered_map<INDEX_TYPE, ExplicitWaveFrontTry2*>& waves) {
			char fname[2048];
			printf("doint round %d, parallelism = %d\n", round, waves.size());
			sprintf(fname, "test_s_%d.txt", round);
			printf("writing temp file %s\n", fname);
			FILE* ftest = fopen(fname, "w");
			for (auto hit = waves.begin(); hit != waves.end(); hit++) {

				ExplicitWaveFrontTry2* h = (*hit).second;


				for (auto it = h->mPriors.begin(); it != h->mPriors.end(); it++) {
					fprintf(ftest, "p %d %d\n", (*hit).first, (*it).first);
				}
				for (auto it = h->mWaiting.begin(); it != h->mWaiting.end(); it++) {
					if (h->IsObstructed((*it))) {
						fprintf(ftest, "o %d %d\n", (*hit).first, (*it));
					}
					else {
						fprintf(ftest, "w %d %d\n", (*hit).first, (*it));
					}

					DEPGRAPH_TYPE::neighbor_iterator nit(h->mGlobalGraph);
					for (nit.begin((*it)); nit.valid(); nit.advance()) {
						INDEX_TYPE neighbor_id = nit.value();
						if (h->ABeforeB((*it), neighbor_id)) {
							fprintf(ftest, "e %d %d\n", (*it), neighbor_id);
						}
						else {
							fprintf(ftest, "e %d %d\n", neighbor_id, (*it));
						}
					}


				}
				for (auto it = h->mReady2Go.begin(); it != h->mReady2Go.end(); it++) {
					fprintf(ftest, "r %d %d\n", (*it), *it);
					DEPGRAPH_TYPE::neighbor_iterator nit(h->mGlobalGraph);
					for (nit.begin((*it)); nit.valid(); nit.advance()) {
						INDEX_TYPE neighbor_id = nit.value();
						if (h->ABeforeB((*it), neighbor_id)) {
							fprintf(ftest, "e %d %d\n", (*it), neighbor_id);
						}
						else {
							fprintf(ftest, "e %d %d\n", neighbor_id, (*it));
						}
					}

				}

			}
			for (INDEX_TYPE i = 0; i < mD->IndexSpaceSize(); i++) {
				if (ExplicitWaveFrontTry2::mProcessed->GetLabel(i) == 1) {
					fprintf(ftest, "d %d %d\n", i, i);
				}
			}
			round++;
			fclose(ftest);
			printf("done\n");

		}

		enum WorkType { EXPAND_ONLY, MERGE_THEN_EXPAND, TREAT, MERGE_THEN_SPLIT };
		struct WorkRequest {
			WorkType work_type;
			INDEX_TYPE a;
			INDEX_TYPE b;
			ExplicitWaveFrontTry2* h1;
			ExplicitWaveFrontTry2* h2;
		};

		stack<WorkRequest> mWorkToDo;
		stack<WorkRequest> mTempWorkToDo;

		void doStuff() {
			if (has_explicit_thred_num) omp_set_num_threads(mNumThreads);
			ExplicitWaveFrontTry2::mProcessed = new DenseLabeling<char>(mD->IndexSpaceSize());
			ExplicitWaveFrontTry2::mProcessed->SetAll(0);

			// create a new wavefront for every "minimum"		
			unordered_map<INDEX_TYPE, ExplicitWaveFrontTry2*> waves;
			unordered_map<INDEX_TYPE, ExplicitWaveFrontTry2*> global_waves;
			std::chrono::steady_clock::time_point g_start_time = std::chrono::steady_clock::now();;
			std::chrono::steady_clock::time_point start_split;

			std::vector<INDEX_TYPE> topo_index_partition;
			INDEX_TYPE tTotalNumToDo = 0;
			int num_threads = mNumThreads;
#pragma omp parallel
			{
#pragma omp single
				{
					if (!has_explicit_thred_num) num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mD->IndexSpaceSize(), num_threads, topo_index_partition);
				}

				int thread_num = omp_get_thread_num();
				INDEX_TYPE tLocalNumToDo = 0;
				DEPGRAPH_TYPE::vertex_iterator vit(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (vit.begin(); vit.valid(); vit.advance()) {
					INDEX_TYPE id = vit.value();
					tLocalNumToDo++;
					bool ismin = true;
					DEPGRAPH_TYPE::neighbor_iterator nit(mD);
					for (nit.begin(id); nit.valid(); nit.advance()) {
						if (mD->Before(nit.value(), id)) {
							ismin = false;
							break;
						}
					}
					if (ismin) {
						//(*ip) = 0;
#pragma omp critical
						{
							int zero = 0;
							ExplicitWaveFrontTry2* w = new ExplicitWaveFrontTry2(mD);
							printf("created w%d\n", w->my_wave_id);
							w->SetSeed(id, zero);
							mWorkToDo.push(WorkRequest{ EXPAND_ONLY, id, 0, w, NULL });
							//waves[id] = w;
						}
					}
				}
#pragma omp critical
				{
					printf("thread %d did %d ids\n", thread_num, tLocalNumToDo);
					tTotalNumToDo += tLocalNumToDo;
				}
			}
			printf("total number of edges: %d\n", tTotalNumToDo);

			stack<WorkRequest> spillover;
#pragma omp parallel
			{
				bool first_time = true;
				int thread_num = omp_get_thread_num();
				bool completely_done = false;
				while (!completely_done) {
					// if this is not the first time doing work, then we start with a massive split
					if (!first_time) {


					}
					first_time = true;


					bool keep_looping = true;
					WorkRequest work;
					while (keep_looping) { // begin while
#pragma omp critical
					{
						if (mWorkToDo.empty()) {
							keep_looping = false;
						}
						else {
							work = mWorkToDo.top();
							mWorkToDo.pop();
						}

					}
						if (!keep_looping) continue;

						if (work.work_type == EXPAND_ONLY) {
							ExplicitWaveFrontTry2* h = work.h1;
#pragma omp critical
							{
								printf("thread %d: expanding w%d\n", thread_num, h->my_wave_id);
							}
							bool cont_expand = true;
							while (cont_expand) {
								INDEX_TYPE id;
								cont_expand = h->FindNextToProcess(id);
								if (!cont_expand) break;
								//tLocalProcCount++;

								//if (waves.count(id) != 0 && waves[id] != h) {
								//	printf("about to process w%d obstruction %d in w%d\n", waves[id]->my_wave_id, id, h->my_wave_id);
								//	break;
								//}
								
								int poop;
								h->SetAsProcessed(id, poop);
							}
							if (h->FrontSize() == 0) {
								// done
#pragma omp critical
							{
								printf("thread %d: deleting w%d\n", thread_num, h->my_wave_id);
							}
								delete h;
							}
							else {
								INDEX_TYPE obsid = h->FindLowestSeed();
#pragma omp critical
								{
									if (waves.count(obsid) > 0) {
										ExplicitWaveFrontTry2* otherh = waves[obsid];
										waves.erase(obsid);
										mWorkToDo.push(WorkRequest{ MERGE_THEN_EXPAND, obsid, obsid, h, otherh });

									}
									else {
										waves[obsid] = h;
									}
								}



							}
						}
						else if (work.work_type == MERGE_THEN_EXPAND) {
#pragma omp critical
						{
							printf("thread %d: merging w%d <== w%d\n", thread_num, work.h1->my_wave_id, work.h2->my_wave_id);
						}
							work.h1->mergeAndConsumeWaveFront(work.h2);
#pragma omp critical
							{
								if (mWorkToDo.size() < num_threads ) {
									keep_looping = false;
									spillover.push(WorkRequest{ EXPAND_ONLY, 0, 0, work.h1, NULL });
								}
								else {
									mWorkToDo.push(WorkRequest{ EXPAND_ONLY, 0, 0, work.h1, NULL });
								}
							}
						}




					} // end while loop for this round
#pragma omp barrier			
					/// COUNT HOW MUCH THERE IS LEFT TO DO
	#pragma omp single
				{
					start_split = std::chrono::steady_clock::now();
					ArrayIndexPartitioner::EvenChunkSplit(mD->IndexSpaceSize(), num_threads, topo_index_partition);
					tTotalNumToDo = 0;
				}

				INDEX_TYPE tLocalNumToDo = 0;
				DEPGRAPH_TYPE::vertex_iterator vit(mD, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (vit.begin(); vit.valid(); vit.advance()) {
					INDEX_TYPE id = vit.value();
					if (ExplicitWaveFrontTry2::mProcessed->GetLabel(id) == 0) tLocalNumToDo++;
				}
#pragma omp critical
				{
					tTotalNumToDo += tLocalNumToDo;
				}
#pragma omp barrier
#pragma omp single
				{
					printf("after round, there are %d left to do\n", tTotalNumToDo);

					while (!spillover.empty()) {
						mWorkToDo.push(spillover.top());
						spillover.pop();
					}
					for (auto it = waves.begin(); it != waves.end(); it++) {
						mWorkToDo.push(WorkRequest{ TREAT, (*it).first, 0, (*it).second, NULL });
					}
					waves.clear();
				}
#pragma omp barrier
				if (tTotalNumToDo == 0)	completely_done = true;			
				
				if (!completely_done) {
					// fix up some work!
					//vector<WorkRequest> linearized;
					vector<pair<INDEX_TYPE, ExplicitWaveFrontTry2*> > splits;

					bool cont = true;
					while (cont) {
						WorkRequest work;
#pragma omp critical
						{
							if (mWorkToDo.empty()) {
								cont = false;
							}
							else {
								work = mWorkToDo.top();
								mWorkToDo.pop();
							}
						}
						// now do the split
						if (work.work_type == EXPAND_ONLY) {
							work.h1->deepSplitWaveFront(splits);
							delete work.h1;
						}
						else if (work.work_type == MERGE_THEN_EXPAND) {
							work.h1->mergeAndConsumeWaveFront(work.h2);
							work.h1->deepSplitWaveFront(splits);
						}
						else if (work.work_type == TREAT){
							work.h1->deepSplitWaveFront(splits);
						} // end if worktype

						// now all threads have local list of split wavefronts.
						for (auto it = splits.begin(); it != splits.end(); it++) {
							INDEX_TYPE id = (*it).first;
							INDEX_TYPE placeholder;
							ExplicitWaveFrontTry2* h = (*it).second;
							if (h->FindNextToProcess(placeholder)) {
								// has one to process, add to work request
#pragma omp critical
								{
									mTempWorkToDo.push(WorkRequest{ EXPAND_ONLY, placeholder, 0, h, NULL });
								}

							}
							else {

#pragma omp critical
								{
									if (waves.count(id) != 0) {
										ExplicitWaveFrontTry2* otherh = waves[id];
										waves.erase(id);
										mTempWorkToDo.push(WorkRequest{ MERGE_THEN_EXPAND, id, id, h, otherh });

									}
									else {
										waves[id] = h;
									}
								}
							}							
						}
						splits.clear();



					} // end while loop to split
					



				} // end work for splitting
#pragma omp barrier
#pragma omp single
				{
					mWorkToDo = mTempWorkToDo;
					std::chrono::steady_clock::time_point end_split = std::chrono::steady_clock::now();
					printf("split took %d ms\n", std::chrono::duration_cast<std::chrono::milliseconds>(end_split - start_split).count());
				}

				} // end global while loop




//#pragma omp critical
//						{
//							std::chrono::steady_clock::time_point my_end_time = std::chrono::steady_clock::now();
//
//							printf("thread %d: done %d ms\n", thread_num, std::chrono::duration_cast<std::chrono::milliseconds>(my_end_time - g_start_time).count());
//						}
			} // end LOOPING parallel section





















			//set<INDEX_TYPE> global_seen;
			//INDEX_TYPE tTotalProcCount = 0;
			//int round = 0;
			//bool bcont = true;
			//while (bcont) {

			//	if (round % 5 == 0) OutputState(round/5, waves);
			//	round++;

			//	

			//	//// do first expansion
			//	int tLocalProcCount = 0;
			//	int randocount = 0;
			//	for (auto it = waves.begin(); it != waves.end(); it++) {
			//		ExplicitWaveFrontTry2* h = (*it).second;
			//		printf("-- expanding w%d\n", h->my_wave_id);
			//		bool cont = true;
			//		while (cont) {
			//			INDEX_TYPE id;
			//			cont = h->FindNextToProcess(id);
			//			if (!cont) break;
			//			tLocalProcCount++;

			//			if (waves.count(id) != 0 && waves[id] != h) {
			//				printf("about to process w%d obstruction %d in w%d\n", waves[id]->my_wave_id, id, h->my_wave_id);
			//				break;
			//			}

			//			h->SetAsProcessed(id, randocount);
			//		}
			//	}

			//	//OutputState(round++, waves);


			//	tTotalProcCount += tLocalProcCount;
			//	printf("did %d of %d...\n", tTotalProcCount, tTotalNumToDo);

			//	if (tTotalProcCount >= tTotalNumToDo) {
			//		printf("did %d of %d, returning\n", tTotalProcCount, tTotalNumToDo);
			//		bcont = false;
			//		return;
			//	}

			//	printf("waves.size() = %d\n", waves.size());
			//	unordered_map<INDEX_TYPE, ExplicitWaveFrontTry2*> mergemap;
			//	for (auto it = waves.begin(); it != waves.end(); it++) {
			//		ExplicitWaveFrontTry2* h = (*it).second;
			//		//printf("-- splitting w%d\n", h->my_wave_id);
			//		vector<INDEX_TYPE> obs;
			//		h->FindObstructions(obs);
			//		bool consumed = false;
			//		for (int i = 0; i < obs.size(); i++) {
			//			if (mergemap.count(obs[i]) != 0) {
			//				printf("merged at %d w%d <== w%d\n", obs[i], mergemap[obs[i]]->my_wave_id, h->my_wave_id);
			//				mergemap[obs[i]]->mergeAndConsumeWaveFront(h);
			//				consumed = true;
			//				i = obs.size();
			//			}
			//		}
			//		if (!consumed) {
			//			for (int i = 0; i < obs.size(); i++) {
			//				mergemap[obs[i]] = h;
			//			}
			//		}
			//	}
			//	waves.clear();
			//	set<ExplicitWaveFrontTry2*> uniqs;
			//	for (auto it = mergemap.begin(); it != mergemap.end(); it++) {
			//		ExplicitWaveFrontTry2* h = (*it).second;
			//		if (uniqs.count(h) != 0) continue;
			//		uniqs.insert(h);
			//		//printf("-- splitting w%d\n", h->my_wave_id);
			//		INDEX_TYPE ob = h->FindLowestSeed();
			//		if (waves.count(ob) != 0) {
			//			printf("SHOSLSKDKDKDKDK %d w%d w%d\n", ob, h->my_wave_id, waves[ob]->my_wave_id);
			//			waves[ob]->mergeAndConsumeWaveFront(h);
			//		}
			//		else {
			//			waves[ob] = h;
			//		}
			//	}
			//	//waves = mergemap;

			//	//vector<pair<INDEX_TYPE, ExplicitWaveFrontTry2*> > splits;
			//	//for (auto it = waves.begin(); it != waves.end(); it++) {

			//	//	ExplicitWaveFrontTry2* h = (*it).second;
			//	//	printf("-- splitting w%d\n", h->my_wave_id);
			//	//	vector<INDEX_TYPE> obs;
			//	//	h->FindObstructions(obs);
			//	//	int ii = splits.size();
			//	//	if (obs.size() > 1) {
			//	//		if (waves.size() > 100) {
			//	//			splits.push_back(pair<INDEX_TYPE, ExplicitWaveFrontTry2*>(h->FindLowestSeed(), h));
			//	//		}
			//	//		else {
			//	//			h->deepSplitWaveFront(splits);
			//	//			for (; ii < splits.size(); ii++) {
			//	//				printf("-- split w%d -> w%d %d\n", h->my_wave_id, splits[ii].second->my_wave_id, splits[ii].first);
			//	//			}
			//	//			delete h;
			//	//		}
			//	//	}
			//	//	else if (obs.size() == 1) {
			//	//		printf("-- splitpassthrough w%d %d\n", h->my_wave_id, obs[0]);

			//	//		splits.push_back(pair<INDEX_TYPE, ExplicitWaveFrontTry2*>(obs[0], h));
			//	//	}
			//	//	else {
			//	//		printf("-- deleting w%d\n", h->my_wave_id);
			//	//		delete h;
			//	//	}
			//	//}
			//	//printf("done split -> %d\n", splits.size());
			//	////OutputState(round++, splits);
			//	////waves.clear();
			//	////if (splits.size() == 0) return;
			//	////printf("merging...\n");
			//	////for (int i = 1; i < splits.size(); i++) {
			//	////	splits[0].second->mergeAndConsumeWaveFront(splits[i].second);
			//	////}
			//	////waves.clear();
			//	////waves[splits[0].first] = splits[0].second;
			//	//unordered_map<INDEX_TYPE, ExplicitWaveFrontTry2*> mergemap;
			//	//for (auto it = splits.begin(); it != splits.end(); it++) {
			//	//	INDEX_TYPE id = (*it).first;
			//	//	ExplicitWaveFrontTry2* h = (*it).second;
			//	//	if (h->FrontSize() == 0) continue;

			//	//	if (mergemap.count(id) > 0) {
			//	//		printf("-- merge %d w%d <== w%d\n", id, mergemap[id]->my_wave_id, h->my_wave_id);

			//	//		mergemap[id]->mergeAndConsumeWaveFront(h);
			//	//	}
			//	//	else {
			//	//		printf("-- adding to merge %d w%d\n", id, h->my_wave_id);
			//	//		mergemap[id] = h;
			//	//	}				//if (h->FindLowestSeed() != id) printf("lowest = %d, seed = %d\n", h->FindLowestSeed(), id);
			//	//	//waves[id] = h;
			//	//}
			//	//waves = mergemap;



			//	printf("done  merge..\n");
			//	//round++;

			//	if (waves.size() == 0) return;


			//}
		}
	};


}
#endif
