/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef GI_EXPERIMENTAL2_H
#define GI_EXPERIMENTAL2_H

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


	class HashingWaveFrontTry2 {
	public:

		HashingWaveFrontTry2(DEPGRAPH_TYPE* global_graph) : mGlobalGraph(global_graph) { my_wave_id = count++; }

		struct PriorElement {
			BYTE_TYPE count_unprocessed_afters;
			PAYLOAD_TYPE payload;

			PriorElement() {
				count_unprocessed_afters = 0;
			}
			PriorElement( PAYLOAD_TYPE p) : payload(p) {
				count_unprocessed_afters = 0;
			}
		};


		// MEMBER VARIABLES

		DEPGRAPH_TYPE* mGlobalGraph;
		static int count;
		int my_wave_id;
		map<INDEX_TYPE, PriorElement > mPriors;
		map<INDEX_TYPE, BYTE_TYPE> mWaiting;
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
		BYTE_TYPE AddEnqueuedIdToPriorCount(INDEX_TYPE my_id, INDEX_TYPE fromid) {
			BYTE_TYPE countlower = 0;
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(my_id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (neighbor_id != fromid && ABeforeB(neighbor_id, my_id)) {
					if (InPriors(neighbor_id)) {
						mPriors[neighbor_id].count_unprocessed_afters++;
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
				AddEnqueuedIdToPriorCount(my_id, fromid);
				mWaiting[my_id] = countlower;
			}
		}

		void UpdateNeighborWaitingCountOrInsert(INDEX_TYPE id, INDEX_TYPE fromid) {
			if (InWaiting(id)) {
				if (mWaiting[id] == 1) {
					mWaiting.erase(id);
					mReady2Go.insert(id);
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

		void UpdatePriorCountOrDelete(INDEX_TYPE id) {
			PriorElement& e = mPriors[id];
			if (e.count_unprocessed_afters == 1) {
				mPriors.erase(id);
				return;
			}
			e.count_unprocessed_afters--;
		}

		void SetAsProcessed(INDEX_TYPE id, PAYLOAD_TYPE& p) {
			RemoveFromReady(id);
			PriorElement ne(p);
			mPriors[id] = ne;
			mProcessed->SetLabel(id,1);
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (ABeforeB(id, neighbor_id)) {
					UpdateNeighborWaitingCountOrInsert(neighbor_id, id);
					mPriors[id].count_unprocessed_afters++;
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
			return true;
		}

		HashingWaveFrontTry2* createFromIncompleteSeed(INDEX_TYPE base_id) {
			HashingWaveFrontTry2* newwave = new HashingWaveFrontTry2(mGlobalGraph);
			stack<INDEX_TYPE> front;
			front.push(base_id);
			while (!front.empty()) {
				INDEX_TYPE front_id = front.top();
				front.pop();
				if (newwave->InWaiting(front_id)) continue;
				newwave->mWaiting[front_id] = 0;
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				for (nit.begin(front_id); nit.valid(); nit.advance()) {
					INDEX_TYPE neighbor_id = nit.value();
					if (ABeforeB(front_id, neighbor_id)) {
						if (InWaiting(neighbor_id)) {
							front.push(neighbor_id);
						}
					}
					else {
						if (!InPriors(neighbor_id)) {
							newwave->mWaiting[front_id]++;
						}
						else {
							if (!newwave->InPriors(neighbor_id)) {
								newwave->mPriors[neighbor_id] = mPriors[neighbor_id];
								newwave->mPriors[neighbor_id].count_unprocessed_afters = 1;
							}
							else {
								newwave->mPriors[neighbor_id].count_unprocessed_afters++; // COME BACK TO THIS MAYBE
							}
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
					if (InWaiting(neighbor_id) || InReady(neighbor_id)) {
						return false;// has_internal = true;
						//break;
					}
					if (!InPriors(neighbor_id)) {
						has_unassigned_lower = true;
					}
				}
			}
			return has_unassigned_lower;
		}

		bool IsProcessed(INDEX_TYPE id) { return mProcessed->GetLabel(id) == 1; }
		bool IsRelativeMin(INDEX_TYPE id) {
			if (IsProcessed(id)) { return false; }
			DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
			bool has_unassigned_lower = false;
			for (nit.begin(id); nit.valid(); nit.advance()) {
				INDEX_TYPE neighbor_id = nit.value();
				if (! IsProcessed(neighbor_id) && ABeforeB(neighbor_id, id)) {
					return false;
				}
			}
			return true;
		}

		void deepSplitWaveFront(vector < pair<INDEX_TYPE, HashingWaveFrontTry2*> > & obstructions) {
			for (auto it = mWaiting.begin(); it != mWaiting.end(); it++) {
				INDEX_TYPE waitid = (*it).first;
				if (IsObstructed(waitid) != IsRelativeMin(waitid)) printf("ERROR: relative mismatch %d %d \n", IsObstructed(waitid) , IsRelativeMin(waitid));
				if (IsObstructed(waitid)) {
					obstructions.push_back(pair<INDEX_TYPE, HashingWaveFrontTry2*>(waitid, createFromIncompleteSeed(waitid)));
				}

			}
		}


		void mergeAndConsumeWaveFront(HashingWaveFrontTry2* other) {
			mPriors.insert(other->mPriors.begin(), other->mPriors.end());
			mWaiting.insert(other->mWaiting.begin(), other->mWaiting.end());
			mReady2Go.insert(other->mReady2Go.begin(), other->mReady2Go.end());
			delete other;
			for (auto it = mPriors.begin(); it != mPriors.end(); it++) {
				(*it).second.count_unprocessed_afters = 0;
			}
			vector<INDEX_TYPE> newreadytogo;
			for (auto waiting_iterator = mWaiting.begin(); waiting_iterator != mWaiting.end(); waiting_iterator++) {
				INDEX_TYPE waiting_id = waiting_iterator->first;
				waiting_iterator->second = 0; // how many this item is waiting on 
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				for (nit.begin(waiting_id); nit.valid(); nit.advance()) {
					INDEX_TYPE nid = nit.value();
					// for everything that is not a before with prior payload computed, count it up
					if (ABeforeB(nid, waiting_id)) {
						if (!InPriors(nid)) {
							waiting_iterator->second++;
						}
					}
				}
				if (waiting_iterator->second == 0) {
					mReady2Go.insert(waiting_id);
					newreadytogo.push_back(waiting_id);
				}
			}

			/// count number of unprocessed later neighbors
			for (auto prior_it = mPriors.begin(); prior_it != mPriors.end(); prior_it++) {
				DEPGRAPH_TYPE::neighbor_iterator nit(mGlobalGraph);
				for (nit.begin(prior_it->first); nit.valid(); nit.advance()) {
					INDEX_TYPE neighbor_id = nit.value();
					if (ABeforeB(prior_it->first, neighbor_id)) {
						//if (InWaiting(neighbor_id)){
							prior_it->second.count_unprocessed_afters++; // should this be NOT PROCESSED
						//}
					}
				}
			}

			for (auto it = newreadytogo.begin(); it != newreadytogo.end(); it++) {
				mWaiting.erase(*it);
			}

		}

		int FrontSize() { return mWaiting.size() + mReady2Go.size(); }

		INDEX_TYPE FindLowestSeed() {
			bool haslowest = false;
			INDEX_TYPE lowest;
			for (auto it = mWaiting.begin(); it != mWaiting.end(); it++) {
				INDEX_TYPE waitid = (*it).first;

				if (IsObstructed(waitid)) {
					if (haslowest) {
						if (ABeforeB(waitid, lowest)) {
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
		int FindObstructions(vector<INDEX_TYPE>& res) {
			for (auto it = mWaiting.begin(); it != mWaiting.end(); it++) {
				INDEX_TYPE waitid = (*it).first;
				if (IsObstructed(waitid)) {
					res.push_back(waitid);
				}
			}
			return res.size();
		}
	};




	class TestingDependencyGraphParallelTraversalTry2 {
	protected:
		DEPGRAPH_TYPE* mD;
		bool has_explicit_thred_num;
		int mNumThreads;
	public:
		TestingDependencyGraphParallelTraversalTry2(DEPGRAPH_TYPE* d) : mD(d), has_explicit_thred_num(false) {


		}
		void SetNumThreads(int i) { mNumThreads = i; }

		void OutputState(int round, unordered_map<INDEX_TYPE, HashingWaveFrontTry2*>& waves) {
			char fname[2048];
			printf("doint round %d, parallelism = %d\n", round, waves.size());
			sprintf(fname, "test_s_%d.txt", round);
			printf("writing temp file %s\n", fname);
			FILE* ftest = fopen(fname, "w");
			for (auto hit = waves.begin(); hit != waves.end(); hit++) {

				HashingWaveFrontTry2* h = (*hit).second;

				
				for (auto it = h->mPriors.begin(); it != h->mPriors.end(); it++) {
					fprintf(ftest, "p %d %d\n", (*hit).first, (*it).first);
				}
				for (auto it = h->mWaiting.begin(); it != h->mWaiting.end(); it++) {
					if (h->IsObstructed((*it).first)) {
						fprintf(ftest, "o %d %d\n", (*hit).first, (*it).first);
					}
					else {
						fprintf(ftest, "w %d %d\n", (*hit).first, (*it).first);
					}

					DEPGRAPH_TYPE::neighbor_iterator nit(h->mGlobalGraph);
					for (nit.begin((*it).first); nit.valid(); nit.advance()) {
						INDEX_TYPE neighbor_id = nit.value();
						if (h->ABeforeB((*it).first, neighbor_id)) {
							fprintf(ftest, "e %d %d\n", (*it).first, neighbor_id);
						}
						else {
							fprintf(ftest, "e %d %d\n", neighbor_id, (*it).first);
						}
					}


				}
				for (auto it = h->mReady2Go.begin(); it != h->mReady2Go.end(); it++) {
					fprintf(ftest, "r %d %d\n", (*it), *it);
				}

			}
			for (INDEX_TYPE i = 0; i < mD->IndexSpaceSize(); i++) {
				if (HashingWaveFrontTry2::mProcessed->GetLabel(i) == 1) {
					fprintf(ftest, "d %d %d\n", i, i);
				}
			}
			round++;
			fclose(ftest);
			printf("done\n");

		}

		void doStuff() {

			HashingWaveFrontTry2::mProcessed = new DenseLabeling<char>(mD->IndexSpaceSize());
			HashingWaveFrontTry2::mProcessed->SetAll(0);

			// create a new wavefront for every "minimum"		
			unordered_map<INDEX_TYPE, HashingWaveFrontTry2*> waves;
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
							HashingWaveFrontTry2* w = new HashingWaveFrontTry2(mD);
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
					HashingWaveFrontTry2* h = (*it).second;
					bool cont = true;
					while (cont) {
						INDEX_TYPE id;
						cont = h->FindNextToProcess(id);
						if (!cont) break;
						tLocalProcCount++;
						h->SetAsProcessed(id, randocount);
					}
				}

				OutputState(round++, waves);


				tTotalProcCount += tLocalProcCount;
				printf("did %d of %d...\n", tTotalProcCount, tTotalNumToDo);

				if (tTotalProcCount >= tTotalNumToDo) {
					printf("did %d of %d, returning\n", tTotalProcCount, tTotalNumToDo);
					bcont = false;
					return;
				}

				vector<pair<INDEX_TYPE, HashingWaveFrontTry2*> > splits;
				for (auto it = waves.begin(); it != waves.end(); it++) {
					
					HashingWaveFrontTry2* h = (*it).second;
					vector<INDEX_TYPE> obs;
					h->FindObstructions(obs);
					if (obs.size() > 1) {
						h->deepSplitWaveFront(splits);
						delete h;
					}
					else if (obs.size() == 1) {
						splits.push_back(pair<INDEX_TYPE, HashingWaveFrontTry2*>(obs[0], h));
					}
					else {
						delete h;
					}
				}
				printf("done split -> %d\n", splits.size());
				waves.clear();
				unordered_map<INDEX_TYPE, HashingWaveFrontTry2*> mergemap;
				for (auto it = splits.begin(); it != splits.end(); it++) {
					INDEX_TYPE id = (*it).first;
					HashingWaveFrontTry2* h = (*it).second;
					if (h->FrontSize() == 0) continue;

					if (mergemap.count(id) > 0) {
						mergemap[id]->mergeAndConsumeWaveFront(h);
					}
					else {
						mergemap[id] = h;
					}				//if (h->FindLowestSeed() != id) printf("lowest = %d, seed = %d\n", h->FindLowestSeed(), id);
					//waves[id] = h;
				}
				waves = mergemap;
				OutputState(round++, waves);


				////vector < pair<INDEX_TYPE, HashingWaveFrontTry2*> > ohG;
				//for (auto it = waves.begin(); it != waves.end(); it++) {
				//	HashingWaveFrontTry2* h = (*it).second;

				//	vector<INDEX_TYPE> obs;
				//	h->FindObstructions(obs);
				//	
				//	//bool hasmerge = false;
				//	if (obs.size() == 0) continue;
				//	INDEX_TYPE lowest = obs[0];
				//	for (int i = 1; i < obs.size(); i++) {
				//		INDEX_TYPE oid = obs[i];
				//		if (h->ABeforeB(oid, lowest)) lowest = oid;
				//	}
				//	if (mergemap.count(lowest) > 0) {
				//		mergemap[lowest]->mergeAndConsumeWaveFront(h);
				//	}
				//	else {
				//		mergemap[lowest] = h;
				//	}

				//	//if (mergemap.count(oid) > 0) {
				//	//	for (int i = 0; i < obs.size(); i++) {
				//	//	INDEX_TYPE oid = obs[i];
				//	//	if (mergemap.count(oid) > 0) {
				//	//		printf("merging %d: %d <== %d\n", oid, mergemap[oid]->my_wave_id, h->my_wave_id);
				//	//		mergemap[oid]->mergeAndConsumeWaveFront(h);
				//	//		break;
				//	//	}
				//	//	printf("adding %d -> %d\n", oid, h->my_wave_id);
				//	//	mergemap[oid] = h;
				//	//}
				//}
				//waves = mergemap;

				//OutputState(round++, waves);


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

				//waves.clear();
				//waves = mergemap;
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


//////	class TestingDependencyGraphTraversal {
//////	protected:
//////		DEPGRAPH_TYPE* mD;
//////	public:
//////		TestingDependencyGraphTraversal(DEPGRAPH_TYPE* d) : mD(d) {
//////
//////
//////		}
//////
//////		void doStuff() {
//////
//////
//////
//////			HashingWaveFront* wave = new HashingWaveFront(mD);
//////
//////			int targetnum = 100; //mD->NumElements() / 1000;
//////#ifdef USECUSTOMALLOCATOR
//////
//////			FSBAllocator<int> int_alloc;
//////#endif			
//////
//////			DEPGRAPH_TYPE::vertex_iterator vit(mD);
//////			for (vit.begin(); vit.valid(); vit.advance()) {
//////				INDEX_TYPE id = vit.value();
//////				bool ismin = true;
//////				DEPGRAPH_TYPE::neighbor_iterator nit(mD);
//////				for (nit.begin(id); nit.valid(); nit.advance()) {
//////					if (mD->Before(nit.value(), id)) {
//////						ismin = false;
//////						break;
//////					}
//////				}
//////				if (ismin) {
//////#ifdef USECUSTOMALLOCATOR
//////					int* ip = int_alloc.allocate(1); // new int;
//////#else
//////					//int* ip = new int;
//////#endif			
//////
//////					//(*ip) = 0;
//////					int zero = 0;
//////					wave->SetSeed(id, zero);
//////				}
//////			}
//////			FILE* fout = fopen("outtest.txt", "w");
//////			fprintf(fout, "ASDFASDF %d %d %d \n", wave->mPriors.size(), wave->mWaiting.size(), wave->mReady2Go.size());
//////
//////			int count = 0;
//////			while (true) {
//////				INDEX_TYPE id;
//////				bool  hasnext = wave->GetNextToProcess(id);
//////				if (!hasnext) return;
//////
//////				DEPGRAPH_TYPE::neighbor_iterator nit(mD);
//////				int mval = 0;
//////				for (nit.begin(id); nit.valid(); nit.advance()) {
//////					if (mD->Before(nit.value(), id)) {
//////						int i = wave->GetPayload(nit.value());
//////
//////						if (i > mval) mval = i;
//////						break;
//////					}
//////				}
//////				mval++;
//////#ifdef USECUSTOMALLOCATOR
//////				int* pmval = int_alloc.allocate(1); //new int;
//////#else
//////				//int* pmval = new int;
//////#endif
//////				//(*pmval) = mval;
//////				count++;
//////				wave->SetAsProcessed(id, mval);
//////				if (count % targetnum == 0) fprintf(fout, "ASDFASDF %d %d %d %d \n", count, wave->mPriors.size(), wave->mWaiting.size(), wave->mReady2Go.size());
//////
//////				if (count == 10000){ //mD->NumElements() / 2) {
//////
//////					vector<HashingWaveFront*> splits;
//////					vector < pair<INDEX_TYPE, HashingWaveFront*> > obstructions;
//////					printf("about to split\n");
//////
//////					wave->deepSplitWaveFront(splits, obstructions);
//////					printf("done\n");
//////
//////					FILE* ftest = fopen("testids.txt", "w");
//////					for (auto it = wave->mPriors.begin(); it != wave->mPriors.end(); it++) {
//////						fprintf(ftest, "p %d\n", (*it).first);
//////					}
//////					for (auto it = wave->mWaiting.begin(); it != wave->mWaiting.end(); it++) {
//////						fprintf(ftest, "w %d\n", (*it).first);
//////					}
//////#ifdef QUEQUE
//////					stack<INDEX_TYPE> ms;// = wave->mReady2Go;
//////					queue<INDEX_TYPE> mq = wave->mReady2Go;
//////#endif
//////#ifdef RANDQ
//////					stack<INDEX_TYPE> ms = wave->mReady2Go.ms;
//////					queue<INDEX_TYPE> mq = wave->mReady2Go.mq;
//////#endif
//////#ifdef STACKQ
//////					stack<INDEX_TYPE> ms = wave->mReady2Go;
//////					queue<INDEX_TYPE> mq;// = wave->mReady2Go;
//////#endif
//////					while (!ms.empty()) {
//////						fprintf(ftest, "r %d\n", ms.top());
//////						ms.pop();
//////					}
//////					while (!mq.empty()) {
//////						fprintf(ftest, "r %d\n", mq.front());
//////						mq.pop();
//////					}
//////					fclose(ftest);
//////					ftest = fopen("testids_s.txt", "w");
//////					for (int i = 0; i < splits.size(); i++) {
//////						//printf("doing %d\n", i);
//////						for (auto it = splits[i]->mPriors.begin(); it != splits[i]->mPriors.end(); it++) {
//////							fprintf(ftest, "p %d %d\n", i, (*it).first);
//////						}
//////						for (auto it = splits[i]->mWaiting.begin(); it != splits[i]->mWaiting.end(); it++) {
//////							fprintf(ftest, "w %d %d\n", i, (*it).first);
//////						}
//////
//////
//////#ifdef QUEQUE
//////						stack<INDEX_TYPE> ms;// = wave->mReady2Go;
//////						queue<INDEX_TYPE> mq = splits[i]->mReady2Go;
//////#endif
//////#ifdef RANDQ
//////						stack<INDEX_TYPE> ms = splits[i]->mReady2Go.ms;
//////						queue<INDEX_TYPE> mq = splits[i]->mReady2Go.mq;
//////#endif
//////#ifdef STACKQ
//////						stack<INDEX_TYPE> ms = splits[i]->mReady2Go;
//////						queue<INDEX_TYPE> mq;// = wave->mReady2Go;
//////#endif
//////						while (!ms.empty()) {
//////							fprintf(ftest, "r %d %d\n", i, ms.top());
//////							ms.pop();
//////						}
//////						while (!mq.empty()) {
//////							fprintf(ftest, "r %d %d\n", i, mq.front());
//////							mq.pop();
//////						}
//////					}
//////					fclose(ftest);
//////				}
//////
//////
//////			}
//////			fclose(fout);
//////
//////
//////		}
//////
//////
//////	};


//////
//////	class TestingDependencyGraphIndependentTraversal {
//////	protected:
//////		DEPGRAPH_TYPE* mD;
//////	public:
//////		TestingDependencyGraphIndependentTraversal(DEPGRAPH_TYPE* d) : mD(d) {
//////
//////
//////		}
//////
//////		void doStuff() {
//////
//////
//////
//////			//HashingWaveFront* wave = new HashingWaveFront(mD);
//////
//////
//////#ifdef USECUSTOMALLOCATOR
//////
//////			FSBAllocator<int> int_alloc;
//////#endif			
//////			map<INDEX_TYPE, HashingWaveFront*> waves;
//////			DEPGRAPH_TYPE::vertex_iterator vit(mD);
//////			for (vit.begin(); vit.valid(); vit.advance()) {
//////				INDEX_TYPE id = vit.value();
//////				bool ismin = true;
//////				DEPGRAPH_TYPE::neighbor_iterator nit(mD);
//////				for (nit.begin(id); nit.valid(); nit.advance()) {
//////					if (mD->Before(nit.value(), id)) {
//////						ismin = false;
//////						break;
//////					}
//////				}
//////				if (ismin) {
//////#ifdef USECUSTOMALLOCATOR
//////					int* ip = int_alloc.allocate(1); // new int;
//////#else
//////					//int* ip = new int;
//////#endif			
//////
//////					//(*ip) = 0;
//////					int zero = 0;
//////					HashingWaveFront* w = new HashingWaveFront(mD);
//////					w->SetSeed(id, zero);
//////					waves[id] = w;
//////				}
//////			}
//////
//////
//////			//// do first expansion
//////			int randocount = 0;
//////			for (auto it = waves.begin(); it != waves.end(); it++) {
//////				HashingWaveFront* h = (*it).second;
//////				bool cont = true;
//////				while (cont) {
//////					INDEX_TYPE id;
//////					cont = h->GetNextToProcess(id);
//////					if (!cont) break;
//////					// this is where work goes for payload
//////					h->SetAsProcessed(id, randocount);
//////					randocount++;
//////
//////				}
//////			}
//////
//////			vector<HashingWaveFront*> h2;
//////			vector < pair<INDEX_TYPE, HashingWaveFront*> > oh2;
//////			for (auto it = waves.begin(); it != waves.end(); it++) {
//////				HashingWaveFront* h = (*it).second;
//////				vector<HashingWaveFront*> th;
//////				vector < pair<INDEX_TYPE, HashingWaveFront*> > oh;
//////				h->deepSplitWaveFront(th, oh);
//////				if (th.size() > 0)  {
//////					printf("whoa should not get here\n");
//////				}
//////				oh2.insert(oh2.end(), oh.begin(), oh.end());
//////
//////			}
//////
//////			FILE* ftest = fopen("testids_s.txt", "w");
//////
//////			for (auto hit = oh2.begin(); hit != oh2.end(); hit++) {
//////
//////				HashingWaveFront* h = (*hit).second;
//////				//printf("doing %d\n", i);
//////				for (auto it = h->mPriors.begin(); it != h->mPriors.end(); it++) {
//////					fprintf(ftest, "p %d %d\n", (*hit).first, (*it).first);
//////				}
//////				for (auto it = h->mWaiting.begin(); it != h->mWaiting.end(); it++) {
//////					fprintf(ftest, "w %d %d\n", (*hit).first, (*it).first);
//////				}
//////				fprintf(ftest, "r %d %d\n", (*hit).first, (*hit).first);
//////				//
//////				//#ifdef QUEQUE
//////				//				stack<INDEX_TYPE> ms;// = wave->mReady2Go;
//////				//				queue<INDEX_TYPE> mq = h->mReady2Go;
//////				//#endif
//////				//#ifdef RANDQ
//////				//				stack<INDEX_TYPE> ms = h->mReady2Go.ms;
//////				//				queue<INDEX_TYPE> mq = h->mReady2Go.mq;
//////				//#endif
//////				//#ifdef STACKQ
//////				//				stack<INDEX_TYPE> ms = h->mReady2Go;
//////				//				queue<INDEX_TYPE> mq;// = wave->mReady2Go;
//////				//#endif
//////				//				while (!ms.empty()) {
//////				//					fprintf(ftest, "r %d %d\n", (*hit).first, ms.top());
//////				//					ms.pop();
//////				//				}
//////				//				while (!mq.empty()) {
//////				//					fprintf(ftest, "r %d %d\n", (*hit).first, mq.front());
//////				//					mq.pop();
//////				//				}
//////			}
//////			fclose(ftest);
//////		}
//////
//////
//////
//////
//////
//////
//////
//////
//////
//////	};
//////
//////
}
#endif