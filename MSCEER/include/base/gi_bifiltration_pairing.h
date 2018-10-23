/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef BIFILTRATION_PAIRING_H
#define BIFILTRATION_PAIRING_H

#include "gi_labeling.h"
#include "gi_discrete_gradient_labeling.h"
//#include "gi_topological_simplicial_complex.h"

namespace GInt {

	template<class MeshType, class MaxVLType1, class MaxVLType2>
	class BifiltrationPairing
	{
	protected:

		MeshType* mMesh;
		MaxVLType1* mLabel1;
		MaxVLType2* mLabel2;
		DenseLabeling<INDEX_TYPE>* mPairs;

		typedef std::pair<INDEX_TYPE, INDEX_TYPE> vpair;
		vpair make_pair(INDEX_TYPE a, INDEX_TYPE b) {
			if (a < b) return vpair(a, b);
			return vpair(b, a);
		}

		struct cell_rec {
			DIM_TYPE dim;
			INDEX_TYPE id;

			cell_rec(DIM_TYPE d, INDEX_TYPE i) : dim(d), id(i) {}
			bool operator<(const cell_rec& a) const {
				if (this->dim < a.dim) return false;
				if (this->dim > a.dim) return true;
				return this->id > a.id;
			}
		};

		struct cell_pairing {
			INDEX_TYPE id;
			int num_missing;
			INDEX_TYPE pair;
			bool paired;
			cell_pairing() {
				//printf("SHOULD NEVER CALL\n");
			}
			cell_pairing(INDEX_TYPE id) : id(id) {
				pair = -1;
				paired = false;
				num_missing = 0;
			}
		};

		inline bool has_id(const std::vector<INDEX_TYPE>& cells, INDEX_TYPE id) {
			for (auto i : cells)
				if (i == id) return true;
			return false;
		}

		//std::queue<INDEX_TYPE> readytogo;

		void decrementCofacets(INDEX_TYPE id, std::map<INDEX_TYPE, cell_pairing>& cp) {
			typename TopologicalSimplicialComplex2d::CofacetsIterator cfit(mMesh);
			for (cfit.begin(id); cfit.valid(); cfit.advance()) {
				INDEX_TYPE cid = cfit.value();
				if (cp.count(cid) == 0) continue;
				cp[cid].num_missing--;
			}
		}

		void HomotopyExpand(std::vector<INDEX_TYPE>& cells, std::map<INDEX_TYPE, cell_pairing>& cp) {
			//std::priority_queue<cell_rec> sorted;

			//printf("in:");
			//for (auto id : cells) printf(" %llu", id);
			//printf("\n");

			std::vector<INDEX_TYPE> cells_of_dim[3];
			for (auto id : cells) {
				cells_of_dim[mMesh->dimension(id)].push_back(id);
				//sorted.push(cell_rec(mMesh->dimension(id), id));
				cell_pairing c(id);
				typename TopologicalSimplicialComplex2d::FacetsIterator fit(mMesh);
				for (fit.begin(id); fit.valid(); fit.advance()) {
					INDEX_TYPE nid = fit.value();
					if (has_id(cells, nid)) c.num_missing++;
				}
				//if (c.num_missing == 1) readytogo.push(c.id);
				cp[id] = c;
			}

			for (int i = 0; i < 3; i++) {
				//while (!sorted.empty()) {



				// logic is we need to process every cell of dimension i
				// until all have been processed, first try to pair
				// if no pairing was successful, make one critical and repeat
				int num_processed = 0;
				int total_to_process = 0;
				for (auto id : cells_of_dim[i]) if (!cp[id].paired) total_to_process++;
				while (num_processed < total_to_process) {

					int start_num_proc = num_processed;
					// try to pair as many as possible
					for (auto id : cells_of_dim[i]) {
						if (cp[id].paired) continue; // already paired
						if (cp[id].num_missing > 0)
							printf("ERROR: should never get here1\n");


						typename TopologicalSimplicialComplex2d::CofacetsIterator cfit(mMesh);
						for (cfit.begin(id); cfit.valid(); cfit.advance()) {
							INDEX_TYPE cfid = cfit.value();
							if (cp.count(cfid) == 0) continue; // not in our lower star
							if (cp[cfid].paired) printf("ERROR: should never get here2\n");
							if (cp[cfid].num_missing == 1) {
								// pair cells
								cp[id].pair = cfid;
								cp[id].paired = true;
								cp[cfid].pair = id;
								cp[cfid].paired = true;
								decrementCofacets(id, cp);
								decrementCofacets(cfid, cp);

								num_processed++;
								break;
							}
						}
					}

					if (start_num_proc == num_processed) {
						// then no more pairs were possible 
						for (auto id : cells_of_dim[i]) {
							if (cp[id].paired) continue; // already paired
							// make one critical and break
							cp[id].pair = id;
							cp[id].paired = true;
							decrementCofacets(id, cp);
							num_processed++;
							break;
						}
					}
				}

			}
			//printf("out:");
			//for (auto c : cp) {
			//	if (c.second.pair == c.second.id) printf(" (%d:%llu)", mMesh->dimension(c.second.id), c.second.pair);
			//	if (mMesh->dimension(c.second.id) < mMesh->dimension(c.second.pair)) 
			//		printf(" (%d:%llu->%d:%llu)", mMesh->dimension(c.second.id), c.second.id, mMesh->dimension(c.second.pair), c.second.pair);
			//}
			//printf("\n");
		}
	public:
		BifiltrationPairing(MeshType* mesh, MaxVLType1* label1, MaxVLType2* label2) : mMesh(mesh), mLabel1(label1), mLabel2(label2) {
			mPairs = new DenseLabeling<INDEX_TYPE>(mesh->numCells());
			mvarray1 = new INDEX_TYPE[mesh->numCells(0)];
			mvarray2 = new INDEX_TYPE[mesh->numCells(0)];
		}
		~BifiltrationPairing() {
			delete mPairs;
			delete[] mvarray1;
			delete[] mvarray2;
		}
		std::set<INDEX_TYPE> mHasCriticalNeighbor;
		INDEX_TYPE* mvarray1;
		INDEX_TYPE* mvarray2;

		struct sorter1{
			MaxVLType1* l;
			bool operator()(INDEX_TYPE a, INDEX_TYPE b) const {
				return l->Before(a, b);
			}
		};
		struct sorter2{
			MaxVLType2* l;
			bool operator()(INDEX_TYPE a, INDEX_TYPE b) const {
				return l->Before(a, b);
			}
		};
		void ComputePairing() {

			printf("computing pairing1\n");
			INDEX_TYPE* tvarray1 = new INDEX_TYPE[mMesh->numCells(0)];
			INDEX_TYPE* tvarray2 = new INDEX_TYPE[mMesh->numCells(0)];
			for (INDEX_TYPE i = 0; i < mMesh->numCells(0); i++) {
				tvarray1[i] = i;
				tvarray2[i] = i;
			}
			printf("computing pairing2\n");

			sorter1 s1; s1.l = mLabel1;
			std::sort(tvarray1, tvarray1 + mMesh->numCells(0), s1);
			printf("computing pairing3\n");
			sorter2 s2; s2.l = mLabel2;
			std::sort(tvarray2, tvarray2 + mMesh->numCells(0), s2);
			printf("computing pairing4\n");

			for (INDEX_TYPE i = 0; i < mMesh->numCells(0); i++) {
				mvarray1[tvarray1[i]] = i;
				mvarray2[tvarray2[i]] = i;
			}
			printf("computing pairing5\n");
			delete[] tvarray1;
			delete[] tvarray2;
			printf("computing pairing6\n");

			typename TopologicalSimplicialComplex2d::DCellsIterator verts(mMesh, 0);
			for (verts.begin(); verts.valid(); verts.advance()){
				INDEX_TYPE vert_GI = verts.value();
				std::map<vpair, std::vector<INDEX_TYPE> > subsets;

				subsets[make_pair(vert_GI, vert_GI)].push_back(vert_GI);
				typename TopologicalSimplicialComplex2d::AdjacentCellsIterator star(mMesh);
				for (star.begin(vert_GI); star.valid(); star.advance()) {
					INDEX_TYPE ncell_GI = star.value();

					// make the sets!
					INDEX_TYPE vid1 = mLabel1->Cell2HighestVertex(ncell_GI);
					if (vid1 != vert_GI) continue; // not in lower star of f1

					INDEX_TYPE vid2 = mLabel2->Cell2HighestVertex(ncell_GI);

					vpair vp = make_pair(vid1, vid2);
					subsets[vp].push_back(ncell_GI);

				}
				//if (subsets.size() != 1) printf("WHOATHAHAHHSHSHSHSH\n");
				for (auto S : subsets) {
					//if (S.second.size() == 1) {
					//	mPairs->SetLabel(S.second[0], S.second[0]);
					//}
					//else if (S.second.size() == 2) {
					//	mPairs->SetLabel(S.second[0], S.second[1]);
					//	mPairs->SetLabel(S.second[1], S.second[0]);
					//}
					//else {
					// do homotopy expand
					std::map<INDEX_TYPE, cell_pairing> cp;
					HomotopyExpand(S.second, cp);
					bool hascritn = false;
					for (auto p : cp) {
						mPairs->SetLabel(p.first, p.second.pair);
						mPairs->SetLabel(p.second.pair, p.first);

						if (p.first == p.second.pair) hascritn = true;

					}
					if (hascritn){
						for (auto p : cp) {
							mHasCriticalNeighbor.insert(p.first);
							mHasCriticalNeighbor.insert(p.second.pair);
						}
					}
					//}
				}

			}
		}

DenseLabeling<INDEX_TYPE>* GetLabeling() { return mPairs; }

	};





}
#endif
