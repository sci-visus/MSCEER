#ifndef BIFILTRATION_PAIRING_H
#define BIFILTRATION_PAIRING_H

#include "gi_labeling.h"
#include "gi_discrete_gradient_labeling.h"
#include "gi_topological_simplicial_complex.h"

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

		struct sorter1 {
			MaxVLType1* l;
			bool operator()(INDEX_TYPE a, INDEX_TYPE b) const {
				return l->Before(a, b);
			}
		};
		struct sorter2 {
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
			for (verts.begin(); verts.valid(); verts.advance()) {
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
					if (hascritn) {
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



	template<class MeshType, class MaxVLType, class GradType>
	class MyRobins
	{
	protected:

		MeshType* mMesh;
		MaxVLType* mLabel1;
		DenseLabeling<char>* mLabel2;
		//DenseLabeling<INDEX_TYPE>* mPairs;
		GradType* mGrad;

		typedef std::pair<INDEX_TYPE, INDEX_TYPE> vpair;
		vpair make_pair(INDEX_TYPE a, INDEX_TYPE b) {
			return vpair(a, b);
			//return vpair(b, a);
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
			typename MeshType::CofacetsIterator cfit(mMesh);
			for (cfit.begin(id); cfit.valid(); cfit.advance()) {
				INDEX_TYPE cid = cfit.value();
				if (cp.count(cid) == 0) continue;
				cp[cid].num_missing--;
			}
		}

		bool is_steeper(INDEX_TYPE edge1, INDEX_TYPE edge2, INDEX_TYPE vid) const {
			INDEX_TYPE vido1;
			INDEX_TYPE vido2;
			typename MeshType::FacetsIterator fit1(mMesh);
			typename MeshType::FacetsIterator fit2(mMesh);
			fit1.begin(edge1);
			fit2.begin(edge2);
			if (vid == fit1.value()) fit1.advance();
			if (vid == fit2.value()) fit2.advance();
			vido1 = fit1.value();
			vido2 = fit2.value();
			return mLabel1->Before(mMesh->VertexNumberFromCellID(vido1), mMesh->VertexNumberFromCellID(vido2));
		}
	public:
		INDEX_TYPE lowest_vertex(INDEX_TYPE cid) {
			typename MeshType::CellVerticesIterator vit(mMesh);
			vit.begin(cid);
			INDEX_TYPE lowest = vit.value();
			vit.advance();
			while (vit.valid()) {
				INDEX_TYPE other = vit.value();
				if (mLabel1->Before(mMesh->VertexNumberFromCellID(other), mMesh->VertexNumberFromCellID(lowest))) lowest = other;
				vit.advance();
			}
			return lowest;
		}
	protected:
		INDEX_TYPE PickLowestCandidate(std::vector<INDEX_TYPE>& cands) {
			if (cands.size() == 1) return cands[0];

			INDEX_TYPE curr_lowest = cands[0];
			INDEX_TYPE lv = lowest_vertex(curr_lowest);
			for (int i = 1; i < cands.size(); i++) {
				INDEX_TYPE olv = lowest_vertex(cands[i]);
				if (mLabel1->Before(mMesh->VertexNumberFromCellID(olv), mMesh->VertexNumberFromCellID(lv))) {
					lv = olv;
					curr_lowest = cands[i];
				}
			}
			return curr_lowest;
		}


		void HomotopyExpand(std::vector<INDEX_TYPE>& cells, std::map<INDEX_TYPE, cell_pairing>& cp) {
			//std::priority_queue<cell_rec> sorted;

			//printf("in:");
			//for (auto id : cells) printf(" %llu", id);
			//printf("\n");

			std::vector<INDEX_TYPE> cells_of_dim[4];
			for (auto id : cells) {
				cells_of_dim[mMesh->dimension(id)].push_back(id);
				//sorted.push(cell_rec(mMesh->dimension(id), id));
				cell_pairing c(id);
				typename MeshType::FacetsIterator fit(mMesh);
				for (fit.begin(id); fit.valid(); fit.advance()) {
					INDEX_TYPE nid = fit.value();
					if (has_id(cells, nid)) c.num_missing++;
				}
				//if (c.num_missing == 1) readytogo.push(c.id);
				cp[id] = c;
			}

			// is there a vertex?
			if (cells_of_dim[0].size() > 0) {

				if (cells_of_dim[0].size() > 1) {
					printf("ERROR: too many vertices %d\n", cells_of_dim[0].size());
					for (auto asdf : cells_of_dim[0]) printf("%llu ", asdf);
					printf("\n");
				}
				INDEX_TYPE vid = cells_of_dim[0][0];
				// if so, find the lowest edge, pair that way
				if (cells_of_dim[1].size() == 0) {
					// make vertex critical
					cp[vid].pair = vid;
					cp[vid].paired = true;
					decrementCofacets(vid, cp);
				}
				else {

					INDEX_TYPE lowest_edge_id;
					if (cells_of_dim[1].size() == 1) {
						// just pair with only option 
						lowest_edge_id = cells_of_dim[1][0];
					}
					else {
						// find minimal edge

						lowest_edge_id = cells_of_dim[1][0];

						for (int i = 1; i < cells_of_dim[1].size(); i++) {
							if (is_steeper(cells_of_dim[1][i], lowest_edge_id, vid)) {
								lowest_edge_id = cells_of_dim[1][i];
							}
						}


					}
					// pair in direction of steepest descent
					cp[vid].pair = lowest_edge_id;
					cp[vid].paired = true;
					cp[lowest_edge_id].pair = vid;
					cp[lowest_edge_id].paired = true;
					decrementCofacets(vid, cp);
					decrementCofacets(lowest_edge_id, cp);

				}
			}



			for (int i = 0; i < 4; i++) {
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

						std::vector<INDEX_TYPE> candidates;
						typename MeshType::CofacetsIterator cfit(mMesh);
						for (cfit.begin(id); cfit.valid(); cfit.advance()) {
							INDEX_TYPE cfid = cfit.value();
							if (cp.count(cfid) == 0) continue; // not in our lower star
							if (cp[cfid].paired) printf("ERROR: should never get here2\n");
							if (cp[cfid].num_missing == 1) {
								// pair cells
								candidates.push_back(cfid);

							}
						}

						//if (candidates.size() > 1) printf("got here candidates: %d\n", candidates.size());
						if (candidates.size() > 0) {
							INDEX_TYPE cfid = PickLowestCandidate(candidates);
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

					if (start_num_proc == num_processed) {
						// then no more pairs were possible 
						std::vector<INDEX_TYPE> candidates;
						for (auto id : cells_of_dim[i]) {
							if (cp[id].paired) continue; // already paired
							// make one critical and break
							candidates.push_back(id);
						}
						INDEX_TYPE id = PickLowestCandidate(candidates);
						//asdf want to make lowest critical!?
						cp[id].pair = id;
						cp[id].paired = true;
						decrementCofacets(id, cp);
						num_processed++;

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
		MyRobins(MeshType* mesh, MaxVLType* label1, DenseLabeling<char>* label2,
			GradType* grad) : mMesh(mesh), mLabel1(label1), mLabel2(label2), mGrad(grad) {
		}
		~MyRobins() {
		}

		void ComputeLowerStar(INDEX_TYPE vert_GI) {
			// for each vertex compute the sets that will make up its lower start with restriction
			std::map<char, std::vector<INDEX_TYPE> > subsets;

			// add the vertex to restriction label set
			//subsets[mLabel2->GetLabel(vert_GI)].push_back(vert_GI);

			// now add all lower star to restriciton sets 
			typename MeshType::AdjacentCellsIterator star(mMesh);
			for (star.begin(vert_GI); star.valid(); star.advance()) {
				INDEX_TYPE ncell_GI = star.value();

				// discard a cell if its highest vertex is NOT the vertex, hence not part of lower star
				INDEX_TYPE vid1 = mLabel1->Cell2HighestVertex(ncell_GI);
				if (vid1 != vert_GI) continue; // not in lower star of f1

				// get the restriction label of th lower start cell
				char vid2 = mLabel2->GetLabel(ncell_GI);
				vid2 = vid2 * mMesh->maxDim() + mMesh->boundaryValue(vid2);
				// add the cell to the right subset
				subsets[vid2].push_back(ncell_GI);

			}

			// now do homotopy expand on each subset!
			for (auto S : subsets) {

				// do homotopy expand
				std::map<INDEX_TYPE, cell_pairing> cp;
				HomotopyExpand(S.second, cp);
				for (auto p : cp) {
					if (p.first == p.second.pair) {
						mGrad->setCritical(p.first, true);
						mGrad->setAssigned(p.first, 1);
					}
					else {
						mGrad->setPair(p.first, p.second.pair);
						mGrad->setPair(p.second.pair, p.first);
						mGrad->setAssigned(p.first, 1);
						mGrad->setAssigned(p.second.pair, 1);
					}
				}

			}

		}
		void ComputePairing() {
			std::vector<INDEX_TYPE> topo_index_partition;
			int num_threads;
#pragma omp parallel
			{
#pragma omp single
				{
					num_threads = omp_get_num_threads();
					ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, topo_index_partition);
				}

				int thread_num = omp_get_thread_num();

				// iterate over all vertices
				typename MeshType::DCellsIterator verts(mMesh, 0, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (verts.begin(); verts.valid(); verts.advance()) {
					INDEX_TYPE vert_GI = verts.value();

					ComputeLowerStar(vert_GI);


				}
			}
		}
		//DenseLabeling<INDEX_TYPE>* GetLabeling() { return mPairs; }

	};

	template<class MeshType, class MaxVLType, class GradType, int LMAXDIM, int LMAX_LABEL>
	class MyRobinsNoalloc
	{
	protected:

		MeshType* mMesh;
		MaxVLType* mLabel1;
		DenseLabeling<char>* mLabel2;
		//DenseLabeling<INDEX_TYPE>* mPairs;
		GradType* mGrad;

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
				pair = -1;
				paired = false;
				num_missing = 0;
			}

		};

		struct myStaticMap {
			INDEX_TYPE idarray[27];
			cell_pairing cparray[27];
			int size;

			void push_back(INDEX_TYPE val) {
				idarray[size] = val;
				++size;
			}

			myStaticMap() : size(0) {}

			void sort() {
				// simple bubble sort since we actually expect in-order things most of the time
				int i, j;
				for (i = 0; i < size - 1; i++)

					// Last i elements are already in place   
					for (j = 0; j < size - i - 1; j++)
						if (idarray[j] > idarray[j + 1]) {
							//swap(&arr[j], &arr[j + 1]);
							INDEX_TYPE temp = idarray[j];
							idarray[j] = idarray[j + 1];
							idarray[j + 1] = temp;
						}
				for (i = 0; i < size; i++) {
					cparray[i].id = idarray[i];
				}
			};
			//#define BINARYFIND
#ifdef BINARYFIND
			// binary search to find index of id, if not found return -1
			int find(INDEX_TYPE id) const {
				int first = 0;                 //first array element
				int last = size - 1;            //last array element
				int	middle;                        //mid point of search

				while (first <= last)
				{
					middle = (first + last) >> 1; //this finds the mid point
					if (idarray[middle] == id) {
						return middle;
					}
					else if (idarray[middle] > id) // if it's in the lower half
					{
						last = middle - 1;
					}
					else {
						first = middle + 1;                 //if it's in the upper half
					}
				}
				return -1;  // not found
			}
#else
			int find(INDEX_TYPE id) const {
				for (int i = 0; i < size; i++) {
					if (idarray[i] == id) return i;
				}
				return -1;
			}
#endif

		};


		inline bool has_id(const small_INDEX_vector& cells, INDEX_TYPE id) {
			for (int i = 0; i < cells.size; i++)
				if (i == id) return true;
			return false;
		}

		//std::queue<INDEX_TYPE> readytogo;

		void decrementCofacets(INDEX_TYPE id, myStaticMap& cp) {
			typename MeshType::CofacetsIterator cfit(mMesh);
			for (cfit.begin(id); cfit.valid(); cfit.advance()) {
				INDEX_TYPE cid = cfit.value();
				int offset = cp.find(cid);
				//if (offset != cp.slowfind(cid)) {
				//	printf("S:KLSDJFSL:DKFJSDKL:FJ\n");
				//}
				if (offset == -1) continue;
				cp.cparray[offset].num_missing--;
			}
		}

		INDEX_TYPE is_steeper(INDEX_TYPE edge1, INDEX_TYPE vertex_number, INDEX_TYPE vid) const {
			INDEX_TYPE vido2;
			typename MeshType::FacetsIterator fit2(mMesh);
			fit2.begin(edge1);
			if (vid == fit2.value()) fit2.advance();
			vido2 = fit2.value();
			INDEX_TYPE r_vertex_number = mMesh->VertexNumberFromCellID(vido2);
			if (mLabel1->Before(vertex_number, r_vertex_number)) return vertex_number;
			return r_vertex_number;
			//return mLabel1->Before(mMesh->VertexNumberFromCellID(vido1), mMesh->VertexNumberFromCellID(vido2));
		}
	public:
		INDEX_TYPE lowest_vertex(INDEX_TYPE cid) {
			return mMesh->VertexNumberFromCellID(mLabel1->Cell2LowestVertex(cid));

			//typename MeshType::CellVerticesIterator vit(mMesh);
			//vit.begin(cid);
			//INDEX_TYPE lowest = vit.value();
			//INDEX_TYPE lowest_vertex_id = mMesh->VertexNumberFromCellID(lowest);
			//vit.advance();
			//while (vit.valid()) {
			//	INDEX_TYPE other = vit.value();
			//	INDEX_TYPE other_vertex_id = mMesh->VertexNumberFromCellID(other);
			//	if (mLabel1->Before(other_vertex_id, lowest_vertex_id)) {
			//		lowest = other;
			//		lowest_vertex_id = other_vertex_id;
			//	}
			//	vit.advance();
			//}
			//return lowest_vertex_id;
		}
	protected:
		INDEX_TYPE PickLowestCandidate(small_INDEX_vector& cands, myStaticMap& cp) {
			if (cands.size == 1) return cands[0];

			INDEX_TYPE curr_lowest_listid = cands[0];
			INDEX_TYPE lv_vid = lowest_vertex(cp.idarray[cands[0]]);
			if (lv_vid < 0) {
				printf("got here\n");
			}

			for (int i = 1; i < cands.size; i++) {
				INDEX_TYPE olv_vid = lowest_vertex(cp.idarray[cands[i]]);
				if (mLabel1->Before(olv_vid, lv_vid)) {
					lv_vid = olv_vid;
					curr_lowest_listid = cands[i];
				}
			}
			return curr_lowest_listid;
		}


		//		//void HomotopyExpand(std::vector<INDEX_TYPE>& cells, std::map<INDEX_TYPE, cell_pairing>& cp) {
		//		void HomotopyExpand(small_INDEX_vector& cells, myStaticMap& cp, std::vector<std::pair<INDEX_TYPE, INDEX_TYPE> >& ordering) {
		//			//std::priority_queue<cell_rec> sorted;
		//
		//			//printf("in:");
		//			//for (auto id : cells) printf(" %llu", id);
		//			//printf("\n");
		//
		//
		//			small_INDEX_vector listid_of_d_cells[4];
		//			for (int i = 0; i < cells.size; i++) {
		//				INDEX_TYPE id = cells[i];
		//				cp.push_back(id);
		//			}
		//#ifdef BINARYFIND
		//			cp.sort();
		//#endif
		//
		//			for (int i = 0; i < cp.size; i++) {
		//				INDEX_TYPE id = cp.idarray[i];
		//				//for (auto id : cells) {
		//				listid_of_d_cells[mMesh->dimension(id)].push_back(i);
		//				//sorted.push(cell_rec(mMesh->dimension(id), id));
		//
		//				typename MeshType::FacetsIterator fit(mMesh);
		//				for (fit.begin(id); fit.valid(); fit.advance()) {
		//					INDEX_TYPE nid = fit.value();
		//					int id2 = cp.find(nid);
		//					//if (id2 != cp.slowfind(nid)) {
		//					//	printf("S:KLSDJFSL:DKFJSDKL:FJ\n");
		//					//}
		//					if (id2 != -1) cp.cparray[i].num_missing++;
		//				}
		//			}
		//
		//
		//			// is there a vertex?
		//			if (listid_of_d_cells[0].size > 0) {
		//
		//				if (listid_of_d_cells[0].size > 1) {
		//					printf("ERROR: too many vertices %d\n", listid_of_d_cells[0].size);
		//					//for (auto asdf : listid_of_d_cells[0]) printf("%llu ", asdf);
		//					printf("\n");
		//				}
		//				INDEX_TYPE vid_listid = listid_of_d_cells[0][0];
		//				INDEX_TYPE vid = cp.idarray[vid_listid];
		//				if (vid_listid == -1) {
		//					printf("ERROR SHOULD NEVER GET HERE\n");
		//				}
		//				// if so, find the lowest edge, pair that way
		//				if (listid_of_d_cells[1].size == 0) {
		//					// make vertex critical
		//					cp.cparray[vid_listid].pair = cp.idarray[vid_listid];
		//					cp.cparray[vid_listid].paired = true;
		//					decrementCofacets(cp.idarray[vid_listid], cp);
		//					ordering.push_back(std::pair<INDEX_TYPE, INDEX_TYPE>(vid, vid));
		//				}
		//				else {
		//
		//					INDEX_TYPE lowest_edge_id;
		//					int lowest_edge_listid;
		//					if (listid_of_d_cells[1].size == 1) {
		//						// just pair with only option 
		//						lowest_edge_listid = listid_of_d_cells[1][0];
		//						lowest_edge_id = cp.idarray[lowest_edge_listid];
		//					}
		//					else {
		//						// find minimal edge
		//
		//						lowest_edge_listid = listid_of_d_cells[1][0];
		//						lowest_edge_id = cp.idarray[lowest_edge_listid];
		//
		//						/* this block to cache the lowest vertex number*/
		//						typename MeshType::FacetsIterator fit1(mMesh);
		//						fit1.begin(lowest_edge_id);
		//						if (vid == fit1.value()) fit1.advance();
		//						INDEX_TYPE lowest_edge_other_vertex_number = mMesh->VertexNumberFromCellID(fit1.value());
		//						/* endblock */
		//
		//						for (int i = 1; i < listid_of_d_cells[1].size; i++) {
		//							int other_listid = listid_of_d_cells[1][i];
		//							INDEX_TYPE other_edge_lower_vertex_number = is_steeper(cp.idarray[other_listid], lowest_edge_other_vertex_number, vid);
		//							if (other_edge_lower_vertex_number != lowest_edge_other_vertex_number) {
		//								lowest_edge_id = cp.idarray[other_listid];
		//								lowest_edge_listid = other_listid;
		//								lowest_edge_other_vertex_number = other_edge_lower_vertex_number;
		//							}
		//						}
		//
		//
		//					}
		//					if (vid_listid == -1) {
		//						printf("ERROR SHOULD NEVER GET HERE2\n");
		//					}
		//					// pair in direction of steepest descent
		//					cp.cparray[vid_listid].pair = lowest_edge_id;
		//					cp.cparray[vid_listid].paired = true;
		//					cp.cparray[lowest_edge_listid].pair = vid;
		//					cp.cparray[lowest_edge_listid].paired = true;
		//					decrementCofacets(vid, cp);
		//					decrementCofacets(lowest_edge_id, cp);
		//					ordering.push_back(std::pair<INDEX_TYPE, INDEX_TYPE>(vid, lowest_edge_id));
		//
		//				}
		//			}
		//
		//
		//
		//			for (int i = 0; i < 4; i++) {
		//				//while (!sorted.empty()) {
		//
		//
		//
		//				// logic is we need to process every cell of dimension i
		//				// until all have been processed, first try to pair
		//				// if no pairing was successful, make one critical and repeat
		//				int num_processed = 0;
		//				int total_to_process = 0;
		//				for (int j = 0; j < listid_of_d_cells[i].size; j++) {
		//					INDEX_TYPE i_cell_listid = listid_of_d_cells[i][j];
		//					if (!cp.cparray[i_cell_listid].paired) total_to_process++;
		//				}
		//				while (num_processed < total_to_process) {
		//
		//					int start_num_proc = num_processed;
		//					// try to pair as many as possible
		//					for (int j = 0; j < listid_of_d_cells[i].size; j++) {
		//						INDEX_TYPE i_cell_listid = listid_of_d_cells[i][j];
		//						if (cp.cparray[i_cell_listid].paired) continue; // already paired
		//						if (cp.cparray[i_cell_listid].num_missing > 0) {
		//							printf("ERROR: should never get here1\n");
		//						}
		//						INDEX_TYPE i_cell = cp.idarray[i_cell_listid];
		//						small_INDEX_vector candidates;
		//						typename MeshType::CofacetsIterator cfit(mMesh);
		//						for (cfit.begin(i_cell); cfit.valid(); cfit.advance()) {
		//							INDEX_TYPE cfid = cfit.value();
		//							int cofacet_listid = cp.find(cfid);
		//							//if (cofacet_listid != cp.slowfind(cfid)) {
		//							//	printf("S:KLSDJFSL:DKFJSDKL:FJ\n");
		//							//}
		//
		//							if (cofacet_listid == -1) continue; // not in our lower star
		//							if (cp.cparray[cofacet_listid].paired) {
		//								printf("ERROR: should never get here2\n");
		//							}
		//							if (cp.cparray[cofacet_listid].num_missing == 1) {
		//								// pair cells
		//								candidates.push_back(cofacet_listid);
		//
		//							}
		//						}
		//
		//						//if (candidates.size() > 1) printf("got here candidates: %d\n", candidates.size());
		//						if (candidates.size > 0) {
		//							INDEX_TYPE lowest_cofacet_listid = PickLowestCandidate(candidates, cp);
		//							INDEX_TYPE cfid = cp.idarray[lowest_cofacet_listid];
		//							cp.cparray[i_cell_listid].pair = cfid;
		//							cp.cparray[i_cell_listid].paired = true;
		//							cp.cparray[lowest_cofacet_listid].pair = i_cell;
		//							cp.cparray[lowest_cofacet_listid].paired = true;
		//							decrementCofacets(i_cell, cp);
		//							decrementCofacets(cfid, cp);
		//							ordering.push_back(std::pair<INDEX_TYPE, INDEX_TYPE>(i_cell, cfid));
		//
		//							num_processed++;
		//							break;
		//						}
		//					}
		//
		//					if (start_num_proc == num_processed) {
		//						// then no more pairs were possible 
		//						small_INDEX_vector candidates;
		//						for (int j = 0; j < listid_of_d_cells[i].size; j++) {
		//							INDEX_TYPE i_cell_listid = listid_of_d_cells[i][j];
		//							if (cp.cparray[i_cell_listid].paired) continue; // already paired
		//							// make one critical and break
		//							candidates.push_back(i_cell_listid);
		//						}
		//						INDEX_TYPE offset2 = PickLowestCandidate(candidates, cp);
		//						INDEX_TYPE id = cp.idarray[offset2];
		//						//asdf want to make lowest critical!?
		//						cp.cparray[offset2].pair = id;
		//						cp.cparray[offset2].paired = true;
		//						decrementCofacets(id, cp);
		//						ordering.push_back(std::pair<INDEX_TYPE, INDEX_TYPE>(id, id));
		//						num_processed++;
		//
		//					}
		//				}
		//
		//			}
		//			//printf("out:");
		//			//for (auto c : cp) {
		//			//	if (c.second.pair == c.second.id) printf(" (%d:%llu)", mMesh->dimension(c.second.id), c.second.pair);
		//			//	if (mMesh->dimension(c.second.id) < mMesh->dimension(c.second.pair)) 
		//			//		printf(" (%d:%llu->%d:%llu)", mMesh->dimension(c.second.id), c.second.id, mMesh->dimension(c.second.pair), c.second.pair);
		//			//}
		//			//printf("\n");
		//		}

				//void HomotopyExpand(std::vector<INDEX_TYPE>& cells, std::map<INDEX_TYPE, cell_pairing>& cp) {
		void HomotopyExpand(small_INDEX_vector& cells, myStaticMap& cp) {
			//std::priority_queue<cell_rec> sorted;

		//printf("in:");
		//for (auto id : cells) printf(" %llu", id);
		//printf("\n");
			vector<Vec2l> coords;
	
			for (int i = 0; i < cells.size; i++) {
				Vec2l v;
				mMesh->cellid2Coords(cells[i], v);
				coords.push_back(v);
			}
			small_INDEX_vector listid_of_d_cells[4];
			for (int i = 0; i < cells.size; i++) {
				INDEX_TYPE id = cells[i];
				cp.push_back(id);
			}
#ifdef BINARYFIND
			cp.sort();
#endif
			for (int i = 0; i < cp.size; i++) {
				INDEX_TYPE id = cp.idarray[i];
				//for (auto id : cells) {
				listid_of_d_cells[mMesh->dimension(id)].push_back(i);
				//sorted.push(cell_rec(mMesh->dimension(id), id));

				typename MeshType::FacetsIterator fit(mMesh);
				for (fit.begin(id); fit.valid(); fit.advance()) {
					INDEX_TYPE nid = fit.value();
					int id2 = cp.find(nid);
					//if (id2 != cp.slowfind(nid)) {
					//	printf("S:KLSDJFSL:DKFJSDKL:FJ\n");
					//}
					if (id2 != -1) cp.cparray[i].num_missing++;
				}
			}


			// is there a vertex?
			if (listid_of_d_cells[0].size > 0) {

				if (listid_of_d_cells[0].size > 1) {
					printf("ERROR: too many vertices %d\n", listid_of_d_cells[0].size);
					//for (auto asdf : listid_of_d_cells[0]) printf("%llu ", asdf);
					printf("\n");
				}
				INDEX_TYPE vid_listid = listid_of_d_cells[0][0];
				INDEX_TYPE vid = cp.idarray[vid_listid];
				if (vid_listid == -1) {
					printf("ERROR SHOULD NEVER GET HERE\n");
				}
				// if so, find the lowest edge, pair that way
				if (listid_of_d_cells[1].size == 0) {
					// make vertex critical
					cp.cparray[vid_listid].pair = cp.idarray[vid_listid];
					cp.cparray[vid_listid].paired = true;
					decrementCofacets(cp.idarray[vid_listid], cp);
				}
				else {

					INDEX_TYPE lowest_edge_id;
					int lowest_edge_listid;
					if (listid_of_d_cells[1].size == 1) {
						// just pair with only option 
						lowest_edge_listid = listid_of_d_cells[1][0];
						lowest_edge_id = cp.idarray[lowest_edge_listid];
					}
					else {
						// find minimal edge

						lowest_edge_listid = listid_of_d_cells[1][0];
						lowest_edge_id = cp.idarray[lowest_edge_listid];

						/* this block to cache the lowest vertex number*/
						typename MeshType::FacetsIterator fit1(mMesh);
						fit1.begin(lowest_edge_id);
						if (vid == fit1.value()) fit1.advance();
						INDEX_TYPE lowest_edge_other_vertex_number = mMesh->VertexNumberFromCellID(fit1.value());
						/* endblock */

						for (int i = 1; i < listid_of_d_cells[1].size; i++) {
							int other_listid = listid_of_d_cells[1][i];
							INDEX_TYPE other_edge_lower_vertex_number = is_steeper(cp.idarray[other_listid], lowest_edge_other_vertex_number, vid);
							if (other_edge_lower_vertex_number != lowest_edge_other_vertex_number) {
								lowest_edge_id = cp.idarray[other_listid];
								lowest_edge_listid = other_listid;
								lowest_edge_other_vertex_number = other_edge_lower_vertex_number;
							}
						}


					}
					if (vid_listid == -1) {
						printf("ERROR SHOULD NEVER GET HERE2\n");
					}
					// pair in direction of steepest descent
					cp.cparray[vid_listid].pair = lowest_edge_id;
					cp.cparray[vid_listid].paired = true;
					cp.cparray[lowest_edge_listid].pair = vid;
					cp.cparray[lowest_edge_listid].paired = true;
					decrementCofacets(vid, cp);
					decrementCofacets(lowest_edge_id, cp);

				}
			}



			for (int i = 0; i < 4; i++) {
				//while (!sorted.empty()) {



				// logic is we need to process every cell of dimension i
				// until all have been processed, first try to pair
				// if no pairing was successful, make one critical and repeat
				int num_processed = 0;
				int total_to_process = 0;
				for (int j = 0; j < listid_of_d_cells[i].size; j++) {
					INDEX_TYPE i_cell_listid = listid_of_d_cells[i][j];
					if (!cp.cparray[i_cell_listid].paired) total_to_process++;
				}
				while (num_processed < total_to_process) {

					int start_num_proc = num_processed;
					// try to pair as many as possible
					for (int j = 0; j < listid_of_d_cells[i].size; j++) {
						INDEX_TYPE i_cell_listid = listid_of_d_cells[i][j];
						if (cp.cparray[i_cell_listid].paired) continue; // already paired
						if (cp.cparray[i_cell_listid].num_missing > 0) {
							printf("ERROR: should never get here1\n");
						}
						INDEX_TYPE i_cell = cp.idarray[i_cell_listid];
						small_INDEX_vector candidates;
						typename MeshType::CofacetsIterator cfit(mMesh);
						for (cfit.begin(i_cell); cfit.valid(); cfit.advance()) {
							INDEX_TYPE cfid = cfit.value();
							int cofacet_listid = cp.find(cfid);
							//if (cofacet_listid != cp.slowfind(cfid)) {
							//	printf("S:KLSDJFSL:DKFJSDKL:FJ\n");
							//}

							if (cofacet_listid == -1) continue; // not in our lower star
							if (cp.cparray[cofacet_listid].paired) {
								printf("ERROR: should never get here2\n");
							}
							if (cp.cparray[cofacet_listid].num_missing == 1) {
								// pair cells
								candidates.push_back(cofacet_listid);

							}
						}

						//if (candidates.size() > 1) printf("got here candidates: %d\n", candidates.size());
						if (candidates.size > 0) {
							INDEX_TYPE lowest_cofacet_listid = PickLowestCandidate(candidates, cp);
							INDEX_TYPE cfid = cp.idarray[lowest_cofacet_listid];
							cp.cparray[i_cell_listid].pair = cfid;
							cp.cparray[i_cell_listid].paired = true;
							cp.cparray[lowest_cofacet_listid].pair = i_cell;
							cp.cparray[lowest_cofacet_listid].paired = true;
							decrementCofacets(i_cell, cp);
							decrementCofacets(cfid, cp);

							num_processed++;
							break;
						}
					}

					if (start_num_proc == num_processed) {
						// then no more pairs were possible 
						small_INDEX_vector candidates;
						for (int j = 0; j < listid_of_d_cells[i].size; j++) {
							INDEX_TYPE i_cell_listid = listid_of_d_cells[i][j];
							if (cp.cparray[i_cell_listid].paired) continue; // already paired
							// make one critical and break
							candidates.push_back(i_cell_listid);
						}
						INDEX_TYPE offset2 = PickLowestCandidate(candidates, cp);
						INDEX_TYPE id = cp.idarray[offset2];
						//asdf want to make lowest critical!?
						cp.cparray[offset2].pair = id;
						cp.cparray[offset2].paired = true;
						decrementCofacets(id, cp);
						num_processed++;

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
		MyRobinsNoalloc(MeshType* mesh, MaxVLType* label1, DenseLabeling<char>* label2,
			GradType* grad) : mMesh(mesh), mLabel1(label1), mLabel2(label2), mGrad(grad) {
		}

		MyRobinsNoalloc(MeshType* mesh, MaxVLType* label1,
			GradType* grad) : mMesh(mesh), mLabel1(label1), mLabel2(NULL), mGrad(grad) {
		}
		~MyRobinsNoalloc() {
		}





		void LookAtLowerStarOrder(INDEX_TYPE vert_GI, std::vector<std::pair<INDEX_TYPE, INDEX_TYPE> >& ordering) {
			// for each vertex compute the sets that will make up its lower start with restriction
			//std::map<char, std::vector<INDEX_TYPE> > subsets;

			small_INDEX_vector subsets_map[LMAXDIM * LMAX_LABEL];

			// add the vertex to restriction label set
			//subsets[mLabel2->GetLabel(vert_GI)].push_back(vert_GI);

			// now add all lower star to restriciton sets 
			typename MeshType::AdjacentCellsIterator star(mMesh);
			for (star.begin(vert_GI); star.valid(); star.advance()) {
				INDEX_TYPE ncell_GI = star.value();

				// discard a cell if its highest vertex is NOT the vertex, hence not part of lower star
				INDEX_TYPE vid1 = mLabel1->Cell2HighestVertex(ncell_GI);
				if (vid1 != vert_GI) continue; // not in lower star of f1

				// get the restriction label of th lower start cell
				char vid2;
				if (mLabel2 == NULL)
					vid2 = 0;
				else
					vid2 = mLabel2->GetLabel(ncell_GI);
				vid2 = vid2 * mMesh->maxDim() + mMesh->boundaryValue(ncell_GI);
				if (vid2 > LMAXDIM * LMAX_LABEL - 1) {
					printf("ERROR: not einough space for label: %d > %d * %d, maxDim=%d, lab=%d, bv=%d\n", vid2,
						LMAXDIM, LMAX_LABEL, mMesh->maxDim(), mLabel2->GetLabel(ncell_GI), mMesh->boundaryValue(ncell_GI));
					vid2 = LMAXDIM * LMAX_LABEL - 1;
				}
				// add the cell to the right subset
				//subsets[vid2].push_back(ncell_GI);
				subsets_map[vid2].push_back(ncell_GI);
			}

			// now do homotopy expand on each subset!
			//for (auto S : subsets) {
			for (int i = 0; i < LMAXDIM * LMAX_LABEL; i++) {
				if (subsets_map[i].size == 0) continue;

				// do homotopy expand
				myStaticMap cp;
				//HomotopyExpand(S.second, cp);
				HomotopyExpand(subsets_map[i], cp, ordering);

				for (int j = 0; j < cp.size; j++) {
					INDEX_TYPE id = cp.idarray[j];
					INDEX_TYPE pair = cp.cparray[j].pair;
					if (id == pair) {
						mGrad->setCritical(id, true);
						mGrad->setAssigned(id, 1);
					}
					else {
						mGrad->setPair(id, pair);
						mGrad->setPair(pair, id);
						mGrad->setAssigned(id, 1);
						mGrad->setAssigned(pair, 1);
					}
				}

			}

		}



		void ComputeLowerStar(INDEX_TYPE vert_GI) {
			// for each vertex compute the sets that will make up its lower start with restriction
			//std::map<char, std::vector<INDEX_TYPE> > subsets;

			small_INDEX_vector subsets_map[LMAXDIM * LMAX_LABEL];

			// add the vertex to restriction label set
			//subsets[mLabel2->GetLabel(vert_GI)].push_back(vert_GI);

			// now add all lower star to restriciton sets 
			typename MeshType::AdjacentCellsIterator star(mMesh);
			for (star.begin(vert_GI); star.valid(); star.advance()) {
				INDEX_TYPE ncell_GI = star.value();

				// discard a cell if its highest vertex is NOT the vertex, hence not part of lower star
				INDEX_TYPE vid1 = mLabel1->Cell2HighestVertex(ncell_GI);
				if (vid1 != vert_GI) continue; // not in lower star of f1

				// get the restriction label of th lower start cell
				char vid2;
				if (mLabel2 == NULL)
					vid2 = 0;
				else
					vid2 = mLabel2->GetLabel(ncell_GI);
				vid2 = vid2 * mMesh->maxDim() + mMesh->boundaryValue(ncell_GI);
				// add the cell to the right subset
				//subsets[vid2].push_back(ncell_GI);
				if (vid2 > LMAXDIM * LMAX_LABEL - 1) {
					printf("ERROR: not einough space for label: %d > %d * %d, maxDim=%d, lab=%d, bv=%d\n", vid2,
						LMAXDIM, LMAX_LABEL, mMesh->maxDim(), mLabel2->GetLabel(ncell_GI), mMesh->boundaryValue(ncell_GI));
					vid2 = LMAXDIM * LMAX_LABEL - 1;
				}
				subsets_map[vid2].push_back(ncell_GI);
			}

			// now do homotopy expand on each subset!
			//for (auto S : subsets) {
			for (int i = 0; i < LMAXDIM * LMAX_LABEL; i++) {
				if (subsets_map[i].size == 0) continue;

				// do homotopy expand
				myStaticMap cp;
				//HomotopyExpand(S.second, cp);
				HomotopyExpand(subsets_map[i], cp);

				for (int j = 0; j < cp.size; j++) {
					INDEX_TYPE id = cp.idarray[j];
					INDEX_TYPE pair = cp.cparray[j].pair;
					if (id == pair) {
						mGrad->setCritical(id, true);
						mGrad->setAssigned(id, 1);
					}
					else {
						mGrad->setPair(id, pair);
						mGrad->setPair(pair, id);
						mGrad->setAssigned(id, 1);
						mGrad->setAssigned(pair, 1);
					}
				}

			}

		}
		void ComputePairing() {
			std::vector<INDEX_TYPE> topo_index_partition;
			int num_threads;
#pragma omp parallel
			{
#pragma omp single
				{
					num_threads = omp_get_num_threads();
					//printf("MyRobinsNoAlloc::ComputePairing() --> num_threads = %d\n", num_threads);
					ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, topo_index_partition);
					//printf("MyRobinsNoAlloc::ComputePairing() --> splitvec numcells = %lld\n", mMesh->numCells());
					//for (auto id : topo_index_partition) {
					//	printf("%lld\n", id);
					//}
				}
#pragma omp barrier
				int thread_num = omp_get_thread_num();
				//printf("MyRobinsNoAlloc::ComputePairing() --> thread num %d doing work\n", thread_num);

				// iterate over all vertices
				typename MeshType::DCellsIterator verts(mMesh, 0, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
				for (verts.begin(); verts.valid(); verts.advance()) {
					INDEX_TYPE vert_GI = verts.value();

					ComputeLowerStar(vert_GI);


				}
				//printf("MyRobinsNoAlloc::ComputePairing() --> thread num %d done\n", thread_num);

			}

			//DenseLabeling<INDEX_TYPE>* GetLabeling() { return mPairs; }
/*
num_threads = 8;
printf("MyRobinsNoAlloc::ComputePairing() --> num_threads = %d\n", num_threads);
ArrayIndexPartitioner::EvenChunkSplit(mMesh->numCells(), num_threads, topo_index_partition);
printf("MyRobinsNoAlloc::ComputePairing() --> splitvec numcells = %lld\n", mMesh->numCells());
for (auto id : topo_index_partition) {
	printf("%lld\n", id);
}

for (int thread_num = 0; thread_num < num_threads; thread_num++) {
	//int thread_num = omp_get_thread_num();
	printf("MyRobinsNoAlloc::ComputePairing() --> thread num %d doing work\n", thread_num);

	// iterate over all vertices
	typename MeshType::DCellsIterator verts(mMesh, 0, topo_index_partition[thread_num], topo_index_partition[thread_num + 1]);
	for (verts.begin(); verts.valid(); verts.advance()) {
		INDEX_TYPE vert_GI = verts.value();

		ComputeLowerStar(vert_GI);


	}
	printf("MyRobinsNoAlloc::ComputePairing() --> thread num %d done\n", thread_num);
}
*/
		}


	

};





}
#endif
