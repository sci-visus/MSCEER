#ifndef GI_CLS_MAP_GRADIENT_H
#define GI_CLS_MAP_GRADIENT_H

// Lower-star gradient with std::map bookkeeping: the earlier, allocating form of the
// rule in gi_cls_gradient.h. Kept for the tools that still use it; it passes the
// restriction label, not the cell id, to the mesh boundary query.

#include "gi_labeling.h"
#include "gi_discrete_gradient_labeling.h"
#include "gi_topological_simplicial_complex.h"
using namespace std;
namespace GInt {

	template<class MeshType, class MaxVLType, class GradType>
	class CLSMapGradient
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
		CLSMapGradient(MeshType* mesh, MaxVLType* label1, DenseLabeling<char>* label2,
			GradType* grad) : mMesh(mesh), mLabel1(label1), mLabel2(label2), mGrad(grad) {
		}
		CLSMapGradient(MeshType* mesh, MaxVLType* label1, GradType* grad) :
			mMesh(mesh), mLabel1(label1), mLabel2(NULL), mGrad(grad) {
		}
		~CLSMapGradient() {
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
				char vid2;
				if (mLabel2 == NULL) {
					vid2 = 0;
				}
				else {
					vid2 = mLabel2->GetLabel(ncell_GI);
				}
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

	// Previous name of CLSMapGradient.
	template<class MeshType, class MaxVLType, class GradType>
	using MyRobins = CLSMapGradient<MeshType, MaxVLType, GradType>;


}
#endif
