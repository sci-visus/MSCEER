/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/


#ifndef UNION_FIND_LABELING_H
#define UNION_FIND_LABELING_H

#include "gi_basic_types.h"
#include "gi_labeling.h"
//#include "gi_experimental.h"
#include "gi_regular_grid_3d.h"


namespace GInt {

	class VolumeConnectedComponents {
	public:

		INDEX_TYPE Find(INDEX_TYPE id) {
			INDEX_TYPE tmp1 = mIDVol->GetLabel(id);
			if (tmp1 == id) return id;
			INDEX_TYPE tmp2 = Find(tmp1);
			mIDVol->SetLabel(id, tmp2);
			return tmp2;	
		}

		void Merge(INDEX_TYPE id1, INDEX_TYPE id2) {
			INDEX_TYPE fid1 = Find(id1);
			INDEX_TYPE fid2 = Find(id2);
			if (fid1 < fid2) {
				mIDVol->SetLabel(fid2, fid1);
			}
			else {
				mIDVol->SetLabel(fid1, fid2);
			}
		}

		void AddVoxel(INDEX_TYPE id) {
			mIDVol->SetLabel(id, id);
		}


		VolumeConnectedComponents(RegularGrid3D* grid) : mGrid(grid), mIDVol(NULL){}

		RegularGrid3D* mGrid;

		DenseLabeling<INDEX_TYPE>* mIDVol;

		void PerformUnionFind(DenseLabeling<char>* maskvol) {
			if (mIDVol != NULL) delete mIDVol;
			mIDVol = new DenseLabeling<INDEX_TYPE>(maskvol->GetNumLabels());
			mIDVol->SetAll(-1);
			for (INDEX_TYPE i = 0; i < maskvol->GetNumLabels(); i++) {
				if (maskvol->GetLabel(i) > 0) {
					AddVoxel(i);
				}
			}

			Vec3l negs[6];
			for (INDEX_TYPE i = 0; i < maskvol->GetNumLabels(); i++) {
				if (maskvol->GetLabel(i) == 0) continue;
				Vec3l coords = mGrid->XYZ3d(i);
				int count = mGrid->GatherExistingNeighborsSameBdry6(coords, negs);
				for (int j = 0; j < count; j++) {
					INDEX_TYPE nid = mGrid->Index3d(negs[j]);
					if (maskvol->GetLabel(nid) == 0) continue;
					Merge(nid, i);
				}
			}
			
			for (INDEX_TYPE i = 0; i < maskvol->GetNumLabels(); i++) {
				if (maskvol->GetLabel(i) == 0) continue;
				Find(i);
			}

		}









	};






	//class UnionFindLabeling : public DenseLabeling < INDEX_TYPE > {
	//protected:
	//	INDEX_TYPE* m_labels;
	//	INDEX_TYPE m_num_labels;
	//public:

	//	UnionFindLabeling(INDEX_TYPE num_labels) :DenseLabeling<INDEX_TYPE>(num_labels) {
	//	}

	//	void SetLabel(INDEX_TYPE id, INDEX_TYPE label) {
	//		m_labels[id] = label;
	//	}

	//	INDEX_TYPE GetLabel(INDEX_TYPE id)  {
	//		return Find(id);
	//	}

	//	INDEX_TYPE Find(INDEX_TYPE id) {
	//		INDEX_TYPE tlabel = this->GetLabel(id);
	//		if (this->GetLabel(id) == id) {
	//			return id;
	//		}
	//		// hopefully we don't have 1-cycles!!!!
	//		//else if (a[a[s]] == s) {
	//		//	return s;
	//		//}
	//		tlabel = Find(tlabel);
	//		this->SetLabel(id, tlabel);
	//		return tlabel;
	//	}

	//	void OrderedUnion(INDEX_TYPE keep, INDEX_TYPE merge) {
	//		SetLabel(merge, keep);
	//	}

	//};
}

#if 0    
	class BubbleSet {
    protected:
        TopoVertexGraph* m_graph;
        unordered_map<INDEX_TYPE, BYTE_TYPE> m_counts;

    public:

        BubbleSet(TopoVertexGraph* graph) : m_graph(graph) {
        }


        void insert(INDEX_TYPE id) {
            m_counts[id] = 0;
            TopoVertexGraph::neighbor_iterator nit(m_graph);
            for (nit.begin(id); nit.valid(); nit.advance()) {
                INDEX_TYPE other_id = nit.value();
                if (m_graph->Before(other_id, id)) {
                    if (m_counts[other_id] == 1) {
                        m_counts.erase(other_id);
                    }
                    else {
                        m_counts[other_id]--;
                    }
                }
                else {
                    m_counts[id]++;
                }


            }


        }

        // the reason this is a "bubbleset" is that
        // the only calls to contains on the boundary of the set work
        // like a bubble - the air dosn't know if its inside or outside the bubble
        bool contains(INDEX_TYPE id) {
            return m_counts.count(id) > 0;
        }

        void merge(BubbleSet& other) {

        }

        // define iterators

    };




    template< class Comparer>
    class FlatRegionExpansion {
    public:

        FlatRegionExpansion(RegularGridTrilinearFunction* func, RegularGrid3D* grid) :
            m_func(func), m_grid(grid) {
            mCompare = new Comparer(func);
        }



    protected:

        //struct bridge {
        //	INDEX_TYPE labela;
        //	INDEX_TYPE labelb;
        //	float value;
        //};
        //
        RegularGrid3D* m_grid;
        Comparer* mCompare;
        UnionFindLabeling* m_destinations;
        RegularGridTrilinearFunction* m_func;

        map<INDEX_TYPE, vector<INDEX_TYPE>> m_merge_map;

        bool IsExtremeVertexIn6Neighborhood(INDEX_TYPE id) const {
            Vec3l t_neighbors[6];
            Vec3l t_coords = m_grid->XYZ3d(id);
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
            int nn = m_grid->GatherExistingNeighborsSameBdry6(m_grid->XYZ3d(tid), neighbors);
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
            Vec3l xyz = m_grid->XYZ3d(startid);
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
                    Enqueue_Later_Neighbors(m_grid->XYZ3d(currid), growing_front, seen);
                }

            }
        }


        void CreateLabeling() {

            m_destinations = new UnionFindLabeling(m_grid->NumElements());
            m_destinations->SetAll(-1);
            const INDEX_TYPE t_num_vertices = m_grid->NumElements();
            std::vector<INDEX_TYPE> extrema;

            // set all potential extrema, so we terminate near them
#pragma omp parallel for
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                if (IsExtremeVertexIn6Neighborhood(i)) {
                    //m_destinations->SetLabel(i, i);
                    vector<INDEX_TYPE> myext;
                    myext.push_back(i);
#pragma omp critical
                    {
                        extrema.push_back(i);
                        m_merge_map[i] = myext;
                    }
                }
                else {
                    //m_destinations->SetLabel(i, -1);
                }

            }

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
        }


    };
}

#endif
#endif
