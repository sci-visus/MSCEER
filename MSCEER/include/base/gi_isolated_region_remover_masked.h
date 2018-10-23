/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef ISOLATED_REGION_REMOVER_MASKED_H
#define ISOLATED_REGION_REMOVER_MASKED_H


#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_labeling.h"
#include "gi_regular_grid_3d.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_adaptive_euler_advector_3d.h"
#include "gi_timing.h"

#include <map>
#include <set>

namespace GInt {
    template< class Comparer>
    class IsolatedRegionRemoverMasked {

    protected:
        const DenseLabeling<int>* m_input;
        std::map<INDEX_TYPE, int> m_inversemask;

        DenseLabeling<int>* m_destinations_unmasked;

        DenseLabeling<INDEX_TYPE>* m_destinations;
        const RegularGrid3D* m_grid;
        RegularGridTrilinearFunction* m_func;

		bool ABeforeB(INDEX_TYPE a, INDEX_TYPE b) const {
			return mCompare->Compare(a, b);
		}
        // if the vertex is a max, return own index
        // if the vertex has a higher neighbor with same label, return that (or highest such)
        // if the vertex does not have higher with same label, return highest id
        INDEX_TYPE GetEarliestVertexIn6Neighborhood(INDEX_TYPE id) const {
            Vec3l t_neighbors[6];
            Vec3l t_coords = m_grid->XYZ3d(id);
            int t_num_neighbors = m_grid->GatherExistingNeighborsSameBdry6(t_coords, t_neighbors);

           INDEX_TYPE overall_earliest_id = id;
 
            for (int i = 0; i < t_num_neighbors; i++) {
                INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
				if (ABeforeB(t_neighbor_vertex, overall_earliest_id)) {
                    overall_earliest_id = t_neighbor_vertex;
                }
     //           if (m_input->GetLabel(id) == m_input->GetLabel(t_neighbor_vertex) &&
					//ABeforeB(t_neighbor_vertex, same_label_earliest_id)) {
     //               same_label_earliest_id = t_neighbor_vertex;
     //               has_samelabel_earlier = true;
     //           }
            }
            //if (is_local_region_earliest) return id;
           // if (has_samelabel_earlier) return same_label_earliest_id;
            return overall_earliest_id;
        }


        INDEX_TYPE PathCompressFind(INDEX_TYPE id) {
            if (m_destinations->GetLabel(id) == id) return id;
            INDEX_TYPE retval = PathCompressFind(m_destinations->GetLabel(id));
			m_destinations->SetLabel(id, retval);
//            INDEX_TYPE idref = m_destinations->operator[](id);
////#pragma omp atomic
//            m_destinations->operator[](id) += retval - idref; // this is stupid - openmp 2.0 supports only binops= for atomics
            return retval;
        }

		void UnionByBefore(INDEX_TYPE set_a, INDEX_TYPE set_b) {
			//if (set_a != set_b) printf("abeforeb = %d, %f %n", ABeforeB(set_a, set_b));
			if (ABeforeB(set_a, set_b)) {
				m_destinations->SetLabel(set_b, set_a);
			}
			else {
				m_destinations->SetLabel(set_a, set_b);
			}
		}



		INDEX_TYPE UnionNeighbors(INDEX_TYPE id)  {
			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->XYZ3d(id);
			int t_num_neighbors = m_grid->GatherExistingNeighborsSameBdry6(t_coords, t_neighbors);

			for (int i = 0; i < t_num_neighbors; i++) {
				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
		
				if (m_input->GetLabel(id) == m_input->GetLabel(t_neighbor_vertex)) {
					INDEX_TYPE set_a = PathCompressFind(id);
					INDEX_TYPE set_b = PathCompressFind(t_neighbor_vertex);
					if (set_a != set_b)	UnionByBefore(set_a, set_b);
				}
			}
			return 0;
		}



        Comparer* mCompare;

    public:

        IsolatedRegionRemoverMasked(RegularGridTrilinearFunction* func, const DenseLabeling<int> *input, const std::vector<INDEX_TYPE> *mask) :
            m_input(input), m_func(func) {
            m_grid = func->GetGrid();
            mCompare = new Comparer(func);

            for(int i = 0; i < mask->size(); i++) {
                m_inversemask.insert( std::pair<INDEX_TYPE, int> (mask->at(i), i) );
            }
        }

        ~IsolatedRegionRemoverMasked(){
            delete mCompare;
            delete m_destinations;
        }

        DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }
        DenseLabeling<int>* GetOutputLabelsUnmasked() { return m_destinations_unmasked; }
        const RegularGrid3D* GetGrid() { return m_grid; }
        RegularGridTrilinearFunction* GetFunction() { return m_func; }
        void ComputeOutput(bool verbose = false) {
			

            ThreadedTimer gtimer(1);
            gtimer.StartGlobal();

            if(verbose){
                printf(" -- Cleaning noise in integration...");
                fflush(stdout);
            }

            const INDEX_TYPE t_num_vertices = m_grid->NumElements();


            std::set<int> incount;
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				incount.insert(m_input->GetLabel(i));
			}
			printf("regioncleaner initially %d region labels\n");

            m_destinations = new DenseLabeling<INDEX_TYPE>(t_num_vertices);


            // set all potential extrema, so we terminate near them
#pragma omp parallel for
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                //INDEX_TYPE nid = GetEarliestVertexIn6Neighborhood(i);
                m_destinations->SetLabel(i, i);
            }

//#pragma omp parallel for

			// identify all connected components
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				UnionNeighbors(i);
            }

			// find all minima
            std::vector<INDEX_TYPE> lmins;
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				INDEX_TYPE set_id = PathCompressFind(i);

				if (i == set_id) {
					lmins.push_back(i);
				}
				
			}
			printf("found %d minima\n", lmins.size());



			for (auto id : lmins) {
				INDEX_TYPE lid = GetEarliestVertexIn6Neighborhood(id);
				if (lid != id) {
					UnionByBefore(id, PathCompressFind(lid));
					//if (PathCompressFind(id) == id) printf("ERROR:asdfasdfasdfasdf\n");
				}
			}

         std::set<INDEX_TYPE> newmins;
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				INDEX_TYPE set_id = PathCompressFind(i);
				newmins.insert(set_id);
				m_destinations->SetLabel(i, set_id);
			}
			printf("after merge = %d minima\n", newmins.size());

            ThreadedTimer ltimer(1);
            ltimer.StartGlobal();

            // create unmasked versions!
//            m_destinations_unmasked = new DenseLabeling<int>(t_num_vertices);
//
////#pragma omp parallel for
//            for(int i = 0; i < t_num_vertices; i++) {
//
//                if (m_input->GetLabel(i) < 0){
//                    m_destinations->SetLabel(i, m_input->GetLabel(i));
//                    continue;
//                }
//                m_destinations_unmasked->SetLabel(i, m_inversemask.at( m_destinations->GetLabel(i) ));
//            }

            ltimer.EndGlobal();
            gtimer.EndGlobal ();
            if (verbose){
                printf(" done!");   gtimer.PrintAll();
                //ltimer.PrintAll();
            }
        }
    };
#if 0
    class PathCompressor {

    protected:
        DenseLabeling<INDEX_TYPE>* m_input;
        DenseLabeling<INDEX_TYPE>* m_destinations;
        const RegularGrid3D* m_grid;




        // if the vertex is a max, return own index
        // if the vertex has a higher neighbor with same label, return that (or highest such)
        // if the vertex does not have higher with same label, return highest id

        INDEX_TYPE PathCompressFind(INDEX_TYPE id) {
            if (m_destinations->GetLabel(id) == id) return id;
            INDEX_TYPE retval = PathCompressFind(m_destinations->GetLabel(id));
//			INDEX_TYPE idref = m_destinations->operator[](id);
//#pragma omp atomic
//			m_destinations->operator[](id) += retval - idref; // this is stupid - openmp 2.0 supports only binops= for atomics
            m_destinations->SetLabel(id, retval);
            return retval;
        }

        //Comparer* mCompare;

    public:
        PathCompressor(RegularGridTrilinearFunction* func, DenseLabeling<INDEX_TYPE>* input) :
            m_input(input) {
            m_grid = func->GetGrid();
        }


        DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }
        const RegularGrid3D* GetGrid() { return m_grid; }
        void ComputeOutput() {


            m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
            const INDEX_TYPE t_num_vertices = m_grid->NumElements();
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                m_destinations->SetLabel(i, m_input->GetLabel(i));
            }

            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
                PathCompressFind(i);
            }

        }



    };

#endif

}
#endif
