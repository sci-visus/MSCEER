/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef ISOLATED_REGION_REMOVER_H
#define ISOLATED_REGION_REMOVER_H


#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_labeling.h"
#include "gi_regular_grid_3d.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_adaptive_euler_advector_3d.h"



namespace GInt {

	template< class Comparer, class GridFuncType>
	class IsolatedCCRegionRemoverNEW {

	protected:
		DenseLabeling<INT_TYPE>* m_input;
		DenseLabeling<INDEX_TYPE>* m_destinations;
		const RegularGrid3D* m_grid;
		GridFuncType* m_func;
		Comparer* mCompare;


		inline bool AComesBeforeB(INDEX_TYPE a, INDEX_TYPE b) const {
			return mCompare->Compare(a, b);
		}

		void MergeByOrder(INDEX_TYPE a_rep, INDEX_TYPE b_rep) {
			// merge by rank
			if (AComesBeforeB(a_rep, b_rep)) {
				m_destinations->SetLabel(b_rep, a_rep);
			}
			else {
				m_destinations->SetLabel(a_rep, b_rep);
			}
		}
		void MergeAToB(INDEX_TYPE a_rep, INDEX_TYPE b_rep) {
				m_destinations->SetLabel(a_rep, b_rep);
		}
		void MergeWithSameNeighbors(INDEX_TYPE id) {
			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->XYZ3d(id);
			int t_num_neighbors = m_grid->GatherExistingNeighborsAll6(t_coords, t_neighbors);
			

			for (int i = 0; i < t_num_neighbors; i++) {
				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);

				if (m_input->GetLabel(id) == m_input->GetLabel(t_neighbor_vertex)) {
					INDEX_TYPE id_rep = PathCompressFind(id);
					INDEX_TYPE neg_rep = PathCompressFind(t_neighbor_vertex);
					MergeByOrder(id_rep, neg_rep);
				}
			}
		}

		int WhichRegion(INDEX_TYPE id) {
			return m_input->GetLabel(PathCompressFind(id));
		}

		// so id is the earliest vertex for some contiguous region
		// now we check if it is interior
		// if not, we check if there is an earlier in a different region
		// if there is, we merge the regions
		void MergeBoundaryExtrema(INDEX_TYPE id)  {
			
			if (PathCompressFind(id) != id) printf("WHOATHERE pathcompressfind got non extrema\n");

			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->XYZ3d(id);
			int t_num_neighbors = m_grid->GatherExistingNeighborsAll6(t_coords, t_neighbors);

			int id_reg = WhichRegion(id);

			bool original_vertex_is_extremum = true;
			INDEX_TYPE t_overall_earliest = id;

			for (int i = 0; i < t_num_neighbors; i++) {
				original_vertex_is_extremum = false;
				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
				if (AComesBeforeB(t_neighbor_vertex, id)) {
					if (WhichRegion(id) == WhichRegion(t_neighbor_vertex)) continue;// this is case when UF walked uphill - but is fine printf("as;dlkfjasd;lfkjas;dl\n");
					//INDEX_TYPE neg_rep = PathCompressFind(t_neighbor_vertex);
					if (AComesBeforeB(t_neighbor_vertex, t_overall_earliest)) {
						t_overall_earliest = t_neighbor_vertex;
					}
				}
			}
			if (original_vertex_is_extremum) return; // we got nothing to do
			// else merge the regions!
			mNumRegionsRemoved++;
			this->MergeAToB(id, PathCompressFind(t_overall_earliest));
		}

		// if the vertex is a max, return own index
		// if the vertex has a higher neighbor with same label, return that (or highest such)
		// if the vertex does not have higher with same label, return highest id
		INDEX_TYPE GetEarliestVertexIn6Neighborhood(INDEX_TYPE id) const {
			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->XYZ3d(id);
			int t_num_neighbors = m_grid->GatherExistingNeighborsAll6(t_coords, t_neighbors);

			bool has_samelabel_earlier = false;
			bool original_vertex_is_extremum = true;
			INDEX_TYPE t_overall_earliest = id;
			INDEX_TYPE t_earliest_same_label = id;

			for (int i = 0; i < t_num_neighbors; i++) {
				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
				if (AComesBeforeB(t_neighbor_vertex, t_overall_earliest)) {
					t_overall_earliest = t_neighbor_vertex;
					original_vertex_is_extremum = false;
				}
				if (m_input->GetLabel(id) == m_input->GetLabel(t_neighbor_vertex) &&
					AComesBeforeB(t_neighbor_vertex, t_earliest_same_label)) {
					t_earliest_same_label = t_neighbor_vertex;
					has_samelabel_earlier = true;
				}
			}
			if (original_vertex_is_extremum) return id;
			if (has_samelabel_earlier) return t_earliest_same_label;
			return t_overall_earliest;
		}


		INDEX_TYPE PathCompressFind(INDEX_TYPE id) {
			if (m_destinations->GetLabel(id) == id) return id;
			INDEX_TYPE retval = PathCompressFind(m_destinations->GetLabel(id));
			m_destinations->SetLabel(id, retval);
			return retval;
		}

		INDEX_TYPE Find(INDEX_TYPE id) const {
			if (m_destinations->GetLabel(id) == id) return id;
			return Find(m_destinations->GetLabel(id));
		}
		int mNumRegionsRemoved;
	public:
		IsolatedCCRegionRemoverNEW(GridFuncType* func, DenseLabeling<INT_TYPE>* input) :
			m_input(input), m_func(func) {
			m_grid = func->GetGrid();
			mCompare = new Comparer(func);
		}
		~IsolatedCCRegionRemoverNEW() {
			delete mCompare;
			delete m_destinations;
		}

		DenseLabeling<INT_TYPE>* GetOutputLabels() { return m_input; }
		const RegularGrid3D* GetGrid() { return m_grid; }
		GridFuncType* GetFunction() { return m_func; }
		void ComputeOutput() {


			m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
			const INDEX_TYPE t_num_vertices = m_grid->NumElements();
			printf("fart\n");

#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				m_destinations->SetLabel(i, i);
			}
//#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				MergeWithSameNeighbors(i);
			}
			std::vector<INDEX_TYPE> tomergeids;
#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				if (Find(i) == i) {
#pragma omp critical
					{
						tomergeids.push_back(i);
					}
				}
			}
			mNumRegionsRemoved = 0;
			for (INDEX_TYPE id : tomergeids) {
				MergeBoundaryExtrema(id);
			}
			// set all potential extrema, so we terminate near them
#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				m_input->SetLabel(i, m_input->GetLabel(Find(i)));
			}

			printf("  -- Removed %d extra regions\n", mNumRegionsRemoved);
		}



	};



	template< class Comparer>
	class IsolatedRegionRemoverNEW {

	protected:
		DenseLabeling<INT_TYPE>* m_input;
		DenseLabeling<INDEX_TYPE>* m_destinations;
		const RegularGrid3D* m_grid;
		RegularGridTrilinearFunction* m_func;
		Comparer* mCompare;


		inline bool AComesBeforeB(INDEX_TYPE a, INDEX_TYPE b) const {
			return mCompare->Compare(a, b);
		}

		// if the vertex is a max, return own index
		// if the vertex has a higher neighbor with same label, return that (or highest such)
		// if the vertex does not have higher with same label, return highest id
		INDEX_TYPE GetEarliestVertexIn6Neighborhood(INDEX_TYPE id) const {
			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->XYZ3d(id);
			int t_num_neighbors = m_grid->GatherExistingNeighborsAll6(t_coords, t_neighbors);

			bool has_samelabel_earlier = false;
			bool original_vertex_is_extremum = true;
			INDEX_TYPE t_overall_earliest = id;
			INDEX_TYPE t_earliest_same_label = id;

			for (int i = 0; i < t_num_neighbors; i++) {
				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
				if (AComesBeforeB(t_neighbor_vertex, t_overall_earliest)) {
					t_overall_earliest = t_neighbor_vertex;
					original_vertex_is_extremum = false;
				}
				if (m_input->GetLabel(id) == m_input->GetLabel(t_neighbor_vertex) &&
					AComesBeforeB(t_neighbor_vertex, t_earliest_same_label)) {
					t_earliest_same_label = t_neighbor_vertex;
					has_samelabel_earlier = true;
				}
			}
			if (original_vertex_is_extremum) return id;
			if (has_samelabel_earlier) return t_earliest_same_label;
			return t_overall_earliest;
		}


		INDEX_TYPE PathCompressFind(INDEX_TYPE id) {
			if (m_destinations->GetLabel(id) == id) return id;
			INDEX_TYPE retval = PathCompressFind(m_destinations->GetLabel(id));
			m_destinations->SetLabel(id, retval); 
			return retval;
		}


	public:
		IsolatedRegionRemoverNEW(RegularGridTrilinearFunction* func, DenseLabeling<INT_TYPE>* input) :
			m_input(input), m_func(func) {
			m_grid = func->GetGrid();
			mCompare = new Comparer(func);
		}
		~IsolatedRegionRemoverNEW(){
			delete mCompare;
			delete m_destinations;
		}

		DenseLabeling<INT_TYPE>* GetOutputLabels() { return m_input; }
		const RegularGrid3D* GetGrid() { return m_grid; }
		RegularGridTrilinearFunction* GetFunction() { return m_func; }
		void ComputeOutput() {


			m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
			const INDEX_TYPE t_num_vertices = m_grid->NumElements();
			printf("fart\n");
			// set all potential extrema, so we terminate near them
#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				INDEX_TYPE nid = GetEarliestVertexIn6Neighborhood(i);
				m_destinations->SetLabel(i, nid);
			}
			printf("here1a\n");
			//#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				m_input->SetLabel(i, m_input->GetLabel(PathCompressFind(i)));
			}

		}



	};




	template< class Comparer>
	class IsolatedRegionRemover {

	protected:
		DenseLabeling<INDEX_TYPE>* m_input;
		DenseLabeling<INDEX_TYPE>* m_destinations;
		const RegularGrid3D* m_grid;
		RegularGridTrilinearFunction* m_func;

	
		inline bool AComesBeforeB(INDEX_TYPE a, INDEX_TYPE b) const {
			return mCompare->Compare(a, b);
		}

		// if the vertex is a max, return own index
		// if the vertex has a higher neighbor with same label, return that (or highest such)
		// if the vertex does not have higher with same label, return highest id
		INDEX_TYPE GetEarliestVertexIn6Neighborhood(INDEX_TYPE id) const {
			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->XYZ3d(id);
			int t_num_neighbors = m_grid->GatherExistingNeighborsSameBdry6(t_coords, t_neighbors);

			bool has_samelabel_earlier = false;
			bool original_vertex_is_extremum = true;
			INDEX_TYPE t_overall_earliest = id;
			INDEX_TYPE t_earliest_same_label = id;

			for (int i = 0; i < t_num_neighbors; i++) {
				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
				if (AComesBeforeB(t_neighbor_vertex, t_overall_earliest)) {
					t_overall_earliest = t_neighbor_vertex;
					original_vertex_is_extremum = false;
				}
				if (m_input->GetLabel(id) == m_input->GetLabel(t_neighbor_vertex) &&
					AComesBeforeB(t_neighbor_vertex, t_earliest_same_label)) {
					t_earliest_same_label = t_neighbor_vertex;
					has_samelabel_earlier = true;
				}
			}
			if (original_vertex_is_extremum) return id;
			if (has_samelabel_earlier) return t_earliest_same_label;
			return t_overall_earliest;
		}


		INDEX_TYPE PathCompressFind(INDEX_TYPE id) {
			if (m_destinations->GetLabel(id) == id) return id;
			INDEX_TYPE retval = PathCompressFind(m_destinations->GetLabel(id));
			INDEX_TYPE idref = m_destinations->operator[](id);
#pragma omp atomic
			m_destinations->operator[](id) += retval - idref; // this is stupid - openmp 2.0 supports only binops= for atomics
			return retval;
		}

		Comparer* mCompare;

	public:
		IsolatedRegionRemover(RegularGridTrilinearFunction* func, DenseLabeling<INDEX_TYPE>* input) :
			m_input(input), m_func(func) {
			m_grid = func->GetGrid();
			mCompare = new Comparer(func);
		}
		~IsolatedRegionRemover(){
			delete mCompare;
			delete m_destinations;
		}

		DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }
		const RegularGrid3D* GetGrid() { return m_grid; }
		RegularGridTrilinearFunction* GetFunction() { return m_func; }
		void ComputeOutput() {

	
			m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
			const INDEX_TYPE t_num_vertices = m_grid->NumElements();

			// set all potential extrema, so we terminate near them
#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				INDEX_TYPE nid = GetEarliestVertexIn6Neighborhood(i);
				m_destinations->SetLabel(i, nid);
			}
//#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				PathCompressFind(i);
			}

	}
	


	};
	
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



}
#endif