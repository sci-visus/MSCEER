/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <map>
#include <set>
#include <queue>
#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_labeling.h"
#include "gi_regular_grid.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_adaptive_euler_advector.h"



namespace GInt {
	template <class GridFuncType>
	class StrictlyNumericIntegrator {

	protected:
		DenseLabeling<INDEX_TYPE>* m_destinations;
		RegularGrid3D* m_grid;
		GridFuncType* m_func;


		Vec3i m_xyz;
		Vec3b m_periodic;
		FLOATTYPE m_error_threshold;
		FLOATTYPE m_gradient_threshold;

		char* m_fname;
		int m_rkindex;

		bool IsHighestVertexIn6Neighborhood(INDEX_TYPE id) const {
			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->Coords(id);
			int t_num_neighbors = m_grid->GatherExistingNeighborsSameBdry6(t_coords, t_neighbors);

			INDEX_TYPE t_current_lowest = id;
			for (int i = 0; i < t_num_neighbors; i++) {
				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
				if (m_func->IsGreater(t_neighbor_vertex, t_current_lowest)) {
					return false;
				}
			}
			return true;
		}


	public:
		StrictlyNumericIntegrator(Vec3i xyz, Vec3b periodic, char* fname, FLOATTYPE error_threshold, FLOATTYPE gradient_threshold) :
			m_xyz(xyz), m_periodic(periodic), m_fname(fname), m_rkindex(1),
			m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold) {
		}
		void SetRKLevel(int level) {
			m_rkindex = level;
		}

		DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }
		RegularGrid3D* GetGrid() { return m_grid; }
		GridFuncType* GetFunction() { return m_func; }

		void BeginIntegration() {

			// SET UP CONTEXT FOR INTEGRATION
			m_grid = new RegularGrid3D(m_xyz, m_periodic);
			m_func = new GridFuncType(m_grid);
			m_func->LoadImageFromFile(m_fname);

			m_func->ComputeGradFromImage(m_rkindex);

			m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
			const INDEX_TYPE t_num_vertices = m_grid->NumElements();
			AdvectionChecker* inside_voxel_advection_checker = new TerminateNearCrits(m_destinations, m_grid);
			//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);

			// set all potential maxima, so we terminate near them
#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				if (IsHighestVertexIn6Neighborhood(i)) {
					m_destinations->SetLabel(i, i);
				}
				else {
					m_destinations->SetLabel(i, -1);
				}
			}

			int t1, t2; t1 = t2 = 0;
#pragma omp parallel 
	{
				AdaptiveEulerAdvector3D<GridFuncType, 1> t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_advection_checker);

#pragma omp for schedule(dynamic)  nowait 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				// early skip if this is already a maximum
				if (m_destinations->GetLabel(i) == i) {
					continue;
				}

				Vec3l t_coords = m_grid->Coords(i); // get the coordinates of the poitn
				if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
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
					// if we terminated or hit a critical point, then we are done
					if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT||
						t_return_code == ADVECTION_EVENT::HIT_EXTREMUM || 
						t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
						t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
						m_destinations->SetLabel(i, t_next_id);
						t_continue = false;
					}
				}


				//if (m_grid->DistToBoundary()) {}
				//if (IsLowestVertexIn6Neighborhood(i)) {
				//	m_destinations->SetLabel(i, i);
				//}
				//else {
				//	m_destinations->SetLabel(i, -1);
				//}
			}


			printf("%d boundary, %d noboundary\n", t1, t2);
			}

	}
	


	};
	



}
#endif
