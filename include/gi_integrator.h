/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

//#ifndef INTEGRATOR_H
//#define INTEGRATOR_H
//
//
//#include "gi_basic_types.h"
//#include "gi_vectors.h"
//#include "gi_labeling.h"
//#include "gi_regular_grid.h"
//#include "gi_regular_grid_trilinear_function.h"
//#include "gi_topological_regular_grid.h"
//#include "gi_discrete_gradient_labeling.h"
//#include "gi_topological_utility_functions.h"
//
//namespace GInt {
//	class Integrator {
//
//	protected:
//		DenseLabeling<INDEX_TYPE>* m_destinations;
//		RegularGrid3D* m_grid;
//		RegularGridTrilinearFunction* m_func;
//		TopologicalRegularGrid3D* m_mesh;
//		DiscreteGradientLabeling* m_dgrad;
//
//		Vec3i m_xyz;
//		Vec3b m_periodic;
//
//		char* m_fname;
//		int m_rkindex;
//
//	public:
//		Integrator(Vec3i xyz, Vec3b periodic, char* fname) : m_xyz(xyz), m_periodic(periodic), m_fname(fname), m_rkindex(1) {
//
//
//		}
//		void SetRKLevel(int level) {
//			m_rkindex = level;
//		}
//
//		void BeginIntegration() {
//
//			// SET UP CONTEXT FOR INTEGRATION
//			m_grid = new RegularGrid3D(m_xyz, m_periodic);
//			m_func = new RegularGridTrilinearFunction(m_grid);
//
//			m_func->LoadImageFromFile(m_fname);
//
//			m_func->ComputeGradFromImage(m_rkindex);
//
//			m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
//			const INDEX_TYPE t_num_vertices = m_grid->NumElements();
//
//			//m_mesh = new TopologicalRegularGrid3D(m_grid);
//
//			//m_dgrad = new DiscreteGradientLabeling(m_mesh);
//			//m_dgrad->ClearAllGradient();
//
//			// set all potential minima
//			std::vector<INDEX_TYPE> maxvec;
//#pragma omp parallel for 
//			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
//				if (IsLowestVertexIn6Neighborhood(i, m_func)) {
//					m_destinations->SetLabel(i, i);
//					//INDEX_TYPE t_vertex_cellid = m_mesh->CellIDFromVertexNumber(i);
//					//m_dgrad->setAssigned(t_vertex_cellid, true);
//					//m_dgrad->setCritical(t_vertex_cellid, true);
//#pragma omp critical 
//					{
//						maxvec.push_back(i);
//					}
//				}
//				else {
//					m_destinations->SetLabel(i, -1);
//				}
//			}
//
//
//
//		}
//
//
//
//	};
//	
//
//
//
//}
//#endif