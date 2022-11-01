/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

//#ifndef NUMERIC_INTEGRATOR_DOWNWARDS_SINGLE_LINE_H
//#define NUMERIC_INTEGRATOR_DOWNWARDS_SINGLE_LINE_H
//
//#include <stdio.h>
//#include <set>
//#include <queue>
//#include "gi_basic_types.h"
//#include "gi_vectors.h"
//#include "gi_labeling.h"
//#include "gi_regular_grid.h"
//#include "gi_regular_grid_trilinear_function.h"
//#include "gi_adaptive_euler_advector.h"
//
//
//
//namespace GInt {
//	class NumericIntegratorDownwardSingleLine {
//
//	protected:
//		
//		//struct bridge {
//		//	INDEX_TYPE labela;
//		//	INDEX_TYPE labelb;
//		//	float value;
//		//};
//		//
//		
//		DenseLabeling<INDEX_TYPE>* m_destinations;
//		DenseLabeling<int>* m_certains;
//		RegularGrid3D* m_grid;
//		RegularGridTrilinearFunction* m_func;
//
//
//		Vec3i m_xyz;
//		Vec3b m_periodic;
//		FLOATTYPE m_error_threshold;
//		FLOATTYPE m_gradient_threshold;
//
//		char* m_fname;
//
//		bool IsLowestVertexIn6Neighborhood(INDEX_TYPE id) const {
//			Vec3l t_neighbors[6];
//			Vec3l t_coords = m_grid->XYZ3d(id);
//			int t_num_neighbors = m_grid->GatherExistingNeighborsSameBdry6(t_coords, t_neighbors);
//
//			INDEX_TYPE t_current_lowest = id;
//			for (int i = 0; i < t_num_neighbors; i++) {
//				INDEX_TYPE t_neighbor_vertex = m_grid->Index3d(t_neighbors[i]);
//				if (m_func->IsGreater(t_current_lowest, t_neighbor_vertex)) {
//					return false;
//				}
//			}
//			return true;
//		}
//		//void Enqueue_Higher_Neighbors(Vec3l xyz, std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > &expansion, std::set<INDEX_TYPE>&seen) {
//		//	INDEX_TYPE tid = m_grid->Index3d(xyz);
//
//		//	Vec3l neighbors[6];
//		//	int nn = m_grid->GatherExistingNeighborsSameBdry6(xyz, neighbors);
//		//	for (int i = 0; i < nn; i++) {
//		//		INDEX_TYPE tneg = m_grid->Index3d(neighbors[i]);
//		//		if (m_certains->GetLabel(tneg) == -1 && m_func->IsGreater(tid, tneg) && seen.count(tneg) == 0) { 
//		//			seen.insert(tneg); 
//		//			expansion.push(std::pair<FLOATTYPE, INDEX_TYPE>(m_func->SampleImage(tneg), tneg)); 
//		//		}
//
//		//	}
//		//}
//
//		//int Inspect_Higher_Certains(INDEX_TYPE tid) {
//		//	INDEX_TYPE tneg;
//		//	int higher_certain = m_certains->GetLabel(tid);
//		//	bool has_higher = false;
//
//		//	Vec3l neighbors[6];
//		//	int nn = m_grid->GatherExistingNeighborsSameBdry6(m_grid->XYZ3d(tid), neighbors);
//		//	for (int i = 0; i < nn; i++) {
//		//		INDEX_TYPE tneg = m_grid->Index3d(neighbors[i]);
//		//		if (m_func->IsGreater(tneg, tid)) {
//		//			if (m_certains->GetLabel(tneg) < 0) return -1; // if a higher one is uncertain, we are uncertain
//		//			if (!has_higher) {
//		//				higher_certain = m_certains->GetLabel(tneg);
//		//				has_higher = true;
//		//			}
//		//			else {
//		//				if (higher_certain != m_certains->GetLabel(tneg)) return -1;
//		//			}
//		//		}
//		//	}
//
//		//	if (!has_higher) {
//		//		printf("ERROR should never get here\n");
//		//		return -1;
//		//	}
//		//	return higher_certain;
//
//		//}
//
//		//void Expand_Lower_Neighborhood(INDEX_TYPE startid) {
//		//	Vec3l xyz = m_grid->XYZ3d(startid);
//		//	std::set<INDEX_TYPE> seen;
//
//		//	INDEX_TYPE tid = startid;
//		//	// the natural ordering using the < operator on pairs will give us the highest
//		//	// element first, simulating region growing from high to low
//		//	std::priority_queue<std::pair<FLOATTYPE, INDEX_TYPE> > growing_front;
//		//	seen.insert(startid);
//		//	Enqueue_Lower_Neighbors(xyz, growing_front, seen);
//
//		//	while (!growing_front.empty()) {
//
//		//		INDEX_TYPE currid = growing_front.top().second;
//		//		growing_front.pop();
//
//		//		int cellvale = Inspect_Higher_Certains(currid);
//		//		// find highers
//
//		//		// cellvalue >=0 indicates that there is certainty here, so lets expand 
//		//		if (cellvale >= 0) {
//		//			m_certains->SetLabel(currid, cellvale);
//		//			m_destinations->SetLabel(currid, startid);
//		//			Enqueue_Lower_Neighbors(m_grid->XYZ3d(currid), growing_front, seen);
//		//		}
//
//		//	}
//		//}
//
//
//
//	public:
//
//
//		NumericIntegratorDownwardSingleLine(RegularGridTrilinearFunction* func, RegularGrid3D* grid, FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int interation_limit) :
//			m_xyz(func->GetGrid()->XYZ()), m_periodic(func->GetGrid()->Periodic()), m_func(func), m_grid(grid),
//			m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold) {
//		}
//
//		DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }
//		RegularGrid3D* GetGrid() { return m_grid; }
//		RegularGridTrilinearFunction* GetFunction() { return m_func; }
//		void BeginIntegration() {
//
//	
//			m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
//			m_certains = new DenseLabeling<int>(m_grid->NumElements());
//			const INDEX_TYPE t_num_vertices = m_grid->NumElements();
//			AdvectionChecker* inside_voxel_advection_checker = new TerminateNearCertain(m_certains, m_grid);
//			//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);
//
//			std::vector<INDEX_TYPE> minima;
//
//			// set all potential maxima, so we terminate near them
//#pragma omp parallel for 
//			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
//				m_certains->SetLabel(i, -1);
//				if (IsLowestVertexIn6Neighborhood(i)) {
//					m_destinations->SetLabel(i, i);
//#pragma omp critical
//					{
//						minima.push_back(i);
//					}
//				}
//				else {
//					m_destinations->SetLabel(i, -1);
//				}
//
//			}
//			printf("found %d maxima!\nExpanding maxima certain regions\n", maxima.size());
//
//			int num_maxima = maxima.size();
//#pragma omp parallel shared(maxima) 
//			{
//#pragma omp for schedule(dynamic)  nowait 
//				for (int m = 0; m < num_maxima; m++) {
//					INDEX_TYPE maximum = maxima[m];
//					m_certains->SetLabel(maximum, m);
//					Expand_Lower_Neighborhood(maximum);
//				}
//			}
//
//
//			int t1, t2; t1 = t2 = 0;
//#pragma omp parallel 
//	{
//			AdaptiveEulerAdvector3D t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_advection_checker);
//			std::vector<INDEX_TYPE> t_path;
//			t_path.reserve(100);
//#pragma omp for schedule(guided)  nowait 
//			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
//				// early skip if this is already a maximum
//				if (m_certains->GetLabel(i) != -1 || m_destinations->GetLabel(i) != -1) {
//					continue;
//				}
//
//				t_path.clear();
//				t_path.push_back(i);
//				Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
//				if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
//				Vec3d t_current_point = t_coords;
//				int t_num_iterations_left = 500;
//				bool t_continue = true;
//
//				while (t_continue) {
//					Vec3d t_next_point;
//					ADVECTION_EVENT t_return_code;
//					if (m_grid->DistToBoundary(t_coords) <= 1) {
//						t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
//						t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
//						t1++;
//					}
//					else {
//						t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
//						t_coords = (t_current_point + 0.5);
//						t2++;
//					}
//					INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
//					t_path.push_back(t_next_id);
//					// if we terminated or hit a critical point, then we are done
//					if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
//						t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
//						t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
//						t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
//						INDEX_TYPE t_dest_label = m_destinations->GetLabel(t_next_id);
//						int t_certain_label = m_certains->GetLabel(t_next_id);
//						for (int j = 0; j < t_path.size(); j++) {
//							m_destinations->SetLabel(t_path[j], t_dest_label);
//							//m_certains->SetLabel(t_path[j], t_certain_label);
//						}
//						t_continue = false;
//					}
//				}
//			}
//#pragma omp for schedule(static)  
//			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
//				// early skip if this is already a maximum
//				if (m_certains->GetLabel(i) != -1) {
//					continue;
//				}
//
//				Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
//				if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
//				Vec3l negs[6];
//
//				int nn = m_grid->GatherExistingNeighborsSameBdry6(t_coords, negs);
//
//				bool needsupdate = false;
//				INDEX_TYPE t_dest_label = m_destinations->GetLabel(i);
//				for (int j = 0; j < nn; j++)  {
//					if (m_destinations->GetLabel(m_grid->Index3d(negs[j])) != t_dest_label) {
//						needsupdate = true;
//						break;
//					}
//				}
//
//				if (needsupdate) continue;
//
//				m_certains->SetLabel(i, 1);
//			}
//#pragma omp for schedule(guided)  nowait 
//				for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
//					// early skip if this is already a maximum
//					if (m_certains->GetLabel(i) != -1) {
//						continue;
//					}
//
//					Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
//					if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
//
//
//
//					Vec3d t_current_point = t_coords;
//					int t_num_iterations_left = 500;
//					bool t_continue = true;
//
//					while (t_continue) {
//						Vec3d t_next_point;
//						ADVECTION_EVENT t_return_code;
//						if (m_grid->DistToBoundary(t_coords) <= 1) {
//							t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
//							t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
//							t1++;
//						}
//						else {
//							t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
//							t_coords = (t_current_point + 0.5);
//							t2++;
//						}
//						INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
//						t_path.push_back(t_next_id);
//						// if we terminated or hit a critical point, then we are done
//						if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
//							t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
//							t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
//							t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
//							m_destinations->SetLabel(i, m_destinations->GetLabel(t_next_id));
//							t_continue = false;
//						}
//					}
//				}
//
//				//if (m_grid->DistToBoundary()) {}
//				//if (IsLowestVertexIn6Neighborhood(i)) {
//				//	m_destinations->SetLabel(i, i);
//				//}
//				//else {
//				//	m_destinations->SetLabel(i, -1);
//				//}
//			
//
//
//			printf("%d boundary, %d noboundary\n", t1, t2);
//			}
//
//	}
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