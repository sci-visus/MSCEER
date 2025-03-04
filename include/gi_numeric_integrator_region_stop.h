#ifndef NUMERIC_INTEGRATOR_REGION_STOP_H
#define NUMERIC_INTEGRATOR_REGION_STOP_H


#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_labeling.h"
#include "gi_regular_grid.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_adaptive_euler_advector.h"



namespace GInt {
	template< class Advector, class Comparer>
	class NumericIntegratorRegionStop {

	protected:

		//struct bridge {
		//	INDEX_TYPE labela;
		//	INDEX_TYPE labelb;
		//	float value;
		//};
		//

		DenseLabeling<INDEX_TYPE>* m_destinations;
		DenseLabeling<int>* m_certains;
		RegularGrid3D* m_grid;
		RegularGridTrilinearFunction* m_func;


		Vec3i m_xyz;
		Vec3b m_periodic;
		FLOATTYPE m_error_threshold;
		FLOATTYPE m_gradient_threshold;
		int m_num_iterations_left;
		char* m_fname;

		//struct cmp{
		//	bool operator()(INDEX_TYPE a, INDEX_TYPE b) {
		//		return mCompare->Compare(a, b);
		//	}
		//};

		bool IsExtremeVertexIn6Neighborhood(INDEX_TYPE id) const {
			Vec3l t_neighbors[6];
			Vec3l t_coords = m_grid->Coords(id);
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
			int nn = m_grid->GatherExistingNeighborsSameBdry6(m_grid->Coords(tid), neighbors);
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
			Vec3l xyz = m_grid->Coords(startid);
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
					Enqueue_Later_Neighbors(m_grid->Coords(currid), growing_front, seen);
				}

			}
		}


		Comparer* mCompare;
		std::unordered_set<INDEX_TYPE> m_extrema;
	public:


		NumericIntegratorRegionStop(RegularGridTrilinearFunction* func, RegularGrid3D* grid, FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int iteration_limit) :
			m_num_iterations_left(iteration_limit), m_xyz(func->GetGrid()->XYZ()), m_periodic(func->GetGrid()->Periodic()), m_func(func), m_grid(grid),
			m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold) {
			mCompare = new Comparer(func);
		}

		DenseLabeling<INDEX_TYPE>* GetOutputLabels() { return m_destinations; }
		RegularGrid3D* GetGrid() { return m_grid; }
		RegularGridTrilinearFunction* GetFunction() { return m_func; }
		std::unordered_set<INDEX_TYPE>& GetExtrema() { return m_extrema; }
		void BeginIntegration() {

			printf("Doing numeric integration for volume assignment\n");

			//m_func->ComputeGradFromImage(m_rkindex);


			m_destinations = new DenseLabeling<INDEX_TYPE>(m_grid->NumElements());
			m_certains = new DenseLabeling<int>(m_grid->NumElements());
			const INDEX_TYPE t_num_vertices = m_grid->NumElements();
			AdvectionChecker* inside_voxel_critical_advection_checker = new TerminateNearCertain(m_certains, m_grid);
			AdvectionChecker* no_check = new NoTermination();//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);

			std::vector<INDEX_TYPE> extrema;
			printf("-- finding extrema to terminate integral lines\n");
			// set all potential extrema, so we terminate near them
#pragma omp parallel for 
			for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
				m_certains->SetLabel(i, -1);
				if (IsExtremeVertexIn6Neighborhood(i)) {
					m_destinations->SetLabel(i, i);
#pragma omp critical
					{
						extrema.push_back(i);
						m_extrema.insert(i);
					}
				}
				else {
					m_destinations->SetLabel(i, -1);
				}

			}
			printf("  -- found %d Extrema!\n-- Expanding Extrema certain regions\n", extrema.size());

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

			printf("  -- expansion done\n-- doing numerical integration\n");
			//m_destinations->OutputToIntFile("dest.raw");
			int t1, t2; t1 = t2 = 0;
#pragma omp parallel 
			{
				Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_critical_advection_checker);

#pragma omp for schedule(guided)  nowait 
				for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {
					// early skip if this is already a maximum
					if (m_certains->GetLabel(i) != -1) {
						continue;
					}
					Vec3l t_coords = m_grid->Coords(i); // get the coordinates of the poitn
					//if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
					Vec3d t_current_point = t_coords;
					int t_num_iterations_left = m_num_iterations_left;

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
						if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
							t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
							t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
							t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
							INDEX_TYPE t_dest_label = m_destinations->GetLabel(t_next_id);

							if (t_dest_label == -1) {
								m_destinations->SetLabel(i, t_next_id);
							}
							else {
								m_destinations->SetLabel(i, t_dest_label);
							}
							t_continue = false;
						}
					}
				}
				printf("-- Done numeric integration of 3d volume!\n");

			}
		}


	};


}
#endif
