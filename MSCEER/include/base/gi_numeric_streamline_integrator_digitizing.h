/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef GI_NUMERIC_STREAMLINE_INTEGRATOR_DIGITIZING_H
#define GI_NUMERIC_STREAMLINE_INTEGRATOR_DIGITIZING_H

#include <set>
#include <queue>
#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_labeling.h"
#include "gi_regular_grid_3d.h"
#include "gi_regular_grid_2d.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_regular_grid_bilinear_function.h"
#include "gi_adaptive_euler_advector_3d.h"
#include "gi_adaptive_euler_advector_2d.h"
//#include "gi_adaptive_euler_advector_reverse.h"

namespace GInt {

	//template< class Advector>
	//class NumericStreamlineIntegrator {


	//	FLOATTYPE m_error_threshold;
	//	FLOATTYPE m_gradient_threshold;
	//	int m_max_steps;

	//	RegularGrid3D* m_grid;
	//	RegularGridTrilinearFunction* m_func;
	//	bool in_bounds(const Vec3d &x) const {
	//		return (0.0 <= x[0] && x[0] < 1.0 &&
	//			0.0 <= x[1] && x[1] < 1.0 &&
	//			0.0 <= x[2] && x[2] < 1.0);
	//	}

	//	AdvectionChecker* m_checker;
	//public:


	//	void SetAdvectionChecker(AdvectionChecker* checker) { m_checker = checker; }
	//	NumericStreamlineIntegrator(RegularGrid3D* grid, RegularGridTrilinearFunction* func,
	//		FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int max_steps) :
	//		m_grid(grid), m_func(func),
	//		m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold), m_max_steps(max_steps) {
	//		m_checker = new NoTermination();
	//	}


	//	GInt::ADVECTION_EVENT IntegrateStreamline(const Vec3d &seed, std::vector<Vec3d> &sline) const {

	//		Vec3l start_coords = (seed + 0.5); // since casting rounds down to integer... this is starting vertex
	//		INDEX_TYPE start_id = m_grid->Index3d(start_coords); // vertex id of start point
	//		int t1, t2; t1 = t2 = 0;
	//		Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, m_checker);
	//		sline.reserve(100);

	//		sline.push_back(seed);

	//		Vec3l t_coords = start_coords; // get the coordinates of the poitn
	//		//if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
	//		Vec3d t_current_point = seed;
	//		int t_num_iterations_left = m_max_steps;
	//		bool t_continue = true;

	//		while (t_continue) {
	//			ADVECTION_EVENT t_return_code;
	//			if (m_grid->DistToBoundary(t_coords) <= 1) {
	//				t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
	//				t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
	//				t1++;
	//			}
	//			else {
	//				t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
	//				t_coords = (t_current_point + 0.5);
	//				t2++;
	//			}


	//			sline.push_back(t_current_point);
	//			//INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);

	//			// if we terminated or hit a critical point, then we are done
	//			if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
	//				t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
	//				t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
	//				t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {


	//				return t_return_code;
	//				t_continue = false;
	//			}
	//		}
	//		return ADVECTION_EVENT::NONE;
	//	}



	//};

	template<class MeshType,  class GridFuncType, class Advector>
	class DigitizingNumericStreamlineIntegrator3dASC {


		FLOATTYPE m_error_threshold;
		FLOATTYPE m_gradient_threshold;
		int m_max_steps;

		RegularGrid3D* m_grid;
		GridFuncType* m_func;
		MeshType* m_mesh;
		DenseLabeling<char>* m_labeling;
		bool in_bounds(const Vec3d &x) const {
			return (0.0 <= x[0] && x[0] < 1.0 &&
				0.0 <= x[1] && x[1] < 1.0 &&
				0.0 <= x[2] && x[2] < 1.0);
		}

		AdvectionChecker* m_checker; 
		int m_target_label;
	public:

		bool already_labeled(INDEX_TYPE id) const { return m_labeling->GetLabel(id) == m_target_label; }
		void set_label(INDEX_TYPE id) { 
			//if (m_labeling->GetLabel(id) == 0) {
			//	printf("setting background to m_target_label, instead of manifold\n");
			//}
			m_labeling->SetLabel(id, m_target_label);
		}

		DenseLabeling<char>* get_output() {
			return m_labeling;
		}
		ADVECTION_EVENT DigitizeSegmentInternal(Vec3d a, Vec3d b, std::vector<INDEX_TYPE>& dline/*, std::set<INDEX_TYPE>& seen_so_far*/) {

			// solve x y z intersects;
			bool has_intersection[3];
			INDEX_TYPE intersection_index[3];
			double intersection_t[3];
			int num_crossed = 0;
			//for (auto i = 0; i < count; i++) {
			//printf("here1\n");
			Vec3d a2b = b - a;
			//a.PrintFloat();
			//b.PrintFloat();
			//printf("here1a\n");
			for (int d = 0; d < 3; d++) {
					
				//printf("here1b\n");
				// if it does not cross plane, don't need further tests
					int voxel_a = floor(a[d]);
					int voxel_b = floor(b[d]);
					if (voxel_a == voxel_b) {
						has_intersection[d] = false;
						continue;
					}
					//printf("here1c\n");
					has_intersection[d] = true;

					int plane = (voxel_a > voxel_b ? voxel_a : voxel_b); // std::max(voxel_a, voxel_b); // the floor rounds down, this picks which one was crossed
					double t = (plane - a[d]) / (b[d] - a[d]); // get parameter value from a --> b
					intersection_t[num_crossed] = t;
				//	printf("here1d\n");

					if (t < 0 || t > 1) {
						printf("dont think should ever get here\n");
						return ADVECTION_EVENT::OTHER;
					}
					// get topological coordinate of the intersecting voxel
					Vec3l facecooords = (a2b * t + a).IntFloor();

					//if (facecooords[d] != plane) printf("WHOATHERE %d != %d\n", facecooords[d], plane);

					facecooords = facecooords * 2;
					//facecooords.PrintInt();
					// now add in 1 offset of face
					for (int td = 0; td < 3; td++) {
						if (d == td) continue;
						facecooords[td] += 1;
					}
					//printf("here1e\n");
					//facecooords.PrintInt();

					INDEX_TYPE faceid = m_mesh->coords2Cellid(facecooords);
					//printf("here1ea\n");
					intersection_index[num_crossed] = faceid;
					num_crossed++;

					//printf("here1f\n");


				}
				//printf("here2 %d\n", num_crossed);

			//}

			// ok now we have the (up to 3) faces that are intersected, 
			// return if none were crossed

			if (num_crossed == 0) return ADVECTION_EVENT::HIT_EXTREMUM;

			//// now re-order any that might be out of order
			//if (num_crossed == 2) {
			//	// swap if out of order
			//	if (intersection_t[0] > intersection_t[1]) {
			//		INDEX_TYPE temp = intersection_index[0];
			//		intersection_index[0] = intersection_index[1];
			//		intersection_index[1] = temp;
			//	}
			//}
			//else if (num_crossed == 3) {
				// sort
				for (int i = 0; i < num_crossed; i++) {
					for (int j = i; j < num_crossed; j++) {
						if (intersection_t[i] > intersection_t[j]) {
							// swap
							INDEX_TYPE temp = intersection_index[i];
							intersection_index[i] = intersection_index[j];
							intersection_index[j] = temp;
							double t = intersection_t[i];
							intersection_t[i] = intersection_t[j];
							intersection_t[j] = t;
						}
					}
				}

			//}
			//else if (num_crossed > 3) {
			//	printf("WHOATHERENELLY< num crossed = %d\n", num_crossed);
			//}

			// get the starting voxel
			INDEX_TYPE start_hex;
			Vec3d start_v = a;
			Vec3l base_v = (start_v.IntFloor() * 2) + Vec3l(1, 1, 1); // get the vertex coordinate, then get the hex
			start_hex = m_mesh->coords2Cellid(base_v);

			for (int i = 0; i < num_crossed; i++) {

				INDEX_TYPE next_hex;
				INDEX_TYPE faceid = intersection_index[i];
				bool has_start = false;

				typename MeshType::CofacetsIterator cfit(this->m_mesh);
				for (cfit.begin(faceid); cfit.valid(); cfit.advance()) {
					INDEX_TYPE hex_id = cfit.value();

					if (hex_id == start_hex) {
						has_start = true;
						continue;
					}

					// now hex_id is next hex in hex->face->hex path
					next_hex = hex_id;
					continue; // continue, since we want to check all hexes to make sure start is there
				}

				if (!has_start) {
					//printf("whoa no continuity\n"); 
				}
				else {
					// do yer thing!

					//if (seen_so_far.count(start_hex) == 0) {
						if (already_labeled(start_hex)) return ADVECTION_EVENT::HIT_PREASSIGNED;
						set_label(start_hex);
						dline.push_back(start_hex);
						//seen_so_far.insert(start_hex);
					//}
					//if (seen_so_far.count(faceid) == 0) {
						if (already_labeled(faceid)) return ADVECTION_EVENT::HIT_PREASSIGNED;
						set_label(faceid);
						dline.push_back(faceid);
						//seen_so_far.insert(faceid);
					//}
					//if (seen_so_far.count(next_hex) == 0) {
					//	if (already_labeled(next_hex)) return ADVECTION_EVENT::HIT_PREASSIGNED;
					//	set_label(next_hex);
					//	dline.push_back(next_hex);
					//	seen_so_far.insert(next_hex);
					//}

					//
					start_hex = next_hex;
				}

			}

			

		}

	public:


		void SetAdvectionChecker(AdvectionChecker* checker) { m_checker = checker; }
		DigitizingNumericStreamlineIntegrator3dASC(RegularGrid3D* grid, GridFuncType* func, MeshType* mesh,
			FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int max_steps) :
			m_grid(grid), m_func(func),
			m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold), m_max_steps(max_steps), m_mesh(mesh) {
			m_labeling = NULL;
			m_checker = NULL;// new TerminateNearAssignedHex<MeshType>(m_labeling, mesh);
			//m_labeling->SetAll(0);
			m_target_label = 1;
		}
		void SetDigitizingTarget(DenseLabeling<char>* other_labeling, int target) {
			m_labeling = other_labeling;
			m_target_label = target;
			m_checker = new TerminateNearAssignedHex<MeshType>(other_labeling, m_mesh, target);
		}

		GInt::ADVECTION_EVENT IntegrateStreamline(const Vec3d &seed, std::vector<Vec3d> &sline, std::vector<INDEX_TYPE>& dline)  {
			if (m_labeling == NULL) {
				printf("Error: no destination set for digitizing\n");
				return ADVECTION_EVENT::OTHER;
			}

			Vec3l start_coords = (seed + 0.5); // since casting rounds down to integer... this is starting vertex
			INDEX_TYPE start_id = m_grid->Index3d(start_coords); // vertex id of start point
			int t1, t2; t1 = t2 = 0;
			Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, m_checker);
			sline.reserve(100);

			sline.push_back(seed);
			//std::set<INDEX_TYPE> seen_so_far;
			Vec3l t_coords = start_coords; // get the coordinates of the poitn
			//if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
			Vec3d t_current_point = seed;
			int t_num_iterations_left = m_max_steps;
			bool t_continue = true;

			while (t_continue) {
				Vec3d t_old_position = t_current_point;
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
					// digitizing here
					DigitizeSegmentInternal(t_old_position, t_current_point, dline/*, seen_so_far*/);
					// end digitizing here

				}
				sline.push_back(t_current_point);


				// if we terminated or hit a critical point, then we are done
				if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
					t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
					t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
					t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {


					return t_return_code;
					t_continue = false;
				}
			}
			return ADVECTION_EVENT::NONE;
		}



	};

	 
	template<class MeshType, class GridFuncType, class Advector>
	class DigitizingNumericStreamlineIntegrator3dDSC {


		FLOATTYPE m_error_threshold;
		FLOATTYPE m_gradient_threshold;
		int m_max_steps;

		RegularGrid3D* m_grid;
		GridFuncType* m_func;
		MeshType* m_mesh;
		DenseLabeling<char>* m_labeling;
		bool in_bounds(const Vec3d &x) const {
			return (0.0 <= x[0] && x[0] < 1.0 &&
				0.0 <= x[1] && x[1] < 1.0 &&
				0.0 <= x[2] && x[2] < 1.0);
		}

		AdvectionChecker* m_checker;
		int m_target_label;
	public:
	
		bool already_labeled(INDEX_TYPE id) const { return m_labeling->GetLabel(id) == m_target_label; }
		void set_label(INDEX_TYPE id) {
			//if (m_labeling->GetLabel(id) == 0) {
			//	printf("setting background to m_target_label, instead of manifold\n");
			//}
			m_labeling->SetLabel(id, m_target_label);
		}

		DenseLabeling<char>* get_output() {
			return m_labeling;
		}
		ADVECTION_EVENT DigitizeSegmentInternal(Vec3d a, Vec3d b, std::vector<INDEX_TYPE>& dline/*, std::set<INDEX_TYPE>& seen_so_far*/) {
  
			// solve x y z intersects;
			bool has_intersection[3];
			INDEX_TYPE intersection_index[3];
			double intersection_t[3];
			int num_crossed = 0;
			//for (auto i = 0; i < count; i++) {
			//printf("here1\n");
			//a.PrintFloat();
			//b.PrintFloat();
			//printf("here1a\n");

			// NOTE THAT WE MOVE THESE SO THAT WHEN WE TAKE FLOOR WE GET WHICH HEX
			// IT IS IN
			a += Vec3d(0.5, 0.5, 0.5); // move so that we can use same test as hex
			b += Vec3d(0.5, 0.5, 0.5); // move so that we can use same test as hex
			Vec3d a2b = b - a;

			for (int d = 0; d < 3; d++) {

				//printf("here1b\n");
				// if it does not cross between voxels, don't need further tests
				int voxel_a = floor(a[d]);
				int voxel_b = floor(b[d]);
				if (voxel_a == voxel_b) {
					has_intersection[d] = false;
					continue;
				}
				//printf("here1c\n");
				has_intersection[d] = true;

				int plane =(voxel_a > voxel_b ? voxel_a : voxel_b); // the floor rounds down, this picks which one was crossed
				double t = (plane - a[d]) / (b[d] - a[d]); // get parameter value from a --> b
				intersection_t[num_crossed] = t;
				//	printf("here1d\n");

				if (t < 0 || t > 1) {
					printf("dont think should ever get here\n");
					return ADVECTION_EVENT::OTHER;
				}
				// get topological coordinate of the intersecting voxel
				Vec3l facecooords = (a2b * t + a).IntFloor();

				if (facecooords[d] != plane) {
					//printf("WHOATHERE %d != %d\n", facecooords[d], plane);
					
				}

				facecooords = facecooords * 2;
				//facecooords.PrintInt();
				// now add in 1 offset of face
				for (int td = 0; td < 3; td++) {
					if (d == td) continue;
					facecooords[td] += 1;
				}
				//printf("here1e\n");
				//facecooords.PrintInt();
				// to go from face to edge, subtract 1
				facecooords += Vec3l(-1, -1, -1);

				INDEX_TYPE faceid = m_mesh->coords2Cellid(facecooords);
				//printf("here1ea\n");
				intersection_index[num_crossed] = faceid;
				num_crossed++;

				//printf("here1f\n");


			}
			//printf("here2 %d\n", num_crossed);

			//}

			// ok now we have the (up to 3) edges that are intersected, 
			// return if none were crossed

			if (num_crossed == 0) return ADVECTION_EVENT::HIT_EXTREMUM;

			//// now re-order any that might be out of order
			//if (num_crossed == 2) {
			//	// swap if out of order
			//	if (intersection_t[0] > intersection_t[1]) {
			//		INDEX_TYPE temp = intersection_index[0];
			//		intersection_index[0] = intersection_index[1];
			//		intersection_index[1] = temp;
			//	}
			//}
			//else if (num_crossed == 3) {
			// sort
			for (int i = 0; i < num_crossed; i++) {
				for (int j = i; j < num_crossed; j++) {
					if (intersection_t[i] > intersection_t[j]) {
						// swap
						INDEX_TYPE temp = intersection_index[i];
						intersection_index[i] = intersection_index[j];
						intersection_index[j] = temp;
						double t = intersection_t[i];
						intersection_t[i] = intersection_t[j];
						intersection_t[j] = t;
					}
				}
			}

			//}
			//else if (num_crossed > 3) {
			//	printf("WHOATHERENELLY< num crossed = %d\n", num_crossed);
			//}

			// get the starting vertex
			INDEX_TYPE start_hex;
			Vec3d start_v = a;
			Vec3l base_v = ((start_v).IntFloor() * 2); // get the vertex coordinate, then get the hex
			start_hex = m_mesh->coords2Cellid(base_v);

			for (int i = 0; i < num_crossed; i++) {

				INDEX_TYPE next_hex;
				INDEX_TYPE faceid = intersection_index[i];
				bool has_start = false;

				typename MeshType::FacetsIterator cfit(this->m_mesh);
				for (cfit.begin(faceid); cfit.valid(); cfit.advance()) {
					INDEX_TYPE hex_id = cfit.value();

					if (hex_id == start_hex) {
						has_start = true;
						continue;
					}

					// now hex_id is next hex in hex->face->hex path
					next_hex = hex_id;
					continue; // continue, since we want to check all hexes to make sure start is there
				}

				if (!has_start) {
					//printf("whoa no continuity\n");
				}
				else {
					// do yer thing!

					//if (seen_so_far.count(start_hex) == 0) {
					if (already_labeled(start_hex)) return ADVECTION_EVENT::HIT_PREASSIGNED;
					set_label(start_hex);
					dline.push_back(start_hex);
					//seen_so_far.insert(start_hex);
					//}
					//if (seen_so_far.count(faceid) == 0) {
					if (already_labeled(faceid)) return ADVECTION_EVENT::HIT_PREASSIGNED;
					set_label(faceid);
					dline.push_back(faceid);
					//seen_so_far.insert(faceid);
					//}
					//if (seen_so_far.count(next_hex) == 0) {
					//	if (already_labeled(next_hex)) return ADVECTION_EVENT::HIT_PREASSIGNED;
					//	set_label(next_hex);
					//	dline.push_back(next_hex);
					//	seen_so_far.insert(next_hex);
					//}

					//
					start_hex = next_hex;
				}

			}



		}

	public:

  
		void SetAdvectionChecker(AdvectionChecker* checker) { m_checker = checker; }
		DigitizingNumericStreamlineIntegrator3dDSC(RegularGrid3D* grid, GridFuncType* func, MeshType* mesh,
			FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int max_steps) :
			m_grid(grid), m_func(func),
			m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold), m_max_steps(max_steps), m_mesh(mesh) {
			m_labeling = NULL;
			m_checker = NULL;// new TerminateNearAssignedHex<MeshType>(m_labeling, mesh);
							 //m_labeling->SetAll(0);
			m_target_label = 1;
		}
		void SetDigitizingTarget(DenseLabeling<char>* other_labeling, int target) {
			m_labeling = other_labeling;
			m_target_label = target;
			m_checker = new TerminateNearAssignedVert<MeshType>(other_labeling, m_mesh, target);
		}

		GInt::ADVECTION_EVENT IntegrateStreamline(const Vec3d &seed, std::vector<Vec3d> &sline, std::vector<INDEX_TYPE>& dline) {
			if (m_labeling == NULL) {
				printf("Error: no destination set for digitizing\n");
				return ADVECTION_EVENT::OTHER;
			}

			Vec3l start_coords = (seed + 0.5); // since casting rounds down to integer... this is starting vertex
			INDEX_TYPE start_id = m_grid->Index3d(start_coords); // vertex id of start point
			int t1, t2; t1 = t2 = 0;
			Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, m_checker);
			sline.reserve(100);

			sline.push_back(seed);
			//std::set<INDEX_TYPE> seen_so_far;
			Vec3l t_coords = start_coords; // get the coordinates of the poitn
										   //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
			Vec3d t_current_point = seed;
			int t_num_iterations_left = m_max_steps;
			bool t_continue = true;

			while (t_continue) {
				Vec3d t_old_position = t_current_point;
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
					// digitizing here
					DigitizeSegmentInternal(t_old_position, t_current_point, dline/*, seen_so_far*/);
					// end digitizing here

				}
				sline.push_back(t_current_point);


				// if we terminated or hit a critical point, then we are done
				if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
					t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
					t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
					t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {


					return t_return_code;
					t_continue = false;
				}
			}
			return ADVECTION_EVENT::NONE;
		}



	};




	template<class MeshType, class GridFuncType, class Advector>
	class DigitizingNumericStreamlineIntegrator2dASC {


		FLOATTYPE m_error_threshold;
		FLOATTYPE m_gradient_threshold;
		int m_max_steps;

		RegularGrid2D* m_grid;
		GridFuncType* m_func;
		MeshType* m_mesh;
		DenseLabeling<char>* m_labeling;
		bool in_bounds(const Vec2d &x) const {
			return (0.0 <= x[0] && x[0] < 1.0 &&
				0.0 <= x[1] && x[1] < 1.0);
		}

		AdvectionChecker2D* m_checker;
		int m_target_label;
	public:

		bool already_labeled(INDEX_TYPE id) const { return m_labeling->GetLabel(id) == m_target_label; }
		void set_label(INDEX_TYPE id) {
			//if (m_labeling->GetLabel(id) == 0) {
			//	printf("setting background to m_target_label, instead of manifold\n");
			//}
			m_labeling->SetLabel(id, m_target_label);
		}

		DenseLabeling<char>* get_output() {
			return m_labeling;
		}
		ADVECTION_EVENT DigitizeSegmentInternal(Vec2d a, Vec2d b, std::vector<INDEX_TYPE>& dline/*, std::set<INDEX_TYPE>& seen_so_far*/) {

			// solve x y z intersects;
			bool has_intersection[2];
			INDEX_TYPE intersection_index[2];
			double intersection_t[2];
			int num_crossed = 0;
			//for (auto i = 0; i < count; i++) {
			//printf("here1\n");
			Vec2d a2b = b - a;
			//a.PrintFloat();
			//b.PrintFloat();
			//printf("here1a\n");
			for (int d = 0; d < 2; d++) {

				//printf("here1b\n");
				// if it does not cross plane, don't need further tests
				int voxel_a = floor(a[d]);
				int voxel_b = floor(b[d]);
				if (voxel_a == voxel_b) {
					has_intersection[d] = false;
					continue;
				}
				//printf("here1c\n");
				has_intersection[d] = true;

				int plane = (voxel_a > voxel_b ? voxel_a : voxel_b); // std::max(voxel_a, voxel_b); // the floor rounds down, this picks which one was crossed
				double t = (plane - a[d]) / (b[d] - a[d]); // get parameter value from a --> b
				intersection_t[num_crossed] = t;
				//	printf("here1d\n");

				if (t < 0 || t > 1) {
					printf("dont think should ever get here\n");
					return ADVECTION_EVENT::OTHER;
				}
				// get topological coordinate of the intersecting voxel
				Vec2l facecooords = (a2b * t + a).IntFloor();

				//if (facecooords[d] != plane) printf("WHOATHERE %d != %d\n", facecooords[d], plane);

				facecooords = facecooords * 2;
				//facecooords.PrintInt();
				// now add in 1 offset of face
				for (int td = 0; td < 2; td++) {
					if (d == td) continue;
					facecooords[td] += 1;
				}
				//printf("here1e\n");
				//facecooords.PrintInt();

				INDEX_TYPE faceid = m_mesh->coords2Cellid(facecooords);
				//printf("here1ea\n");
				intersection_index[num_crossed] = faceid;
				num_crossed++;

				//printf("here1f\n");


			}
			//printf("here2 %d\n", num_crossed);

			//}

			// ok now we have the (up to 3) faces that are intersected, 
			// return if none were crossed

			if (num_crossed == 0) return ADVECTION_EVENT::HIT_EXTREMUM;

			//// now re-order any that might be out of order
			//if (num_crossed == 2) {
			//	// swap if out of order
			//	if (intersection_t[0] > intersection_t[1]) {
			//		INDEX_TYPE temp = intersection_index[0];
			//		intersection_index[0] = intersection_index[1];
			//		intersection_index[1] = temp;
			//	}
			//}
			//else if (num_crossed == 3) {
			// sort
			for (int i = 0; i < num_crossed; i++) {
				for (int j = i; j < num_crossed; j++) {
					if (intersection_t[i] > intersection_t[j]) {
						// swap
						INDEX_TYPE temp = intersection_index[i];
						intersection_index[i] = intersection_index[j];
						intersection_index[j] = temp;
						double t = intersection_t[i];
						intersection_t[i] = intersection_t[j];
						intersection_t[j] = t;
					}
				}
			}

			//}
			//else if (num_crossed > 3) {
			//	printf("WHOATHERENELLY< num crossed = %d\n", num_crossed);
			//}

			// get the starting voxel
			INDEX_TYPE start_hex;
			Vec2d start_v = a;
			Vec2l base_v = (start_v.IntFloor() * 2) + Vec2l(1, 1); // get the vertex coordinate, then get the hex
			start_hex = m_mesh->coords2Cellid(base_v);

			for (int i = 0; i < num_crossed; i++) {

				INDEX_TYPE next_hex;
				INDEX_TYPE faceid = intersection_index[i];
				bool has_start = false;

				typename MeshType::CofacetsIterator cfit(this->m_mesh);
				for (cfit.begin(faceid); cfit.valid(); cfit.advance()) {
					INDEX_TYPE hex_id = cfit.value();

					if (hex_id == start_hex) {
						has_start = true;
						continue;
					}

					// now hex_id is next hex in hex->face->hex path
					next_hex = hex_id;
					continue; // continue, since we want to check all hexes to make sure start is there
				}

				if (!has_start) {
					//printf("whoa no continuity\n"); 
				}
				else {
					// do yer thing!

					//if (seen_so_far.count(start_hex) == 0) {
					if (already_labeled(start_hex)) return ADVECTION_EVENT::HIT_PREASSIGNED;
					set_label(start_hex);
					dline.push_back(start_hex);
					//seen_so_far.insert(start_hex);
					//}
					//if (seen_so_far.count(faceid) == 0) {
					if (already_labeled(faceid)) return ADVECTION_EVENT::HIT_PREASSIGNED;
					set_label(faceid);
					dline.push_back(faceid);
					//seen_so_far.insert(faceid);
					//}
					//if (seen_so_far.count(next_hex) == 0) {
					//	if (already_labeled(next_hex)) return ADVECTION_EVENT::HIT_PREASSIGNED;
					//	set_label(next_hex);
					//	dline.push_back(next_hex);
					//	seen_so_far.insert(next_hex);
					//}

					//
					start_hex = next_hex;
				}

			}



		}

	public:


		void SetAdvectionChecker(AdvectionChecker* checker) { m_checker = checker; }
		DigitizingNumericStreamlineIntegrator2dASC(RegularGrid2D* grid, GridFuncType* func, MeshType* mesh,
			FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int max_steps) :
			m_grid(grid), m_func(func),
			m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold), m_max_steps(max_steps), m_mesh(mesh) {
			m_labeling = NULL;
			m_checker = NULL;// new TerminateNearAssignedHex<MeshType>(m_labeling, mesh);
			//m_labeling->SetAll(0);
			m_target_label = 1;
		}
		void SetDigitizingTarget(DenseLabeling<char>* other_labeling, int target) {
			m_labeling = other_labeling;
			m_target_label = target;
			m_checker = new TerminateNearAssignedQuad2D<MeshType>(other_labeling, m_mesh, target);
		}

		GInt::ADVECTION_EVENT IntegrateStreamline(const Vec2d &seed, std::vector<Vec2d> &sline, std::vector<INDEX_TYPE>& dline)  {
			if (m_labeling == NULL) {
				printf("Error: no destination set for digitizing\n");
				return ADVECTION_EVENT::OTHER;
			}

			Vec2l start_coords = (seed + 0.5); // since casting rounds down to integer... this is starting vertex
			INDEX_TYPE start_id = m_grid->Index2d(start_coords); // vertex id of start point
			int t1, t2; t1 = t2 = 0;
			Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, m_checker);
			sline.reserve(100);

			sline.push_back(seed);
			//std::set<INDEX_TYPE> seen_so_far;
			Vec2l t_coords = start_coords; // get the coordinates of the poitn
			//if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
			Vec2d t_current_point = seed;
			int t_num_iterations_left = m_max_steps;
			bool t_continue = true;

			while (t_continue) {
				Vec2d t_old_position = t_current_point;
				Vec2d t_next_point;
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
					// digitizing here
					DigitizeSegmentInternal(t_old_position, t_current_point, dline/*, seen_so_far*/);
					// end digitizing here

				}
				sline.push_back(t_current_point);


				// if we terminated or hit a critical point, then we are done
				if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
					t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
					t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
					t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {


					return t_return_code;
					t_continue = false;
				}
			}
			return ADVECTION_EVENT::NONE;
		}



	};

	template<class MeshType, class GridFuncType, class Advector>
	class DigitizingNumericStreamlineIntegrator2dDSC {


		FLOATTYPE m_error_threshold;
		FLOATTYPE m_gradient_threshold;
		int m_max_steps;

		RegularGrid2D* m_grid;
		GridFuncType* m_func;
		MeshType* m_mesh;
		DenseLabeling<char>* m_labeling;
		bool in_bounds(const Vec2d &x) const {
			return (0.0 <= x[0] && x[0] < 1.0 &&
				0.0 <= x[1] && x[1] < 1.0);
		}

		AdvectionChecker2D* m_checker;
		int m_target_label;
	public:

		bool already_labeled(INDEX_TYPE id) const { return m_labeling->GetLabel(id) == m_target_label; }
		void set_label(INDEX_TYPE id) {
			//if (m_labeling->GetLabel(id) == 0) {
			//	printf("setting background to m_target_label, instead of manifold\n");
			//}
			m_labeling->SetLabel(id, m_target_label);
		}

		DenseLabeling<char>* get_output() {
			return m_labeling;
		}
		ADVECTION_EVENT DigitizeSegmentInternal(Vec2d a, Vec2d b, std::vector<INDEX_TYPE>& dline/*, std::set<INDEX_TYPE>& seen_so_far*/) {

			// solve x y z intersects;
			bool has_intersection[2];
			INDEX_TYPE intersection_index[2];
			double intersection_t[2];
			int num_crossed = 0;
			//for (auto i = 0; i < count; i++) {
			//printf("here1\n");
			//a.PrintFloat();
			//b.PrintFloat();
			//printf("here1a\n");

			// NOTE THAT WE MOVE THESE SO THAT WHEN WE TAKE FLOOR WE GET WHICH HEX
			// IT IS IN
			a += Vec2d(0.5, 0.5); // move so that we can use same test as hex
			b += Vec2d(0.5, 0.5); // move so that we can use same test as hex
			Vec2d a2b = b - a;

			for (int d = 0; d < 2; d++) {

				//printf("here1b\n");
				// if it does not cross between voxels, don't need further tests
				int voxel_a = floor(a[d]);
				int voxel_b = floor(b[d]);
				if (voxel_a == voxel_b) {
					has_intersection[d] = false;
					continue;
				}
				//printf("here1c\n");
				has_intersection[d] = true;

				int plane = (voxel_a > voxel_b ? voxel_a : voxel_b); // the floor rounds down, this picks which one was crossed
				double t = (plane - a[d]) / (b[d] - a[d]); // get parameter value from a --> b
				intersection_t[num_crossed] = t;
				//	printf("here1d\n");

				if (t < 0 || t > 1) {
					printf("dont think should ever get here\n");
					return ADVECTION_EVENT::OTHER;
				}
				// get topological coordinate of the intersecting voxel
				Vec2l facecooords = (a2b * t + a).IntFloor();

				if (facecooords[d] != plane) {
					//printf("WHOATHERE %d != %d\n", facecooords[d], plane);

				}

				facecooords = facecooords * 2;
				//facecooords.PrintInt();
				// now add in 1 offset of face
				for (int td = 0; td < 2; td++) {
					if (d == td) continue;
					facecooords[td] += 1;
				}
				//printf("here1e\n");
				//facecooords.PrintInt();
				// to go from face to edge, subtract 1
				facecooords += Vec2l(-1, -1);

				INDEX_TYPE faceid = m_mesh->coords2Cellid(facecooords);
				//printf("here1ea\n");
				intersection_index[num_crossed] = faceid;
				num_crossed++;

				//printf("here1f\n");


			}
			//printf("here2 %d\n", num_crossed);

			//}

			// ok now we have the (up to 3) edges that are intersected, 
			// return if none were crossed

			if (num_crossed == 0) return ADVECTION_EVENT::HIT_EXTREMUM;

			//// now re-order any that might be out of order
			//if (num_crossed == 2) {
			//	// swap if out of order
			//	if (intersection_t[0] > intersection_t[1]) {
			//		INDEX_TYPE temp = intersection_index[0];
			//		intersection_index[0] = intersection_index[1];
			//		intersection_index[1] = temp;
			//	}
			//}
			//else if (num_crossed == 3) {
			// sort
			for (int i = 0; i < num_crossed; i++) {
				for (int j = i; j < num_crossed; j++) {
					if (intersection_t[i] > intersection_t[j]) {
						// swap
						INDEX_TYPE temp = intersection_index[i];
						intersection_index[i] = intersection_index[j];
						intersection_index[j] = temp;
						double t = intersection_t[i];
						intersection_t[i] = intersection_t[j];
						intersection_t[j] = t;
					}
				}
			}

			//}
			//else if (num_crossed > 3) {
			//	printf("WHOATHERENELLY< num crossed = %d\n", num_crossed);
			//}

			// get the starting vertex
			INDEX_TYPE start_hex;
			Vec2d start_v = a;
			Vec2l base_v = ((start_v).IntFloor() * 2); // get the vertex coordinate, then get the hex
			start_hex = m_mesh->coords2Cellid(base_v);

			for (int i = 0; i < num_crossed; i++) {

				INDEX_TYPE next_hex;
				INDEX_TYPE faceid = intersection_index[i];
				bool has_start = false;

				typename MeshType::FacetsIterator cfit(this->m_mesh);
				for (cfit.begin(faceid); cfit.valid(); cfit.advance()) {
					INDEX_TYPE hex_id = cfit.value();

					if (hex_id == start_hex) {
						has_start = true;
						continue;
					}

					// now hex_id is next hex in hex->face->hex path
					next_hex = hex_id;
					continue; // continue, since we want to check all hexes to make sure start is there
				}

				if (!has_start) {
					//printf("whoa no continuity\n");
				}
				else {
					// do yer thing!

					//if (seen_so_far.count(start_hex) == 0) {
					if (already_labeled(start_hex)) return ADVECTION_EVENT::HIT_PREASSIGNED;
					set_label(start_hex);
					dline.push_back(start_hex);
					//seen_so_far.insert(start_hex);
					//}
					//if (seen_so_far.count(faceid) == 0) {
					if (already_labeled(faceid)) return ADVECTION_EVENT::HIT_PREASSIGNED;
					set_label(faceid);
					dline.push_back(faceid);
					//seen_so_far.insert(faceid);
					//}
					//if (seen_so_far.count(next_hex) == 0) {
					//	if (already_labeled(next_hex)) return ADVECTION_EVENT::HIT_PREASSIGNED;
					//	set_label(next_hex);
					//	dline.push_back(next_hex);
					//	seen_so_far.insert(next_hex);
					//}

					//
					start_hex = next_hex;
				}

			}



		}

	public:


		void SetAdvectionChecker(AdvectionChecker2D* checker) { m_checker = checker; }
		DigitizingNumericStreamlineIntegrator2dDSC(RegularGrid2D* grid, GridFuncType* func, MeshType* mesh,
			FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int max_steps) :
			m_grid(grid), m_func(func),
			m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold), m_max_steps(max_steps), m_mesh(mesh) {
			m_labeling = NULL;
			m_checker = NULL;// new TerminateNearAssignedHex<MeshType>(m_labeling, mesh);
			//m_labeling->SetAll(0);
			m_target_label = 1;
		}
		void SetDigitizingTarget(DenseLabeling<char>* other_labeling, int target) {
			m_labeling = other_labeling;
			m_target_label = target;
			m_checker = new TerminateNearAssignedVert2D<MeshType>(other_labeling, m_mesh, target);
		}

		GInt::ADVECTION_EVENT IntegrateStreamline(const Vec2d &seed, std::vector<Vec2d> &sline, std::vector<INDEX_TYPE>& dline) {
			if (m_labeling == NULL) {
				printf("Error: no destination set for digitizing\n");
				return ADVECTION_EVENT::OTHER;
			}

			Vec2l start_coords = (seed + 0.5); // since casting rounds down to integer... this is starting vertex
			INDEX_TYPE start_id = m_grid->Index2d(start_coords); // vertex id of start point
			int t1, t2; t1 = t2 = 0;
			Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, m_checker);
			sline.reserve(100);

			sline.push_back(seed);
			//std::set<INDEX_TYPE> seen_so_far;
			Vec2l t_coords = start_coords; // get the coordinates of the poitn
			//if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
			Vec2d t_current_point = seed;
			int t_num_iterations_left = m_max_steps;
			bool t_continue = true;

			while (t_continue) {
				Vec2d t_old_position = t_current_point;
				Vec2d t_next_point;
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
					// digitizing here
					DigitizeSegmentInternal(t_old_position, t_current_point, dline/*, seen_so_far*/);
					// end digitizing here

				}
				sline.push_back(t_current_point);


				// if we terminated or hit a critical point, then we are done
				if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
					t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
					t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
					t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {


					return t_return_code;
					t_continue = false;
				}
			}
			return ADVECTION_EVENT::NONE;
		}



	};


}
#endif
