/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef TOPOLOGICAL_REGULAR_GRID_2D_H
#define TOPOLOGICAL_REGULAR_GRID_2D_H


#include <unordered_map>
#include <set>

#include "base/gi_basic_types.h"
#include "base/gi_vectors.h"
#include "base/gi_regular_grid_2d.h"



namespace GInt {


	class TopologicalRegularGrid2D {
	public:
		void display() const {

			for (int i = 0; i < 8; i++) {
				printf(" m_cell_vertex_offsets for celltype = %d : (", i);

				/*for(int j = 1; j < m_cell_vertex_offsets[i][0]; j++) {
				printf(" %d ", m_cell_vertex_offsets[i][j]);
				}
				/*/for (int j = 1; j < m_cell_vertex_offsets[i][0]; j++) {
					printf(" %d [%d %d] ", m_cell_vertex_offsets[i][j],
						m_cell_vertex_is_in_positive_axis_directions[i][j][0],
						m_cell_vertex_is_in_positive_axis_directions[i][j][1]);
				}

				printf(")\n");
			}

			printf("\n m_adjacency_offsets = (");
			for (int i = 0; i < 4; i++) {
				printf(" %d ", m_adjacency_offsets[i]);
			}
			printf(")\n");

			printf("\n m_wraparound_adjacency_offsets = (");
			for (int i = 0; i < 4; i++) {
				printf(" %d ", m_wraparound_adjacency_offsets[i]);
			}
			printf(")\n");
		}

		INDEX_TYPE get8NeighborOffset(int pos) const {
			return m_adjacent_cell_offsets[pos];
		}
	protected:

		RegularGrid2D* m_base_grid;
		INDEX_TYPE m_num_cells;			// total number of cells (vertices, edges, quads, hex)
		Vec2l		m_base_grid_xy;	// the number of vertices in each dimension of the base grid
		Vec2l		m_mesh_xy;			// the number of *CELLS* in each dimension of our topological mesh
		Vec2b		m_mesh_periodic;
		Vec2l		m_mesh_not_periodic;// whether or not each axis is periodic 0 if yes, 1 if no
		INDEX_TYPE m_num_dcells[4];		// store the number of vertices, edges, quads, hexahedra in our mesh
		INDEX_TYPE m_adjacency_offsets[4]; // precompute and store the integer offset to the 4 adjacent cellst
		INDEX_TYPE m_wraparound_adjacency_offsets[4];
		INDEX_TYPE m_facet_positions[4][5];
		DIM_TYPE		m_facet_along_which_axis[4][5]; // slight misnomer - actually, the result of this / 2 (integer math) is the axis
		INDEX_TYPE m_cofacet_positions[4][5];
		DIM_TYPE		m_cofacet_along_which_axis[4][5];
		std::unordered_map<INDEX_TYPE, BYTE_TYPE> m_offset_2_direction_map;
		std::unordered_map<INDEX_TYPE, BYTE_TYPE> m_offset_2_position_vertices_map;

		INDEX_TYPE m_adjacent_cell_offsets[9];

		INDEX_TYPE m_cell_vertex_offsets[4][5];                 // there are 8 kinds of cells, and each can have up to 8 vertices

		INDEX_TYPE m_wrapping_offsets[2];                       // new.. created by Harsh.. to handle wrapping boundarues
		Vec2b m_cell_vertex_is_in_positive_axis_directions[4][5]; // records for each element of m_cell_vertex_offsets if that vertex is in positive x,y,z directions (used for boundary of periodic meshes)

		INT_TYPE Extent(DIM_TYPE val) const {
			return m_mesh_xy[val];
		}

		int GetCellTypeAndOrientation(Vec2l p) const {
			return (p[0] % 2) + 2 * (p[1] % 2);
		}

		inline INDEX_TYPE dot(const Vec2i &w, const Vec2l &v) {
			return w[0] * v[0] + w[1] * v[1];
		}


		void InitializeValues() {
			// set the dimesions of the topological mesh - double the number of vertices in each axis,
			// but subtract one from each non-periodic axis
			m_mesh_xy = m_base_grid_xy * 2;
			for (int i = 0; i < 2; i++){
				if (m_base_grid->Periodic()[i])
					m_mesh_not_periodic[i] = 0;
				else
					m_mesh_not_periodic[i] = 1;
			}
			m_mesh_xy -= m_mesh_not_periodic;
			m_num_cells = m_mesh_xy[0] * m_mesh_xy[1];

			m_mesh_periodic = m_base_grid->Periodic();

			/*printf("input (%d,%d,%d), new dim (%d,%d,%d), num_cells = %'d\n",
			m_base_grid_xyz[0], m_base_grid_xyz[1], m_base_grid_xyz[2],
			m_mesh_xyz[0], m_mesh_xyz[1], m_mesh_xyz[2],
			m_num_cells);*/


			//numbers of cells of each dimension
			// number of vertices is easy 
			m_num_dcells[0] = m_base_grid_xy[0] * m_base_grid_xy[1];

			// number of edges depends on periodicity - compute each for edges pointing in xdirection, y, and z
			INDEX_TYPE t_num_1cells_aligned_in_x = (m_mesh_xy[0] / 2)*m_base_grid_xy[1];
			INDEX_TYPE t_num_1cells_aligned_in_y = m_base_grid_xy[0] * (m_mesh_xy[1] / 2);
			m_num_dcells[1] = t_num_1cells_aligned_in_x + t_num_1cells_aligned_in_y;


			// and for 2-cells
			m_num_dcells[2] = (m_mesh_xy[0] / 2)*(m_mesh_xy[1] / 2);

			//printf("num_dcells (%d,%d,%d,%d)\n", m_num_dcells[0], m_num_dcells[1], m_num_dcells[2], m_num_dcells[3] );


			// useful vector to compute many properties needed below!
			const Vec2l meshoffset(1, m_mesh_xy[0]);     // offsets in two dimensions

			// harsh created a new wrapping offset in three dimensions
			m_wrapping_offsets[0] = -m_mesh_xy[0];
			m_wrapping_offsets[1] = -m_mesh_xy[0] * m_mesh_xy[1];

			//offset_list
			m_adjacency_offsets[0] = dot(Vec2i(1, 0), meshoffset);           // +x
			m_adjacency_offsets[1] = dot(Vec2i(-1, 0), meshoffset);           // -x
			m_adjacency_offsets[2] = dot(Vec2i(0, 1), meshoffset);           // +y
			m_adjacency_offsets[3] = dot(Vec2i(0, -1), meshoffset);           // -y


			/*m_adjacency_offsets[0] = (1) + (0) * m_mesh_xyz[0] + (0) * m_mesh_xyz[0] * m_mesh_xyz[1];
			m_adjacency_offsets[1] = (-1) + (0) * m_mesh_xyz[0] + (0) * m_mesh_xyz[0] * m_mesh_xyz[1];
			m_adjacency_offsets[2] = (0) + (1) * m_mesh_xyz[0] + (0) * m_mesh_xyz[0] * m_mesh_xyz[1];
			m_adjacency_offsets[3] = (0) + (-1) * m_mesh_xyz[0] + (0) * m_mesh_xyz[0] * m_mesh_xyz[1];
			m_adjacency_offsets[4] = (0) + (0) * m_mesh_xyz[0] + (1) * m_mesh_xyz[0] * m_mesh_xyz[1];
			m_adjacency_offsets[5] = (0) + (0) * m_mesh_xyz[0] + (-1) * m_mesh_xyz[0] * m_mesh_xyz[1];*/

			// we will also compute the offsets that wrap around the periodic dimensions - the nonperiodic ones will not be used
			m_wraparound_adjacency_offsets[0] = dot(Vec2i(-(m_mesh_xy[0] - 1), 0), meshoffset); // im at X-1, want offset getting me to 0
			m_wraparound_adjacency_offsets[1] = dot(Vec2i((m_mesh_xy[0] - 1), 0), meshoffset); // im at 0, want offset getting me to X-1
			m_wraparound_adjacency_offsets[2] = dot(Vec2i(0, -(m_mesh_xy[1] - 1)), meshoffset);
			m_wraparound_adjacency_offsets[3] = dot(Vec2i(0, (m_mesh_xy[1] - 1)), meshoffset);


			/*m_wraparound_adjacency_offsets[0] = (-1 * (m_mesh_xyz[0] - 1)) + (0) * m_mesh_xyz[0] + (0) * m_mesh_xyz[0] * m_mesh_xyz[1]; // im at X-1, want offset getting me to 0
			m_wraparound_adjacency_offsets[1] = (1 * (m_mesh_xyz[0] - 1)) + (0) * m_mesh_xyz[0] + (0) * m_mesh_xyz[0] * m_mesh_xyz[1]; // im at 0, want offset getting me to X-1
			m_wraparound_adjacency_offsets[2] = (0) + (-1 * (m_mesh_xyz[1] - 1)) * m_mesh_xyz[0] + (0) * m_mesh_xyz[0] * m_mesh_xyz[1];
			m_wraparound_adjacency_offsets[3] = (0) + (1 * (m_mesh_xyz[1] - 1)) * m_mesh_xyz[0] + (0) * m_mesh_xyz[0] * m_mesh_xyz[1];
			m_wraparound_adjacency_offsets[4] = (0) + (0) * m_mesh_xyz[0] + (-1 * (m_mesh_xyz[2] - 1)) * m_mesh_xyz[0] * m_mesh_xyz[1];
			m_wraparound_adjacency_offsets[5] = (0) + (0) * m_mesh_xyz[0] + (1 * (m_mesh_xyz[2] - 1)) * m_mesh_xyz[0] * m_mesh_xyz[1];*/

			// fill in all adjacent cells
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					INDEX_TYPE tid = (i + 1) + 3 * (j + 1);
					m_adjacent_cell_offsets[tid] = dot(Vec2i(i, j), meshoffset);
					//i + j *  m_mesh_xyz[0] + k * m_mesh_xyz[0] * m_mesh_xyz[1];
				}
			}


			// use this to fill in adjacency lists
			//printf(" fill adjacency lists...\n");
			for (int j = 0; j < 2; j++) {
				for (int i = 0; i < 2; i++) {

					int celltype = i + 2 * j;
					int place = 1;

					for (int ii = 0; ii <= i % 2; ii++) {
						for (int jj = 0; jj <= j % 2; jj++) {


							int lindex = (ii * 2 + 1 - i) +
								(jj * 2 + 1 - j) * 3;

							Vec2b t(ii == 1, jj == 1);

							//t[0] = ii == 1;
							//t[1] = jj == 1;
							//t[2] = kk == 1;
							//m_cell_vertex_is_in_positive_axis_directions[i + 2 * j + 4 * k][place] = t;
							//m_cell_vertex_offsets[i + 2*j + 4*k][place++] = m_adjacent_cell_offsets[lindex];

							m_cell_vertex_is_in_positive_axis_directions[celltype][place] = t;
							m_cell_vertex_offsets[celltype][place] = m_adjacent_cell_offsets[lindex];

							place++;

						}
					}
					m_cell_vertex_offsets[celltype][0] = place;             // size
				}
			}



			// create map from index diference to "direction" -- position in m_adjacent_cell_offsets list
			for (int i = 0; i < 9; i++) {
				m_offset_2_position_vertices_map[m_adjacent_cell_offsets[i]] = i;
				// do nothing for the i=j=k=0 point
			}

			//// fill in all adjacent cells
			//for (int i = -1; i <= 1; i++) {
			//	for (int j = -1; j <= 1; j++) {
			//		for (int k = -1; k <= 1; k++) {
			//			m_adjacent_cell_offsets[(i + 1) + 3 * (j + 1) + 9 * (k + 1)] =
			//				i + j *  m_mesh_xyz[0] + k * m_mesh_xyz[0] * m_mesh_xyz[1];
			//		}
			//	}
			//}

			// go through each cell type - 8 cell types, characterized by even or odd
			//  coordinates. e.g. all even means it's a vertex. if only the x is odd, 
			// then its an edge aligned with x axis, if all three odd, then a 3-cell
			// each type of cell can have facets ONLY alond X, Y, Z directions. so 
			// simply check each coordinate - if it's odd, add the offsets to the list
			// in those directions, and record which coordinate it corresponds to
			// in mFacetAlongAxisList
			// the cell type and orientation is as follows:
			// id = 0 - vertex
			// id = 1 - x-edge
			// id = 2 - y-edge
			// id = 3 - xy-face

			// for each facet list, the first element is always the number of facets, then a sequence of 
			// that many offsets

			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {

					int m_cell_type_and_orientation = GetCellTypeAndOrientation(Vec2i(i, j));
					DIM_TYPE dim = i + j;

					int place = 1;
					if (i % 2 == 1) {
						m_facet_along_which_axis[m_cell_type_and_orientation][place] = 0;
						m_facet_positions[m_cell_type_and_orientation][place++] = m_adjacency_offsets[0];
						m_facet_along_which_axis[m_cell_type_and_orientation][place] = 1;
						m_facet_positions[m_cell_type_and_orientation][place++] = m_adjacency_offsets[1];

					}
					if (j % 2 == 1) {
						m_facet_along_which_axis[m_cell_type_and_orientation][place] = 2;
						m_facet_positions[m_cell_type_and_orientation][place++] = m_adjacency_offsets[2];
						m_facet_along_which_axis[m_cell_type_and_orientation][place] = 3;
						m_facet_positions[m_cell_type_and_orientation][place++] = m_adjacency_offsets[3];
					}

					m_facet_positions[m_cell_type_and_orientation][0] = place - 1;
				}
			}


			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {

					int m_cell_type_and_orientation = GetCellTypeAndOrientation(Vec2i(i, j));
					DIM_TYPE dim = i + j;

					int place = 1;
					if (i % 2 == 0) {
						m_cofacet_along_which_axis[m_cell_type_and_orientation][place] = 0;
						m_cofacet_positions[m_cell_type_and_orientation][place++] = m_adjacency_offsets[0];
						m_cofacet_along_which_axis[m_cell_type_and_orientation][place] = 1;
						m_cofacet_positions[m_cell_type_and_orientation][place++] = m_adjacency_offsets[1];
					}
					if (j % 2 == 0) {
						m_cofacet_along_which_axis[m_cell_type_and_orientation][place] = 2;
						m_cofacet_positions[m_cell_type_and_orientation][place++] = m_adjacency_offsets[2];
						m_cofacet_along_which_axis[m_cell_type_and_orientation][place] = 3;
						m_cofacet_positions[m_cell_type_and_orientation][place++] = m_adjacency_offsets[3];
					}

					m_cofacet_positions[m_cell_type_and_orientation][0] = place - 1;
				}
			}


			for (BYTE_TYPE i = 0; i < 4; i++) {
				m_offset_2_direction_map[m_adjacency_offsets[i]] = i;
				// we have a potential bug on the next line, but quick fix for now
				// basically the offsets willnot be unique if any of the extents are = 0,1,2
				// if wraparound offsets are equal to offset position
				// at least run a check and spit out error message
				if (m_offset_2_direction_map.count(m_wraparound_adjacency_offsets[i]) != 0) {
					printf("WARNING: offset2position map already has %ll, old pos = %d, new pos = %d\n",
						m_wraparound_adjacency_offsets[i], m_offset_2_direction_map[m_wraparound_adjacency_offsets[i]], i);
				}
				m_offset_2_direction_map[m_wraparound_adjacency_offsets[i]] = i; // potential bug! 
			}
		}


		//BOUNDARY_TYPE BlockBoundaryValue(Vec2l coords) {
		//	return (BOUNDARY_TYPE)((coords[0] == 0) || (coords[0] == m_mesh_xyz[0] - 1)) +
		//		((coords[1] == 0) || (coords[1] == m_mesh_xyz[1] - 1)) +
		//		((coords[2] == 0) || (coords[2] == m_mesh_xyz[2] - 1));
		//}

		bool IsPeriodic(int axis) const {
			return m_base_grid->Periodic()[axis];
		}

	public:
		class AllCellsIterator {
		protected:
			INDEX_TYPE m_pos;
			const INDEX_TYPE m_start;
			const INDEX_TYPE m_end;

		public:
			AllCellsIterator(TopologicalRegularGrid2D* m) :
				m_pos(0), m_start(0), m_end(m->numCells()) {}

			AllCellsIterator(TopologicalRegularGrid2D* m, INDEX_TYPE start, INDEX_TYPE end) :
				m_pos(start), m_start(start), m_end(end) {}

			void begin() {
				m_pos = m_start;
			}
			void advance() {
				m_pos++;
			}
			bool valid() const {
				return m_pos < m_end;
			}
			INDEX_TYPE value() const {
				return m_pos;
			}
		};

		class DCellsIterator {
		protected:
			INDEX_TYPE m_pos;
			const INDEX_TYPE m_start;
			const INDEX_TYPE m_end;
			TopologicalRegularGrid2D* const m_mesh;
			const DIM_TYPE mDim;
		public:
			DCellsIterator(TopologicalRegularGrid2D* mesh, DIM_TYPE dim) :
				m_pos(0), m_start(0), m_end(mesh->numCells()), mDim(dim), m_mesh(mesh) {}
			DCellsIterator(TopologicalRegularGrid2D* mesh, DIM_TYPE dim, INDEX_TYPE start, INDEX_TYPE end) :
				m_pos(start), m_end(end), m_start(start), mDim(dim), m_mesh(mesh) {}

			void begin() {
				m_pos = m_start;
				if (m_mesh->dimension(m_pos) != mDim) advance();
			}
			void advance() {
				m_pos++;
				while (valid() && m_mesh->dimension(m_pos) != mDim) m_pos++;
			}
			bool valid() const {
				return m_pos < m_end;
			}
			INDEX_TYPE value() const {
				return m_pos;
			}
		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class FacetsIterator {
		protected:
			int m_pos, m_end;
			INDEX_TYPE m_base_index;
			bool m_use_wraparound;
			Vec2l m_coords;
			const TopologicalRegularGrid2D* const m_mesh;
			BOUNDARY_TYPE m_block_boundary_value; // whether we are on boundary of block, as laid out in memory
			int m_cell_type_and_orientation; // where to look in the facets list for a cell of this type (dimesion and orientation) - based on order we assigned things in 

			void advance_boundary() {
				// increment position
				m_pos++;
				// if invalid, return
				if (m_pos > m_end) return;

				// now keep advancing while we are valid and are not pointing INSIDE global mesh
				while (m_pos <= m_end) {
					// THe logic of this block is as follows:
					// - future calls to value() of this iterator simply add the offset of
					//   m_pos to the m_base_index cell index
					// - we need to increment m_pos such that the offset represents a step in the
					//   grid that is actually inside the mesh
					// - we have stored offsets in the list such that ODD elements of the list
					//   point in the positive direction on an axis, and even ones point in a 
					//   negative direction. furthermore the DIRECTION we are looking can be 
					//   derived from m_pos - X axis is first 2 elements, then 2 for y, then 2 
					//   for z. 
					// - so if we are even, check if we are at coordinate 0 of the axis m_pos/2
					// - if we are on a boundary of the block, do a further check to see if that
					//   direction is periodic, if it is, m_pos is a valid direction, but set the 
					//   bool that will tell us to use the wraparound list in future value() calls
					// - otherwise the direction is truly out of the domain, and keep looking for a 
					//   m_pos that works
					DIM_TYPE tOffsetAxis = m_mesh->m_facet_along_which_axis[m_cell_type_and_orientation][m_pos] / 2;
					if (m_pos % 2 == 0) { // then we are in the NEGATIVE direction
						if (m_coords[tOffsetAxis] != 0) {
							m_use_wraparound = false;
							return;
						}
						else {
							if (m_mesh->IsPeriodic(tOffsetAxis)) { // use wraparound next
								m_use_wraparound = true;
								return;
							}
						}
					}
					else { //  therefore (m_pos % 2 == 1) MUST BE TRUE so we are looking in positive direction
						if (m_coords[tOffsetAxis] != m_mesh->Extent(tOffsetAxis) - 1) {
							m_use_wraparound = false;
							return;
						}
						else {
							if (m_mesh->IsPeriodic(tOffsetAxis)) { // use wraparound next
								m_use_wraparound = true;
								return;
							}
						}
					}
					m_pos++;
				}
			}

		private:
			void init() {

				m_block_boundary_value = m_mesh->MemoryBlockBoundaryValue(m_coords); // get global boundary information
				m_cell_type_and_orientation = m_mesh->GetCellTypeAndOrientation(m_coords); // get cell type and orientation
				m_end = m_mesh->m_facet_positions[m_cell_type_and_orientation][0]; // set the number of cells

				if (m_block_boundary_value == 0) {
					// our position is simply the first element in the list
					m_use_wraparound = false;
					m_pos = 1;
				}
				else {
					// the first element might not be in bounds, so advance boundary
					m_pos = 0;
					advance_boundary();
				}
			}

		public:
			FacetsIterator(const TopologicalRegularGrid2D *const mesh) : m_mesh(mesh) {}


			void begin(Vec2l const &coords) {

				// what happens for all facet iterators. set the base, get the coordinates of the cell, 
				// use those coordinates to determine if boundary
				// then use the coordinates to figure out which facet list we are on - 0 - 8 for each possible cell:
				// 0 = vertex, 1 = edge along x, 2 = edge along y, 3 = xy face, etc...
				m_coords = coords;
				m_base_index = m_mesh->coords2Cellid(coords); // record teh base

				init();
			}
			void begin(INDEX_TYPE const &cellid) {

				// what happens for all facet iterators. set the base, get the coordinates of the cell, 
				// use those coordinates to determine if boundary
				// then use the coordinates to figure out which facet list we are on - 0 - 8 for each possible cell:
				// 0 = vertex, 1 = edge along x, 2 = edge along y, 3 = xy face, etc...
				m_base_index = cellid; // record teh base
				m_mesh->cellid2Coords(cellid, m_coords); // get integer coordinates

				init();

			}
			void advance() {
				if (m_block_boundary_value == 0) {
					m_pos++;
				}
				else {
					advance_boundary();
				}
			}
			bool valid() const {
				return m_pos <= m_end;
			}
			INDEX_TYPE value() const {
				if (m_use_wraparound) {
					// pick out which axis and which direction to use in the wraparound offset list
					return m_base_index + m_mesh->m_wraparound_adjacency_offsets[m_mesh->m_facet_along_which_axis[m_cell_type_and_orientation][m_pos]];
				}
				else {
					return m_base_index + m_mesh->m_facet_positions[m_cell_type_and_orientation][m_pos];
				}
			}
		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class CofacetsIterator {
		protected:
			int m_pos, m_end;
			INDEX_TYPE m_base_index;
			bool m_use_wraparound;
			Vec2l m_coords;
			const TopologicalRegularGrid2D* const m_mesh;
			BOUNDARY_TYPE m_block_boundary_value; // whether we are on boundary of block, as laid out in memory
			int m_cell_type_and_orientation; // where to look in the cofacets list for a cell of this type (dimesion and orientation) - based on order we assigned things in 

			void advance_boundary() {
				// increment position
				m_pos++;
				// if invalid, return
				if (m_pos > m_end) return;

				// now keep advancing while we are valid and are not pointing INSIDE global mesh
				while (m_pos <= m_end) {
					// THe logic of this block is symmetric to that of facets iterator - look there
					DIM_TYPE tOffsetAxis = m_mesh->m_cofacet_along_which_axis[m_cell_type_and_orientation][m_pos] / 2;
					if (m_pos % 2 == 0) { // then we are in the NEGATIVE direction
						if (m_coords[tOffsetAxis] != 0) {
							m_use_wraparound = false;
							return;
						}
						else {
							if (m_mesh->IsPeriodic(tOffsetAxis)) { // use wraparound next
								m_use_wraparound = true;
								return;
							}
						}
					}
					else  { // so (m_pos % 2 == 1) is true so we are looking in positive direction
						if (m_coords[tOffsetAxis] != m_mesh->Extent(tOffsetAxis) - 1) {
							m_use_wraparound = false;
							return;
						}
						else {
							if (m_mesh->IsPeriodic(tOffsetAxis)) { // use wraparound next
								m_use_wraparound = true;
								return;
							}
						}
					}
					m_pos++;
				}
			}

		private:
			void init(){
				m_block_boundary_value = m_mesh->MemoryBlockBoundaryValue(m_coords); // get global boundary information
				m_cell_type_and_orientation = m_mesh->GetCellTypeAndOrientation(m_coords); // get cell type and orientation
				m_end = m_mesh->m_cofacet_positions[m_cell_type_and_orientation][0]; // set the number of cells

				if (m_block_boundary_value == 0) {
					// our position is simply the first element in the list
					m_use_wraparound = false;
					m_pos = 1;
				}
				else {
					// the first element might not be in bounds, so advance boundary
					m_pos = 0;
					advance_boundary();
				}
			}

		public:
			CofacetsIterator(const TopologicalRegularGrid2D *const m) : m_mesh(m) {}

			void begin(Vec2l const &coords) {

				// what happens for all facet iterators. set the base, get the coordinates of the cell, 
				// use those coordinates to determine if boundary
				// then use the coordinates to figure out which facet list we are on - 0 - 8 for each possible cell:
				// 0 = vertex, 1 = edge along x, 2 = edge along y, 3 = xy face, etc...
				m_coords = coords; // record teh base
				m_base_index = m_mesh->coords2Cellid(m_coords); // get integer coordinates

				init();

			}
			void begin(INDEX_TYPE const &cellid) {

				// what happens for all facet iterators. set the base, get the coordinates of the cell, 
				// use those coordinates to determine if boundary
				// then use the coordinates to figure out which facet list we are on - 0 - 8 for each possible cell:
				// 0 = vertex, 1 = edge along x, 2 = edge along y, 3 = xy face, etc...
				m_base_index = cellid; // record teh base
				m_mesh->cellid2Coords(cellid, m_coords); // get integer coordinates

				init();

			}
			void advance() {
				if (m_block_boundary_value == 0) {
					m_pos++;
				}
				else {
					advance_boundary();
				}
			}
			bool valid() const {
				return m_pos <= m_end;
			}
			INDEX_TYPE value() const {
				if (m_use_wraparound) {
					return m_base_index + m_mesh->m_wraparound_adjacency_offsets[m_mesh->m_cofacet_along_which_axis[m_cell_type_and_orientation][m_pos]];
				}
				else {
					return m_base_index + m_mesh->m_cofacet_positions[m_cell_type_and_orientation][m_pos];
				}
			}
		};


		//// also needs boundary!!! WILL USE IF rather than virutal function
		class AdjacentCellsIterator {
		protected:
			int m_pos, m_end;
			INDEX_TYPE m_base_index;
			INDEX_TYPE m_boundary_offset_index;
			Vec2l m_coords;
			const TopologicalRegularGrid2D* const m_mesh;
			BOUNDARY_TYPE m_block_boundary_value; // whether we are on boundary of block, as laid out in memory

			void advance_boundary() {
				// increment position
				m_pos++;
				// if invalid, return
				if (m_pos > m_end) return;

				// now keep advancing while we are valid and are not pointing INSIDE global mesh
				while (m_pos <= m_end) {
					// m_pos iterates over 27 neighborhood - key insight is that the neighbor index is separable
					// into the x, y, and z offsets from our base, so we compute the the offsets in each direction
					// separately. we add up those offsetsThen any call to value will simply sum up those offsets.
					//
					// offsets are calculated as follows: 
					// - tpos can be analyzed to find if the current position is -1, 0, or 1 along each axis from
					//   our current point. 
					// - for each axis (x,y,z) check, if cell at tpos is at -1 from offset, are we at coordinate 0? if so,
					//   if the mesh is periodic, get the wraparound offset, else this is out of bounds, and continue;
					int tpos = m_pos - 1;
					Vec2l m_offset_xy; // use offset for x,y,z?
					// do X -----------------------------------------------------------------
					if (tpos % 3 == 0) { // then we are in the NEGATIVE x direction
						if (m_coords[0] == 0) {
							if (m_mesh->IsPeriodic(0)) {
								m_offset_xy[0] = m_mesh->m_wraparound_adjacency_offsets[1];
							}
							else  {
								m_pos++;
								continue;
							}
						}
						else {
							m_offset_xy[0] = m_mesh->m_adjacency_offsets[1]; // just normal offset
						}
					}
					else if (tpos % 3 == 2) {
						if (m_coords[0] == m_mesh->Extent(0) - 1) {
							if (m_mesh->IsPeriodic(0)) {
								m_offset_xy[0] = m_mesh->m_wraparound_adjacency_offsets[0];
							}
							else  {
								m_pos++;
								continue;
							}
						}
						else {
							m_offset_xy[0] = m_mesh->m_adjacency_offsets[0]; // just normal offset
						}
					}
					else {
						m_offset_xy[0] = 0;
					}

					// do Y -----------------------------------------------------------------
					if ((tpos / 3) % 3 == 0) { // then we are in the NEGATIVE y direction
						if (m_coords[1] == 0) {
							if (m_mesh->IsPeriodic(1)) {
								m_offset_xy[1] = m_mesh->m_wraparound_adjacency_offsets[3];
							}
							else  {
								m_pos++;
								continue;
							}
						}
						else {
							m_offset_xy[1] = m_mesh->m_adjacency_offsets[3]; // just normal offset
						}
					}
					else if ((tpos / 3) % 3 == 2) {
						if (m_coords[1] == m_mesh->Extent(1) - 1) {
							if (m_mesh->IsPeriodic(1)) {
								m_offset_xy[1] = m_mesh->m_wraparound_adjacency_offsets[2];
							}
							else  {
								m_pos++;
								continue;
							}
						}
						else {
							m_offset_xy[1] = m_mesh->m_adjacency_offsets[2]; // just normal offset
						}
					}
					else {
						m_offset_xy[1] = 0;
					}

					// if we get this far then the point is in the block and we just add up the offsets
					m_boundary_offset_index = m_offset_xy[0] + m_offset_xy[1] ;
					return;

				}
			}
		private:
			void init() {
				m_block_boundary_value = m_mesh->MemoryBlockBoundaryValue(m_coords); // get global boundary information
				//if (m_block_boundary_value) m_domain_boundary_value = m_mesh->boundaryValue(m_coords); // also check logical
				m_end = 9; // set the number of cells

				// not boundary so we don't need to skip elements
				if (m_block_boundary_value == 0) {
					// our position is simply the first element in the list

					m_pos = 1;
				}
				else {
					// the first element might not be in bounds, so advance boundary
					m_pos = 0;
					advance_boundary();
				}
			}

			friend class TopologicalRegularGrid2D;
			bool _nearMemoryBoundary() const {
				return m_block_boundary_value != 0;
			}
			int _posOfIterator() const {
				return m_pos - 1;
			}

		public:
			AdjacentCellsIterator(const TopologicalRegularGrid2D *const mesh) : m_mesh(mesh) {}

			// this encodes both orientation AND the position in a number <= 72 (i.e. 8 * 9)
			void begin(INDEX_TYPE const &cellid) {

				// what happens for all facet iterators. set the base, get the coordinates of the cell, 
				// use those coordinates to determine if boundary
				// then use the coordinates to figure out which facet list we are on - 0 - 8 for each possible cell:
				// 0 = vertex, 1 = edge along x, 2 = edge along y, 3 = xy face, etc...
				m_base_index = cellid; // record teh base
				m_mesh->cellid2Coords(cellid, m_coords); // get integer coordinates

				init();

			}
			void begin(Vec2l const &coords) {

				// what happens for all facet iterators. set the base, get the coordinates of the cell, 
				// use those coordinates to determine if boundary
				// then use the coordinates to figure out which facet list we are on - 0 - 8 for each possible cell:
				// 0 = vertex, 1 = edge along x, 2 = edge along y, 3 = xy face, etc...
				m_coords = coords;
				m_base_index = m_mesh->coords2Cellid(m_coords);; // record teh base
				// get integer coordinates

				init();

			}
			void advance() {
				if (m_block_boundary_value == 0) {
					m_pos++;
				}
				else {
					advance_boundary();
				}
			}

			bool valid() const {
				return m_pos <= m_end;
			}

			INDEX_TYPE value() const {
				// BROKEN FOR PERIODIC
				if (m_block_boundary_value == 0) {
					return m_base_index + m_mesh->m_adjacent_cell_offsets[m_pos - 1];
				}
				else{
					return m_base_index + m_boundary_offset_index;
				}
				//}
			}

			INDEX_TYPE value(int pos) {
				if (pos == m_pos - 1) return value();
				m_pos = pos;
				this->advance();
				return value();

			}


		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class CellVerticesIterator {

		protected:
			INDEX_TYPE m_base_index;
			Vec2l m_coords;
			int m_cell_type_and_orientation; // where to look in the facets list for a cell of this type (dimesion and orientation) - based on order we assigned things in

			int m_pos, m_end;
			const TopologicalRegularGrid2D* const m_mesh;

			Vec2b m_wrap_dim;
			bool m_wrap_needed;

			//BOUNDARY_TYPE m_block_boundary_value; // whether we are on boundary of block, as laid out in memory -
			// this actually only matters if we are on the extent boundary, and the mesh is periodic
			// - if the mesh is not periodic then ALL cells will have ALL their vertices in the expected
			//   locations, if it IS periodic, then some cells will have to wrap around to get their vertices
			// 
			//BOUNDARY_TYPE m_domain_boundary_value; // whether we are on a boundary of the domain, i.e. a non-periodic extent of the domain


			void init() {
				m_cell_type_and_orientation = m_mesh->GetCellTypeAndOrientation(m_coords);
				m_end = m_mesh->m_cell_vertex_offsets[m_cell_type_and_orientation][0]; // set the number of cells 

				for (int d = 0; d < 2; d++){
					m_wrap_dim[d] = (m_coords[d] == m_mesh->Extent(d) - 1) && m_mesh->IsPeriodic(d);
				}
				m_wrap_needed = m_wrap_dim[0] || m_wrap_dim[1] ;
				m_pos = 1;

				/*
				printf("\n CellVerticesIterator initialized...\n");
				printf(" base_idx = %d, coords = (%d %d %d), end = %d\n wrap_needed = %d, wrap_dims = (%d %d %d), cell_type_and_orientation = %d \n",
				m_base_index, m_coords[0], m_coords[1], m_coords[2], m_end,
				m_wrap_needed, m_wrap_dim[0], m_wrap_dim[1], m_wrap_dim[2], m_cell_type_and_orientation);
				*/
			}

		public:
			CellVerticesIterator(const TopologicalRegularGrid2D *const mesh) : m_mesh(mesh) {}

			void begin(INDEX_TYPE const &cellid) {
				//printf("CellVerticesIterator::begin(%d)\n", cellid);
				m_base_index = cellid; // record teh base
				m_mesh->cellid2Coords(cellid, m_coords); // get integer coordinates
				init();

			}
			void begin(Vec2l const &coords) {
				//printf("CellVerticesIterator::begin((%d,%d,%d))\n", coords[0], coords[1], coords[2]);;
				m_coords = coords;
				m_base_index = m_mesh->coords2Cellid(m_coords);; // record teh base
				init();
			}
			void advance() {
				m_pos++;
			}

			bool valid() const {
				return m_pos < m_end;
			}

			bool _nearMemoryBoundary() const {
				return m_wrap_needed;
			}
			// this encodes both orientation AND the position in a number <= 72 (i.e. 8 * 9)
			int _posOfIterator() const {
				return m_cell_type_and_orientation + 8 * m_pos;
			}

			INDEX_TYPE value() const {

				return value(m_pos);

				/*
				// BROKEN FOR PERIODIC
				// ALSO BROKEN FOR COMPLEMENT MESHES because of boundary
				if (!m_wrap_needed) {
				return m_base_index + m_mesh->m_cell_vertex_offsets[m_cell_type_and_orientation][m_pos];

				}
				else {
				return boundary_value();
				}*/
			}

			// added by Harsh to directly get a certain vertex
			// fixed for periodic boundary on 01.18.2017
			INDEX_TYPE value(int pos) const {

				INDEX_TYPE voffset = m_mesh->m_cell_vertex_offsets[m_cell_type_and_orientation][pos];

				if (m_wrap_needed) {

					const Vec2b &is_wrapping = m_mesh->m_cell_vertex_is_in_positive_axis_directions[m_cell_type_and_orientation][pos];

					for (int i = 0; i < 2; i++) {

						if (m_wrap_dim[i] && is_wrapping[i]) {
							voffset += m_mesh->m_wrapping_offsets[i];
						}
					}
				}
				return m_base_index + voffset;

				/*
				// BROKEN FOR PERIODIC
				// ALSO BROKEN FOR COMPLEMENT MESHES because of boundary
				//TODO: Harsh used only periodic(0). this assumes either all dimensions are periodic or none
				if (!m_wrap_needed) {
				return m_base_index + m_mesh->m_cell_vertex_offsets[m_cell_type_and_orientation][pos];

				}
				else {
				return boundary_value();
				}*/
			}
		};




	public:

		TopologicalRegularGrid2D(RegularGrid2D* base_grid) :
			m_base_grid(base_grid) {
			m_base_grid_xy = base_grid->XY();
			this->InitializeValues();

			printf(" -- Created TopologicalRegularGrid2D [%lld %lld] = %lld cells. dcells = (%lld, %lld, %lld)\n",
				m_mesh_xy[0], m_mesh_xy[1], m_mesh_xy[2], m_num_cells,
				m_num_dcells[0], m_num_dcells[1], m_num_dcells[2]);
		}


		virtual  ~TopologicalRegularGrid2D() {
			printf("delete: TopologicalRegularGrid2D \n");
		}



		// returns an integer from 0-5 which is the "direction"
		// to look (positive x, negative x, positive y, negative y, ...)
		// to get from base_cell to other_cell
		// undefined behavior if other_cell is not a facet or cofacet of 
		// base_cell
		BYTE_TYPE Compress6NeighborOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell) {
			INDEX_TYPE t_offset = other_cell - base_cell;
			return m_offset_2_direction_map[t_offset];
		}
		INDEX_TYPE UncompressByteTo6NeighborOffset(INDEX_TYPE base_cell, BYTE_TYPE direction) {
			Vec2l t_coords;
			cellid2Coords(base_cell, t_coords);
			if (MemoryBlockBoundaryValue(t_coords)) {
				if (direction % 2 == 1) {
					// we are looking in negative direction
					if (t_coords[direction / 2] == 0) {
						return base_cell + m_wraparound_adjacency_offsets[direction];
					}
					return base_cell + m_adjacency_offsets[direction];
				}
				else {
					// we are looking in positive direction
					if (t_coords[direction / 2] == m_mesh_xy[direction / 2] - 1){
						return base_cell + m_wraparound_adjacency_offsets[direction];
					}
					return base_cell + m_adjacency_offsets[direction];
				}
			}
			else {
				return base_cell + m_adjacency_offsets[direction];
			}

		}

		// these return values from 0-27 - an offset computed by the "27" version 
		// MUST be turned into an id with the matching "27" version
		//BYTE_TYPE CompressVertexOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell, CellVerticesIterator& cell_vertices_iter) const {
		//	if (cell_vertices_iter._nearMemoryBoundary())
		//		return cell_vertices_iter._posOfIterator() + 128;

		//	INDEX_TYPE t_offset = other_cell - base_cell;
		//	return (*m_offset_2_position_vertices_map.find(t_offset)).second;
		//}
		BYTE_TYPE CompressVertexOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell) const {

			INDEX_TYPE t_offset = other_cell - base_cell;
			//if (m_offset_2_position_vertices_map.count(t_offset) == 0) printf("SHOSDFKL:JSDL:FDJKSL:DKF\\n");
			return (*m_offset_2_position_vertices_map.find(t_offset)).second;
		}



		//INDEX_TYPE UncompressByteToVertexOffset(INDEX_TYPE base_cell, BYTE_TYPE direction) const {
		//	if (direction >= 128) {
		//		CellVerticesIterator vit(this);
		//		vit.begin(base_cell);
		//		return vit.value((direction - 128) / 8); // we encoded the position as orientation + pos * 8 + 128 - so to get back the pos, we subtract 128, the (begin) set the orientation, and now we have the pos
		//	}
		//	else {
		//		return base_cell + m_adjacent_cell_offsets[direction];
		//	}
		//}
		INDEX_TYPE UncompressByteToVertexOffset(INDEX_TYPE base_cell, BYTE_TYPE direction) const {
			return base_cell + m_adjacent_cell_offsets[direction];
		}

		BYTE_TYPE CompressAdjacentOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell, AdjacentCellsIterator& cell_vertices_iter) const {
			if (cell_vertices_iter._nearMemoryBoundary())
				return cell_vertices_iter._posOfIterator() + 128;

			INDEX_TYPE t_offset = other_cell - base_cell;
			return (*m_offset_2_position_vertices_map.find(t_offset)).second;
		}

		INDEX_TYPE UncompressByteToAdjacentOffset(INDEX_TYPE base_cell, BYTE_TYPE direction) const {
			if (direction >= 128) {
				AdjacentCellsIterator vit(this);
				vit.begin(base_cell);
				return vit.value((direction - 128)); // we encoded the position as orientation + pos * 8 + 128 - so to get back the pos, we subtract 128, the (begin) set the orientation, and now we have the pos
			}
			else {
				return base_cell + m_adjacent_cell_offsets[direction];
			}
		}

		// General mesh information
		INDEX_TYPE numCellsAxis(int i) const {
			return m_mesh_xy[i];
		}

		INDEX_TYPE numCells() const {
			return m_num_cells;
		}

		INDEX_TYPE numCells(DIM_TYPE dim) const {
			return m_num_dcells[dim];
		}

		DIM_TYPE maxDim() const {
			return 2;
		}
		inline const Vec2l& XY() const { return this->m_mesh_xy; }

		// Queries regarding individual cells
		inline void cellid2Coords(INDEX_TYPE cellid, Vec2l &coords) const {
			coords[0] = cellid % m_mesh_xy[0];
			coords[1] = (cellid / m_mesh_xy[0]) % m_mesh_xy[1];
		}

		INDEX_TYPE coords2Cellid(const Vec2l &coords) const {
			return ((INDEX_TYPE)coords[0])
				+ ((INDEX_TYPE)coords[1]) * m_mesh_xy[0];
		}

		BOUNDARY_TYPE MemoryBlockBoundaryValue(const Vec2l& coords) const {
			return (BOUNDARY_TYPE)((coords[0] == 0 || coords[0] == m_mesh_xy[0] - 1) +
				(coords[1] == 0 || coords[1] == m_mesh_xy[1] - 1)) ;
		}

		BOUNDARY_TYPE boundaryValue(const Vec2l& coords) const {

			return (BOUNDARY_TYPE)(
				(!m_mesh_periodic[0] && ((coords[0] == 0) || (coords[0] == m_mesh_xy[0] - 1))) +
				(!m_mesh_periodic[1] && ((coords[1] == 0) || (coords[1] == m_mesh_xy[1] - 1))) 
				);

		}
		BOUNDARY_TYPE boundaryValue(INDEX_TYPE cellid) const {
			Vec2l coords;
			cellid2Coords(cellid, coords);
			return boundaryValue(coords);
		}



		// need a boundary value for coords
		// and a dimension for coords
		DIM_TYPE dimension(const Vec2l& coords) const {
			return (coords[0] % 2 + coords[1] % 2);
		}
		DIM_TYPE dimension(const INDEX_TYPE& cellid) const {
			Vec2l coords;
			cellid2Coords(cellid, coords);
			return dimension(coords);
		}




		void centroid(INDEX_TYPE cellid, Vec2d &coords) const {
			Vec2l icoords;
			cellid2Coords(cellid, icoords);
			coords = icoords;
		}


		INDEX_TYPE VertexNumberFromCellID(const INDEX_TYPE cellid) const {
			Vec2l icoords;
			cellid2Coords(cellid, icoords);
			return icoords[0] / 2 + (icoords[1] / 2) * m_base_grid_xy[0] ;
		}

		INDEX_TYPE CellIDFromVertexNumber(const INDEX_TYPE vertex_number) const {
			//Vec2l icoords;
			//this->m_base_grid->XYZ3d(vertex_number);
			Vec2l icoords = this->m_base_grid->XY2d(vertex_number);
			return coords2Cellid(icoords * 2);
		}

	};



}


#endif
