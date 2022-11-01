/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef TOPOLOGICAL_COMPLEMENT_GRID_H
#define TOPOLOGICAL_COMPLEMENT_GRID_H


#include <map>
#include <set>
#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_regular_grid.h"



namespace GInt {


	// TopologicalRegularGrid3D is a specialized class to work with the connectivity of a
	// regular 3d grid - with quantities derived from an input regular 3d grid. 
	// This class handles no storage on its own - it simply
	// provides a mechanism for iterating over cells, iterating over facets and cofacets of cells
	// and answering questions about a cell, such as its dimension, whether or not it sits on a 
	// boundary, what its coordinates are, etc. 
	//
	// This class uses a cell index numbering to make neighborhood queries fast - simply reduces
	// to adding a pre-computed offset to the current cell id. E.g. For a nonperiodic grid with X*Y*Z values,
	// the indices generated by this class uses a grid of size (2X-1)*(2Y-1)*(2Z-1), basically
	// representing every cell (vertex, edge, quad, hex) with a unique id, such that facets/cofacets
	// of a cell id can be computed by simply adding an offset to the  cell id. 

	// NOTE: we use unsigned long longs for the axes, since we do not want any arithmetic we do with them
	// to be accidentally truncated to 32-bit integers.
	template <class BaseMeshType> 
	class TopologicalComplementGrid {
	protected:
		BaseMeshType* mBaseMesh;
	
	public:
            void display() const {

				mBaseMesh->display();
            }


	protected:

		INT_TYPE Extent(DIM_TYPE val) const {
			return mBaseMesh->Extent(val);
		}

        int GetCellTypeAndOrientation(Vec3l p) const {
			return mBaseMesh->GetCellTypeAndOrientation(p);
		}

  /*      inline INDEX_TYPE dot(const Vec3i &w, const Vec3l &v) {
            return w[0]*v[0] + w[1]*v[1] + w[2]*v[2];
        }*/

        bool IsPeriodic(int axis) const {
			return mBaseMesh->IsPeriodic(axis);
		}

	public:
		class AllCellsIterator {
		protected:
			BaseMeshType::AllCellsIterator mit;
		public:
			AllCellsIterator(TopologicalComplementGrid* m) 	{
				mit = BaseMeshType::AllCellsIterator(m->mBaseMesh);
			}

			AllCellsIterator(BaseMeshType* m, INDEX_TYPE start, INDEX_TYPE end)  {
				mit = BaseMeshType::AllCellsIterator(m->mBaseMesh, start, end);
			}

			void begin() {
				mit.begin();
			}
			void advance() {
				mit.advance();
			}
			bool valid() const {
				return mit.valid();
			}
			INDEX_TYPE value() const {
				return mit.value();
			}
		};

		class DCellsIterator {
		protected:
			BaseMeshType::DCellsIterator mit;
		public:
			DCellsIterator(BaseMeshType* mesh, DIM_TYPE dim) {
				mit = BaseMeshType::DCellsIterator(mesh, mesh->maxDim() - dim);
			}
			DCellsIterator(BaseMeshType* mesh, DIM_TYPE dim, INDEX_TYPE start, INDEX_TYPE end) {
				mit = BaseMeshType::DCellsIterator(mesh, mesh->maxDim() - dim, start, end);
			}

			void begin() {
				mit.begin();
			}
			void advance() {
				mit.advance();
			}
			bool valid() const {
				return mit.valid();
			}
			INDEX_TYPE value() const {
				return mit.value();
			}
		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class FacetsIterator {
		protected:
			BaseMeshType::CofacetsIterator mit;

		public:
			FacetsIterator(const BaseMeshType *const mesh) {
				mit = BaseMeshType::CofacetsIterator(mesh);
			}


			void begin(Vec3l const &coords) {
				mit.begin(coords);
			}
			void begin(INDEX_TYPE const &cellid) {
				mit.begin(cellid);
			}
			void advance() {
				mit.advance();
			}
			bool valid() const {
				return mit.valid()
			}
			INDEX_TYPE value() const {
				return mit.value();
			}
		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class CofacetsIterator {
		protected:
			BaseMeshType::FacetsIterator mit;
		public:
			CofacetsIterator(const BaseMeshType *const m) : m_mesh(m) {
						mit = BaseMeshType::FacetsIterator(mesh); 
			}

			void begin(Vec3l const &coords) {

				// what happens for all facet iterators. set the base, get the coordinates of the cell, 
				// use those coordinates to determine if boundary
				// then use the coordinates to figure out which facet list we are on - 0 - 8 for each possible cell:
				// 0 = vertex, 1 = edge along x, 2 = edge along y, 3 = xy face, etc...
				m_coords = coords; // record teh base
				m_base_index = m_mesh->coords2Cellid(m_coords); // get integer coordinates

                init();

			}

			void begin(Vec3l const &coords) {
				mit.begin(coords);
			}
			void begin(INDEX_TYPE const &cellid) {
				mit.begin(cellid);
			}
			void advance() {
				mit.advance();
			}
			bool valid() const {
				return mit.valid()
			}
			INDEX_TYPE value() const {
				return mit.value();
			}
		};


		//// also needs boundary!!! WILL USE IF rather than virutal function
		class AdjacentCellsIterator {
		protected:
			BaseMeshType::AdjacentCellsIterator mit;

		public:
			AdjacentCellsIterator(const BaseMeshType *const mesh) {
				mit = BaseMeshType::AdjacentCellsIterator(mesh);
			}

			void begin(Vec3l const &coords) {
				mit.begin(coords);
			}
			void begin(INDEX_TYPE const &cellid) {
				mit.begin(cellid);
			}
			void advance() {
				mit.advance();
			}
			bool valid() const {
				return mit.valid()
			}
			INDEX_TYPE value() const {
				return mit.value();
			}

			INDEX_TYPE value(int pos) {
				return mit.value(pos);
			}
		};

		// CHECK THIS FOR COMPLEMENT
		//// also needs boundary!!! WILL USE IF rather than virutal function
		class CellVerticesIterator {

		protected:
            INDEX_TYPE m_base_index;
            Vec3l m_coords;
            int m_cell_type_and_orientation; // where to look in the facets list for a cell of this type (dimesion and orientation) - based on order we assigned things in

            int m_pos, m_end;
            const TopologicalRegularGrid3D* const m_mesh;

            Vec3b m_wrap_dim;
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

                for(int d = 0; d < 3; d++){
                    m_wrap_dim[d] = (m_coords[d] == m_mesh->Extent(d)-1) && m_mesh->IsPeriodic(d);
                }
                m_wrap_needed = m_wrap_dim[0] || m_wrap_dim[1] || m_wrap_dim[2];
                m_pos = 1;

                /*
                printf("\n CellVerticesIterator initialized...\n");
                printf(" base_idx = %d, coords = (%d %d %d), end = %d\n wrap_needed = %d, wrap_dims = (%d %d %d), cell_type_and_orientation = %d \n",
                       m_base_index, m_coords[0], m_coords[1], m_coords[2], m_end,
                       m_wrap_needed, m_wrap_dim[0], m_wrap_dim[1], m_wrap_dim[2], m_cell_type_and_orientation);
                */
            }

#if 0
            INDEX_TYPE wraparound_index() {

            }


            INDEX_TYPE boundary_value() const {

                INDEX_TYPE tempval = m_base_index;
                for (int i = 0; i < 3; i++) {
                    if (m_coords[i] == m_mesh->Extent(i) - 1
                            && m_mesh->m_cell_vertex_is_in_positive_axis_directions[m_cell_type_and_orientation][m_pos][i]) {
                        tempval += m_mesh->m_wraparound_adjacency_offsets[i * 2];
                    }
                    else {
                        tempval += m_mesh->m_adjacency_offsets[i * 2];
                    }
                }
                return tempval;
            }
#endif
		public:
			CellVerticesIterator(const TopologicalRegularGrid3D *const mesh) : m_mesh(mesh) {}

			void begin(INDEX_TYPE const &cellid) {
                //printf("CellVerticesIterator::begin(%d)\n", cellid);
				m_base_index = cellid; // record teh base
				m_mesh->cellid2Coords(cellid, m_coords); // get integer coordinates
				init();

			}
			void begin(Vec3l const &coords) {
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

                if(m_wrap_needed) {

                    const Vec3b &is_wrapping = m_mesh->m_cell_vertex_is_in_positive_axis_directions[m_cell_type_and_orientation][pos];

                    for(int i = 0; i < 3; i++) {

                        if(m_wrap_dim[i] && is_wrapping[i]) {
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

		TopologicalComplementGrid(BaseMeshType* base_grid) : mBaseMesh(base_grid) {}
           

		virtual  ~TopologicalComplementGrid() {
			printf("delete: TopologicalComplementGrid \n");
		}

		BYTE_TYPE Compress6NeighborOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell) {
			return mBaseMesh->Compress6NeighborOffsetToByte(base_cell, other_cell);
		}
		INDEX_TYPE UncompressByteTo6NeighborOffset(INDEX_TYPE base_cell, BYTE_TYPE direction) {
			return mBaseMesh->UncompressByteTo6NeighborOffset(base_cell, direction);
		}

		
		// CHECK THIS FOR COMPLEMENT
		BYTE_TYPE CompressVertexOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell, CellVerticesIterator& cell_vertices_iter) const {
			if (cell_vertices_iter._nearMemoryBoundary())
				return cell_vertices_iter._posOfIterator() + 128;
			
			INDEX_TYPE t_offset = other_cell - base_cell;
				return (*m_offset_2_position_vertices_map.find(t_offset)).second;
		}
		


		// CHECK THIS FOR COMPLEMENT
		INDEX_TYPE UncompressByteToVertexOffset(INDEX_TYPE base_cell, BYTE_TYPE direction) const {
			if (direction >= 128) {
				CellVerticesIterator vit(this);
				vit.begin(base_cell);
				return vit.value((direction - 128) / 8); // we encoded the position as orientation + pos * 8 + 128 - so to get back the pos, we subtract 128, the (begin) set the orientation, and now we have the pos
			}
			else {
				return base_cell + m_adjacent_cell_offsets[direction];
			}
		}
		
		// CHECK THIS FOR COMPLEMENT
		BYTE_TYPE CompressAdjacentOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell, AdjacentCellsIterator& cell_vertices_iter) const {
			if (cell_vertices_iter._nearMemoryBoundary())
				return cell_vertices_iter._posOfIterator() + 128;

			INDEX_TYPE t_offset = other_cell - base_cell;
			return (*m_offset_2_position_vertices_map.find(t_offset)).second;
		}

		// CHECK THIS FOR COMPLEMENT
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
			return mBaseMesh->numCellsAxis(i);
		}

		INDEX_TYPE numCells() const {
			return mBaseMesh->numCells();
		}

		INDEX_TYPE numCells(DIM_TYPE dim) const {
			return mBaseMesh->numCells(mBaseMesh->maxDim() - dim);
		}

		DIM_TYPE maxDim() const {
			return mBaseMesh->maxDim();
		}


		// Queries regarding individual cells
		inline void cellid2Coords(INDEX_TYPE cellid, Vec3l &coords) const {
			mBaseMesh->cellid2Coords(cellid, coords);
		}

        INDEX_TYPE coords2Cellid(const Vec3l &coords) const {
			return mBaseMesh->coords2Cellid(coords);
		}

        BOUNDARY_TYPE MemoryBlockBoundaryValue(const Vec3l& coords) const {
			return mBaseMesh->MemoryBlockBoundaryValue(coords);
		}

		BOUNDARY_TYPE boundaryValue(const Vec3l& coords) const {
			return mBaseMesh->boundaryValue(coords);
		}
		BOUNDARY_TYPE boundaryValue(INDEX_TYPE cellid) const {
			return mBaseMesh->boundaryValue(cellid);
		}



		// need a boundary value for coords
		// and a dimension for coords
		DIM_TYPE dimension(const Vec3l& coords) const {
			return mBaseMesh->maxDim() - mBaseMesh->dimension(coords);
		}
		DIM_TYPE dimension(const INDEX_TYPE& cellid) const {
			return mBaseMesh->maxDim() - mBaseMesh->dimension(cellid);
		}




		void centroid(INDEX_TYPE cellid, Vec3d &coords) const {
			mBaseMesh->centroid(cellid, coords);
		}

		// CHECK THIS FOR COMPLEMENT
		INDEX_TYPE VertexNumberFromCellID(const INDEX_TYPE cellid) const {
			Vec3l icoords;
			cellid2Coords(cellid, icoords);
			return icoords[0] / 2 + (icoords[1] / 2) * m_base_grid_xyz[0] + (icoords[2] / 2) * m_base_grid_xyz[0] * m_base_grid_xyz[1];
		}
		// CHECK THIS FOR COMPLEMENT
		INDEX_TYPE CellIDFromVertexNumber(const INDEX_TYPE vertex_number) const {
            //Vec3l icoords;
            //this->m_base_grid->XYZ3d(vertex_number);
            Vec3l icoords = this->m_base_grid->XYZ3d(vertex_number);
			return coords2Cellid(icoords*2);
		}
		
 
	};
}


#endif
