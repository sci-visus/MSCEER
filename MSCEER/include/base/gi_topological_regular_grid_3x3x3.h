/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef TOPOLOGICAL_REGULAR_GRID_3x3x3_H
#define TOPOLOGICAL_REGULAR_GRID_3x3x3_H


#include <unordered_map>
#include <set>

#include "base/gi_basic_types.h"
#include "base/gi_vectors.h"




namespace GInt {





	class Explicit3x3x3SmallRegularGrid {
	public:
		void display() const {

		}

		INDEX_TYPE get27NeighborOffset(int pos) const {
			return m_adjacent_cell_offsets[pos];
		}
	protected:
		static const int VertNumFromCellidArray[125];
		static const int DimArray[125];
		static const int m_adjacent_cell_offsets[27];
	public:
		class AllCellsIterator {
		protected:


		public:
			AllCellsIterator(Explicit3x3x3SmallRegularGrid* m) {}

			AllCellsIterator(Explicit3x3x3SmallRegularGrid* m, INDEX_TYPE start, INDEX_TYPE end) {}

			void begin() {
				printf("NOT SUPPORTED\n");
			}
			void advance() {
				printf("NOT SUPPORTED\n");
			}
			bool valid() const {
				return false;
			}
			INDEX_TYPE value() const {
				return 0;
			}
		};

		class DCellsIterator {
		protected:
		public:
			DCellsIterator(Explicit3x3x3SmallRegularGrid* mesh, DIM_TYPE dim) {}
			DCellsIterator(Explicit3x3x3SmallRegularGrid* mesh, DIM_TYPE dim, INDEX_TYPE start, INDEX_TYPE end)  {}

			void begin() {

			}
			void advance() {
			}
			bool valid() const {
				return false;
			}
			INDEX_TYPE value() const {
				return 0;
			}
		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class FacetsIterator {
		protected:
			static const int FacetsIdList[125][7];
			INDEX_TYPE m_base_index;
			int m_pos;
			int m_end;

		public:
			FacetsIterator(const Explicit3x3x3SmallRegularGrid *const mesh) {}


			void begin(Vec3l const &coords) {

				printf("Explicit3x3x3SmallRegularGrid::FacetsIterator::begin(Vec3l const &coords) NOT SUPPORTED\n");
			}
			void begin(INDEX_TYPE const &cellid) {

				// what happens for all facet iterators. set the base, get the coordinates of the cell, 
				// use those coordinates to determine if boundary
				// then use the coordinates to figure out which facet list we are on - 0 - 8 for each possible cell:
				// 0 = vertex, 1 = edge along x, 2 = edge along y, 3 = xy face, etc...
				m_base_index = cellid; // record teh base
				m_pos = 1;
				m_end = 1 + FacetsIdList[cellid][0];
			}
			void advance() {
				m_pos++;
			}
			bool valid() const {
				return m_pos <= m_end;
			}
			INDEX_TYPE value() const {
				return FacetsIdList[m_base_index][m_pos];
			}
		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class CofacetsIterator {
		protected:
			static const int CoFacetsIdList[125][7];
			INDEX_TYPE m_base_index;
			int m_pos;
			int m_end;

		public:
			CofacetsIterator(const Explicit3x3x3SmallRegularGrid *const mesh) {}


			void begin(Vec3l const &coords) {

				printf("Explicit3x3x3SmallRegularGrid::CofacetsIterator::begin(Vec3l const &coords) NOT SUPPORTED\n");
			}
			void begin(INDEX_TYPE const &cellid) {

				// what happens for all facet iterators. set the base, get the coordinates of the cell, 
				// use those coordinates to determine if boundary
				// then use the coordinates to figure out which facet list we are on - 0 - 8 for each possible cell:
				// 0 = vertex, 1 = edge along x, 2 = edge along y, 3 = xy face, etc...
				m_base_index = cellid; // record teh base
				m_pos = 1;
				m_end = 1 + CoFacetsIdList[cellid][0];
			}
			void advance() {
				m_pos++;
			}
			bool valid() const {
				return m_pos <= m_end;
			}
			INDEX_TYPE value() const {
				return CoFacetsIdList[m_base_index][m_pos];
			}
		};


		//// also needs boundary!!! WILL USE IF rather than virutal function
		class AdjacentCellsIterator {
		protected:
		protected:
			

			INDEX_TYPE m_base_index;
			int m_pos;
			const Explicit3x3x3SmallRegularGrid * mMesh;

		public:
			AdjacentCellsIterator(const Explicit3x3x3SmallRegularGrid *const mesh) : mMesh(mesh) {}


			void begin(Vec3l const &coords) {

				printf("Explicit3x3x3SmallRegularGrid::AdjacentCellsIterator::begin(Vec3l const &coords) NOT SUPPORTED\n");
			}
			void begin(INDEX_TYPE const &cellid) {
				m_base_index = cellid; // record teh base
				m_pos = 0;			
			}
			void advance() {
				m_pos++;
			}
			bool valid() const {
				return m_pos <27;
			}
			INDEX_TYPE value() const {
				return mMesh->get27NeighborOffset(m_pos) + m_base_index;
			}
		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class CellVerticesIterator {

		protected:
		public:
			CellVerticesIterator(const Explicit3x3x3SmallRegularGrid *const mesh)  {}

			void begin(INDEX_TYPE const &cellid) {

			}
			void begin(Vec3l const &coords) {

			}
			void advance() {

			}

			bool valid() const {
				return false;
			}

			INDEX_TYPE value() const {

				return -1;
			}

			// added by Harsh to directly get a certain vertex
			// fixed for periodic boundary on 01.18.2017
			INDEX_TYPE value(int pos) const {

				return -1;
			}
		};

	public:

		Explicit3x3x3SmallRegularGrid() {			
		}


		virtual  ~Explicit3x3x3SmallRegularGrid() {
			printf("delete: Explicit3x3x3SmallRegularGrid \n");
		}



		//// returns an integer from 0-5 which is the "direction"
		//// to look (positive x, negative x, positive y, negative y, ...)
		//// to get from base_cell to other_cell
		//// undefined behavior if other_cell is not a facet or cofacet of 
		//// base_cell
		//BYTE_TYPE Compress6NeighborOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell) {
		//	printf("NOT SUPPOERTED\n");
		//	return 0;
		//}
		//INDEX_TYPE UncompressByteTo6NeighborOffset(INDEX_TYPE base_cell, BYTE_TYPE direction) {

		//}

		//// these return values from 0-27 - an offset computed by the "27" version 
		//// MUST be turned into an id with the matching "27" version
		////BYTE_TYPE CompressVertexOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell, CellVerticesIterator& cell_vertices_iter) const {
		////	if (cell_vertices_iter._nearMemoryBoundary())
		////		return cell_vertices_iter._posOfIterator() + 128;

		////	INDEX_TYPE t_offset = other_cell - base_cell;
		////	return (*m_offset_2_position_vertices_map.find(t_offset)).second;
		////}
		//BYTE_TYPE CompressVertexOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell) const {

		//	INDEX_TYPE t_offset = other_cell - base_cell;
		//	//if (m_offset_2_position_vertices_map.count(t_offset) == 0) printf("SHOSDFKL:JSDL:FDJKSL:DKF\\n");
		//	return (*m_offset_2_position_vertices_map.find(t_offset)).second;
		//}



		////INDEX_TYPE UncompressByteToVertexOffset(INDEX_TYPE base_cell, BYTE_TYPE direction) const {
		////	if (direction >= 128) {
		////		CellVerticesIterator vit(this);
		////		vit.begin(base_cell);
		////		return vit.value((direction - 128) / 8); // we encoded the position as orientation + pos * 8 + 128 - so to get back the pos, we subtract 128, the (begin) set the orientation, and now we have the pos
		////	}
		////	else {
		////		return base_cell + m_adjacent_cell_offsets[direction];
		////	}
		////}
		INDEX_TYPE UncompressByteToVertexOffset(INDEX_TYPE base_cell, BYTE_TYPE direction) const {
			return base_cell + m_adjacent_cell_offsets[direction];
		}

		//BYTE_TYPE CompressAdjacentOffsetToByte(INDEX_TYPE base_cell, INDEX_TYPE other_cell, AdjacentCellsIterator& cell_vertices_iter) const {
		//	if (cell_vertices_iter._nearMemoryBoundary())
		//		return cell_vertices_iter._posOfIterator() + 128;

		//	INDEX_TYPE t_offset = other_cell - base_cell;
		//	return (*m_offset_2_position_vertices_map.find(t_offset)).second;
		//}

		//INDEX_TYPE UncompressByteToAdjacentOffset(INDEX_TYPE base_cell, BYTE_TYPE direction) const {
		//	if (direction >= 128) {
		//		AdjacentCellsIterator vit(this);
		//		vit.begin(base_cell);
		//		return vit.value((direction - 128)); // we encoded the position as orientation + pos * 8 + 128 - so to get back the pos, we subtract 128, the (begin) set the orientation, and now we have the pos
		//	}
		//	else {
		//		return base_cell + m_adjacent_cell_offsets[direction];
		//	}
		//}

		// General mesh information
		INDEX_TYPE numCellsAxis(int i) const {
			return 5;
		}

		INDEX_TYPE numCells() const {
			return 125;
		}

		//INDEX_TYPE numCells(DIM_TYPE dim) const {
		//	return m_num_dcells[dim];
		//}

		DIM_TYPE maxDim() const {
			return 3;
		}
		//inline const Vec3l& XYZ() const { return this->m_mesh_xyz; }

		//// Queries regarding individual cells
		//inline void cellid2Coords(INDEX_TYPE cellid, Vec3l &coords) const {
		//	coords[0] = cellid % m_mesh_xyz[0];
		//	coords[1] = (cellid / m_mesh_xyz[0]) % m_mesh_xyz[1];
		//	coords[2] = (cellid / (m_mesh_xyz[0] * m_mesh_xyz[1]));
		//}

		INDEX_TYPE coords2Cellid(const Vec3l &coords) const {
			return ((INDEX_TYPE)coords[0])
				+ ((INDEX_TYPE)coords[1]) * 5
				+ ((INDEX_TYPE)coords[2]) * 25;
		}

		//BOUNDARY_TYPE MemoryBlockBoundaryValue(const Vec3l& coords) const {
		//	return (BOUNDARY_TYPE)((coords[0] == 0 || coords[0] == m_mesh_xyz[0] - 1) +
		//		(coords[1] == 0 || coords[1] == m_mesh_xyz[1] - 1) +
		//		(coords[2] == 0 || coords[2] == m_mesh_xyz[2] - 1)
		//		);

		//	/*return (BOUNDARY_TYPE)((coords[0] == 0) || (coords[0] == m_mesh_xyz[0] - 1)) +
		//	((coords[1] == 0) || (coords[1] == m_mesh_xyz[1] - 1)) +
		//	((coords[2] == 0) || (coords[2] == m_mesh_xyz[2] - 1));*/
		//}

		BOUNDARY_TYPE boundaryValue(const Vec3l& coords) const {
			return 0;
		}
		BOUNDARY_TYPE boundaryValue(INDEX_TYPE cellid) const {
			return 0;
		}



		//// need a boundary value for coords
		//// and a dimension for coords
		//DIM_TYPE dimension(const Vec3l& coords) const {
		//	return (coords[0] % 2 + coords[1] % 2 + coords[2] % 2);
		//}
		DIM_TYPE dimension(const INDEX_TYPE& cellid) const {
			return DimArray[cellid];
		}




		//void centroid(INDEX_TYPE cellid, Vec3d &coords) const {
		//	Vec3l icoords;
		//	cellid2Coords(cellid, icoords);
		//	coords = icoords;
		//}

		//INDEX_TYPE Extent(DIM_TYPE i) const {
		//return mExtentList[i*2];
		//}


		INDEX_TYPE VertexNumberFromCellID(const INDEX_TYPE cellid) const {
			return VertNumFromCellidArray[cellid];
		}

		//INDEX_TYPE CellIDFromVertexNumber(const INDEX_TYPE vertex_number) const {
		//	//Vec3l icoords;
		//	//this->m_base_grid->XYZ3d(vertex_number);
		//	Vec3l icoords = this->m_base_grid->XYZ3d(vertex_number);
		//	return coords2Cellid(icoords * 2);
		//}

		//int num_faces(INDEX_TYPE cellid) const {

		//    int cnt = 0;
		//    FacetsIterator fiter (this);
		//    for(fiter.begin (cellid); fiter.valid (); fiter.advance ()) {
		//        cnt++;
		//    }
		//    return cnt;
		//}

		//// check if a is a coface of b
		//bool is_cofacet_of(INDEX_TYPE a, INDEX_TYPE b)  const {

		//    CofacetsIterator fiter (this);
		//    for(fiter.begin (b); fiter.valid (); fiter.advance ()) {
		//        if(a == fiter.value ()) {
		//            return true;
		//        }
		//    }
		//    return false;
		//}


	};


}


#endif
