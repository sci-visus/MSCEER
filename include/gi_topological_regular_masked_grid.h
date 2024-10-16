#ifndef TOPOLOGICAL_REGULAR_MASKED_GRID_H
#define TOPOLOGICAL_REGULAR_MASKED_GRID_H


#include <map>
#include <set>
#include <unordered_set>
#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_regular_grid.h"
#include "gi_topological_regular_grid.h"
#include "gi_labeling.h"

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
	class CellTester {
	public:
		virtual bool TestCell(INDEX_TYPE id) {
			return true;
		}
	};

	class CellTesterDefaultTrue : public CellTester {
	public:
		virtual bool TestCell(INDEX_TYPE id) {
			return true;
		}
	};
	class CellTesterLaberInput : public CellTester {
	protected:
		DenseLabeling<char>* m_labeling;
	public:
		void SetLabeling(DenseLabeling<char>* inp) { m_labeling = inp; }
		virtual bool TestCell(INDEX_TYPE id) {
			return m_labeling->GetLabel(id) == 1;
		}
	};

	class CellTesterSparse : public CellTester {
	protected:
		std::unordered_set<INDEX_TYPE> m_cells;
	public:
		void AddCell(INDEX_TYPE id) { m_cells.insert(id); }
		void Reset() { m_cells.clear(); }
		virtual bool TestCell(INDEX_TYPE id) {
			return m_cells.count(id);
		}
	};

	
	template<typename LABEL_TYPE>
	class CellTesterSparseLabeled : public CellTester {
	protected:
		SparseLabeling<LABEL_TYPE>* m_labeling;
	public:
		CellTesterSparseLabeled(SparseLabeling<LABEL_TYPE>* labels) : m_labeling(labels) {}
		//void SetUnderlyingLabeling(SparseLabeling<LABEL_TYPE>* labels) {
		//	m_labeling = labels;
		//}
		virtual bool TestCell(INDEX_TYPE id) {
			return m_labeling->Has(id);
		}
	};

	class TopologicalRegularMaskedGrid  : virtual public TopologicalRegularGrid3D {
	protected:
		CellTester* m_tester;
	
	public:
		bool InMesh(INDEX_TYPE id) const {
			return m_tester->TestCell(id);
		}
	public:
		class AllCellsIterator : public TopologicalRegularGrid3D::AllCellsIterator {
		protected:
			const TopologicalRegularMaskedGrid* const m_mesh;
		public:
			AllCellsIterator(TopologicalRegularMaskedGrid* m) :
				TopologicalRegularGrid3D::AllCellsIterator(m), m_mesh(m) {}

			AllCellsIterator(TopologicalRegularMaskedGrid* m, INDEX_TYPE start, INDEX_TYPE end) :
				TopologicalRegularGrid3D::AllCellsIterator(m, start, end), m_mesh(m) {}

			void begin() {
				TopologicalRegularGrid3D::AllCellsIterator::begin();
				while (TopologicalRegularGrid3D::AllCellsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::AllCellsIterator::value())) 
					TopologicalRegularGrid3D::AllCellsIterator::advance();
			}
			void advance() {
				TopologicalRegularGrid3D::AllCellsIterator::advance();
				while (TopologicalRegularGrid3D::AllCellsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::AllCellsIterator::value()))
					TopologicalRegularGrid3D::AllCellsIterator::advance(); 
			}
			bool valid() const {
				return TopologicalRegularGrid3D::AllCellsIterator::valid();
			}
			INDEX_TYPE value() const {
				return TopologicalRegularGrid3D::AllCellsIterator::value();
			}
		};

		class DCellsIterator : public TopologicalRegularGrid3D::DCellsIterator {
		protected:
			const TopologicalRegularMaskedGrid* const m_mesh;
		public:
			DCellsIterator(TopologicalRegularMaskedGrid* mesh, DIM_TYPE dim) :
				TopologicalRegularGrid3D::DCellsIterator(mesh, dim), m_mesh(mesh) {}
			DCellsIterator(TopologicalRegularMaskedGrid* mesh, DIM_TYPE dim, INDEX_TYPE start, INDEX_TYPE end) :
				TopologicalRegularGrid3D::DCellsIterator(mesh, dim, start, end), m_mesh(mesh) {}

			void begin() {
				TopologicalRegularGrid3D::DCellsIterator::begin();
				while (TopologicalRegularGrid3D::DCellsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::DCellsIterator::value()))
					TopologicalRegularGrid3D::DCellsIterator::advance();
			}
			void advance() {
				TopologicalRegularGrid3D::DCellsIterator::advance();
				while (TopologicalRegularGrid3D::DCellsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::DCellsIterator::value()))
					TopologicalRegularGrid3D::DCellsIterator::advance();
			}
			bool valid() const {
				return TopologicalRegularGrid3D::DCellsIterator::valid();
			}
			INDEX_TYPE value() const {
				return TopologicalRegularGrid3D::DCellsIterator::value();
			}
		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class FacetsIterator : public TopologicalRegularGrid3D::FacetsIterator {
		protected:
            const TopologicalRegularMaskedGrid* const m_mesh;


 

		public:
            FacetsIterator(const TopologicalRegularMaskedGrid *const mesh) :
				TopologicalRegularGrid3D::FacetsIterator(mesh), m_mesh(mesh) {}


			void begin(Vec3l const &coords) {
				TopologicalRegularGrid3D::FacetsIterator::begin(coords);
				while (TopologicalRegularGrid3D::FacetsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::FacetsIterator::value()))
					TopologicalRegularGrid3D::FacetsIterator::advance();
			}
			void begin(INDEX_TYPE const &cellid) {
				TopologicalRegularGrid3D::FacetsIterator::begin(cellid);
				while (TopologicalRegularGrid3D::FacetsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::FacetsIterator::value()))
					TopologicalRegularGrid3D::FacetsIterator::advance();
			}
			void advance() {
				TopologicalRegularGrid3D::FacetsIterator::advance();
				while (TopologicalRegularGrid3D::FacetsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::FacetsIterator::value()))
					TopologicalRegularGrid3D::FacetsIterator::advance();
			}
			bool valid() const {
				return TopologicalRegularGrid3D::FacetsIterator::valid();
			}
			INDEX_TYPE value() const {
				return TopologicalRegularGrid3D::FacetsIterator::value();
			}
		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class CofacetsIterator : public TopologicalRegularGrid3D::CofacetsIterator {
		public:
            const TopologicalRegularMaskedGrid* const m_mesh;
			CofacetsIterator(TopologicalRegularMaskedGrid* mesh) :
			TopologicalRegularGrid3D::CofacetsIterator(mesh), m_mesh(mesh) {}


			void begin(Vec3l const &coords) {
				TopologicalRegularGrid3D::CofacetsIterator::begin(coords);
				while (TopologicalRegularGrid3D::CofacetsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::CofacetsIterator::value()))
					TopologicalRegularGrid3D::CofacetsIterator::advance();
			}
			void begin(INDEX_TYPE const &cellid) {
				TopologicalRegularGrid3D::CofacetsIterator::begin(cellid);
				while (TopologicalRegularGrid3D::CofacetsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::CofacetsIterator::value()))
					TopologicalRegularGrid3D::CofacetsIterator::advance();
			}
			void advance() {
				TopologicalRegularGrid3D::CofacetsIterator::advance();
				while (TopologicalRegularGrid3D::CofacetsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::CofacetsIterator::value()))
					TopologicalRegularGrid3D::CofacetsIterator::advance();
			}
			bool valid() const {
				return TopologicalRegularGrid3D::CofacetsIterator::valid();
			}
			INDEX_TYPE value() const {
				return TopologicalRegularGrid3D::CofacetsIterator::value();
			}
		};


		//// also needs boundary!!! WILL USE IF rather than virutal function
		class AdjacentCellsIterator : public TopologicalRegularGrid3D::AdjacentCellsIterator {
		protected:
			const TopologicalRegularMaskedGrid* const m_mesh;

		public:
			AdjacentCellsIterator(const TopologicalRegularMaskedGrid *const mesh) : 
			TopologicalRegularGrid3D::AdjacentCellsIterator(mesh), m_mesh(mesh) {}


			void begin(Vec3l const &coords) {
				TopologicalRegularGrid3D::AdjacentCellsIterator::begin(coords);
				while (TopologicalRegularGrid3D::AdjacentCellsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::AdjacentCellsIterator::value()))
					TopologicalRegularGrid3D::AdjacentCellsIterator::advance();
			}
			void begin(INDEX_TYPE const &cellid) {
				TopologicalRegularGrid3D::AdjacentCellsIterator::begin(cellid);
				while (TopologicalRegularGrid3D::AdjacentCellsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::AdjacentCellsIterator::value()))
					TopologicalRegularGrid3D::AdjacentCellsIterator::advance();
			}
			void advance() {
				TopologicalRegularGrid3D::AdjacentCellsIterator::advance();
				while (TopologicalRegularGrid3D::AdjacentCellsIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::AdjacentCellsIterator::value()))
					TopologicalRegularGrid3D::AdjacentCellsIterator::advance();
			}
			bool valid() const {
				return TopologicalRegularGrid3D::AdjacentCellsIterator::valid();
			}
			INDEX_TYPE value() const {
				return TopologicalRegularGrid3D::AdjacentCellsIterator::value();
			}
		};

		//// also needs boundary!!! WILL USE IF rather than virutal function
		class CellVerticesIterator : public TopologicalRegularGrid3D::CellVerticesIterator {

		protected:
          const TopologicalRegularMaskedGrid* const m_mesh;
		public:
			CellVerticesIterator(const TopologicalRegularMaskedGrid *const mesh) :
				TopologicalRegularGrid3D::CellVerticesIterator(mesh), m_mesh(mesh) {}


			void begin(Vec3l const &coords) {
				TopologicalRegularGrid3D::CellVerticesIterator::begin(coords);
				while (TopologicalRegularGrid3D::CellVerticesIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::CellVerticesIterator::value()))
					TopologicalRegularGrid3D::CellVerticesIterator::advance();
			}
			void begin(INDEX_TYPE const &cellid) {
				TopologicalRegularGrid3D::CellVerticesIterator::begin(cellid);
				while (TopologicalRegularGrid3D::CellVerticesIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::CellVerticesIterator::value()))
					TopologicalRegularGrid3D::CellVerticesIterator::advance();
			}
			void advance() {
				TopologicalRegularGrid3D::CellVerticesIterator::advance();
				while (TopologicalRegularGrid3D::CellVerticesIterator::valid() && !
					m_mesh->InMesh(TopologicalRegularGrid3D::CellVerticesIterator::value()))
					TopologicalRegularGrid3D::CellVerticesIterator::advance();
			}
			bool valid() const {
				return TopologicalRegularGrid3D::CellVerticesIterator::valid();
			}
			INDEX_TYPE value() const {
				return TopologicalRegularGrid3D::CellVerticesIterator::value();
			}
		};





	public:

		TopologicalRegularMaskedGrid(RegularGrid3D* base_grid): TopologicalRegularGrid3D(base_grid) 
		{
		}

		void SetTester(CellTester* tester) { m_tester = tester; }

		virtual  ~TopologicalRegularMaskedGrid() {
			printf("delete: TopologicalRegularMaskedGrid \n");
		}



	};

}


#endif