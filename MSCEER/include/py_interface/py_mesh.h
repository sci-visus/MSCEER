#ifndef PY_MESH_H
#define PY_MESH_H

#include <vector>
#include <set>
#include <queue>
#include <time.h>

#include "base/gi_basic_types.h"
#include "base/gi_labeling.h"
#include "base/gi_discrete_gradient_labeling.h"
#include "base/gi_topological_regular_grid_2d.h"
#include "base/gi_regular_grid_2d.h"
#include "base/gi_topological_regular_grid_3d.h"
#include "base/gi_regular_grid_3d.h"

#include "py_interface/py_basic_types.h"

// wrapper class for meshes
class Mesh {

public:

	// global information about mesh
	virtual  int meshDimension() { return 0; }
	virtual IndexType numCells() { return 0; }
	virtual IndexType numDCells(int dim) { return 0; }

	// information about individual cells
	virtual int boundaryValue(IndexType cellid) { return 0; }
	virtual int dimension(IndexType cellid) { return 0; }
	virtual Coord Centroid(IndexType cellid) {
		return{ 0, 0, 1 };
	}

	// neighborhoods
	virtual CellList Facets(IndexType cellid) {
		return{ 0 };
	}
	virtual CellList Cofacets(IndexType cellid) {
		return{ 0 };
	}

};

class RegularGrid2D : public Mesh {
protected:
	GInt::RegularGrid2D* mGrid;
	GInt::TopologicalRegularGrid2D* mMesh;

public:
	RegularGrid2D(int x, int y) {
		mGrid = new GInt::RegularGrid2D(GInt::Vec2l(x, y), GInt::Vec2b(false, false));
		mMesh = new GInt::TopologicalRegularGrid2D(mGrid);
	}

	// global information about mesh
	virtual  int meshDimension() { return 2; }
	virtual IndexType numCells() { return mMesh->numCells(); }
	virtual IndexType numDCells(int dim) { return mMesh->numCells(dim); }

	// information about individual cells
	virtual int boundaryValue(IndexType cellid) { return mMesh->boundaryValue(cellid); }
	virtual int dimension(IndexType cellid) { return mMesh->dimension(cellid); }
	virtual Coord Centroid(IndexType cellid) {
		GInt::Vec2d coords;
		mMesh->centroid(cellid, coords);
		Coord ret(2); ret[0] = coords[0]; ret[1] = coords[1];
		return ret;
	}

	// neighborhoods
	virtual CellList Facets(IndexType cellid) {
		CellList ret;
		GInt::TopologicalRegularGrid2D::FacetsIterator fit(mMesh);
		for (fit.begin(cellid); fit.valid(); fit.advance()) {
			ret.push_back(fit.value());
		}
		return ret;
	}
	virtual CellList Cofacets(IndexType cellid) {
		CellList ret;
		GInt::TopologicalRegularGrid2D::CofacetsIterator fit(mMesh);
		for (fit.begin(cellid); fit.valid(); fit.advance()) {
			ret.push_back(fit.value());
		}
		return ret;
	}
};

class RegularGrid3D : public Mesh {
protected:
	GInt::RegularGrid3D* mGrid;
	GInt::TopologicalRegularGrid3D* mMesh;

public:
	RegularGrid3D(int x, int y, int z) {
		mGrid = new GInt::RegularGrid3D(GInt::Vec3l(x, y, z), GInt::Vec3b(false, false, false));
		mMesh = new GInt::TopologicalRegularGrid3D(mGrid);
	}

	// global information about mesh
	virtual  int meshDimension() { return 3; }
	virtual IndexType numCells() { return mMesh->numCells(); }
	virtual IndexType numDCells(int dim) { return mMesh->numCells(dim); }

	// information about individual cells
	virtual int boundaryValue(IndexType cellid) { return mMesh->boundaryValue(cellid); }
	virtual int dimension(IndexType cellid) { return mMesh->dimension(cellid); }
	virtual Coord Centroid(IndexType cellid) {
		GInt::Vec3d coords;
		mMesh->centroid(cellid, coords);
		Coord ret(3); ret[0] = coords[0]; ret[1] = coords[1]; ret[2] = coords[2];
		return ret;
	}

	// neighborhoods
	virtual CellList Facets(IndexType cellid) {
		CellList ret;
		GInt::TopologicalRegularGrid3D::FacetsIterator fit(mMesh);
		for (fit.begin(cellid); fit.valid(); fit.advance()) {
			ret.push_back(fit.value());
		}
		return ret;
	}
	virtual CellList Cofacets(IndexType cellid) {
		CellList ret;
		GInt::TopologicalRegularGrid3D::CofacetsIterator fit(mMesh);
		for (fit.begin(cellid); fit.valid(); fit.advance()) {
			ret.push_back(fit.value());
		}
		return ret;
	}
};

#endif