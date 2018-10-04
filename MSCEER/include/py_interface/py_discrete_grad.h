#ifndef PY_DISCRETE_GRAD_H
#define PY_DISCRETE_GRAD_H


#include "base/gi_basic_types.h"
#include "base/gi_labeling.h"
#include "base/gi_discrete_gradient_labeling.h"
#include "base/gi_topological_regular_grid_2d.h"
#include "base/gi_regular_grid_2d.h"
#include "base/gi_topological_regular_grid_3d.h"
#include "base/gi_regular_grid_3d.h"

#include "py_interface/py_basic_types.h"

// wrapper class for meshes
class DiscreteGrad {

public:
	
	virtual bool isCritical(IndexType id) { return false; }
	virtual bool isHead(IndexType id) { return false; }
	virtual bool isTail(IndexType id) { return false; }
	virtual IndexType pair(IndexType id) { return 0; }


};


class GridDiscreteGrad2D : public DiscreteGrad {
protected:
	GInt::TopologicalRegularGrid2D* mMesh;
	GInt::DiscreteGradientLabeling<GInt::TopologicalRegularGrid2D>* mGrad;

public:
	virtual bool isCritical(IndexType id) { return mGrad->getCritical(id); }
	virtual IndexType pair(IndexType id) { return mGrad->getPair(id); }
	virtual bool isHead(IndexType id) { 
		if (isCritical(id)) return false;
		IndexType tpair = pair(id);
		if (mMesh->dimension(id) > mMesh->dimension(tpair)) return true;
		return false;
	}
	virtual bool isTail(IndexType id) { 
		if (isCritical(id)) return false;
		return !isHead(id);
	}
};

class GridDiscreteGrad3D : public DiscreteGrad {
protected:
	GInt::TopologicalRegularGrid3D* mMesh;
	GInt::DiscreteGradientLabeling<GInt::TopologicalRegularGrid3D>* mGrad;

public:
	virtual bool isCritical(IndexType id) { return mGrad->getCritical(id); }
	virtual IndexType pair(IndexType id) { return mGrad->getPair(id); }
	virtual bool isHead(IndexType id) {
		if (isCritical(id)) return false;
		IndexType tpair = pair(id);
		if (mMesh->dimension(id) > mMesh->dimension(tpair)) return true;
		return false;
	}
	virtual bool isTail(IndexType id) {
		if (isCritical(id)) return false;
		return !isHead(id);
	}
};

#endif