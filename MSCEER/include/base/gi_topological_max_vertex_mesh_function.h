/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef TOPOLOGICAL_MAX_VERTEX_MESH_FUNCTION
#define TOPOLOGICAL_MAX_VERTEX_MESH_FUNCTION

#include <algorithm>
#include <omp.h>
#include "gi_basic_types.h"
//#include "gi_topological_regular_grid.h"
#include "gi_array_index_partition.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_max_vertex_labeling.h"

namespace GInt {
	//
	//typedef TopologicalRegularGrid3D Mesh;
	//typedef float dtype;

	//typedef MaximumVertexLabeling<Mesh, dtype> MaxVLabelingType;
	//typedef RegularGridTrilinearFunction GridFuncType;
	template<class Mesh, class MaxVLabelingType, class GridFuncType, class dtype>
	class TopologicalMaxVertexMeshFunction {
	protected:
		Mesh* mMesh;
		MaxVLabelingType* mMaxVs;
		GridFuncType* mGridFunc;
	public:


		typedef dtype DType;


	public:
		TopologicalMaxVertexMeshFunction()  {
			mMaxVs = NULL;
			mMesh = NULL;
		}

		~TopologicalMaxVertexMeshFunction()  {
		}

		void setMeshAndFuncAndMaxVLabeling(Mesh* m, GridFuncType* func, MaxVLabelingType* maxv) {
			mMesh = m;
			mMaxVs = maxv;
			mGridFunc = func;
		}

		dtype cellValue(INDEX_TYPE cellid) const {
			return mGridFunc->SampleImage(mMesh->VertexNumberFromCellID(mMaxVs->Cell2HighestVertex(cellid)));
		}

        const Mesh* mesh() const {
            return this->mMesh;
        }
        bool lessThan(INDEX_TYPE a, INDEX_TYPE b) const {
            dtype av = cellValue(a);
            dtype bv = cellValue(b);
            if (av < bv) return true;
            if (bv < av) return false;
			INDEX_TYPE avert = mMaxVs->Cell2HighestVertex(a);
			INDEX_TYPE bvert = mMaxVs->Cell2HighestVertex(b);
			if (avert < bvert) return true;
			if (bvert < avert) return false;
            return a < b;
        }
        bool greaterThan(INDEX_TYPE a, INDEX_TYPE b) const {
            return lessThan (b, a);
        }


	};
}
#endif
