#ifndef GI_LAZY_MAXIMUM_VERTEX_LABELING_H
#define GI_LAZY_MAXIMUM_VERTEX_LABELING_H

#include <assert.h>
#include <algorithm>
#include <unordered_map>

#include "gi_basic_types.h"
#include "gi_timing.h"
#include "gi_regular_grid.h"
#include "gi_topological_regular_grid.h"
#include "gi_labeling.h"
#include "gi_array_index_partition.h"

namespace GInt {

	// the purpose of this class is applicaitons where we need the values at a
	// sparse set of cells. This does no caching. better to use non-lazy if number
	// of queries will be on the order of the mesh elements
	template <class MeshType, class GridFuncType>
	class LazyMaximumVertexLabeling {
	protected:
		MeshType* mMesh;
		GridFuncType* mGridFunc;
	public:

		LazyMaximumVertexLabeling(MeshType* mesh, GridFuncType* gridFunc) :
			mMesh(mesh), mGridFunc(gridFunc) {
		}

		DIM_TYPE idxOf_maxVertex(INDEX_TYPE cellid) const {

			typename MeshType::CellVerticesIterator cviter(mMesh);
			cviter.begin(cellid);

			INDEX_TYPE maxV = cviter.value();
			DIM_TYPE test_dir = mMesh->CompressVertexOffsetToByte(cellid, maxV);
#ifdef DEBUGME
			//INDEX_TYPE slice_position_list = mMesh->UncompressByteToVertexOffset(cellid, test_dir);
			//if (slice_position_list != maxV) printf("WHOATHERENELLY!\n");
#endif

			INDEX_TYPE maxVgid = mMesh->VertexNumberFromCellID(maxV);
			cviter.advance();
			for (; cviter.valid(); cviter.advance()) {
				INDEX_TYPE otherid = cviter.value();
				INDEX_TYPE othergid = mMesh->VertexNumberFromCellID(otherid);
				if (mGridFunc->IsGreater(othergid, maxVgid)) {
					maxV = otherid;
					maxVgid = othergid;
					test_dir = mMesh->CompressVertexOffsetToByte(cellid, maxV);

#ifdef DEBUGME
					//slice_position_list = mMesh->UncompressByteToVertexOffset(cellid, test_dir);
					//if (slice_position_list != maxV) printf("WHOATHERENELLY2!\n");
#endif
				}
			}
			return test_dir;
		}

		INDEX_TYPE Cell2HighestVertex(INDEX_TYPE cellid) {
			return mMesh->UncompressByteToVertexOffset(cellid, idxOf_maxVertex(cellid));
		}

		//bool Before(INDEX_TYPE a, INDEX_TYPE b) const {
		//	return mGridFunc->IsGreater(b, a);
		//}

		void ComputeOutput() {
			// do nothing! we are lazy!!
		}
	};

}

#endif
