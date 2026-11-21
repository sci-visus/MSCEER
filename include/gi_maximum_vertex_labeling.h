#ifndef GI_MAXIMUM_VERTEX_LABELING_H
#define GI_MAXIMUM_VERTEX_LABELING_H

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

	template <class MeshType, class GridFuncType>
	class MaximumVertexLabeling {
	protected:
		DenseLabeling<DIM_TYPE>* mMaxVLabel;
		
		MeshType* mMesh;
		GridFuncType* mGridFunc;
	public:

		MaximumVertexLabeling(MeshType* mesh, GridFuncType* gridFunc) :
			mMesh(mesh), mGridFunc(gridFunc) {
			mMaxVLabel = NULL;
		}

		~MaximumVertexLabeling() {
			if (mMaxVLabel != NULL) delete mMaxVLabel;
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
			return mMesh->UncompressByteToVertexOffset(cellid, mMaxVLabel->GetLabel(cellid));
		}

		bool Before(INDEX_TYPE a, INDEX_TYPE b) const {
			return mGridFunc->IsGreater(b, a);
		}

		void ComputeOutput() {

			ThreadedTimer timer(1);
			timer.StartGlobal();

			printf(" -- Creating maxV_labeling ...");
			fflush(stdout);

			if (mMaxVLabel != NULL) delete mMaxVLabel;

			INDEX_TYPE num_cells = mMesh->numCells();
			mMaxVLabel = new DenseLabeling<DIM_TYPE>(num_cells);
			//mMaxVLabel->SetAll(0); // not needed since this is initialized to first vertex anyway
			//return;
#pragma omp parallel
			{
				int num_threads = omp_get_num_threads();
				int thread_num = omp_get_thread_num();

				std::vector<INDEX_TYPE> partition;
				ArrayIndexPartitioner::EvenChunkSplit(num_cells, num_threads, partition);

				for (INDEX_TYPE cellid = partition[thread_num]; cellid < partition[thread_num + 1]; cellid++) {
					mMaxVLabel->SetLabel(cellid, idxOf_maxVertex(cellid));
				}
			}

			//MeshType::AllCellsIterator ait(mMesh);
			//for (ait.begin(); ait.valid(); ait.advance()){
			//	mMaxVLabel->SetLabel(ait.value(), idxOf_maxVertex(ait.value()));
			//}

			timer.EndGlobal();
			printf(" Done! ");
			timer.PrintAll();
		}

 


	};

}

#endif
