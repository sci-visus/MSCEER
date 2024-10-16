#ifndef MAXIMUM_VERTEX_LABELING_H
#define MAXIMUM_VERTEX_LABELING_H

#include <assert.h>
#include <algorithm>
#include <unordered_map> 

#include "gi_basic_types.h"
#include "gi_regular_grid.h"
#include "gi_topological_regular_grid.h"
#include "gi_labeling.h"
#include "gi_array_index_partition.h"

namespace GInt {
	

	// -------------------------------------------------------------------------------
	// comparison of cells based on the max vertex

	template <class MeshType, class GridFuncType>
	class SparseMaximumVertexLabeling {
	protected:
		SparseLabeling<DIM_TYPE>* mMaxVLabel;
		MeshType* mMesh;
		GridFuncType* mGridFunc;
	public:

		SparseMaximumVertexLabeling(MeshType* mesh, GridFuncType* gridFunc) :
			mMesh(mesh), mGridFunc(gridFunc) {
			mMaxVLabel = NULL;
		}

		~SparseMaximumVertexLabeling() {
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
			mMaxVLabel = new SparseLabeling<DIM_TYPE>();
			//mMaxVLabel->SetAll(0); // not needed since this is initialized to first vertex anyway
			//return;
#pragma omp parallel
			{
				int num_threads = omp_get_num_threads();
				int thread_num = omp_get_thread_num();
				std::unordered_map<INDEX_TYPE, DIM_TYPE> local_map;
				std::vector<INDEX_TYPE> partition;
				ArrayIndexPartitioner::EvenChunkSplit(num_cells, num_threads, partition);

				for (INDEX_TYPE cellid = partition[thread_num]; cellid < partition[thread_num + 1]; cellid++) {
					
					if (mMesh->Has(cellid)) {
						local_map[cellid] = idxOf_maxVertex(cellid);
					}
				}
#pragma omp critical
				{
					for (auto& p : local_map) {
						mMaxVLabel->SetLabel(p.first, p.second);
					}
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


	template <class MeshType, class GridFuncType>
	class UncachedMaximumVertexLabeling {
	protected:

		MeshType* mMesh;
		GridFuncType* mGridFunc;
	public:

		UncachedMaximumVertexLabeling(MeshType* mesh, GridFuncType* gridFunc) :
			mMesh(mesh), mGridFunc(gridFunc) {
		}

		~UncachedMaximumVertexLabeling() {
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

		bool Before(INDEX_TYPE a, INDEX_TYPE b) const {
			return mGridFunc->IsGreater(b, a);
		}

		void ComputeOutput() {

			ThreadedTimer timer(1);
			timer.StartGlobal();

			printf(" -- Creating maxV_labeling ...");
			fflush(stdout);

			timer.EndGlobal();
			printf(" Done! ");
			timer.PrintAll();
		}




	};

//	template <class MeshType, class GridFuncType>
//	class RegularGridMaximumVertexLabeling {
//	protected:
//		//DenseLabeling<DIM_TYPE>* mMaxVLabel;
//		//INDEX_TYPE* tmp_ids;
//		BYTE_TYPE* byte_tmp_max_ids;
//		BYTE_TYPE* byte_tmp_min_ids;
//		const RegularGrid3D* mGrid;
//		TopologicalRegularGrid3D* mMesh;
//		RegularGridTrilinearFunction* mGridFunc;
//	public:
//
//		RegularGridMaximumVertexLabeling(MeshType* mesh, GridFuncType* gridFunc) :
//			mMesh(mesh), mGridFunc(gridFunc) {
//			//mMaxVLabel = NULL;
//			mGrid = mGridFunc->GetGrid();
//		}
//
//		~RegularGridMaximumVertexLabeling() {
//			//if (mMaxVLabel != NULL) delete mMaxVLabel;
//		}
//
////		DIM_TYPE idxOf_maxVertex(INDEX_TYPE cellid) const {
////
////			typename MeshType::CellVerticesIterator cviter(mMesh);
////			cviter.begin(cellid);
////
////			INDEX_TYPE maxV = cviter.value();
////			DIM_TYPE test_dir = mMesh->CompressVertexOffsetToByte(cellid, maxV, cviter);
////#ifdef DEBUGME
////			//INDEX_TYPE slice_position_list = mMesh->UncompressByteToVertexOffset(cellid, test_dir);
////			//if (slice_position_list != maxV) printf("WHOATHERENELLY!\n");
////#endif
////
////			INDEX_TYPE maxVgid = mMesh->VertexNumberFromCellID(maxV);
////			cviter.advance();
////			for (; cviter.valid(); cviter.advance()) {
////				INDEX_TYPE otherid = cviter.value();
////				INDEX_TYPE othergid = mMesh->VertexNumberFromCellID(otherid);
////				if (mGridFunc->IsGreater(othergid, maxVgid)) {
////					maxV = otherid;
////					maxVgid = othergid;
////					test_dir = mMesh->CompressVertexOffsetToByte(cellid, maxV, cviter);
////
////#ifdef DEBUGME
////					//slice_position_list = mMesh->UncompressByteToVertexOffset(cellid, test_dir);
////					//if (slice_position_list != maxV) printf("WHOATHERENELLY2!\n");
////#endif
////				}
////			}
////			return test_dir;
////		}
//
//		INDEX_TYPE Cell2HighestVertex(INDEX_TYPE cellid) {
//			return mMesh->UncompressByteToVertexOffset(cellid, byte_tmp_max_ids[cellid]);
//			//return mMesh->UncompressByteToVertexOffset(cellid, mMaxVLabel->GetLabel(cellid));
//		}
//
//		bool Before(INDEX_TYPE a, INDEX_TYPE b) const {
//			return mGridFunc->IsGreater(b, a);
//		}
//
//		struct ValInd {
//			float val;
//			INDEX_TYPE ind;
//		};
//
//
//		void do_fastest_line(int X, const INDEX_TYPE base_offset, float* values, ValInd* o_row) {
//			// X = dim of line
//			// float* values
//			//out = float*, int*
//			int dX = 2 * X - 1;
//			
//			// do 0th index
//			float tmpold = values[0];
//			o_row[0].ind = 0 + base_offset;
//			o_row[0].val = tmpold;
//
//			// now do rest!
//			for (int i = 1; i < X; i ++) {
//				// fisrt get tmp next value and set
//				int di_v = i * 2; // mesh index of vertex
//				float tmpnew = values[i];
//				o_row[di_v].ind = base_offset + di_v;
//				o_row[di_v].val = tmpnew;
//				
//				int di_e = i * 2 - 1; // mesh index of edge
//				if (tmpold > tmpnew) {
//					o_row[di_e].val = tmpold;
//					o_row[di_e].ind = base_offset + di_v - 2; /// check this
//				}
//				else {
//					o_row[di_e].val = tmpnew;
//					o_row[di_e].ind = base_offset + di_v; /// check this
//				}
//
//
//
//				//
//				tmpold = tmpnew;
//			}
//			
//		}
//
//		void do_parallel_lines(int X, 
//			const ValInd* row0, 
//			const ValInd* row1,
//			ValInd* row_o ) {
//			
//			int dX = 2 * X - 1;
//
//			// now do rest!
//			for (int i = 0; i < dX; i++) {
//				// fisrt get tmp next value and set
//				float v0 = row0[i].val;
//				float v1 = row1[i].val;
//				
//				if (v0 > v1) {
//					row_o[i].val = v0;
//					row_o[i].ind = row0[i].ind;
//					continue;
//				}
//				else if (v1 > v0) {
//					row_o[i].val = v1;
//					row_o[i].ind = row1[i].ind;
//					continue;
//				}
//
//				// need index compare			
//				INDEX_TYPE i0 = row0[i].ind;
//				INDEX_TYPE i1 = row1[i].ind;
//				if (i0 > i1) {
//					row_o[i].val = v0;
//					row_o[i].ind = i0;
//				}
//				else {
//					row_o[i].val = v1;
//					row_o[i].ind = i1;
//				}
//			}
//		}
//
//
//		void write_stick(int dX, const INDEX_TYPE base_offset, ValInd* stick) {
//
//			for (INDEX_TYPE i = 0; i < dX; i++) {
//				INDEX_TYPE base_id = i + base_offset;
//				byte_tmp_max_ids[base_id] = mMesh->CompressVertexOffsetToByte(base_id, stick[i].ind);
//			}
//
//		}
//
//		void ComputeOutput() {
//
//			ThreadedTimer timer(1);
//			timer.StartGlobal();
//
//			printf(" -- Creating maxV_labeling ...");
//			fflush(stdout);
//
//			byte_tmp_max_ids = new BYTE_TYPE[mMesh->numCells()];
//			byte_tmp_min_ids = new BYTE_TYPE[mMesh->numCells()];
//			const Vec3l xyz = mGrid->XYZ();
//			const Vec3l dXYZ = mMesh->XYZ();
//			const INDEX_TYPE NUM_X = dXYZ[0];
//			const INDEX_TYPE NUM_Y = dXYZ[1];
//			const INDEX_TYPE SLAB_SIZE = NUM_X * NUM_Y;
//
//
//			//			INDEX_TYPE num_cells = mMesh->numCells();
//#pragma omp parallel
//			{
//				int num_threads = omp_get_num_threads();
//				int thread_num = omp_get_thread_num();
//				
//				// local storage
//				ValInd* slab = new ValInd[SLAB_SIZE * 3];
//				
//				//float* slab_vals = new FLOATTYPE[SLAB_SIZE * 3];
//				//INDEX_TYPE* slab_inds = new INDEX_TYPE[SLAB_SIZE * 3];
//
//				std::vector<INDEX_TYPE> partition;
//				ArrayIndexPartitioner::EvenChunkSplit(xyz[2], num_threads, partition);
//							
//				for (auto v : partition) printf("p %llu\n", v);
//				const INDEX_TYPE START_K = partition[thread_num];
//				INDEX_TYPE localend = (thread_num == num_threads - 1 ? 
//					partition[thread_num + 1] - 1:
//					partition[thread_num + 1]);
//				const INDEX_TYPE END_K = localend;
//
//				ValInd* old_slab = slab;
//				ValInd* mid_slab = &(slab[SLAB_SIZE * 1]);
//				ValInd* new_slab = &(slab[SLAB_SIZE * 2]);
//
//				for (INDEX_TYPE k = START_K; k <= END_K; k++) {
//
//				//float* old_slab_vals = slab_vals;
//					//float* mid_slab_vals; &(slab_vals[SLAB_SIZE * 1);
//					//float* new_slab_vals = &(slab_vals[SLAB_SIZE * 2]);
//					//INDEX_TYPE* old_slab_inds = slab_inds;
//					//INDEX_TYPE* mid_slab_inds; &(slab_inds[SLAB_SIZE * 1);
//					//INDEX_TYPE* new_slab_inds = &(slab_inds[SLAB_SIZE * 2]);
//
//					for (INDEX_TYPE j = 0; j < xyz[1]; j++) {
//						INDEX_TYPE stick_start_j0k0 = j * 2 * NUM_X; // get the row in the slab
//						
//						// J*2, K*2
//						// Do each vertex-edge-vertex-....-edge-vertex stick in the slab
//						INDEX_TYPE base_cell_id = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2));
//						INDEX_TYPE id_data = mGrid->Index3d(Vec3l(0, j, k));
//						do_fastest_line(xyz[0], base_cell_id,
//							&(mGridFunc->GetImage()[id_data]),
//							&(new_slab[stick_start_j0k0]));
//						write_stick(NUM_X, base_cell_id, &(new_slab[stick_start_j0k0]));
//
//						// J*2-1, K*2 = mid of J*2-2, K*2
//						// Do the edge-face-edge-...-face-edge stick in the slab
//						INDEX_TYPE stick_start_j2k0 = (j * 2 - 2) * NUM_X; // get the row in the slab
//						INDEX_TYPE stick_start_j1k0 = (j * 2 - 1) * NUM_X; // get the row in the slab
//						INDEX_TYPE base_cell_j1_id = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2));
//						
//						if (j > 0) {
//							// EDGE BETWEEN v(x, y:y-1, z)
//							do_parallel_lines(xyz[0],
//								&(new_slab[stick_start_j0k0]),
//								&(new_slab[stick_start_j2k0]),
//								&(new_slab[stick_start_j1k0]));
//							write_stick(NUM_X, base_cell_j1_id, &(new_slab[stick_start_j1k0]));
//						}
//
//						// NOW FILL IN BETWEEN SLABS
//						
//						// J*2, K*2-1 = mid of J*2, K*2-2
//						// Do the edge-face-edge-...-face-edge stick in between slabs
//						if (k > START_K) {
//							INDEX_TYPE base_cell_k1_id = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 1));
//							do_parallel_lines(xyz[0],
//								&(new_slab[stick_start_j0k0]),
//								&(old_slab[stick_start_j0k0]),
//								&(mid_slab[stick_start_j0k0]));
//							write_stick(NUM_X, base_cell_k1_id, &(mid_slab[stick_start_j0k0]));
//
//						}
//
//						// J*2-1, K*2-1 = mid of J*2-1, K*2-2
//						// Do the face-quad-face-....-quad-face stick between the slabs
//						if (k > START_K && j > 0) {
//							// face BETWEEN v(x, y:y-1, z:z-1)
//
//							do_parallel_lines(xyz[0],
//								&(new_slab[stick_start_j1k0]),
//								&(old_slab[stick_start_j1k0]),
//								&(mid_slab[stick_start_j1k0]));
//							INDEX_TYPE base_mid_id = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2 - 1));
//							write_stick(NUM_X, base_mid_id, &(mid_slab[stick_start_j1k0]));
//						}
//
//
//
//					}
//
//						// swap new and old slabs
//						ValInd* tmp = old_slab;
//						old_slab = new_slab;
//						new_slab = tmp;				
//						//do_fastest_line(xyz[0], id_cell,
//						//	&(mGridFunc->GetImage()[id_data]),
//						//	&(tmp_vals[id_cell]),
//						//	&(tmp_ids[id_cell]));
//
//						//if (j > 0) {
//						//	// EDGE BETWEEN v(x, y:y-1, z)
//						//	INDEX_TYPE id_cell_last = mMesh->coords2Cellid(Vec3l(0, j * 2 - 2, k * 2));
//						//	INDEX_TYPE id_cell_mid = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2));
//						//	do_parallel_lines(xyz[0],
//						//		&(tmp_vals[id_cell]),
//						//		&(tmp_ids[id_cell]),
//						//		&(tmp_vals[id_cell_last]),
//						//		&(tmp_ids[id_cell_last]),
//						//		&(tmp_vals[id_cell_mid]),
//						//		&(tmp_ids[id_cell_mid]));
//						//}
//						//if (k > START_K) {
//						//	// EDGE BETWEEN v(x, y, z:z-1)
//						//	INDEX_TYPE id_cell_last = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 2));
//						//	INDEX_TYPE id_cell_mid = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 1));
//						//	do_parallel_lines(xyz[0],
//						//		&(tmp_vals[id_cell]),
//						//		&(tmp_ids[id_cell]),
//						//		&(tmp_vals[id_cell_last]),
//						//		&(tmp_ids[id_cell_last]),
//						//		&(tmp_vals[id_cell_mid]),
//						//		&(tmp_ids[id_cell_mid]));
//						//}
//						//if (k > START_K && j > 0) {
//						//	// face BETWEEN v(x, y:y-1, z:z-1)
//						//	INDEX_TYPE id_cell_e0 = mMesh->coords2Cellid(Vec3l(0, j * 2 - 2, k * 2 - 1));
//						//	INDEX_TYPE id_cell_e1 = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 1));
//						//	INDEX_TYPE id_cell_mid = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2 - 1));
//						//	do_parallel_lines(xyz[0],
//						//		&(tmp_vals[id_cell_e0]),
//						//		&(tmp_ids[id_cell_e0]),
//						//		&(tmp_vals[id_cell_e1]),
//						//		&(tmp_ids[id_cell_e1]),
//						//		&(tmp_vals[id_cell_mid]),
//						//		&(tmp_ids[id_cell_mid]));
//						//}
//
//				}
//
//				//INDEX_TYPE num_cells = mMesh->numCells();
//				//for (INDEX_TYPE i = 0; i < num_cells; i++) {
//				//	byte_tmp_max_ids[i] = mMesh->CompressVertexOffsetToByte(i, tmp_ids[i]); // THIS COULD BE IMPROVED WITH BETTER HASH
//				//}
//
//
//				//			INDEX_TYPE num_cells = mMesh->numCells();
//				//			mMaxVLabel = new DenseLabeling<DIM_TYPE>(num_cells);
//				//			//mMaxVLabel->SetAll(0); // not needed since this is initialized to first vertex anyway
//				//			//return;
//				//#pragma omp parallel
//				//			{
//				//				int num_threads = omp_get_num_threads();
//				//				int thread_num = omp_get_thread_num();
//				//
//				//				std::vector<INDEX_TYPE> partition;
//				//				ArrayIndexPartitioner::EvenChunkSplit(num_cells, num_threads, partition);
//				//
//				//				for (INDEX_TYPE cellid = partition[thread_num]; cellid < partition[thread_num + 1]; cellid++) {
//				//					mMaxVLabel->SetLabel(cellid, idxOf_maxVertex(cellid));
//				//				}
//				//			}
//				//
//				//			//MeshType::AllCellsIterator ait(mMesh);
//				//			//for (ait.begin(); ait.valid(); ait.advance()){
//				//			//	mMaxVLabel->SetLabel(ait.value(), idxOf_maxVertex(ait.value()));
//				//			//}
//
//				delete[] slab;
//			}
//			timer.EndGlobal();
//			printf(" Done! ");
//			timer.PrintAll();
//		}
//
//
//
//
//	};
//

	template <class MeshType, class GridFuncType>
	class RegularGridMaxMinVertexLabeling3D {
	protected:

		BYTE_TYPE* byte_tmp_max_ids;
		BYTE_TYPE* byte_tmp_min_ids;
		const RegularGrid3D* mGrid;
		MeshType* mMesh;
		GridFuncType* mGridFunc;
	public:

		RegularGridMaxMinVertexLabeling3D(MeshType* mesh, GridFuncType* gridFunc) :
			mMesh(mesh), mGridFunc(gridFunc) {
			byte_tmp_max_ids = NULL;
			byte_tmp_min_ids = NULL;
			mGrid = mGridFunc->GetGrid();
		}

		~RegularGridMaxMinVertexLabeling3D() {
			if (byte_tmp_max_ids != NULL) delete[] byte_tmp_max_ids;
			if (byte_tmp_min_ids != NULL) delete[] byte_tmp_min_ids;
		}

		void HACK_init() {
			byte_tmp_max_ids = new BYTE_TYPE[mMesh->numCells()];
			byte_tmp_min_ids = new BYTE_TYPE[mMesh->numCells()];
		}

		BYTE_TYPE GetUncompressedMaxVal(INDEX_TYPE cellid) const {
			return byte_tmp_max_ids[cellid];
		}
		BYTE_TYPE GetUncompressedMinVal(INDEX_TYPE cellid) const {
			return byte_tmp_min_ids[cellid];
		}
		void SetUncompressedMaxVal(INDEX_TYPE cellid, BYTE_TYPE val) {
			byte_tmp_max_ids[cellid] = val;
		}
		void SetUncompressedMinVal(INDEX_TYPE cellid, BYTE_TYPE val)  {
			byte_tmp_min_ids[cellid] = val;
		}
		INDEX_TYPE Cell2HighestVertex(INDEX_TYPE cellid) {
			return mMesh->UncompressByteToVertexOffset(cellid, byte_tmp_max_ids[cellid]);
		}
		INDEX_TYPE Cell2LowestVertex(INDEX_TYPE cellid) {
			return mMesh->UncompressByteToVertexOffset(cellid, byte_tmp_min_ids[cellid]);
		}
		bool Before(INDEX_TYPE a, INDEX_TYPE b) const {
			return mGridFunc->IsGreater(b, a);
		}

		struct ValInd {
			float val;
			INDEX_TYPE ind;

			ValInd() {}
			ValInd(float v, INDEX_TYPE i) : val(v), ind(i) {}
		};
		

		void do_fastest_line(int X, const INDEX_TYPE base_offset, float* values, ValInd* o_max_row, ValInd* o_min_row) {
			// X = dim of line
			// float* values
			//out = float*, int*
			int dX = 2 * X - 1;

			// do 0th index
			//ValInd old(values[0], base_offset);
			float tmpold = values[0];
			o_max_row[0].ind = 0 + base_offset;
			o_max_row[0].val = tmpold;
			o_min_row[0].ind = 0 + base_offset;
			o_min_row[0].val = tmpold;

			// now do rest!
			for (int i = 1; i < X; i++) {
				// fisrt get tmp next value and set
				int di_v = i * 2; // mesh index of vertex
				float tmpnew = values[i];
				o_max_row[di_v].ind = base_offset + di_v;
				o_max_row[di_v].val = tmpnew;
				o_min_row[di_v].ind = base_offset + di_v;
				o_min_row[di_v].val = tmpnew;

				int di_e = i * 2 - 1; // mesh index of edge
				if (tmpold > tmpnew) {
					o_max_row[di_e].val = tmpold;
					o_max_row[di_e].ind = base_offset + di_v - 2; /// check this
					o_min_row[di_e].val = tmpnew;
					o_min_row[di_e].ind = base_offset + di_v; /// check this
				}
				else  { // tmpold < tmpnew || di_v-2 < di_v <- 2nd part is always true
					o_max_row[di_e].val = tmpnew;
					o_max_row[di_e].ind = base_offset + di_v; /// check this
					o_min_row[di_e].val = tmpold;
					o_min_row[di_e].ind = base_offset + di_v - 2; /// check this
				}


				//
				tmpold = tmpnew;
			}

		}

		void do_parallel_lines_max(int X,
			const ValInd* row0,
			const ValInd* row1,
			ValInd* row_o) {

			int dX = 2 * X - 1;

			// now do rest!
			for (int i = 0; i < dX; i++) {
				// fisrt get tmp next value and set
				float v0 = row0[i].val;
				float v1 = row1[i].val;

				if (v0 > v1) {
					row_o[i].val = v0;
					row_o[i].ind = row0[i].ind;
					continue;
				}
				else if (v1 > v0) {
					row_o[i].val = v1;
					row_o[i].ind = row1[i].ind;
					continue;
				}

				// need index compare			
				INDEX_TYPE i0 = row0[i].ind;
				INDEX_TYPE i1 = row1[i].ind;
				if (i0 > i1) {
					row_o[i].val = v0;
					row_o[i].ind = i0;
				}
				else {
					row_o[i].val = v1;
					row_o[i].ind = i1;
				}
			}
		}
		void do_parallel_lines_min(int X,
			const ValInd* row0,
			const ValInd* row1,
		ValInd* row_o) {

			int dX = 2 * X - 1;

			// now do rest!
			for (int i = 0; i < dX; i++) {
				// fisrt get tmp next value and set
				float v0 = row0[i].val;
				float v1 = row1[i].val;

				if (v0 < v1) {
					row_o[i].val = v0;
					row_o[i].ind = row0[i].ind;
					continue;
				}
				else if (v1 < v0) {
					row_o[i].val = v1;
					row_o[i].ind = row1[i].ind;
					continue;
				}

				// need index compare			
				INDEX_TYPE i0 = row0[i].ind;
				INDEX_TYPE i1 = row1[i].ind;
				if (i0 < i1) {
					row_o[i].val = v0;
					row_o[i].ind = i0;
				}
				else {
					row_o[i].val = v1;
					row_o[i].ind = i1;
				}
			}
		}

		void write_max_stick(int dX, const INDEX_TYPE base_offset, ValInd* stick) {

			for (INDEX_TYPE i = 0; i < dX; i++) {
				INDEX_TYPE base_id = i + base_offset;
				byte_tmp_max_ids[base_id] = mMesh->CompressVertexOffsetToByte(base_id, stick[i].ind);
			}

		}
		void write_min_stick(int dX, const INDEX_TYPE base_offset, ValInd* stick) {

			for (INDEX_TYPE i = 0; i < dX; i++) {
				INDEX_TYPE base_id = i + base_offset;
				byte_tmp_min_ids[base_id] = mMesh->CompressVertexOffsetToByte(base_id, stick[i].ind);
			}

		}

		void allocate() {
			byte_tmp_max_ids = new BYTE_TYPE[mMesh->numCells()];
			byte_tmp_min_ids = new BYTE_TYPE[mMesh->numCells()];
		}

		void ComputeOutputSerial() {
			const Vec3l xyz = mGrid->XYZ();
			const Vec3l dXYZ = mMesh->XYZ();
			const INDEX_TYPE NUM_X = dXYZ[0];
			const INDEX_TYPE NUM_Y = dXYZ[1];
			const INDEX_TYPE SLAB_SIZE = NUM_X * NUM_Y;


			//			INDEX_TYPE num_cells = mMesh->numCells();
			{
				// local storage
				ValInd* slab = new ValInd[SLAB_SIZE * 3];
				ValInd* slab_min = new ValInd[SLAB_SIZE * 3];

				//float* slab_vals = new FLOATTYPE[SLAB_SIZE * 3];
				//INDEX_TYPE* slab_inds = new INDEX_TYPE[SLAB_SIZE * 3];

				int num_threads = 1;
				int thread_num = 0;
				std::vector<INDEX_TYPE> partition;
				ArrayIndexPartitioner::EvenChunkSplit(xyz[2], num_threads, partition);

				//for (auto v : partition) printf("p %llu\n", v);
				const INDEX_TYPE START_K = partition[thread_num];
				INDEX_TYPE localend = (thread_num == num_threads - 1 ?
					partition[thread_num + 1] - 1 :
					partition[thread_num + 1]);
				const INDEX_TYPE END_K = localend;

				ValInd* old_slab = slab;
				ValInd* mid_slab = &(slab[SLAB_SIZE * 1]);
				ValInd* new_slab = &(slab[SLAB_SIZE * 2]);
			
				ValInd* old_slab_min = slab_min;
				ValInd* mid_slab_min = &(slab_min[SLAB_SIZE * 1]);
				ValInd* new_slab_min = &(slab_min[SLAB_SIZE * 2]);

				for (INDEX_TYPE k = START_K; k <= END_K; k++) {

					//float* old_slab_vals = slab_vals;
					//float* mid_slab_vals; &(slab_vals[SLAB_SIZE * 1);
					//float* new_slab_vals = &(slab_vals[SLAB_SIZE * 2]);
					//INDEX_TYPE* old_slab_inds = slab_inds;
					//INDEX_TYPE* mid_slab_inds; &(slab_inds[SLAB_SIZE * 1);
					//INDEX_TYPE* new_slab_inds = &(slab_inds[SLAB_SIZE * 2]);

					for (INDEX_TYPE j = 0; j < xyz[1]; j++) {
						INDEX_TYPE stick_start_j0k0 = j * 2 * NUM_X; // get the row in the slab

																		// J*2, K*2
																		// Do each vertex-edge-vertex-....-edge-vertex stick in the slab
						INDEX_TYPE base_cell_id = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2));
						INDEX_TYPE id_data = mGrid->Index3d(Vec3l(0, j, k));
						do_fastest_line(xyz[0], base_cell_id,
							&(mGridFunc->GetImage()[id_data]),
							&(new_slab[stick_start_j0k0]),
							&(new_slab_min[stick_start_j0k0]));
						write_max_stick(NUM_X, base_cell_id, &(new_slab[stick_start_j0k0]));
						write_min_stick(NUM_X, base_cell_id, &(new_slab_min[stick_start_j0k0]));

						// J*2-1, K*2 = mid of J*2-2, K*2
						// Do the edge-face-edge-...-face-edge stick in the slab
						INDEX_TYPE stick_start_j2k0 = (j * 2 - 2) * NUM_X; // get the row in the slab
						INDEX_TYPE stick_start_j1k0 = (j * 2 - 1) * NUM_X; // get the row in the slab
						INDEX_TYPE base_cell_j1_id = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2));

						if (j > 0) {
							// EDGE BETWEEN v(x, y:y-1, z)
							do_parallel_lines_max(xyz[0],
								&(new_slab[stick_start_j0k0]),
								&(new_slab[stick_start_j2k0]),
								&(new_slab[stick_start_j1k0]));
							write_max_stick(NUM_X, base_cell_j1_id, &(new_slab[stick_start_j1k0]));
							do_parallel_lines_min(xyz[0],
								&(new_slab_min[stick_start_j0k0]),
								&(new_slab_min[stick_start_j2k0]),
								&(new_slab_min[stick_start_j1k0]));
							write_min_stick(NUM_X, base_cell_j1_id, &(new_slab_min[stick_start_j1k0]));
						}

						// NOW FILL IN BETWEEN SLABS

						// J*2, K*2-1 = mid of J*2, K*2-2
						// Do the edge-face-edge-...-face-edge stick in between slabs
						if (k > START_K) {
							INDEX_TYPE base_cell_k1_id = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 1));
							do_parallel_lines_max(xyz[0],
								&(new_slab[stick_start_j0k0]),
								&(old_slab[stick_start_j0k0]),
								&(mid_slab[stick_start_j0k0]));
							write_max_stick(NUM_X, base_cell_k1_id, &(mid_slab[stick_start_j0k0]));
							do_parallel_lines_min(xyz[0],
								&(new_slab_min[stick_start_j0k0]),
								&(old_slab_min[stick_start_j0k0]),
								&(mid_slab_min[stick_start_j0k0]));
							write_min_stick(NUM_X, base_cell_k1_id, &(mid_slab_min[stick_start_j0k0]));

						}

						// J*2-1, K*2-1 = mid of J*2-1, K*2-2
						// Do the face-quad-face-....-quad-face stick between the slabs
						if (k > START_K && j > 0) {
							// face BETWEEN v(x, y:y-1, z:z-1)

							INDEX_TYPE base_mid_id = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2 - 1));
							do_parallel_lines_max(xyz[0],
								&(new_slab[stick_start_j1k0]),
								&(old_slab[stick_start_j1k0]),
								&(mid_slab[stick_start_j1k0]));
							write_max_stick(NUM_X, base_mid_id, &(mid_slab[stick_start_j1k0]));
							do_parallel_lines_min(xyz[0],
								&(new_slab_min[stick_start_j1k0]),
								&(old_slab_min[stick_start_j1k0]),
								&(mid_slab_min[stick_start_j1k0]));
							write_min_stick(NUM_X, base_mid_id, &(mid_slab_min[stick_start_j1k0]));
						}



					}

					// swap new and old slabs
					ValInd* tmp = old_slab;
					old_slab = new_slab;
					new_slab = tmp;
					tmp = old_slab_min;
					old_slab_min = new_slab_min;
					new_slab_min = tmp;
					//do_fastest_line(xyz[0], id_cell,
					//	&(mGridFunc->GetImage()[id_data]),
					//	&(tmp_vals[id_cell]),
					//	&(tmp_ids[id_cell]));

					//if (j > 0) {
					//	// EDGE BETWEEN v(x, y:y-1, z)
					//	INDEX_TYPE id_cell_last = mMesh->coords2Cellid(Vec3l(0, j * 2 - 2, k * 2));
					//	INDEX_TYPE id_cell_mid = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2));
					//	do_parallel_lines(xyz[0],
					//		&(tmp_vals[id_cell]),
					//		&(tmp_ids[id_cell]),
					//		&(tmp_vals[id_cell_last]),
					//		&(tmp_ids[id_cell_last]),
					//		&(tmp_vals[id_cell_mid]),
					//		&(tmp_ids[id_cell_mid]));
					//}
					//if (k > START_K) {
					//	// EDGE BETWEEN v(x, y, z:z-1)
					//	INDEX_TYPE id_cell_last = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 2));
					//	INDEX_TYPE id_cell_mid = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 1));
					//	do_parallel_lines(xyz[0],
					//		&(tmp_vals[id_cell]),
					//		&(tmp_ids[id_cell]),
					//		&(tmp_vals[id_cell_last]),
					//		&(tmp_ids[id_cell_last]),
					//		&(tmp_vals[id_cell_mid]),
					//		&(tmp_ids[id_cell_mid]));
					//}
					//if (k > START_K && j > 0) {
					//	// face BETWEEN v(x, y:y-1, z:z-1)
					//	INDEX_TYPE id_cell_e0 = mMesh->coords2Cellid(Vec3l(0, j * 2 - 2, k * 2 - 1));
					//	INDEX_TYPE id_cell_e1 = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 1));
					//	INDEX_TYPE id_cell_mid = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2 - 1));
					//	do_parallel_lines(xyz[0],
					//		&(tmp_vals[id_cell_e0]),
					//		&(tmp_ids[id_cell_e0]),
					//		&(tmp_vals[id_cell_e1]),
					//		&(tmp_ids[id_cell_e1]),
					//		&(tmp_vals[id_cell_mid]),
					//		&(tmp_ids[id_cell_mid]));
					//}

				}

				//INDEX_TYPE num_cells = mMesh->numCells();
				//for (INDEX_TYPE i = 0; i < num_cells; i++) {
				//	byte_tmp_max_ids[i] = mMesh->CompressVertexOffsetToByte(i, tmp_ids[i]); // THIS COULD BE IMPROVED WITH BETTER HASH
				//}


				//			INDEX_TYPE num_cells = mMesh->numCells();
				//			mMaxVLabel = new DenseLabeling<DIM_TYPE>(num_cells);
				//			//mMaxVLabel->SetAll(0); // not needed since this is initialized to first vertex anyway
				//			//return;
				//#pragma omp parallel
				//			{
				//				int num_threads = omp_get_num_threads();
				//				int thread_num = omp_get_thread_num();
				//
				//				std::vector<INDEX_TYPE> partition;
				//				ArrayIndexPartitioner::EvenChunkSplit(num_cells, num_threads, partition);
				//
				//				for (INDEX_TYPE cellid = partition[thread_num]; cellid < partition[thread_num + 1]; cellid++) {
				//					mMaxVLabel->SetLabel(cellid, idxOf_maxVertex(cellid));
				//				}
				//			}
				//
				//			//MeshType::AllCellsIterator ait(mMesh);
				//			//for (ait.begin(); ait.valid(); ait.advance()){
				//			//	mMaxVLabel->SetLabel(ait.value(), idxOf_maxVertex(ait.value()));
				//			//}

				delete[] slab;
				delete[] slab_min;
			}
		}

		void ComputeOutput() {

			ThreadedTimer timer(1);
			timer.StartGlobal();

			printf(" -- Creating maxV_labeling ...");
			fflush(stdout);

			byte_tmp_max_ids = new BYTE_TYPE[mMesh->numCells()];
			byte_tmp_min_ids = new BYTE_TYPE[mMesh->numCells()];
			const Vec3l xyz = mGrid->XYZ();
			const Vec3l dXYZ = mMesh->XYZ();
			const INDEX_TYPE NUM_X = dXYZ[0];
			const INDEX_TYPE NUM_Y = dXYZ[1];
			const INDEX_TYPE SLAB_SIZE = NUM_X * NUM_Y;


			//			INDEX_TYPE num_cells = mMesh->numCells();
#pragma omp parallel
			{
				// NOTE(8/5/2019): since we use thread per slice, if there is too few slices we may need to leave some threads idle
				int num_threads = std::min(omp_get_num_threads(), (int)xyz[2]);
				int thread_num = omp_get_thread_num();
				if (thread_num < num_threads) {
					// local storage
					ValInd* slab = new ValInd[SLAB_SIZE * 3];
					ValInd* slab_min = new ValInd[SLAB_SIZE * 3];

					//float* slab_vals = new FLOATTYPE[SLAB_SIZE * 3];
					//INDEX_TYPE* slab_inds = new INDEX_TYPE[SLAB_SIZE * 3];

					std::vector<INDEX_TYPE> partition;
					ArrayIndexPartitioner::EvenChunkSplit(xyz[2], num_threads, partition);

					//for (auto v : partition) printf("p %llu\n", v);
					const INDEX_TYPE START_K = partition[thread_num];
					INDEX_TYPE localend = (thread_num == num_threads - 1 ?
						partition[thread_num + 1] - 1 :
						partition[thread_num + 1]);
					const INDEX_TYPE END_K = localend;

					ValInd* old_slab = slab;
					ValInd* mid_slab = &(slab[SLAB_SIZE * 1]);
					ValInd* new_slab = &(slab[SLAB_SIZE * 2]);
				
					ValInd* old_slab_min = slab_min;
					ValInd* mid_slab_min = &(slab_min[SLAB_SIZE * 1]);
					ValInd* new_slab_min = &(slab_min[SLAB_SIZE * 2]);

					for (INDEX_TYPE k = START_K; k <= END_K; k++) {

						//float* old_slab_vals = slab_vals;
						//float* mid_slab_vals; &(slab_vals[SLAB_SIZE * 1);
						//float* new_slab_vals = &(slab_vals[SLAB_SIZE * 2]);
						//INDEX_TYPE* old_slab_inds = slab_inds;
						//INDEX_TYPE* mid_slab_inds; &(slab_inds[SLAB_SIZE * 1);
						//INDEX_TYPE* new_slab_inds = &(slab_inds[SLAB_SIZE * 2]);

						for (INDEX_TYPE j = 0; j < xyz[1]; j++) {
							INDEX_TYPE stick_start_j0k0 = j * 2 * NUM_X; // get the row in the slab

																			// J*2, K*2
																			// Do each vertex-edge-vertex-....-edge-vertex stick in the slab
							INDEX_TYPE base_cell_id = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2));
							INDEX_TYPE id_data = mGrid->Index3d(Vec3l(0, j, k));
							do_fastest_line(xyz[0], base_cell_id,
								&(mGridFunc->GetImage()[id_data]),
								&(new_slab[stick_start_j0k0]),
								&(new_slab_min[stick_start_j0k0]));
							write_max_stick(NUM_X, base_cell_id, &(new_slab[stick_start_j0k0]));
							write_min_stick(NUM_X, base_cell_id, &(new_slab_min[stick_start_j0k0]));

							// J*2-1, K*2 = mid of J*2-2, K*2
							// Do the edge-face-edge-...-face-edge stick in the slab
							INDEX_TYPE stick_start_j2k0 = (j * 2 - 2) * NUM_X; // get the row in the slab
							INDEX_TYPE stick_start_j1k0 = (j * 2 - 1) * NUM_X; // get the row in the slab
							INDEX_TYPE base_cell_j1_id = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2));

							if (j > 0) {
								// EDGE BETWEEN v(x, y:y-1, z)
								do_parallel_lines_max(xyz[0],
									&(new_slab[stick_start_j0k0]),
									&(new_slab[stick_start_j2k0]),
									&(new_slab[stick_start_j1k0]));
								write_max_stick(NUM_X, base_cell_j1_id, &(new_slab[stick_start_j1k0]));
								do_parallel_lines_min(xyz[0],
									&(new_slab_min[stick_start_j0k0]),
									&(new_slab_min[stick_start_j2k0]),
									&(new_slab_min[stick_start_j1k0]));
								write_min_stick(NUM_X, base_cell_j1_id, &(new_slab_min[stick_start_j1k0]));
							}

							// NOW FILL IN BETWEEN SLABS

							// J*2, K*2-1 = mid of J*2, K*2-2
							// Do the edge-face-edge-...-face-edge stick in between slabs
							if (k > START_K) {
								INDEX_TYPE base_cell_k1_id = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 1));
								do_parallel_lines_max(xyz[0],
									&(new_slab[stick_start_j0k0]),
									&(old_slab[stick_start_j0k0]),
									&(mid_slab[stick_start_j0k0]));
								write_max_stick(NUM_X, base_cell_k1_id, &(mid_slab[stick_start_j0k0]));
								do_parallel_lines_min(xyz[0],
									&(new_slab_min[stick_start_j0k0]),
									&(old_slab_min[stick_start_j0k0]),
									&(mid_slab_min[stick_start_j0k0]));
								write_min_stick(NUM_X, base_cell_k1_id, &(mid_slab_min[stick_start_j0k0]));

							}

							// J*2-1, K*2-1 = mid of J*2-1, K*2-2
							// Do the face-quad-face-....-quad-face stick between the slabs
							if (k > START_K && j > 0) {
								// face BETWEEN v(x, y:y-1, z:z-1)

								INDEX_TYPE base_mid_id = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2 - 1));
								do_parallel_lines_max(xyz[0],
									&(new_slab[stick_start_j1k0]),
									&(old_slab[stick_start_j1k0]),
									&(mid_slab[stick_start_j1k0]));
								write_max_stick(NUM_X, base_mid_id, &(mid_slab[stick_start_j1k0]));
								do_parallel_lines_min(xyz[0],
									&(new_slab_min[stick_start_j1k0]),
									&(old_slab_min[stick_start_j1k0]),
									&(mid_slab_min[stick_start_j1k0]));
								write_min_stick(NUM_X, base_mid_id, &(mid_slab_min[stick_start_j1k0]));
							}



						}

						// swap new and old slabs
						ValInd* tmp = old_slab;
						old_slab = new_slab;
						new_slab = tmp;
						tmp = old_slab_min;
						old_slab_min = new_slab_min;
						new_slab_min = tmp;
						//do_fastest_line(xyz[0], id_cell,
						//	&(mGridFunc->GetImage()[id_data]),
						//	&(tmp_vals[id_cell]),
						//	&(tmp_ids[id_cell]));

						//if (j > 0) {
						//	// EDGE BETWEEN v(x, y:y-1, z)
						//	INDEX_TYPE id_cell_last = mMesh->coords2Cellid(Vec3l(0, j * 2 - 2, k * 2));
						//	INDEX_TYPE id_cell_mid = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2));
						//	do_parallel_lines(xyz[0],
						//		&(tmp_vals[id_cell]),
						//		&(tmp_ids[id_cell]),
						//		&(tmp_vals[id_cell_last]),
						//		&(tmp_ids[id_cell_last]),
						//		&(tmp_vals[id_cell_mid]),
						//		&(tmp_ids[id_cell_mid]));
						//}
						//if (k > START_K) {
						//	// EDGE BETWEEN v(x, y, z:z-1)
						//	INDEX_TYPE id_cell_last = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 2));
						//	INDEX_TYPE id_cell_mid = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 1));
						//	do_parallel_lines(xyz[0],
						//		&(tmp_vals[id_cell]),
						//		&(tmp_ids[id_cell]),
						//		&(tmp_vals[id_cell_last]),
						//		&(tmp_ids[id_cell_last]),
						//		&(tmp_vals[id_cell_mid]),
						//		&(tmp_ids[id_cell_mid]));
						//}
						//if (k > START_K && j > 0) {
						//	// face BETWEEN v(x, y:y-1, z:z-1)
						//	INDEX_TYPE id_cell_e0 = mMesh->coords2Cellid(Vec3l(0, j * 2 - 2, k * 2 - 1));
						//	INDEX_TYPE id_cell_e1 = mMesh->coords2Cellid(Vec3l(0, j * 2, k * 2 - 1));
						//	INDEX_TYPE id_cell_mid = mMesh->coords2Cellid(Vec3l(0, j * 2 - 1, k * 2 - 1));
						//	do_parallel_lines(xyz[0],
						//		&(tmp_vals[id_cell_e0]),
						//		&(tmp_ids[id_cell_e0]),
						//		&(tmp_vals[id_cell_e1]),
						//		&(tmp_ids[id_cell_e1]),
						//		&(tmp_vals[id_cell_mid]),
						//		&(tmp_ids[id_cell_mid]));
						//}

					}

					//INDEX_TYPE num_cells = mMesh->numCells();
					//for (INDEX_TYPE i = 0; i < num_cells; i++) {
					//	byte_tmp_max_ids[i] = mMesh->CompressVertexOffsetToByte(i, tmp_ids[i]); // THIS COULD BE IMPROVED WITH BETTER HASH
					//}


					//			INDEX_TYPE num_cells = mMesh->numCells();
					//			mMaxVLabel = new DenseLabeling<DIM_TYPE>(num_cells);
					//			//mMaxVLabel->SetAll(0); // not needed since this is initialized to first vertex anyway
					//			//return;
					//#pragma omp parallel
					//			{
					//				int num_threads = omp_get_num_threads();
					//				int thread_num = omp_get_thread_num();
					//
					//				std::vector<INDEX_TYPE> partition;
					//				ArrayIndexPartitioner::EvenChunkSplit(num_cells, num_threads, partition);
					//
					//				for (INDEX_TYPE cellid = partition[thread_num]; cellid < partition[thread_num + 1]; cellid++) {
					//					mMaxVLabel->SetLabel(cellid, idxOf_maxVertex(cellid));
					//				}
					//			}
					//
					//			//MeshType::AllCellsIterator ait(mMesh);
					//			//for (ait.begin(); ait.valid(); ait.advance()){
					//			//	mMaxVLabel->SetLabel(ait.value(), idxOf_maxVertex(ait.value()));
					//			//}

					delete[] slab;
					delete[] slab_min;
				}
			}
			timer.EndGlobal();
			printf(" Done! ");
			timer.PrintAll();
		}




	};

	template <class MeshType, class GridFuncType>
	class RegularGridMaxMinVertexLabeling2D {
	protected:

		BYTE_TYPE* byte_tmp_max_ids;
		BYTE_TYPE* byte_tmp_min_ids;
		const RegularGrid2D* mGrid;
		MeshType* mMesh;
		GridFuncType* mGridFunc;
	public:

		RegularGridMaxMinVertexLabeling2D(MeshType* mesh, GridFuncType* gridFunc) :
			mMesh(mesh), mGridFunc(gridFunc) {
			byte_tmp_max_ids = NULL;
			byte_tmp_min_ids = NULL;
			mGrid = mGridFunc->GetGrid();
		}

		~RegularGridMaxMinVertexLabeling2D() {
			if (byte_tmp_max_ids != NULL) delete[] byte_tmp_max_ids;
			if (byte_tmp_min_ids != NULL) delete[] byte_tmp_min_ids;
		}

		void HACK_init() {
			byte_tmp_max_ids = new BYTE_TYPE[mMesh->numCells()];
			byte_tmp_min_ids = new BYTE_TYPE[mMesh->numCells()];
		}

		BYTE_TYPE GetUncompressedMaxVal(INDEX_TYPE cellid) const {
			return byte_tmp_max_ids[cellid];
		}
		BYTE_TYPE GetUncompressedMinVal(INDEX_TYPE cellid) const {
			return byte_tmp_min_ids[cellid];
		}
		void SetUncompressedMaxVal(INDEX_TYPE cellid, BYTE_TYPE val) {
			byte_tmp_max_ids[cellid] = val;
		}
		void SetUncompressedMinVal(INDEX_TYPE cellid, BYTE_TYPE val)  {
			byte_tmp_min_ids[cellid] = val;
		}
		INDEX_TYPE Cell2HighestVertex(INDEX_TYPE cellid) {
			return mMesh->UncompressByteToVertexOffset(cellid, byte_tmp_max_ids[cellid]);
		}
		INDEX_TYPE Cell2LowestVertex(INDEX_TYPE cellid) {
			auto val = mMesh->UncompressByteToVertexOffset(cellid, byte_tmp_min_ids[cellid]);
#ifdef DEBUG_ALL
			if (val < 0) {
				printf("whoatherenelly\n");
			}
#endif
			return val;
		}
		bool Before(INDEX_TYPE a, INDEX_TYPE b) const {
			return mGridFunc->IsGreater(b, a);
		}

		struct ValInd {
			float val;
			INDEX_TYPE ind;

			ValInd() {}
			ValInd(float v, INDEX_TYPE i) : val(v), ind(i) {}
		};


		void do_fastest_line(int X, const INDEX_TYPE base_offset, float* values, ValInd* o_max_row, ValInd* o_min_row) {
			// X = dim of line
			// float* values
			//out = float*, int*
			int dX = 2 * X - 1;

			// do 0th index
			//ValInd old(values[0], base_offset);
			float tmpold = values[0];
			o_max_row[0].ind = 0 + base_offset;
			o_max_row[0].val = tmpold;
			o_min_row[0].ind = 0 + base_offset;
			o_min_row[0].val = tmpold;

			// now do rest!
			for (int i = 1; i < X; i++) {
				// fisrt get tmp next value and set
				int di_v = i * 2; // mesh index of vertex
				float tmpnew = values[i];
				o_max_row[di_v].ind = base_offset + di_v;
				o_max_row[di_v].val = tmpnew;
				o_min_row[di_v].ind = base_offset + di_v;
				o_min_row[di_v].val = tmpnew;

				int di_e = i * 2 - 1; // mesh index of edge
				if (tmpold > tmpnew) {
					o_max_row[di_e].val = tmpold;
					o_max_row[di_e].ind = base_offset + di_v - 2; /// check this
					o_min_row[di_e].val = tmpnew;
					o_min_row[di_e].ind = base_offset + di_v; /// check this
				}
				else  { // tmpold < tmpnew || di_v-2 < di_v <- 2nd part is always true
					o_max_row[di_e].val = tmpnew;
					o_max_row[di_e].ind = base_offset + di_v; /// check this
					o_min_row[di_e].val = tmpold;
					o_min_row[di_e].ind = base_offset + di_v - 2; /// check this
				}


				//
				tmpold = tmpnew;
			}

		}

		void do_parallel_lines_max(int X,
			const ValInd* row0,
			const ValInd* row1,
			ValInd* row_o) {

			int dX = 2 * X - 1;

			// now do rest!
			for (int i = 0; i < dX; i++) {
				// fisrt get tmp next value and set
				float v0 = row0[i].val;
				float v1 = row1[i].val;

				if (v0 > v1) {
					row_o[i].val = v0;
					row_o[i].ind = row0[i].ind;
					continue;
				}
				else if (v1 > v0) {
					row_o[i].val = v1;
					row_o[i].ind = row1[i].ind;
					continue;
				}

				// need index compare			
				INDEX_TYPE i0 = row0[i].ind;
				INDEX_TYPE i1 = row1[i].ind;
				if (i0 > i1) {
					row_o[i].val = v0;
					row_o[i].ind = i0;
				}
				else {
					row_o[i].val = v1;
					row_o[i].ind = i1;
				}
			}
		}
		void do_parallel_lines_min(int X,
			const ValInd* row0,
			const ValInd* row1,
			ValInd* row_o) {

			int dX = 2 * X - 1;

			// now do rest!
			for (int i = 0; i < dX; i++) {
				// fisrt get tmp next value and set
				float v0 = row0[i].val;
				float v1 = row1[i].val;

				if (v0 < v1) {
					row_o[i].val = v0;
					row_o[i].ind = row0[i].ind;
					continue;
				}
				else if (v1 < v0) {
					row_o[i].val = v1;
					row_o[i].ind = row1[i].ind;
					continue;
				}

				// need index compare			
				INDEX_TYPE i0 = row0[i].ind;
				INDEX_TYPE i1 = row1[i].ind;
				if (i0 < i1) {
					row_o[i].val = v0;
					row_o[i].ind = i0;
				}
				else {
					row_o[i].val = v1;
					row_o[i].ind = i1;
				}
			}
		}

		void write_max_stick(int dX, const INDEX_TYPE base_offset, ValInd* stick) {

			for (INDEX_TYPE i = 0; i < dX; i++) {
				INDEX_TYPE base_id = i + base_offset;
				byte_tmp_max_ids[base_id] = mMesh->CompressVertexOffsetToByte(base_id, stick[i].ind);
			}

		}
		void write_min_stick(int dX, const INDEX_TYPE base_offset, ValInd* stick) {

			for (INDEX_TYPE i = 0; i < dX; i++) {
				INDEX_TYPE base_id = i + base_offset;
				byte_tmp_min_ids[base_id] = mMesh->CompressVertexOffsetToByte(base_id, stick[i].ind);
			}

		}
		void ComputeOutput() {

			ThreadedTimer timer(1);
			timer.StartGlobal();

			printf(" -- Creating maxV_labeling ...");
			fflush(stdout);

			byte_tmp_max_ids = new BYTE_TYPE[mMesh->numCells()];
			byte_tmp_min_ids = new BYTE_TYPE[mMesh->numCells()];
#ifdef DEBUG_ALL
			for (int i = 0; i < mMesh->numCells(); i++) byte_tmp_max_ids[i] = byte_tmp_min_ids[i] = 213;
#endif
			const Vec2l xy = mGrid->XY();
			const Vec2l dXY = mMesh->XY();
			const INDEX_TYPE NUM_X = dXY[0];
			const INDEX_TYPE NUM_Y = dXY[1];
			const INDEX_TYPE SLAB_SIZE = NUM_X ;


			//			INDEX_TYPE num_cells = mMesh->numCells();
#pragma omp parallel
			{
				int num_threads = omp_get_num_threads();
				int thread_num = omp_get_thread_num();
				assert(xy[1] >= num_threads);

				// local storage
				ValInd* slab = new ValInd[SLAB_SIZE * 3];
				ValInd* slab_min = new ValInd[SLAB_SIZE * 3];

				//float* slab_vals = new FLOATTYPE[SLAB_SIZE * 3];
				//INDEX_TYPE* slab_inds = new INDEX_TYPE[SLAB_SIZE * 3];

				std::vector<INDEX_TYPE> partition;
				ArrayIndexPartitioner::EvenChunkSplit(xy[1], num_threads, partition);

				//for (auto v : partition) printf("p %llu\n", v);
				const INDEX_TYPE START_J = partition[thread_num];
				INDEX_TYPE localend = (thread_num == num_threads - 1 ?
					partition[thread_num + 1] - 1 :
					partition[thread_num + 1]);
				const INDEX_TYPE END_J = localend;

				ValInd* old_slab = slab;
				ValInd* mid_slab = &(slab[SLAB_SIZE * 1]);
				ValInd* new_slab = &(slab[SLAB_SIZE * 2]);

				ValInd* old_slab_min = slab_min;
				ValInd* mid_slab_min = &(slab_min[SLAB_SIZE * 1]);
				ValInd* new_slab_min = &(slab_min[SLAB_SIZE * 2]);

				for (INDEX_TYPE j = START_J; j <= END_J; j++) {



					// J*2
					// Do each vertex-edge-vertex-....-edge-vertex stick in the slab
					INDEX_TYPE base_cell_id = mMesh->coords2Cellid(Vec2l(0, j * 2));
					INDEX_TYPE id_data = mGrid->Index2d(Vec2l(0, j));
					do_fastest_line(xy[0], base_cell_id,
						&(mGridFunc->GetImage()[id_data]),
						new_slab, new_slab_min);
					write_max_stick(NUM_X, base_cell_id, new_slab);
					write_min_stick(NUM_X, base_cell_id, new_slab_min);

					if (j > START_J) {
						// face BETWEEN v(x, y:y-1, z:z-1)

						INDEX_TYPE base_mid_id = mMesh->coords2Cellid(Vec2l(0, j * 2 - 1));
						do_parallel_lines_max(xy[0],
							new_slab, old_slab, mid_slab);
						write_max_stick(NUM_X, base_mid_id, mid_slab);
						do_parallel_lines_min(xy[0],
							new_slab_min, old_slab_min, mid_slab_min);
						write_min_stick(NUM_X, base_mid_id, mid_slab_min);

					}





					// swap new and old slabs
					ValInd* tmp = old_slab;
					old_slab = new_slab;
					new_slab = tmp;
					tmp = old_slab_min;
					old_slab_min = new_slab_min;
					new_slab_min = tmp;


				}


				delete[] slab;
				delete[] slab_min;
			}
			timer.EndGlobal();
			printf(" Done! ");
			timer.PrintAll();
#ifdef DEBUG_ALL
			for (INDEX_TYPE i = 0; i < mMesh->numCells(); i++) {
				Vec2l coord;
				mMesh->cellid2Coords(i, coord);
				auto dim = mMesh->dimension(coord);
				if (byte_tmp_max_ids[i] > 9) {
					auto v = byte_tmp_max_ids[i];
					printf("asdfasdf;asdf\n");
				}
				if (byte_tmp_min_ids[i] > 9) {
					auto v = byte_tmp_min_ids[i];
					printf("asdfasdf;asdf2\n");

				}
			}
#endif

		}




	};

}

#endif
