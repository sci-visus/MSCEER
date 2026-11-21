#ifndef GI_REGULAR_GRID_MAXMIN_VERTEX_LABELING_2D_H
#define GI_REGULAR_GRID_MAXMIN_VERTEX_LABELING_2D_H

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

		// ---- lazy (unmaterialized) mode -----------------------------------
		// ComputeOutput() is what allocates byte_tmp_*_ids, so a caller that
		// never calls it leaves them NULL and every query below answers on
		// demand instead. That is the right trade when the labeling is read at
		// a SPARSE set of cells (the mesh function's cellValue over critical
		// cells is ~4% of the mesh) rather than per cell, as CLSGradient
		// reads it -- the same reasoning as LazyMaximumVertexLabeling above,
		// but kept inside this class so the choice stays a RUNTIME one: the
		// gradient override that makes CLS unnecessary can fail, and then
		// the arrays are needed after all.
		//
		// The eager byte is a slot index that UncompressByteToVertexOffset
		// immediately turns back into a vertex cell id, so returning that id
		// directly skips a hash probe per improvement
		// (CompressVertexOffsetToByte is an unordered_map lookup) without
		// changing the answer.
		INDEX_TYPE computeExtremeVertex(INDEX_TYPE cellid, bool highest) const {
			typename MeshType::CellVerticesIterator cviter(mMesh);
			cviter.begin(cellid);
			INDEX_TYPE bestCell = cviter.value();
			INDEX_TYPE bestGid = mMesh->VertexNumberFromCellID(bestCell);
			cviter.advance();
			for (; cviter.valid(); cviter.advance()) {
				const INDEX_TYPE otherCell = cviter.value();
				const INDEX_TYPE otherGid = mMesh->VertexNumberFromCellID(otherCell);
				const bool better = highest ? mGridFunc->IsGreater(otherGid, bestGid)
					: mGridFunc->IsGreater(bestGid, otherGid);
				if (better) { bestCell = otherCell; bestGid = otherGid; }
			}
			return bestCell;
		}

		bool IsMaterialized() const { return byte_tmp_max_ids != NULL; }

		BYTE_TYPE GetUncompressedMaxVal(INDEX_TYPE cellid) const {
			if (byte_tmp_max_ids == NULL)
				return mMesh->CompressVertexOffsetToByte(cellid, computeExtremeVertex(cellid, true));
			return byte_tmp_max_ids[cellid];
		}
		BYTE_TYPE GetUncompressedMinVal(INDEX_TYPE cellid) const {
			if (byte_tmp_min_ids == NULL)
				return mMesh->CompressVertexOffsetToByte(cellid, computeExtremeVertex(cellid, false));
			return byte_tmp_min_ids[cellid];
		}
		void SetUncompressedMaxVal(INDEX_TYPE cellid, BYTE_TYPE val) {
			byte_tmp_max_ids[cellid] = val;
		}
		void SetUncompressedMinVal(INDEX_TYPE cellid, BYTE_TYPE val)  {
			byte_tmp_min_ids[cellid] = val;
		}
		INDEX_TYPE Cell2HighestVertex(INDEX_TYPE cellid) {
			if (byte_tmp_max_ids == NULL) return computeExtremeVertex(cellid, true);
			return mMesh->UncompressByteToVertexOffset(cellid, byte_tmp_max_ids[cellid]);
		}
		INDEX_TYPE Cell2LowestVertex(INDEX_TYPE cellid) {
			if (byte_tmp_min_ids == NULL) return computeExtremeVertex(cellid, false);
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
				// NOTE: since we use thread per row, if there are too few rows we leave excess threads idle
				int num_threads = std::min(omp_get_num_threads(), (int)xy[1]);
				int thread_num = omp_get_thread_num();
				if (thread_num < num_threads) {

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
					//write_max_stick(NUM_X, base_cell_id, new_slab);
					//write_min_stick(NUM_X, base_cell_id, new_slab_min);

					if (j > START_J) {
						// face BETWEEN v(x, y:y-1, z:z-1)

						//INDEX_TYPE base_mid_id = mMesh->coords2Cellid(Vec2l(0, j * 2 - 1));
						do_parallel_lines_max(xy[0],
							new_slab, old_slab, mid_slab);
						//write_max_stick(NUM_X, base_mid_id, mid_slab);
						do_parallel_lines_min(xy[0],
							new_slab_min, old_slab_min, mid_slab_min);
						//write_min_stick(NUM_X, base_mid_id, mid_slab_min);

					}

#pragma omp critical
					{
						write_max_stick(NUM_X, base_cell_id, new_slab);
						write_min_stick(NUM_X, base_cell_id, new_slab_min);
						if (j > START_J) {
							INDEX_TYPE base_mid_id = mMesh->coords2Cellid(Vec2l(0, j * 2 - 1));
							write_max_stick(NUM_X, base_mid_id, mid_slab);
							write_min_stick(NUM_X, base_mid_id, mid_slab_min);
						}
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
