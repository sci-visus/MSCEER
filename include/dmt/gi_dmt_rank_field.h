#ifndef GI_DMT_RANK_FIELD_H
#define GI_DMT_RANK_FIELD_H

// A scalar field given as a total order of its vertices.
//
// RankField3D holds one uint32 rank per vertex (RegularGrid3D order, x fastest)
// and compares integers, so it is exact for any grid size; a float image is
// not once ranks pass 2^24.
//
// RankMaxMinVertexLabeling3D is a drop-in for RegularGridMaxMinVertexLabeling3D
// (gi_max_vertex_labeling.h) driven by a RankField3D instead of a float image:
// it provides the interface that CLSGradient and CLSSlidingGradient
// use (Before, Cell2HighestVertex, Cell2LowestVertex, GetUncompressedMaxVal,
// GetUncompressedMinVal, ComputeOutput), so MSCEER's own gradient rules can run
// on exact ranks too.

#include "gi_basic_types.h"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <vector>

namespace GInt {
	namespace dmt {

		class RankField3D {
		public:
			int64_t X = 0, Y = 0, Z = 0;
			std::vector<uint32_t> ranks;    // ranks[x + X*(y + Y*z)]

			RankField3D() {}
			RankField3D(int64_t x, int64_t y, int64_t z) : X(x), Y(y), Z(z), ranks((size_t)(x * y * z)) {}

			int64_t NumVertices() const { return X * Y * Z; }
			bool IsGreater(INDEX_TYPE a, INDEX_TYPE b) const { return ranks[(size_t)a] > ranks[(size_t)b]; }

			// Ranks from arbitrary values given in vertex-id order (x fastest).
			// Ties are broken by `tie_order`: tie_order[i] is the position of
			// vertex i in the tie-breaking order (pass nullptr to break ties by
			// vertex id, which is MSCEER's IsGreater convention).
			template <class T>
			static RankField3D FromValues(const T* values, int64_t x, int64_t y, int64_t z,
				const int64_t* tie_order = nullptr) {
				RankField3D f(x, y, z);
				const int64_t n = f.NumVertices();
				std::vector<int64_t> idx((size_t)n);
				std::iota(idx.begin(), idx.end(), (int64_t)0);
				std::stable_sort(idx.begin(), idx.end(), [&](int64_t a, int64_t b) {
					if (values[a] < values[b]) return true;
					if (values[b] < values[a]) return false;
					const int64_t ta = tie_order ? tie_order[a] : a;
					const int64_t tb = tie_order ? tie_order[b] : b;
					return ta < tb;
					});
				for (int64_t r = 0; r < n; r++) f.ranks[(size_t)idx[(size_t)r]] = (uint32_t)r;
				return f;
			}

			// The first pair of equal values (vertex ids a < b) or false if all
			// values are distinct.
			template <class T>
			static bool FirstDuplicate(const T* values, int64_t n, int64_t& a, int64_t& b) {
				std::vector<int64_t> idx((size_t)n);
				std::iota(idx.begin(), idx.end(), (int64_t)0);
				std::stable_sort(idx.begin(), idx.end(), [&](int64_t i, int64_t j) { return values[i] < values[j]; });
				bool found = false;
				for (int64_t r = 1; r < n; r++) {
					const int64_t i = idx[(size_t)r - 1], j = idx[(size_t)r];
					if (!(values[i] < values[j]) && !(values[j] < values[i])) {
						const int64_t lo = i < j ? i : j, hi = i < j ? j : i;
						if (!found || lo < a || (lo == a && hi < b)) { a = lo; b = hi; found = true; }
					}
				}
				return found;
			}
		};

		template <class MeshType>
		class RankMaxMinVertexLabeling3D {
		protected:
			MeshType* mMesh;
			const RankField3D* mField;
			std::vector<BYTE_TYPE> mMax, mMin;

		public:
			RankMaxMinVertexLabeling3D(MeshType* mesh, const RankField3D* field) : mMesh(mesh), mField(field) {}

			void ComputeOutput() {
				const INDEX_TYPE n = mMesh->numCells();
				mMax.assign((size_t)n, 0);
				mMin.assign((size_t)n, 0);
#pragma omp parallel for schedule(static)
				for (long long c = 0; c < (long long)n; c++) {
					typename MeshType::CellVerticesIterator vit(mMesh);
					INDEX_TYPE best_hi = -1, best_lo = -1;
					for (vit.begin(c); vit.valid(); vit.advance()) {
						const INDEX_TYPE vc = vit.value();
						const INDEX_TYPE vn = mMesh->VertexNumberFromCellID(vc);
						if (best_hi < 0 || mField->IsGreater(vn, mMesh->VertexNumberFromCellID(best_hi))) best_hi = vc;
						if (best_lo < 0 || mField->IsGreater(mMesh->VertexNumberFromCellID(best_lo), vn)) best_lo = vc;
					}
					mMax[(size_t)c] = mMesh->CompressVertexOffsetToByte(c, best_hi);
					mMin[(size_t)c] = mMesh->CompressVertexOffsetToByte(c, best_lo);
				}
			}

			BYTE_TYPE GetUncompressedMaxVal(INDEX_TYPE cellid) const { return mMax[(size_t)cellid]; }
			BYTE_TYPE GetUncompressedMinVal(INDEX_TYPE cellid) const { return mMin[(size_t)cellid]; }
			void SetUncompressedMaxVal(INDEX_TYPE cellid, BYTE_TYPE v) { mMax[(size_t)cellid] = v; }
			void SetUncompressedMinVal(INDEX_TYPE cellid, BYTE_TYPE v) { mMin[(size_t)cellid] = v; }
			INDEX_TYPE Cell2HighestVertex(INDEX_TYPE cellid) const {
				return mMesh->UncompressByteToVertexOffset(cellid, mMax[(size_t)cellid]);
			}
			INDEX_TYPE Cell2LowestVertex(INDEX_TYPE cellid) const {
				return mMesh->UncompressByteToVertexOffset(cellid, mMin[(size_t)cellid]);
			}
			// a, b are vertex numbers: a comes before b in the order.
			bool Before(INDEX_TYPE a, INDEX_TYPE b) const { return mField->IsGreater(b, a); }
		};

	} // namespace dmt
} // namespace GInt

#endif
