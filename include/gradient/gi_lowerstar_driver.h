#ifndef GI_LOWERSTAR_DRIVER_H
#define GI_LOWERSTAR_DRIVER_H

// Host drivers shared by the lower-star gradient rules built on the star
// scaffold of dmt_core/dmt_lowerstar.h. A rule is a type with
//
//   // per-star mask handed to Star as bnd (0 when the rule does not use one)
//   static uint32_t Boundary(const ::dmt::Star27Tables& T, const ::dmt::Cubical3D& G,
//                            int64_t x, int64_t y, int64_t z);
//   static void Star(const ::dmt::Star27Tables& T, const uint32_t key[27],
//                    uint32_t vvalid, uint32_t bnd, ::dmt::StarResult& R,
//                    bool aligned_first);
//
// The rules themselves are in gi_pls_gradient.h and gi_pls_br_gradient.h.
//
// Two drivers with identical results:
//   LowerStarGradient3D<Rule>       flat uint32 ranks in, GradBitfield bytes
//                                   out; no MSCEER mesh is built.
//   LowerStarPQGradient<..., Rule>  the same rule on an MSCEER mesh, max/min
//                                   vertex labeling and gradient labeling
//                                   (ComputeLowerStar per vertex and an OpenMP
//                                   ComputePairing), for MSCEER pipelines.
//
// Every cell lies in exactly one lower star (that of its highest vertex), so
// stars are processed independently and write disjoint cells.

#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "dmt_core/dmt_cubical3d.h"
#include "dmt_core/dmt_hd.h"
#include "dmt_core/dmt_lowerstar.h"
#include "dmt_core/dmt_star27_tables.h"

#include <cstdint>
#include <cstring>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GInt {
	namespace dmt {

		struct LowerStarStats {
			int64_t critical[4] = { 0, 0, 0, 0 };
			int64_t orphans = 0;
			int64_t unassigned = 0;   // cells no star wrote; 0 on success
		};

		inline int ResolveThreads(int num_threads) {
#ifdef _OPENMP
			return num_threads > 0 ? num_threads : omp_get_max_threads();
#else
			(void)num_threads;
			return 1;
#endif
		}

		// Writes the result of one star into per-cell arrays.
		inline void WriteStar(const ::dmt::Cubical3D& G, const ::dmt::Star27Tables& T,
			int64_t x, int64_t y, int64_t z, const ::dmt::StarResult& R,
			uint8_t* out_grad, uint8_t* out_step, uint8_t* out_part) {
			for (int s = 0; s < 27; s++) {
				if (R.pair[s] == -2) continue;
				const int64_t id = G.Id(2 * x + T.off[s][0], 2 * y + T.off[s][1], 2 * z + T.off[s][2]);
				const int code = R.pair[s] < 0 ? ::dmt::kPairCritical : ::dmt::StarDir(T, s, R.pair[s]);
				out_grad[id] = ::dmt::MakeGradByte((uint8_t)code);
				if (out_step) out_step[id] = R.step[s];
				if (out_part) out_part[id] = R.part[s];
			}
		}

		// ranks: one per vertex, x fastest, pairwise distinct (only their order
		// matters). out_grad: (2X-1)(2Y-1)(2Z-1) bytes, x fastest. out_step and
		// out_part are optional (nullptr) per-cell journal arrays. Returns true
		// iff every cell was assigned.
		template <class Rule>
		inline bool LowerStarGradient3D(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
			bool aligned_first, uint8_t* out_grad, uint8_t* out_step, uint8_t* out_part,
			LowerStarStats* stats, int num_threads = 0) {
			const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(X, Y, Z);
			::dmt::Star27Tables T;
			::dmt::BuildStar27TablesCanonical(T);
			const int64_t nc = G.NumCells();
			std::memset(out_grad, 0, (size_t)nc);
			const long long nv = (long long)G.NumVertices();
			long long orphans = 0;
			const int nt = ResolveThreads(num_threads);
#pragma omp parallel for schedule(dynamic, 256) reduction(+:orphans) num_threads(nt)
			for (long long v = 0; v < nv; v++) {
				const int64_t x = v % X, y = (v / X) % Y, z = v / (X * Y);
				uint32_t key[27], vvalid;
				::dmt::GatherStarKeys(G, ranks, x, y, z, key, vvalid);
				::dmt::StarResult R;
				Rule::Star(T, key, vvalid, Rule::Boundary(T, G, x, y, z), R, aligned_first);
				orphans += R.n_orphans;
				WriteStar(G, T, x, y, z, R, out_grad, out_step, out_part);
			}
			LowerStarStats st;
			st.orphans = orphans;
			for (int64_t c = 0; c < nc; c++) {
				if (!::dmt::GradAssigned(out_grad[c])) { st.unassigned++; continue; }
				if (::dmt::GradCritical(out_grad[c])) st.critical[G.Dim(c)]++;
			}
			if (stats) *stats = st;
			return st.unassigned == 0;
		}

		// ---- template driver for MSCEER pipelines ----

		template <class MeshType, class MaxVLType, class GradType, class Rule>
		class LowerStarPQGradient {
		protected:
			MeshType* mMesh;
			MaxVLType* mMaxVL;
			GradType* mGrad;
			uint8_t* mStep;
			uint8_t* mPart;
			bool mAlignedFirst;
			::dmt::Star27Tables mT;
			Vec3l mCells;        // cells per axis
			::dmt::Cubical3D mG;
			int64_t mOrphans = 0;

		public:
			LowerStarPQGradient(MeshType* mesh, MaxVLType* maxvl, GradType* grad,
				uint8_t* step_out = nullptr, uint8_t* part_out = nullptr, bool aligned_first = false)
				: mMesh(mesh), mMaxVL(maxvl), mGrad(grad), mStep(step_out), mPart(part_out),
				mAlignedFirst(aligned_first) {
				::dmt::BuildStar27TablesCanonical(mT);
				mCells = mMesh->XYZ();
				mG = ::dmt::Cubical3D::Make((mCells[0] + 1) / 2, (mCells[1] + 1) / 2, (mCells[2] + 1) / 2);
			}

			int64_t Orphans() const { return mOrphans; }
			bool AlignedFirst() const { return mAlignedFirst; }

			// Lower star of the vertex whose cell id is vert_cell. Returns the
			// number of orphans in it (0 on a lower star).
			int ComputeLowerStar(INDEX_TYPE vert_cell) {
				Vec3l vc;
				mMesh->cellid2Coords(vert_cell, vc);
				// Local order of the 27 vertices around vc via the labeling's
				// Before(): key[j] = number of valid neighbours that come first.
				INDEX_TYPE vnum[27];
				uint32_t vvalid = 0;
				for (int dz = -1; dz <= 1; dz++)
					for (int dy = -1; dy <= 1; dy++)
						for (int dx = -1; dx <= 1; dx++) {
							const int j = ::dmt::StarSlot(dx, dy, dz);
							const long long q[3] = { vc[0] + 2 * dx, vc[1] + 2 * dy, vc[2] + 2 * dz };
							if (q[0] < 0 || q[1] < 0 || q[2] < 0 || q[0] >= mCells[0] || q[1] >= mCells[1] || q[2] >= mCells[2]) {
								vnum[j] = -1;
								continue;
							}
							vnum[j] = mMesh->VertexNumberFromCellID(mMesh->coords2Cellid(Vec3l(q[0], q[1], q[2])));
							vvalid |= 1u << j;
						}
				uint32_t key[27];
				for (int j = 0; j < 27; j++) {
					key[j] = 0;
					if (!((vvalid >> j) & 1u)) continue;
					for (int i = 0; i < 27; i++)
						if (i != j && ((vvalid >> i) & 1u) && mMaxVL->Before(vnum[i], vnum[j])) key[j]++;
				}
				::dmt::StarResult R;
				const uint32_t bnd = Rule::Boundary(mT, mG, vc[0] / 2, vc[1] / 2, vc[2] / 2);
				Rule::Star(mT, key, vvalid, bnd, R, mAlignedFirst);
				for (int s = 0; s < 27; s++) {
					if (R.pair[s] == -2) continue;
					const INDEX_TYPE id = mMesh->coords2Cellid(Vec3l(vc[0] + mT.off[s][0], vc[1] + mT.off[s][1], vc[2] + mT.off[s][2]));
					if (R.pair[s] < 0) {
						mGrad->setCritical(id, true);
					}
					else {
						const int p = R.pair[s];
						const INDEX_TYPE pid = mMesh->coords2Cellid(Vec3l(vc[0] + mT.off[p][0], vc[1] + mT.off[p][1], vc[2] + mT.off[p][2]));
						mGrad->setPair(id, pid);
					}
					mGrad->setAssigned(id, true);
					if (mStep) mStep[id] = R.step[s];
					if (mPart) mPart[id] = R.part[s];
				}
				return R.n_orphans;
			}

			void ComputePairing(int num_threads = 0) {
				const long long X = (mCells[0] + 1) / 2, Y = (mCells[1] + 1) / 2, Z = (mCells[2] + 1) / 2;
				const long long nv = X * Y * Z;
				long long orphans = 0;
				const int nt = ResolveThreads(num_threads);
#pragma omp parallel for schedule(dynamic, 256) reduction(+:orphans) num_threads(nt)
				for (long long v = 0; v < nv; v++) {
					const long long x = v % X, y = (v / X) % Y, z = v / (X * Y);
					orphans += ComputeLowerStar(mMesh->coords2Cellid(Vec3l(2 * x, 2 * y, 2 * z)));
				}
				mOrphans = orphans;
			}
		};

	} // namespace dmt
} // namespace GInt

#endif
