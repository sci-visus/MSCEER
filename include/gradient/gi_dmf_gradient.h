#ifndef GI_DMT_DMF_GRADIENT_H
#define GI_DMT_DMF_GRADIENT_H

// Host driver for the Shivashankar-Natarajan discrete gradient
// (ASSIGNGRADIENT, optionally followed by ASSIGNGRADIENT2), using the
// per-cell functions of dmt_core/dmt_dmf.h.
//
//   Algorithm 1: each cell a of dimension < dmax is paired with the lowest-key
//                cofacet b whose highest facet is a.
//   Algorithm 2: each cell a still unpaired is paired with the lowest-key
//                cofacet b whose second-highest facet is a.
// A claim is refused (and counted as guard1 / guard2) when a or b is already
// paired; no other candidate is tried. Unpaired cells are critical.
//
// The published order visits cells by (dimension, key). Because a cofacet's
// highest and second-highest facets are unique, two cells of one dimension
// never claim the same cofacet, and a claim in dimension d only reads cells
// that dimension d-1 may have written. Each algorithm therefore runs as one
// parallel pass per dimension, d = 0, 1, ..., with exactly the sequential
// result; `sequential_reference` runs the literal (dimension, key)-ordered
// loop instead and is kept to check that.

#include "gi_basic_types.h"
#include "dmt_core/dmt_cubical3d.h"
#include "dmt_core/dmt_dmf.h"
#include "dmt_core/dmt_hd.h"
#include "gi_lowerstar_driver.h"   // ResolveThreads

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <vector>

namespace GInt {
	namespace dmt {

		struct DMFStats {
			int64_t guard1 = 0, guard2 = 0, alg1 = 0, alg2 = 0;
			int64_t critical[4] = { 0, 0, 0, 0 };
		};

		struct DMFOptions {
			bool second_pass = true;
			bool sequential_reference = false;
			int num_threads = 0;
		};

		namespace detail {

			inline int MaxDim(const ::dmt::Cubical3D& G) {
				int d = 0;
				for (int i = 0; i < 3; i++)
					if (G.n[i] > 1) d++;
				return d;
			}

			// (dimension, key) order of all cells.
			inline std::vector<int64_t> DMFOrder(const ::dmt::Cubical3D& G, const uint32_t* ranks) {
				const int64_t nc = G.NumCells();
				std::vector<int64_t> ids((size_t)nc);
				std::iota(ids.begin(), ids.end(), (int64_t)0);
				std::vector<uint32_t> keys((size_t)nc * 8);
				std::vector<uint8_t> dims((size_t)nc);
#pragma omp parallel for schedule(static)
				for (long long c = 0; c < (long long)nc; c++) {
					dims[(size_t)c] = (uint8_t)G.Dim(c);
					uint32_t k[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
					::dmt::DMFKey(G, ranks, c, k);
					for (int i = 0; i < 8; i++) keys[(size_t)c * 8 + i] = k[i];
				}
				std::sort(ids.begin(), ids.end(), [&](int64_t a, int64_t b) {
					if (dims[(size_t)a] != dims[(size_t)b]) return dims[(size_t)a] < dims[(size_t)b];
					return ::dmt::DMFKeyCmp(&keys[(size_t)a * 8], &keys[(size_t)b * 8], 1 << dims[(size_t)a]) < 0;
					});
				return ids;
			}

		} // namespace detail

		// ranks: one per vertex, x fastest, pairwise distinct. out_grad: one
		// GradBitfield byte per cell. out_stage (optional): 1 = paired by
		// Algorithm 1, 2 = by Algorithm 2, 0 = critical, set on both cells of a
		// pair. out_order (optional): each cell's position in the (dimension,
		// key) order, from which the published decision order follows.
		inline bool ComputeDMFGradient3D(const uint32_t* ranks, int64_t X, int64_t Y, int64_t Z,
			const DMFOptions& opt, uint8_t* out_grad, uint8_t* out_stage, int64_t* out_order, DMFStats* stats) {
			const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(X, Y, Z);
			const int64_t nc = G.NumCells();
			const int dmax = detail::MaxDim(G);
			const int nt = ResolveThreads(opt.num_threads);
			const uint8_t kUnpaired = 0;
			std::memset(out_grad, kUnpaired, (size_t)nc);
			if (out_stage) std::memset(out_stage, 0, (size_t)nc);

			// g0 / hi2 direction of every cell, (g0 | hi2 << 4); 0x77 for a vertex.
			std::vector<uint8_t> top((size_t)nc);
#pragma omp parallel for schedule(static) num_threads(nt)
			for (long long c = 0; c < (long long)nc; c++) {
				uint8_t g0 = 7, h2 = 7;
				::dmt::DMFTopFacets(G, ranks, c, g0, h2);
				top[(size_t)c] = (uint8_t)(g0 | (h2 << 4));
			}

			auto paired = [&](int64_t c) { return out_grad[c] != kUnpaired; };
			auto do_pair = [&](int64_t a, int64_t b, uint8_t stage) {
				out_grad[a] = ::dmt::MakeGradByte((uint8_t)G.DirBetween(a, b));
				out_grad[b] = ::dmt::MakeGradByte((uint8_t)G.DirBetween(b, a));
				if (out_stage) { out_stage[a] = stage; out_stage[b] = stage; }
			};

			long long g1 = 0, g2 = 0, a1 = 0, a2 = 0;
			std::vector<int64_t> order;
			if (opt.sequential_reference || out_order) order = detail::DMFOrder(G, ranks);

			if (opt.sequential_reference) {
				for (int alg = 1; alg <= (opt.second_pass ? 2 : 1); alg++) {
					for (int64_t a : order) {
						if (G.Dim(a) >= dmax) continue;
						if (alg == 2 && paired(a)) continue;
						const int64_t b = ::dmt::DMFLowestSelectingCofacet(G, ranks, top.data(), a, alg - 1);
						if (b < 0) continue;
						if (paired(a) || paired(b)) { (alg == 1 ? g1 : g2)++; continue; }
						do_pair(a, b, (uint8_t)alg);
						(alg == 1 ? a1 : a2)++;
					}
				}
			}
			else {
				// Cells of each dimension, in id order (the result does not depend on it).
				std::vector<std::vector<int64_t>> bydim(4);
				for (int64_t c = 0; c < nc; c++) bydim[(size_t)G.Dim(c)].push_back(c);
				for (int alg = 1; alg <= (opt.second_pass ? 2 : 1); alg++) {
					for (int d = 0; d < dmax; d++) {
						const std::vector<int64_t>& cells = bydim[(size_t)d];
						const long long n = (long long)cells.size();
						long long gg = 0, aa = 0;
#pragma omp parallel for schedule(static) reduction(+:gg, aa) num_threads(nt)
						for (long long i = 0; i < n; i++) {
							const int64_t a = cells[(size_t)i];
							if (alg == 2 && paired(a)) continue;
							const int64_t b = ::dmt::DMFLowestSelectingCofacet(G, ranks, top.data(), a, alg - 1);
							if (b < 0) continue;
							if (paired(a) || paired(b)) { gg++; continue; }
							do_pair(a, b, (uint8_t)alg);
							aa++;
						}
						if (alg == 1) { g1 += gg; a1 += aa; }
						else { g2 += gg; a2 += aa; }
					}
				}
			}

			DMFStats st;
			st.guard1 = g1; st.guard2 = g2; st.alg1 = a1; st.alg2 = a2;
			for (int64_t c = 0; c < nc; c++) {
				if (!paired(c)) {
					out_grad[c] = ::dmt::MakeGradByte(::dmt::kPairCritical);
					st.critical[G.Dim(c)]++;
				}
			}
			if (out_order)
				for (int64_t i = 0; i < (int64_t)order.size(); i++) out_order[order[(size_t)i]] = i;
			if (stats) *stats = st;
			return true;
		}

		// The same rule in the shape of the MSCEER gradient builders: global
		// ranks come from the max/min vertex labeling's Before(), then the raw
		// driver runs and the result is written into the gradient labeling.
		template <class MeshType, class MaxVLType, class GradType>
		class DMFGradient {
		protected:
			MeshType* mMesh;
			MaxVLType* mMaxVL;
			GradType* mGrad;
			DMFStats mStats;

		public:
			DMFGradient(MeshType* mesh, MaxVLType* maxvl, GradType* grad) : mMesh(mesh), mMaxVL(maxvl), mGrad(grad) {}

			const DMFStats& Stats() const { return mStats; }

			void ComputePairing(bool second_pass = true, int num_threads = 0) {
				const Vec3l xyz = mMesh->XYZ();
				const int64_t X = (xyz[0] + 1) / 2, Y = (xyz[1] + 1) / 2, Z = (xyz[2] + 1) / 2;
				const int64_t nv = X * Y * Z;
				std::vector<int64_t> vs((size_t)nv);
				std::iota(vs.begin(), vs.end(), (int64_t)0);
				std::sort(vs.begin(), vs.end(), [&](int64_t a, int64_t b) { return mMaxVL->Before(a, b); });
				std::vector<uint32_t> ranks((size_t)nv);
				for (int64_t r = 0; r < nv; r++) ranks[(size_t)vs[(size_t)r]] = (uint32_t)r;
				const ::dmt::Cubical3D G = ::dmt::Cubical3D::Make(X, Y, Z);
				std::vector<uint8_t> bytes((size_t)G.NumCells());
				DMFOptions opt;
				opt.second_pass = second_pass;
				opt.num_threads = num_threads;
				ComputeDMFGradient3D(ranks.data(), X, Y, Z, opt, bytes.data(), nullptr, nullptr, &mStats);
				for (int64_t c = 0; c < G.NumCells(); c++) {
					int64_t q[3];
					G.Coords(c, q);
					const INDEX_TYPE id = mMesh->coords2Cellid(Vec3l(q[0], q[1], q[2]));
					const uint8_t code = ::dmt::GradPairCode(bytes[(size_t)c]);
					if (code == ::dmt::kPairCritical) {
						mGrad->setCritical(id, true);
					}
					else {
						int64_t p[3];
						G.Coords(G.Neighbor(c, code), p);
						mGrad->setPair(id, mMesh->coords2Cellid(Vec3l(p[0], p[1], p[2])));
					}
					mGrad->setAssigned(id, true);
				}
			}
		};

	} // namespace dmt
} // namespace GInt

#endif
