#ifndef GI_DMT_GRADIENT_VIEW_H
#define GI_DMT_GRADIENT_VIEW_H

// A non-owning view of a discrete gradient on the cubical complex of a 3D
// vertex grid: one GradBitfield byte per cell (dmt_core/dmt_hd.h), cell ids as
// in dmt::Cubical3D. Every analysis routine in gi_dmt_*.h takes one of these,
// so it runs the same on a gradient computed here and on one computed by any
// MSCEER pipeline (MakeGradientView bridges a DiscreteGradientLabeling).

#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "dmt_core/dmt_cubical3d.h"
#include "dmt_core/dmt_hd.h"

#include <cstdint>
#include <vector>

namespace GInt {
	namespace dmt {

		struct GradientView3D {
			::dmt::Cubical3D G;
			uint8_t* g = nullptr;

			GradientView3D() {}
			GradientView3D(const ::dmt::Cubical3D& grid, uint8_t* bytes) : G(grid), g(bytes) {}

			bool IsCritical(int64_t id) const { return ::dmt::GradCritical(g[id]); }
			// Partner cell, or -1 for a critical cell.
			int64_t Partner(int64_t id) const {
				const uint8_t code = ::dmt::GradPairCode(g[id]);
				if (code == ::dmt::kPairCritical) return -1;
				return G.Neighbor(id, code);
			}
			// Paired with a cofacet (the cell is the tail of an arrow).
			bool PairedUp(int64_t id) const {
				const int64_t p = Partner(id);
				return p >= 0 && G.Dim(p) == G.Dim(id) + 1;
			}
			bool PairedDown(int64_t id) const {
				const int64_t p = Partner(id);
				return p >= 0 && G.Dim(p) + 1 == G.Dim(id);
			}
			// Pair two adjacent cells (writes both bytes). Returns false and
			// writes nothing if they are not Hasse-adjacent.
			bool SetPair(int64_t a, int64_t b) {
				const int da = G.DirBetween(a, b), db = G.DirBetween(b, a);
				if (da > 5 || db > 5) return false;
				g[a] = ::dmt::MakeGradByte((uint8_t)da);
				g[b] = ::dmt::MakeGradByte((uint8_t)db);
				return true;
			}
			// One side only, for callers that maintain consistency themselves.
			bool SetPartnerOneSided(int64_t a, int64_t b) {
				const int da = G.DirBetween(a, b);
				if (da > 5) return false;
				g[a] = ::dmt::MakeGradByte((uint8_t)da);
				return true;
			}
			void SetCritical(int64_t a) { g[a] = ::dmt::MakeGradByte(::dmt::kPairCritical); }

			// Critical cells of dimension d (or all with d < 0), in id order.
			std::vector<int64_t> Critical(int d = -1) const {
				std::vector<int64_t> out;
				const int64_t n = G.NumCells();
				for (int64_t c = 0; c < n; c++)
					if (IsCritical(c) && (d < 0 || G.Dim(c) == d)) out.push_back(c);
				return out;
			}
		};

		// View over the bytes of a DiscreteGradientLabeling on a non-periodic
		// TopologicalRegularGrid3D. The view aliases the labeling's storage.
		template <class MeshType, class GradType>
		GradientView3D MakeGradientView(MeshType* mesh, GradType* grad) {
			const Vec3l xyz = mesh->XYZ();   // cells per axis
			::dmt::Cubical3D G = ::dmt::Cubical3D::Make((xyz[0] + 1) / 2, (xyz[1] + 1) / 2, (xyz[2] + 1) / 2);
			return GradientView3D(G, reinterpret_cast<uint8_t*>(grad->m_dgrad->LabelArray()));
		}

	} // namespace dmt
} // namespace GInt

#endif
