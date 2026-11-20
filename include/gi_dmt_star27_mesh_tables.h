#ifndef GI_DMT_STAR27_MESH_TABLES_H
#define GI_DMT_STAR27_MESH_TABLES_H

// Star27Tables generated from a real mesh's iterators (FacetsIterator,
// CofacetsIterator, CellVerticesIterator), so the device-safe core and the
// mesh share one definition of the star's topology. Mirrors
// BuildDgrad2DTablesFromMesh in gpu_dgrad/dgrad2d_tables.h.

#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "dmt_core/dmt_star27_tables.h"

namespace GInt {
	namespace dmt {

		// centre must be the cell id of a vertex with all 26 neighbours inside
		// the mesh. Iteration order of the mesh is kept in fac/cof/vert; compare
		// with TablesEquivalent, which ignores order.
		template <class MeshType>
		bool BuildStar27TablesFromMesh(MeshType* mesh, INDEX_TYPE centre, ::dmt::Star27Tables& T) {
			Vec3l cc;
			mesh->cellid2Coords(centre, cc);
			auto slot_of = [&](INDEX_TYPE id, int& s) {
				Vec3l q;
				mesh->cellid2Coords(id, q);
				const long long d[3] = { q[0] - cc[0], q[1] - cc[1], q[2] - cc[2] };
				for (int i = 0; i < 3; i++)
					if (d[i] < -1 || d[i] > 1) return false;
				s = ::dmt::StarSlot((int)d[0], (int)d[1], (int)d[2]);
				return true;
			};
			for (int dz = -1; dz <= 1; dz++)
				for (int dy = -1; dy <= 1; dy++)
					for (int dx = -1; dx <= 1; dx++) {
						const int s = ::dmt::StarSlot(dx, dy, dz);
						const INDEX_TYPE id = mesh->coords2Cellid(Vec3l(cc[0] + dx, cc[1] + dy, cc[2] + dz));
						T.off[s][0] = (int8_t)dx;
						T.off[s][1] = (int8_t)dy;
						T.off[s][2] = (int8_t)dz;
						T.dim[s] = (uint8_t)mesh->dimension(id);

						T.nfac[s] = 0;
						typename MeshType::FacetsIterator fit(mesh);
						for (fit.begin(id); fit.valid(); fit.advance()) {
							int t;
							if (!slot_of(fit.value(), t)) continue;
							if (T.nfac[s] >= 3) return false;
							T.fac[s][T.nfac[s]++] = (uint8_t)t;
						}
						T.ncof[s] = 0;
						typename MeshType::CofacetsIterator cit(mesh);
						for (cit.begin(id); cit.valid(); cit.advance()) {
							int t;
							if (!slot_of(cit.value(), t)) continue;
							if (T.ncof[s] >= 6) return false;
							T.cof[s][T.ncof[s]++] = (uint8_t)t;
						}
						T.nvert[s] = 0;
						typename MeshType::CellVerticesIterator vit(mesh);
						for (vit.begin(id); vit.valid(); vit.advance()) {
							Vec3l q;
							mesh->cellid2Coords(vit.value(), q);
							const long long d[3] = { (q[0] - cc[0]) / 2, (q[1] - cc[1]) / 2, (q[2] - cc[2]) / 2 };
							if (T.nvert[s] >= 8) return false;
							T.vert[s][T.nvert[s]++] = (uint8_t)::dmt::StarSlot((int)d[0], (int)d[1], (int)d[2]);
						}
					}
			::dmt::FinishStar27Masks(T);
			return true;
		}

		// Same topology, ignoring the order of the fac/cof/vert lists.
		inline bool TablesEquivalent(const ::dmt::Star27Tables& A, const ::dmt::Star27Tables& B) {
			for (int s = 0; s < 27; s++) {
				for (int i = 0; i < 3; i++)
					if (A.off[s][i] != B.off[s][i]) return false;
				if (A.dim[s] != B.dim[s] || A.nfac[s] != B.nfac[s] || A.ncof[s] != B.ncof[s] ||
					A.nvert[s] != B.nvert[s] || A.facmask[s] != B.facmask[s] ||
					A.cofmask[s] != B.cofmask[s] || A.vertmask[s] != B.vertmask[s])
					return false;
			}
			for (int d = 0; d < 4; d++)
				if (A.dimmask[d] != B.dimmask[d]) return false;
			return true;
		}

	} // namespace dmt
} // namespace GInt

#endif
