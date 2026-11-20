#ifndef GI_DMT_ACYCLIC_H
#define GI_DMT_ACYCLIC_H

// Validity checks of a discrete gradient.
//
// IsAcyclic: no closed V-path in any dimension. Arrows run tail -> g for each
//   tail paired up with head and each facet g != tail of head that is itself
//   paired up; a Kahn sort over those arrows must consume every tail.
// CreatesCycleAt: after a local change that re-matched some cells, does any
//   new cycle exist? Every new arrow leaves a cell that the change re-matched
//   upwards, so a depth-first search from each of those, following arrows,
//   must not return to its start. Sound for any change; local in cost.
// CheckConsistency: every non-critical cell's partner points back and differs
//   in dimension by one.

#include "gi_dmt_gradient_view.h"

#include <cstdint>
#include <unordered_set>
#include <vector>

namespace GInt {
	namespace dmt {

		// Tails reachable in one V-path step from `tail` (at most 6).
		inline int VPathSuccessors(const GradientView3D& V, int64_t tail, int64_t out[6]) {
			const int64_t head = V.Partner(tail);
			int64_t f[6];
			const int nf = V.G.Facets(head, f);
			int k = 0;
			for (int i = 0; i < nf; i++)
				if (f[i] != tail && V.PairedUp(f[i])) out[k++] = f[i];
			return k;
		}

		inline bool IsAcyclic(const GradientView3D& V) {
			const int64_t n = V.G.NumCells();
			std::vector<int32_t> indeg((size_t)n, 0);
			std::vector<uint8_t> tail((size_t)n, 0);
			int64_t ntails = 0;
			for (int64_t c = 0; c < n; c++) {
				if (!V.PairedUp(c)) continue;
				tail[(size_t)c] = 1;
				ntails++;
			}
			for (int64_t c = 0; c < n; c++) {
				if (!tail[(size_t)c]) continue;
				int64_t s[6];
				const int k = VPathSuccessors(V, c, s);
				for (int i = 0; i < k; i++) indeg[(size_t)s[i]]++;
			}
			std::vector<int64_t> q;
			for (int64_t c = 0; c < n; c++)
				if (tail[(size_t)c] && indeg[(size_t)c] == 0) q.push_back(c);
			int64_t seen = 0;
			while (!q.empty()) {
				const int64_t u = q.back();
				q.pop_back();
				seen++;
				int64_t s[6];
				const int k = VPathSuccessors(V, u, s);
				for (int i = 0; i < k; i++)
					if (--indeg[(size_t)s[i]] == 0) q.push_back(s[i]);
			}
			return seen == ntails;
		}

		inline bool CreatesCycleAt(const GradientView3D& V, const int64_t* cells, int ncells) {
			std::vector<int64_t> stack;
			std::unordered_set<int64_t> seen;
			for (int i = 0; i < ncells; i++) {
				const int64_t start = cells[i];
				if (!V.PairedUp(start)) continue;
				stack.clear();
				seen.clear();
				int64_t s[6];
				const int k = VPathSuccessors(V, start, s);
				for (int j = 0; j < k; j++) stack.push_back(s[j]);
				while (!stack.empty()) {
					const int64_t x = stack.back();
					stack.pop_back();
					if (x == start) return true;
					if (!seen.insert(x).second) continue;
					const int kk = VPathSuccessors(V, x, s);
					for (int j = 0; j < kk; j++) stack.push_back(s[j]);
				}
			}
			return false;
		}

		// Number of inconsistent cells (0 = a valid matching).
		inline int64_t CheckConsistency(const GradientView3D& V) {
			const int64_t n = V.G.NumCells();
			int64_t bad = 0;
			for (int64_t c = 0; c < n; c++) {
				if (!::dmt::GradAssigned(V.g[c])) { bad++; continue; }
				if (V.IsCritical(c)) continue;
				const int64_t p = V.Partner(c);
				int64_t q[3];
				V.G.Coords(p, q);
				if (!V.G.InDomain(q[0], q[1], q[2]) || V.IsCritical(p) || V.Partner(p) != c) { bad++; continue; }
				const int dd = V.G.Dim(p) - V.G.Dim(c);
				if (dd != 1 && dd != -1) bad++;
			}
			return bad;
		}

	} // namespace dmt
} // namespace GInt

#endif
