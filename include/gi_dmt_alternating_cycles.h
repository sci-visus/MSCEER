#ifndef GI_DMT_ALTERNATING_CYCLES_H
#define GI_DMT_ALTERNATING_CYCLES_H

// Changing a discrete gradient without changing its critical cells.
//
// Two matchings on the same Hasse diagram leave the same cells unmatched
// exactly when their symmetric difference is a disjoint union of alternating
// cycles, so every critical-set-preserving change is a set of cycle flips.
// A cycle
//     x_1, m(x_1), x_2, m(x_2), ..., x_k, m(x_k), back to x_1
// has matched edges (x_i, m(x_i)) and unmatched Hasse edges (m(x_i), x_{i+1});
// flipping it pairs x_{i+1} with m(x_i). All 2k cells must be distinct.
// A cycle confined to one dimension stratum would be a closed V-path of the
// gradient it came from, so every cycle found here crosses strata. The
// shortest (k = 2) are the 4-cycle moves:
//     b paired down to v, a paired up to q, b a facet of q, v a facet of a
//     -> b paired up to q, a paired down to v.
// A flip can create a closed V-path; CreatesCycleAt (gi_dmt_acyclic.h),
// seeded with the re-matched cells, decides that soundly and locally.

#include "gi_dmt_acyclic.h"
#include "gi_dmt_gradient_view.h"

#include <algorithm>
#include <cstdint>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace GInt {
	namespace dmt {

		struct Move { int64_t b, v, a, q; };

		// 4-cycle moves with a ranging over `cells` (all cells if empty and
		// `all` is set), in that order.
		inline void MovesFrom(const GradientView3D& V, const std::vector<int64_t>& cells, std::vector<Move>& out) {
			for (int64_t a : cells) {
				if (!V.PairedUp(a)) continue;
				const int64_t q = V.Partner(a);
				int64_t fa[6], fq[6];
				const int na = V.G.Facets(a, fa);
				const int nq = V.G.Facets(q, fq);
				for (int i = 0; i < nq; i++) {
					const int64_t b = fq[i];
					if (b == a || !V.PairedDown(b)) continue;
					const int64_t v = V.Partner(b);
					for (int j = 0; j < na; j++)
						if (fa[j] == v) { out.push_back(Move{ b, v, a, q }); break; }
				}
			}
		}

		inline std::vector<Move> Moves(const GradientView3D& V) {
			std::vector<int64_t> all((size_t)V.G.NumCells());
			for (int64_t c = 0; c < V.G.NumCells(); c++) all[(size_t)c] = c;
			std::vector<Move> out;
			MovesFrom(V, all, out);
			return out;
		}

		inline std::vector<Move> LocalMoves(const GradientView3D& V, const std::vector<int64_t>& zone) {
			std::vector<Move> out;
			MovesFrom(V, zone, out);
			return out;
		}

		// Moves whose four cells are all >= margin from the domain boundary.
		inline std::vector<Move> InteriorMoves(const GradientView3D& V, int64_t margin, const std::vector<int64_t>* zone) {
			std::vector<Move> mv = zone ? LocalMoves(V, *zone) : Moves(V);
			std::vector<Move> out;
			for (const Move& m : mv)
				if (V.G.Bdist(m.b) >= margin && V.G.Bdist(m.v) >= margin && V.G.Bdist(m.a) >= margin && V.G.Bdist(m.q) >= margin)
					out.push_back(m);
			return out;
		}

		inline void ApplyMove(GradientView3D& V, const Move& m) {
			V.SetPair(m.b, m.q);
			V.SetPair(m.a, m.v);
		}
		inline void UndoMove(GradientView3D& V, const Move& m) {
			V.SetPair(m.b, m.v);
			V.SetPair(m.a, m.q);
		}
		inline bool CreatesCycle(const GradientView3D& V, const Move& m) {
			const int64_t seeds[2] = { m.b, m.v };
			return CreatesCycleAt(V, seeds, 2);
		}

		// ---- alternating cycles of any length ----

		struct CycleToken {
			std::vector<int64_t> xs, ms;
		};

		// Hasse neighbours: facets then cofacets.
		inline int HasseNbrs(const ::dmt::Cubical3D& G, int64_t c, int64_t out[12]) {
			int64_t f[6], g[6];
			const int nf = G.Facets(c, f), ng = G.Cofacets(c, g);
			int k = 0;
			for (int i = 0; i < nf; i++) out[k++] = f[i];
			for (int i = 0; i < ng; i++) out[k++] = g[i];
			return k;
		}

		// A zone: the cells a search may touch, as a list (for sampling) and a set.
		struct Zone {
			std::vector<int64_t> list;
			std::unordered_set<int64_t> set;
			void Assign(const std::unordered_set<int64_t>& s) {
				set = s;
				list.assign(s.begin(), s.end());
				std::sort(list.begin(), list.end());
			}
			bool Has(int64_t c) const { return set.count(c) != 0; }
		};

		// One random alternating cycle with k in [kmin, kmax] inside the zone,
		// or an empty vector (rotate.find_cycle, with this generator).
		template <class Rng>
		std::vector<int64_t> FindCycle(const GradientView3D& V, const Zone& zone, Rng& rng, int kmin, int kmax, int tries = 60) {
			std::vector<int64_t> xs;
			if (zone.list.empty()) return xs;
			for (int t = 0; t < tries; t++) {
				const int64_t x0 = zone.list[(size_t)(rng() % zone.list.size())];
				const int64_t m0 = V.Partner(x0);
				if (m0 < 0) continue;
				xs.assign(1, x0);
				std::unordered_set<int64_t> used = { x0, m0 };
				while ((int)xs.size() <= kmax) {
					const int64_t m = V.Partner(xs.back());
					int64_t nb[12];
					const int n = HasseNbrs(V.G, m, nb);
					for (int i = n - 1; i > 0; i--) {             // Fisher-Yates
						const int j = (int)(rng() % (uint64_t)(i + 1));
						std::swap(nb[i], nb[j]);
					}
					if ((int)xs.size() >= kmin)
						for (int i = 0; i < n; i++)
							if (nb[i] == x0) return xs;
					int64_t step = -1;
					for (int i = 0; i < n; i++) {
						const int64_t y = nb[i];
						if (used.count(y)) continue;
						const int64_t my = V.Partner(y);
						if (my < 0 || used.count(my) || !zone.Has(y) || !zone.Has(my)) continue;
						step = y;
						break;
					}
					if (step < 0) break;
					xs.push_back(step);
					used.insert(step);
					used.insert(V.Partner(step));
				}
			}
			xs.clear();
			return xs;
		}

		// A simple alternating cycle whose flip is a valid matching?
		inline bool CycleOk(const GradientView3D& V, const std::vector<int64_t>& xs) {
			const size_t k = xs.size();
			std::vector<int64_t> ms(k);
			for (size_t i = 0; i < k; i++) {
				ms[i] = V.Partner(xs[i]);
				if (ms[i] < 0) return false;
			}
			std::unordered_set<int64_t> cells(xs.begin(), xs.end());
			cells.insert(ms.begin(), ms.end());
			if (cells.size() != 2 * k) return false;
			for (size_t i = 0; i < k; i++) {
				const int64_t a = xs[(i + 1) % k], b = ms[i];
				const int da = V.G.Dim(a), db = V.G.Dim(b);
				if (da - db != 1 && db - da != 1) return false;
				const int64_t lo = da < db ? a : b, hi = da < db ? b : a;
				int64_t f[6];
				const int nf = V.G.Facets(hi, f);
				bool found = false;
				for (int j = 0; j < nf; j++)
					if (f[j] == lo) { found = true; break; }
				if (!found) return false;
			}
			return true;
		}

		inline CycleToken ApplyCycle(GradientView3D& V, const std::vector<int64_t>& xs) {
			CycleToken tok;
			tok.xs = xs;
			const size_t k = xs.size();
			tok.ms.resize(k);
			for (size_t i = 0; i < k; i++) tok.ms[i] = V.Partner(xs[i]);
			for (size_t i = 0; i < k; i++) V.SetPair(xs[(i + 1) % k], tok.ms[i]);
			return tok;
		}

		inline void UndoCycle(GradientView3D& V, const CycleToken& tok) {
			for (size_t i = 0; i < tok.xs.size(); i++) V.SetPair(tok.xs[i], tok.ms[i]);
		}

		inline bool CycleCreatesClosedVPath(const GradientView3D& V, const CycleToken& tok) {
			std::vector<int64_t> cells(tok.xs);
			cells.insert(cells.end(), tok.ms.begin(), tok.ms.end());
			return CreatesCycleAt(V, cells.data(), (int)cells.size());
		}

		// ---- neighbourhoods of a connection ----

		// Cells on some gradient path from a facet of tau to sigma (any p),
		// plus tau and sigma (rotate.pathcells).
		inline std::unordered_set<int64_t> PathCellsPair(const GradientView3D& V, int64_t tau, int64_t sigma) {
			const int d = V.G.Dim(sigma);
			std::unordered_map<int64_t, char> reach;   // 0 grey/false, 1 true
			auto children = [&](int64_t e, int64_t out[6]) -> int {
				const int64_t q = V.Partner(e);
				if (q < 0 || V.G.Dim(q) != d + 1) return 0;
				int64_t f[6];
				const int nf = V.G.Facets(q, f);
				int k = 0;
				for (int i = 0; i < nf; i++)
					if (f[i] != e) out[k++] = f[i];
				return k;
			};
			// post-order evaluation of reach(e)
			auto eval = [&](int64_t start) {
				if (reach.count(start)) return;
				struct Fr { int64_t e; bool done; };
				std::vector<Fr> st = { { start, false } };
				while (!st.empty()) {
					const Fr fr = st.back();
					st.pop_back();
					if (fr.done) {
						int64_t k[6];
						const int nk = children(fr.e, k);
						char ok = 0;
						for (int i = 0; i < nk; i++)
							if (reach[k[i]] == 1) ok = 1;
						reach[fr.e] = ok;
						continue;
					}
					if (reach.count(fr.e)) continue;
					if (fr.e == sigma) { reach[fr.e] = 1; continue; }
					reach[fr.e] = 0;                                  // grey reads as false
					st.push_back(Fr{ fr.e, true });
					int64_t k[6];
					const int nk = children(fr.e, k);
					for (int i = 0; i < nk; i++)
						if (!reach.count(k[i])) st.push_back(Fr{ k[i], false });
				}
			};
			std::unordered_set<int64_t> cells, seen;
			int64_t f[6];
			const int nf = V.G.Facets(tau, f);
			std::vector<int64_t> stack(f, f + nf);
			while (!stack.empty()) {
				const int64_t e = stack.back();
				stack.pop_back();
				if (!seen.insert(e).second) continue;
				eval(e);
				if (reach[e] != 1) continue;
				cells.insert(e);
				const int64_t q = V.Partner(e);
				if (q >= 0 && V.G.Dim(q) == d + 1) {
					cells.insert(q);
					int64_t g[6];
					const int ng = V.G.Facets(q, g);
					for (int i = 0; i < ng; i++)
						if (g[i] != e) stack.push_back(g[i]);
				}
			}
			cells.insert(tau);
			cells.insert(sigma);
			return cells;
		}

		inline std::unordered_set<int64_t> Neighbourhood(const ::dmt::Cubical3D& G, const std::unordered_set<int64_t>& seed, int hops) {
			std::unordered_set<int64_t> cur(seed);
			for (int h = 0; h < hops; h++) {
				std::unordered_set<int64_t> nxt(cur);
				for (int64_t c : cur) {
					int64_t nb[12];
					const int n = HasseNbrs(G, c, nb);
					for (int i = 0; i < n; i++) nxt.insert(nb[i]);
				}
				cur.swap(nxt);
			}
			return cur;
		}

	} // namespace dmt
} // namespace GInt

#endif
