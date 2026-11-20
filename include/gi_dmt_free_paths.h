#ifndef GI_DMT_FREE_PATHS_H
#define GI_DMT_FREE_PATHS_H

// Descending (2,1) gradient paths between arbitrary endpoints, counted per
// facet of the source: two paths leaving a source through different facets
// split at the source itself, two leaving through the same facet and ending
// at the same sink needed a split and a merge on the way.
//
// All of these are views on VPathCounter (unsigned): the paths below a quad q
// entered through its paired edge e are exactly the paths from e. Cells are
// visited in tuple order and facets in facet order, and "best" keeps the first
// maximum, so the reported witnesses are deterministic.
//
//   FacetSinkPaths   {facet of src: {sink: paths}}, sinks given, src any 2-cell,
//                    optional blocked cells the walk may not enter;
//   FacetPaths       the same from a critical 2-cell to the critical 1-cells;
//   MultAt           the largest per-facet count at a critical 2-cell;
//   MultAnywhere     the largest per-facet count from ANY 2-cell (critical or
//                    not) to a critical 1-cell, with its witness;
//   Multiplicities   every critical 2-cell's per-facet counts >= minmult,
//                    best first, plus a histogram of all counts;
//   PathCells        the cells on the paths s2 -> facet -> ... -> s1;
//   FaceToFace       interior quads under the top face down to the bottom face
//                    (or, `interior`, to the 1-cells one layer in, never
//                    entering a boundary cell);
//   CensusBest       busiest source quad under the top face against the edges
//                    at the bottom-middle vertex: (paths, largest single sink).

#include "gi_dmt_gradient_view.h"
#include "gi_dmt_vpath_counter.h"

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace GInt {
	namespace dmt {

		struct FacetRow {
			int64_t facet;
			std::vector<SinkEntry> row;   // unsigned counts in .u
		};

		// {facet: row} from src; a blocked facet, or src's own paired facet, is
		// omitted; a facet that leads nowhere has an empty row.
		inline std::vector<FacetRow> FacetSinkPaths(VPathCounter& C, int64_t src) {
			std::vector<FacetRow> out;
			const GradientView3D& V = C.View();
			int64_t f[6];
			const int nf = V.G.Facets(src, f);
			const int64_t sp = V.Partner(src);
			for (int i = 0; i < nf; i++) {
				if (f[i] == sp) continue;
				if (C.IsBlocked(f[i])) continue;
				FacetRow r;
				r.facet = f[i];
				r.row = C.PathsFrom(f[i]);
				out.push_back(std::move(r));
			}
			return out;
		}

		inline std::vector<int64_t> SortedTuple(const ::dmt::Cubical3D& G, std::vector<int64_t> v) {
			std::sort(v.begin(), v.end(), [&](int64_t a, int64_t b) { return G.TupleLess(a, b); });
			return v;
		}

		inline CounterConfig CriticalSinks(const GradientView3D& V, int p) {
			CounterConfig cfg;
			cfg.sinks = V.Critical(p);
			cfg.signed_counts = false;
			return cfg;
		}

		inline std::vector<FacetRow> FacetPaths(const GradientView3D& V, int64_t s2) {
			VPathCounter C(V, CriticalSinks(V, 1));
			return FacetSinkPaths(C, s2);
		}

		// -1 if s2 is not critical.
		inline int64_t MultAt(const GradientView3D& V, int64_t s2, ::dmt::U128& best) {
			best = ::dmt::U128();
			if (!V.IsCritical(s2)) return -1;
			for (const FacetRow& fr : FacetPaths(V, s2))
				for (const SinkEntry& e : fr.row)
					if (::dmt::Less(best, e.u)) best = e.u;
			return 0;
		}

		struct Witness { int64_t q = -1, facet = -1, sink = -1; };

		inline ::dmt::U128 MultAnywhere(const GradientView3D& V, bool interior_only, Witness& where) {
			VPathCounter C(V, CriticalSinks(V, 1));
			::dmt::U128 best;
			where = Witness();
			const ::dmt::Cubical3D& G = V.G;
			for (int64_t x = 0; x < G.c[0]; x++)
				for (int64_t y = 0; y < G.c[1]; y++)
					for (int64_t z = 0; z < G.c[2]; z++) {
						const int64_t q = G.Id(x, y, z);
						if (G.Dim(q) != 2 || (interior_only && G.OnBoundary(q))) continue;
						std::vector<FacetRow> fp;
						try { fp = FacetSinkPaths(C, q); }
						catch (const ClosedVPath&) { continue; }
						for (const FacetRow& fr : fp)
							for (const SinkEntry& e : fr.row)
								if (::dmt::Less(best, e.u)) { best = e.u; where.q = q; where.facet = fr.facet; where.sink = e.sink; }
					}
			return best;
		}

		struct Multiplicity { ::dmt::U128 n; int64_t s2, facet, s1; };

		// Sorted by (-n, s2, facet, s1) in tuple order. hist counts every
		// (s2, facet, s1) triple reached through a paired quad, n >= 1.
		inline std::vector<Multiplicity> Multiplicities(const GradientView3D& V, bool interior_only, uint64_t minmult,
			std::vector<std::pair<::dmt::U128, int64_t>>* hist) {
			VPathCounter C(V, CriticalSinks(V, 1));
			std::vector<Multiplicity> out;
			std::unordered_map<uint64_t, int64_t> h64;
			std::vector<std::pair<::dmt::U128, int64_t>> hbig;
			const ::dmt::Cubical3D& G = V.G;
			for (int64_t s2 : V.Critical(2)) {
				if (interior_only && G.OnBoundary(s2)) continue;
				int64_t f[6];
				const int nf = G.Facets(s2, f);
				for (int i = 0; i < nf; i++) {
					const int64_t e = f[i];
					if (C.IsSink(e) || !V.PairedUp(e)) continue;
					const std::vector<SinkEntry>* row = nullptr;
					try { row = &C.PathsFrom(e); }
					catch (const ClosedVPath&) { continue; }
					for (const SinkEntry& se : *row) {
						if (hist) {
							if (se.u.hi == 0) h64[se.u.lo]++;
							else {
								bool found = false;
								for (auto& kv : hbig)
									if (kv.first.lo == se.u.lo && kv.first.hi == se.u.hi) { kv.second++; found = true; break; }
								if (!found) hbig.push_back({ se.u, 1 });
							}
						}
						if (se.u.hi != 0 || se.u.lo >= minmult) out.push_back(Multiplicity{ se.u, s2, e, se.sink });
					}
				}
			}
			std::sort(out.begin(), out.end(), [&](const Multiplicity& a, const Multiplicity& b) {
				if (a.n.hi != b.n.hi || a.n.lo != b.n.lo) return ::dmt::Less(b.n, a.n);
				if (a.s2 != b.s2) return G.TupleLess(a.s2, b.s2);
				if (a.facet != b.facet) return G.TupleLess(a.facet, b.facet);
				return G.TupleLess(a.s1, b.s1);
				});
			if (hist) {
				hist->clear();
				for (auto& kv : h64) hist->push_back({ ::dmt::U128(kv.first), kv.second });
				for (auto& kv : hbig) hist->push_back(kv);
			}
			return out;
		}

		// Cells on the descending paths from s2 through `facet` to s1.
		inline std::vector<int64_t> PathCells(const GradientView3D& V, int64_t s2, int64_t facet, int64_t s1) {
			std::unordered_set<int64_t> keep = { s2, facet, s1 };
			std::unordered_map<int64_t, char> state;   // 1 grey, 2 reaches, 3 does not
			const int64_t q0 = V.Partner(facet);
			if (q0 >= 0 && V.G.Dim(q0) == 2) {
				struct Fr { int64_t q; bool done; };
				std::vector<Fr> st = { { q0, false } };
				while (!st.empty()) {
					const Fr fr = st.back();
					st.pop_back();
					const int64_t q = fr.q;
					const int64_t in = V.Partner(q);
					int64_t f[6];
					const int nf = V.G.Facets(q, f);
					if (fr.done) {
						bool hit = false;
						for (int i = 0; i < nf; i++) {
							const int64_t e = f[i];
							if (e == in) continue;
							if (e == s1) { hit = true; keep.insert(e); continue; }
							const int64_t p = V.Partner(e);
							if (p >= 0 && V.G.Dim(p) == 2) {
								auto it = state.find(p);
								if (it != state.end() && it->second == 2) { hit = true; keep.insert(e); }
							}
						}
						state[q] = hit ? 2 : 3;
						if (hit) { keep.insert(q); keep.insert(in); }
						continue;
					}
					if (state.count(q)) continue;   // grey reads as "does not reach", as the reference
					state[q] = 1;
					st.push_back(Fr{ q, true });
					for (int i = 0; i < nf; i++) {
						const int64_t e = f[i];
						if (e == in || e == s1) continue;
						const int64_t p = V.Partner(e);
						if (p >= 0 && V.G.Dim(p) == 2 && !state.count(p)) st.push_back(Fr{ p, false });
					}
				}
			}
			return std::vector<int64_t>(keep.begin(), keep.end());
		}

		// ---- endpoint selectors of the lid probes ----

		// Interior quads with a vertex in the top z layer, tuple order.
		inline std::vector<int64_t> TopSources(const ::dmt::Cubical3D& G) {
			std::vector<int64_t> out;
			for (int64_t c = 0; c < G.NumCells(); c++) {
				if (G.Dim(c) != 2 || G.OnBoundary(c)) continue;
				int64_t q[3];
				G.Coords(c, q);
				const int64_t zmax = (q[2] + 1) / 2;   // highest vertex z of the cell
				if (zmax == G.n[2] - 1) out.push_back(c);
			}
			return SortedTuple(G, out);
		}

		// Edges at the middle vertex of the bottom layer.
		inline std::vector<int64_t> BottomMiddleSinks(const ::dmt::Cubical3D& G) {
			const int64_t mx = 2 * ((G.n[0] - 1) / 2), my = 2 * ((G.n[1] - 1) / 2);
			std::vector<int64_t> out;
			for (int64_t c = 0; c < G.NumCells(); c++) {
				if (G.Dim(c) != 1) continue;
				int64_t q[3];
				G.Coords(c, q);
				if (q[2] > 1) continue;
				bool has = true;
				for (int i = 0; i < 3 && has; i++) {
					const int64_t want = i == 0 ? mx : (i == 1 ? my : 0);
					if (q[i] & 1) has = (q[i] - 1 == want || q[i] + 1 == want);
					else has = q[i] == want;
				}
				if (has) out.push_back(c);
			}
			return SortedTuple(G, out);
		}

		// 1-cells in a face of the domain along `axis`.
		inline std::vector<int64_t> FaceSinks(const ::dmt::Cubical3D& G, int axis, bool at_top) {
			const int64_t end = at_top ? G.Top(axis) : 0;
			std::vector<int64_t> out;
			for (int64_t c = 0; c < G.NumCells(); c++) {
				if (G.Dim(c) != 1) continue;
				int64_t q[3];
				G.Coords(c, q);
				if (q[axis] == end) out.push_back(c);
			}
			return SortedTuple(G, out);
		}

		// Interior 1-cells in the vertex plane `layer` steps in from the low face.
		inline std::vector<int64_t> InteriorSinks(const ::dmt::Cubical3D& G, int layer, int axis) {
			std::vector<int64_t> out;
			for (int64_t c = 0; c < G.NumCells(); c++) {
				if (G.Dim(c) != 1 || G.OnBoundary(c)) continue;
				int64_t q[3];
				G.Coords(c, q);
				if (q[axis] == 2 * layer) out.push_back(c);
			}
			return SortedTuple(G, out);
		}

		struct FaceToFaceResult { ::dmt::U128 best_mult, best_total; int64_t n_sources = 0; };

		inline FaceToFaceResult FaceToFace(const GradientView3D& V, bool interior, int layer) {
			const ::dmt::Cubical3D& G = V.G;
			CounterConfig cfg;
			cfg.signed_counts = false;
			cfg.sinks = interior ? InteriorSinks(G, layer, 2) : FaceSinks(G, 2, false);
			if (interior)
				for (int64_t c = 0; c < G.NumCells(); c++)
					if (G.OnBoundary(c)) cfg.blocked.push_back(c);
			VPathCounter C(V, cfg);
			FaceToFaceResult R;
			const std::vector<int64_t> srcs = TopSources(G);
			R.n_sources = (int64_t)srcs.size();
			for (int64_t s : srcs) {
				std::vector<FacetRow> fp;
				try { fp = FacetSinkPaths(C, s); }
				catch (const ClosedVPath&) { continue; }
				::dmt::U128 tot;
				for (const FacetRow& fr : fp)
					for (const SinkEntry& e : fr.row) {
						tot = ::dmt::Add(tot, e.u);
						if (::dmt::Less(R.best_mult, e.u)) R.best_mult = e.u;
					}
				if (::dmt::Less(R.best_total, tot)) R.best_total = tot;
			}
			return R;
		}

		// freepaths.best: the source with the most paths (first maximum) to the
		// edges at the bottom-middle vertex, and its largest single sink.
		inline std::pair<::dmt::U128, ::dmt::U128> CensusBest(const GradientView3D& V) {
			const ::dmt::Cubical3D& G = V.G;
			CounterConfig cfg;
			cfg.signed_counts = false;
			cfg.sinks = BottomMiddleSinks(G);
			VPathCounter C(V, cfg);
			::dmt::U128 best_total, best_sink;
			bool any = false;
			for (int64_t s : TopSources(G)) {
				const std::vector<FacetRow> fp = FacetSinkPaths(C, s);
				std::unordered_map<int64_t, ::dmt::U128> bys;
				::dmt::U128 tot;
				for (const FacetRow& fr : fp)
					for (const SinkEntry& e : fr.row) {
						tot = ::dmt::Add(tot, e.u);
						bys[e.sink] = ::dmt::Add(bys[e.sink], e.u);
					}
				if (!any || ::dmt::Less(best_total, tot)) {
					::dmt::U128 m;
					for (auto& kv : bys)
						if (::dmt::Less(m, kv.second)) m = kv.second;
					best_total = tot;
					best_sink = m;
					any = true;
				}
			}
			return { best_total, best_sink };
		}

	} // namespace dmt
} // namespace GInt

#endif
