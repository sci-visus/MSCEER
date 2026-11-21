#ifndef GI_DMT_CYCLE_HUNTS_H
#define GI_DMT_CYCLE_HUNTS_H

// Greedy searches over critical-set-preserving flips (gi_dmt_alternating_cycles.h)
// that try to make gradient-path connections unique, and the measures they use.
//
//   Score          (|Gamma|, coefficient) of one connection tau -> sigma;
//   Strangled      critical (p+1, p) connections with more than one path
//                  (and, with require_unit, a coefficient of +-1);
//   GlobalCost     (how many such connections, total excess paths);
//   Hunt           4-cycle moves anywhere, descending on one pair's |Gamma|;
//   HuntLocal      4-cycle moves near the pair's path cells;
//   HuntCycles     alternating cycles of k in [kmin, kmax] near the path cells,
//                  optionally only cells >= margin from the domain boundary;
//   HuntAll        4-cycle moves anywhere, descending on GlobalCost;
//   MinimaRows /   which minima each 1-saddle reaches, and how many 1-saddles a
//   HandleSlides   change of gradient re-wired.
//
// A flip is kept when the objective does not get worse (path count nonzero
// for the pair searches) and undone otherwise; flips that create a closed
// V-path are undone before scoring. The random generator is std::mt19937_64,
// so results are reproducible for a seed but are not those of another
// generator.

#include "gi_dmt_alternating_cycles.h"
#include "gi_dmt_morse_complex.h"
#include "gi_dmt_vpath_counter.h"

#include <cstdint>
#include <map>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace GInt {
	namespace dmt {

		struct Connection { int64_t tau, sigma; ::dmt::U128 u; ::dmt::S128 coeff; };

		inline bool IsUnit(const ::dmt::S128& c) {
			return c.FitsI64() && (c.ToI64() == 1 || c.ToI64() == -1);
		}

		// Connections of the critical (p+1)-cells, tuple order of tau, then row order.
		inline std::vector<Connection> Strangled(const GradientView3D& V, bool require_unit, int p = 1) {
			const MorseMatrix M = MorseBoundary(V, p);
			std::vector<Connection> out;
			for (size_t i = 0; i < M.taus.size(); i++)
				for (int64_t k = M.rowptr[i]; k < M.rowptr[i + 1]; k++) {
					const MorseEntry& e = M.entries[(size_t)k];
					const bool multi = e.u.hi != 0 || e.u.lo > 1;
					if (multi && (!require_unit || IsUnit(e.coeff)))
						out.push_back(Connection{ M.taus[i], e.sigma, e.u, e.coeff });
				}
			return out;
		}

		struct Cost {
			int64_t count = 0;
			::dmt::U128 excess;
			bool operator<(const Cost& o) const { return count != o.count ? count < o.count : ::dmt::Less(excess, o.excess); }
			bool operator<=(const Cost& o) const { return !(o < *this); }
		};

		inline Cost GlobalCost(const GradientView3D& V, bool require_unit, int p = 1) {
			Cost c;
			for (const Connection& s : Strangled(V, require_unit, p)) {
				c.count++;
				::dmt::U128 m = s.u;
				// m - 1, m >= 2
				if (m.lo == 0) { m.hi--; m.lo = ~0ull; }
				else m.lo--;
				c.excess = ::dmt::Add(c.excess, m);
			}
			return c;
		}

		struct HuntOptions {
			int kmin = 2, kmax = 4, hops = 3;
			int64_t rounds = 3000, plateau = 600;
			int64_t margin = -1;          // < 0: no restriction
			uint64_t seed = 0;
			bool verify = false;          // re-check the critical set after every kept flip
		};

		struct HuntResult {
			PairScore best;
			int64_t nflip = 0, tried = 0, accepted = 0, moves_available = 0;
			std::map<int, int64_t> hist;  // kept flips by cycle length k
			Cost best_cost;               // HuntAll only
			bool critical_set_preserved = true;
		};

		namespace detail {

			inline std::vector<uint8_t> CriticalMask(const GradientView3D& V) {
				std::vector<uint8_t> m((size_t)V.G.NumCells());
				for (int64_t c = 0; c < V.G.NumCells(); c++) m[(size_t)c] = V.IsCritical(c) ? 1 : 0;
				return m;
			}

			inline void PairZone(const GradientView3D& V, int64_t tau, int64_t sigma, const HuntOptions& o, Zone& z) {
				std::unordered_set<int64_t> s = Neighbourhood(V.G, PathCellsPair(V, tau, sigma), o.hops);
				if (o.margin >= 0) {
					std::unordered_set<int64_t> t;
					for (int64_t c : s)
						if (V.G.Bdist(c) >= o.margin) t.insert(c);
					s.swap(t);
				}
				z.Assign(s);
			}

			inline bool Accept(const PairScore& nw, const PairScore& cur) {
				return !nw.paths.IsZero() && !::dmt::Less(cur.paths, nw.paths);
			}
			inline bool Better(const PairScore& nw, const PairScore& best) { return ::dmt::Less(nw.paths, best.paths); }
			inline bool IsOne(const PairScore& s) { return s.paths.hi == 0 && s.paths.lo == 1; }

		} // namespace detail

		// rotate.hunt: 4-cycle moves anywhere.
		inline HuntResult Hunt(GradientView3D& V, int64_t tau, int64_t sigma, const HuntOptions& o) {
			std::mt19937_64 rng(o.seed);
			PairScorer S(V.G);
			const std::vector<uint8_t> base = o.verify ? detail::CriticalMask(V) : std::vector<uint8_t>();
			HuntResult R;
			PairScore cur = S.Score(V, tau, sigma);
			R.best = cur;
			for (int64_t it = 0; it < o.rounds; it++) {
				const std::vector<Move> mv = Moves(V);
				R.moves_available = (int64_t)mv.size();
				if (mv.empty()) break;
				const Move m = mv[(size_t)(rng() % mv.size())];
				ApplyMove(V, m);
				R.tried++;
				if (CreatesCycle(V, m)) { UndoMove(V, m); continue; }
				if (o.verify && detail::CriticalMask(V) != base) R.critical_set_preserved = false;
				const PairScore nw = S.Score(V, tau, sigma);
				if (detail::Accept(nw, cur)) {
					cur = nw;
					R.accepted++;
					if (detail::Better(nw, R.best)) R.best = nw;
					if (detail::IsOne(nw)) break;
				}
				else UndoMove(V, m);
			}
			return R;
		}

		// rotate.hunt_local: 4-cycle moves near the pair's path cells.
		inline HuntResult HuntLocal(GradientView3D& V, int64_t tau, int64_t sigma, const HuntOptions& o) {
			std::mt19937_64 rng(o.seed);
			PairScorer S(V.G);
			const std::vector<uint8_t> base = o.verify ? detail::CriticalMask(V) : std::vector<uint8_t>();
			HuntResult R;
			HuntOptions oz = o;
			oz.margin = -1;
			Zone zone;
			detail::PairZone(V, tau, sigma, oz, zone);
			PairScore cur = S.Score(V, tau, sigma);
			R.best = cur;
			int64_t flat = 0;
			for (int64_t it = 0; it < o.rounds; it++) {
				const std::vector<Move> mv = LocalMoves(V, zone.list);
				if (mv.empty()) break;
				const Move m = mv[(size_t)(rng() % mv.size())];
				ApplyMove(V, m);
				if (CreatesCycle(V, m)) { UndoMove(V, m); continue; }
				if (o.verify && detail::CriticalMask(V) != base) R.critical_set_preserved = false;
				const PairScore nw = S.Score(V, tau, sigma);
				if (detail::Accept(nw, cur)) {
					R.nflip++;
					flat = detail::Better(nw, cur) ? 0 : flat + 1;
					cur = nw;
					if (detail::Better(nw, R.best)) {
						R.best = nw;
						detail::PairZone(V, tau, sigma, oz, zone);
					}
					if (detail::IsOne(nw)) break;
				}
				else {
					UndoMove(V, m);
					flat++;
				}
				if (flat > o.plateau) break;
			}
			return R;
		}

		// rotate.hunt_cycles: alternating cycles of k in [kmin, kmax].
		inline HuntResult HuntCycles(GradientView3D& V, int64_t tau, int64_t sigma, const HuntOptions& o) {
			std::mt19937_64 rng(o.seed);
			PairScorer S(V.G);
			const std::vector<uint8_t> base = o.verify ? detail::CriticalMask(V) : std::vector<uint8_t>();
			HuntResult R;
			Zone zone;
			detail::PairZone(V, tau, sigma, o, zone);
			PairScore cur = S.Score(V, tau, sigma);
			R.best = cur;
			int64_t flat = 0;
			for (int64_t it = 0; it < o.rounds; it++) {
				const std::vector<int64_t> xs = FindCycle(V, zone, rng, o.kmin, o.kmax);
				if (xs.empty()) break;
				if (!CycleOk(V, xs)) continue;
				const CycleToken tok = ApplyCycle(V, xs);
				if (CycleCreatesClosedVPath(V, tok)) { UndoCycle(V, tok); continue; }
				if (o.verify && detail::CriticalMask(V) != base) R.critical_set_preserved = false;
				const PairScore nw = S.Score(V, tau, sigma);
				if (detail::Accept(nw, cur)) {
					R.nflip++;
					R.hist[(int)xs.size()]++;
					flat = detail::Better(nw, cur) ? 0 : flat + 1;
					cur = nw;
					if (detail::Better(nw, R.best)) {
						R.best = nw;
						detail::PairZone(V, tau, sigma, o, zone);
					}
					if (detail::IsOne(nw)) break;
				}
				else {
					UndoCycle(V, tok);
					flat++;
				}
				if (flat > o.plateau) break;
			}
			return R;
		}

		// rotate.hunt_all: 4-cycle moves anywhere, descending on GlobalCost.
		inline HuntResult HuntAll(GradientView3D& V, const HuntOptions& o, bool require_unit) {
			std::mt19937_64 rng(o.seed);
			const std::vector<uint8_t> base = o.verify ? detail::CriticalMask(V) : std::vector<uint8_t>();
			HuntResult R;
			Cost cur = GlobalCost(V, require_unit);
			R.best_cost = cur;
			for (int64_t it = 0; it < o.rounds; it++) {
				const std::vector<Move> mv = Moves(V);
				if (mv.empty()) break;
				const Move m = mv[(size_t)(rng() % mv.size())];
				ApplyMove(V, m);
				if (CreatesCycle(V, m)) { UndoMove(V, m); continue; }
				if (o.verify && detail::CriticalMask(V) != base) R.critical_set_preserved = false;
				const Cost nw = GlobalCost(V, require_unit);
				if (nw <= cur) {
					cur = nw;
					R.accepted++;
					if (nw < R.best_cost) R.best_cost = nw;
					if (nw.count == 0) break;
				}
				else UndoMove(V, m);
			}
			return R;
		}

		// 1-saddle -> {minimum: coefficient}, the join structure.
		inline std::unordered_map<int64_t, std::vector<std::pair<int64_t, ::dmt::S128>>> MinimaRows(const GradientView3D& V) {
			const MorseMatrix M = MorseBoundary(V, 0);
			std::unordered_map<int64_t, std::vector<std::pair<int64_t, ::dmt::S128>>> out;
			for (size_t i = 0; i < M.taus.size(); i++) {
				auto& r = out[M.taus[i]];
				for (int64_t k = M.rowptr[i]; k < M.rowptr[i + 1]; k++)
					r.push_back({ M.entries[(size_t)k].sigma, M.entries[(size_t)k].coeff });
			}
			return out;
		}

		struct HandleSlideCount { int64_t saddles = 0, moved = 0, sign_changed = 0; };

		// Of the 1-saddles of `before`: how many reach a different set of minima
		// in `after`, and how many the same set with a changed coefficient.
		inline HandleSlideCount HandleSlides(const GradientView3D& before, const GradientView3D& after) {
			auto m0 = MinimaRows(before), m1 = MinimaRows(after);
			HandleSlideCount h;
			for (auto& kv : m0) {
				h.saddles++;
				std::map<int64_t, std::pair<uint64_t, uint64_t>> a, b;
				for (auto& e : kv.second) a[e.first] = { e.second.lo, e.second.hi };
				auto it = m1.find(kv.first);
				if (it != m1.end())
					for (auto& e : it->second) b[e.first] = { e.second.lo, e.second.hi };
				bool same_keys = a.size() == b.size();
				if (same_keys)
					for (auto& e : a)
						if (!b.count(e.first)) { same_keys = false; break; }
				if (!same_keys) h.moved++;
				else if (a != b) h.sign_changed++;
			}
			return h;
		}

	} // namespace dmt
} // namespace GInt

#endif
