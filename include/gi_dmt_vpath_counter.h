#ifndef GI_DMT_VPATH_COUNTER_H
#define GI_DMT_VPATH_COUNTER_H

// Counting gradient paths (V-paths) of a discrete gradient, unsigned and
// signed (Forman, "Morse theory for cell complexes", Adv. Math. 1998, Defs.
// 8.4 and 8.6).
//
// A gradient path between p-cells steps sigma_i -> sigma_{i+1}, where
// sigma_i is paired with the (p+1)-cell V(sigma_i) and sigma_{i+1} != sigma_i
// is a facet of V(sigma_i). Its multiplicity is the product over its steps of
//     -<d V(sigma_i), sigma_i> * <d V(sigma_i), sigma_{i+1}>,
// so it is +1 or -1. For each p-cell e the counter computes, for every sink
// cell k reachable from e, the number of gradient paths e -> k (U128) and the
// sum of their multiplicities (S128). Paths end at the first sink they reach.
//
// Rows are kept in first-touch order: children in facet order (axes
// ascending, -1 before +1), each child's entries in its own order. That order
// is a property of the gradient, not of the evaluation schedule, so the serial
// DFS and the level-parallel evaluation produce identical rows.
//
// A closed V-path throws ClosedVPath with a cell on the cycle; a count that
// does not fit in 128 bits throws PathCountOverflow.

#include "gi_dmt_gradient_view.h"
#include "dmt_core/dmt_count128.h"
#include "gi_dmt_lowerstar_gradient.h"   // ResolveThreads

#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GInt {
	namespace dmt {

		class ClosedVPath : public std::runtime_error {
		public:
			int64_t cell;
			explicit ClosedVPath(int64_t c) : std::runtime_error("closed V-path at cell " + std::to_string(c)), cell(c) {}
		};

		class PathCountOverflow : public std::overflow_error {
		public:
			PathCountOverflow() : std::overflow_error("gradient-path count exceeds 128 bits") {}
		};

		struct SinkEntry {
			int64_t sink;
			::dmt::U128 u;   // number of paths
			::dmt::S128 s;   // sum of multiplicities
		};

		struct CounterConfig {
			std::vector<int64_t> sinks;     // paths end at the first sink reached
			std::vector<int64_t> blocked;   // cells a path may not enter
			bool signed_counts = true;
		};

		class VPathCounter {
		public:
			VPathCounter(const GradientView3D& view, const CounterConfig& cfg) : V(view), signed_(cfg.signed_counts) {
				node_of.assign((size_t)V.G.NumCells(), kWhite);
				for (size_t i = 0; i < cfg.sinks.size(); i++) {
					const int64_t s = cfg.sinks[i];
					if (sink_index.count(s)) continue;
					const int32_t k = (int32_t)sink_index.size();
					sink_index[s] = k;
					sink_cells.push_back(s);
				}
				for (int64_t b : cfg.blocked) blocked.insert(b);
			}

			int64_t NumSinks() const { return (int64_t)sink_cells.size(); }
			bool IsSink(int64_t c) const { return sink_index.count(c) != 0; }
			bool IsBlocked(int64_t c) const { return !blocked.empty() && blocked.count(c) != 0; }
			const GradientView3D& View() const { return V; }
			size_t MemoBytes() const {
				size_t b = node_of.size() * sizeof(int32_t);
				for (auto& r : rows) b += r.capacity() * sizeof(SinkEntry);
				return b;
			}

			// Paths from p-cell e to the sinks; = morse_signs.grad_paths.
			const std::vector<SinkEntry>& PathsFrom(int64_t e) {
				if (node_of[(size_t)e] >= 0) return rows[(size_t)node_of[(size_t)e]];
				std::vector<Frame> stack;
				stack.push_back(Frame{ e, false });
				Scratch sc(NumSinks());
				while (!stack.empty()) {
					const Frame fr = stack.back();
					stack.pop_back();
					if (fr.done) {
						Finish(fr.cell, sc);
						continue;
					}
					const int32_t st = node_of[(size_t)fr.cell];
					if (st >= 0) continue;
					int64_t kids[6];
					const int nk = Expand(fr.cell, kids);
					if (nk < 0) continue;   // a leaf, already stored
					node_of[(size_t)fr.cell] = kGrey;
					stack.push_back(Frame{ fr.cell, true });
					for (int i = 0; i < nk; i++) {
						const int32_t ks = node_of[(size_t)kids[i]];
						if (ks == kWhite) stack.push_back(Frame{ kids[i], false });
						else if (ks == kGrey) throw ClosedVPath(kids[i]);
					}
				}
				return rows[(size_t)node_of[(size_t)e]];
			}

			// Fill the rows of every cell reachable from `starts`, one level
			// (longest distance to a leaf) at a time, in parallel.
			void PrecomputeLevels(const std::vector<int64_t>& starts, int num_threads = 0) {
				// 1. reachable nodes and their children
				std::vector<int64_t> nodes;
				std::unordered_map<int64_t, int32_t> local;
				std::vector<std::vector<int32_t>> kids;
				std::vector<int64_t> stack(starts.begin(), starts.end());
				auto touch = [&](int64_t c) -> int32_t {
					auto it = local.find(c);
					if (it != local.end()) return it->second;
					const int32_t id = (int32_t)nodes.size();
					local[c] = id;
					nodes.push_back(c);
					kids.emplace_back();
					return id;
				};
				for (int64_t s : starts) touch(s);
				std::vector<char> expanded;
				while (!stack.empty()) {
					const int64_t c = stack.back();
					stack.pop_back();
					const int32_t id = touch(c);
					if ((size_t)id < expanded.size() && expanded[(size_t)id]) continue;
					if (expanded.size() < nodes.size()) expanded.resize(nodes.size(), 0);
					expanded[(size_t)id] = 1;
					if (node_of[(size_t)c] >= 0) continue;     // already computed earlier
					int64_t k[6];
					const int nk = Children(c, k);
					for (int i = 0; i < nk; i++) {
						const int32_t kid = touch(k[i]);
						if (expanded.size() < nodes.size()) expanded.resize(nodes.size(), 0);
						kids[(size_t)id].push_back(kid);
						if (!expanded[(size_t)kid]) stack.push_back(k[i]);
					}
				}
				// 2. Kahn on reversed edges: a node is ready once its children are done
				const size_t n = nodes.size();
				std::vector<int32_t> pending(n, 0);
				std::vector<std::vector<int32_t>> parents(n);
				for (size_t i = 0; i < n; i++) {
					if (node_of[(size_t)nodes[i]] >= 0) continue;
					for (int32_t k : kids[i]) {
						if (node_of[(size_t)nodes[(size_t)k]] >= 0) continue;
						pending[i]++;
						parents[(size_t)k].push_back((int32_t)i);
					}
				}
				std::vector<int32_t> wave;
				for (size_t i = 0; i < n; i++)
					if (node_of[(size_t)nodes[i]] < 0 && pending[i] == 0) wave.push_back((int32_t)i);
				// rows are written into slots reserved up front, so threads never resize `rows`
				size_t done = 0;
				const int nt = ResolveThreads(num_threads);
				for (size_t i = 0; i < n; i++)
					if (node_of[(size_t)nodes[i]] < 0) node_of[(size_t)nodes[i]] = kReserved;
				std::vector<int32_t> slot(n, -1);
				for (size_t i = 0; i < n; i++) {
					if (node_of[(size_t)nodes[i]] != kReserved) continue;
					slot[i] = (int32_t)rows.size();
					rows.emplace_back();
				}
				size_t todo = 0;
				for (size_t i = 0; i < n; i++)
					if (slot[i] >= 0) todo++;
				while (!wave.empty()) {
					const long long w = (long long)wave.size();
#pragma omp parallel num_threads(nt)
					{
						Scratch sc(NumSinks());
#pragma omp for schedule(dynamic, 64)
						for (long long i = 0; i < w; i++) {
							const int32_t id = wave[(size_t)i];
							ComputeRowInto(nodes[(size_t)id], rows[(size_t)slot[(size_t)id]], sc, &slot, &local);
						}
					}
					for (int32_t id : wave) node_of[(size_t)nodes[(size_t)id]] = slot[(size_t)id];
					done += wave.size();
					std::vector<int32_t> next;
					for (int32_t id : wave)
						for (int32_t p : parents[(size_t)id])
							if (--pending[(size_t)p] == 0) next.push_back(p);
					wave.swap(next);
				}
				if (done != todo) {
					for (size_t i = 0; i < n; i++)
						if (slot[i] >= 0 && node_of[(size_t)nodes[i]] == kReserved) {
							// leave the counter consistent, then report
							for (size_t j = 0; j < n; j++)
								if (node_of[(size_t)nodes[j]] == kReserved) node_of[(size_t)nodes[j]] = kWhite;
							throw ClosedVPath(nodes[i]);
						}
				}
			}

		private:
			static constexpr int32_t kWhite = -1, kGrey = -2, kReserved = -3;
			struct Frame { int64_t cell; bool done; };

			// Merge buffer keyed by sink index, reset after each row.
			struct Scratch {
				std::vector<int32_t> pos;
				std::vector<int32_t> touched;
				explicit Scratch(int64_t nsinks) : pos((size_t)nsinks, -1) {}
			};

			GradientView3D V;
			bool signed_;
			std::unordered_map<int64_t, int32_t> sink_index;
			std::vector<int64_t> sink_cells;
			std::unordered_set<int64_t> blocked;
			std::vector<int32_t> node_of;                 // per cell: row index, or white/grey/reserved
			std::vector<std::vector<SinkEntry>> rows;

			// Children of an internal node: facets of its partner other than
			// itself (blocked cells excluded). Returns -1 for a leaf.
			int Children(int64_t c, int64_t out[6]) const {
				if (IsSink(c) || IsBlocked(c) || !V.PairedUp(c)) return -1;
				const int64_t q = V.Partner(c);
				if (IsBlocked(q)) return -1;
				int64_t f[6];
				const int nf = V.G.Facets(q, f);
				int k = 0;
				for (int i = 0; i < nf; i++)
					if (f[i] != c) out[k++] = f[i];
				return k;
			}

			// Store the row of a leaf and return -1, or return the children.
			int Expand(int64_t c, int64_t out[6]) {
				const int nk = Children(c, out);
				if (nk >= 0) return nk;
				node_of[(size_t)c] = (int32_t)rows.size();
				rows.emplace_back();
				if (IsSink(c) && !IsBlocked(c)) {
					SinkEntry e;
					e.sink = c;
					e.u = ::dmt::U128(1);
					e.s = ::dmt::S128::FromInt(1);
					rows.back().push_back(e);
				}
				return -1;
			}

			void Finish(int64_t c, Scratch& sc) {
				std::vector<SinkEntry> row;
				ComputeRowInto(c, row, sc, nullptr, nullptr);
				node_of[(size_t)c] = (int32_t)rows.size();
				rows.push_back(std::move(row));
			}

			// Row of c from its children's rows. With `slot`/`local` (the
			// level-parallel path) children may also be found by local id.
			void ComputeRowInto(int64_t c, std::vector<SinkEntry>& row, Scratch& sc,
				const std::vector<int32_t>* slot, const std::unordered_map<int64_t, int32_t>* local) {
				row.clear();
				int64_t kids[6];
				const int nk = Children(c, kids);
				if (nk < 0) {
					if (IsSink(c) && !IsBlocked(c)) {
						SinkEntry e;
						e.sink = c;
						e.u = ::dmt::U128(1);
						e.s = ::dmt::S128::FromInt(1);
						row.push_back(e);
					}
					return;
				}
				const int64_t q = V.Partner(c);
				const int ie = V.G.Incidence(q, c);
				for (int i = 0; i < nk; i++) {
					const int64_t k = kids[i];
					int32_t r = node_of[(size_t)k];
					if (r < 0 && slot && local) {
						auto it = local->find(k);
						if (it != local->end()) r = (*slot)[(size_t)it->second];
					}
					if (r < 0) throw ClosedVPath(k);
					const int fac = -ie * V.G.Incidence(q, k);
					for (const SinkEntry& ent : rows[(size_t)r]) {
						const int32_t si = sink_index.find(ent.sink)->second;
						int32_t& p = sc.pos[(size_t)si];
						if (p < 0) {
							p = (int32_t)row.size();
							sc.touched.push_back(si);
							SinkEntry z;
							z.sink = ent.sink;
							row.push_back(z);
						}
						SinkEntry& t = row[(size_t)p];
						t.u = ::dmt::Add(t.u, ent.u);
						if (signed_) t.s = ::dmt::Add(t.s, ::dmt::MulSign(ent.s, fac));
						if (t.u.overflow || t.s.overflow) {
							for (int32_t x : sc.touched) sc.pos[(size_t)x] = -1;
							sc.touched.clear();
							throw PathCountOverflow();
						}
					}
				}
				for (int32_t x : sc.touched) sc.pos[(size_t)x] = -1;
				sc.touched.clear();
			}
		};

		// One row of the Morse boundary: for each sink sigma reached from the
		// facets of tau, the number of paths |Gamma|, the coefficient
		// <d~ tau, sigma> = sum over facets e of <d tau, e> * (signed paths
		// e -> sigma), and the largest number of paths through a single facet.
		struct MorseEntry {
			int64_t sigma;
			::dmt::U128 u;
			::dmt::S128 coeff;
			::dmt::U128 best;
		};

		inline void MorseRow(VPathCounter& C, int64_t tau, std::vector<MorseEntry>& out) {
			out.clear();
			std::unordered_map<int64_t, size_t> pos;
			const GradientView3D& V = C.View();
			int64_t f[6];
			const int nf = V.G.Facets(tau, f);
			for (int i = 0; i < nf; i++) {
				const std::vector<SinkEntry>& r = C.PathsFrom(f[i]);
				const int ie = V.G.Incidence(tau, f[i]);
				for (const SinkEntry& e : r) {
					auto it = pos.find(e.sink);
					if (it == pos.end()) {
						it = pos.emplace(e.sink, out.size()).first;
						MorseEntry z;
						z.sigma = e.sink;
						out.push_back(z);
					}
					MorseEntry& t = out[it->second];
					t.u = ::dmt::Add(t.u, e.u);
					t.coeff = ::dmt::Add(t.coeff, ::dmt::MulSign(e.s, ie));
					if (::dmt::Less(t.best, e.u)) t.best = e.u;
					if (t.u.overflow || t.coeff.overflow) throw PathCountOverflow();
				}
			}
			// Python keeps entries with a nonzero path count only.
			size_t w = 0;
			for (size_t i = 0; i < out.size(); i++)
				if (!out[i].u.IsZero()) out[w++] = out[i];
			out.resize(w);
		}

		// |Gamma| and the coefficient for one pair (tau, sigma): a scalar DP
		// toward sigma alone. Exact: with every critical p-cell a sink, the other
		// sinks only add entries for themselves, and with sigma the only sink a
		// critical cell other than sigma contributes nothing -- either way sigma's
		// entry is the same. Scratch arrays are epoch-stamped and reused.
		struct PairScore {
			::dmt::U128 paths;
			::dmt::S128 coeff;
		};

		class PairScorer {
		public:
			explicit PairScorer(const ::dmt::Cubical3D& G)
				: epoch_((size_t)G.NumCells(), 0), state_((size_t)G.NumCells(), 0),
				u_((size_t)G.NumCells()), s_((size_t)G.NumCells()) {}

			PairScore Score(const GradientView3D& V, int64_t tau, int64_t sigma) {
				if (++cur_ == 0) { std::fill(epoch_.begin(), epoch_.end(), 0u); cur_ = 1; }
				PairScore r;
				if (!V.IsCritical(sigma)) return r;   // only critical cells are sinks
				int64_t f[6];
				const int nf = V.G.Facets(tau, f);
				for (int i = 0; i < nf; i++) {
					Eval(V, f[i], sigma);
					const int ie = V.G.Incidence(tau, f[i]);
					r.paths = ::dmt::Add(r.paths, u_[(size_t)f[i]]);
					r.coeff = ::dmt::Add(r.coeff, ::dmt::MulSign(s_[(size_t)f[i]], ie));
				}
				if (r.paths.overflow || r.coeff.overflow) throw PathCountOverflow();
				return r;
			}

		private:
			std::vector<uint32_t> epoch_;
			std::vector<uint8_t> state_;     // 1 grey, 2 done (valid when epoch matches)
			std::vector<::dmt::U128> u_;
			std::vector<::dmt::S128> s_;
			uint32_t cur_ = 0;
			struct Frame { int64_t cell; bool done; };

			bool Done(int64_t c) const { return epoch_[(size_t)c] == cur_ && state_[(size_t)c] == 2; }

			void Eval(const GradientView3D& V, int64_t start, int64_t sigma) {
				if (Done(start)) return;
				std::vector<Frame> stack;
				stack.push_back(Frame{ start, false });
				while (!stack.empty()) {
					const Frame fr = stack.back();
					stack.pop_back();
					const size_t c = (size_t)fr.cell;
					if (fr.done) {
						const int64_t q = V.Partner(fr.cell);
						const int ie = V.G.Incidence(q, fr.cell);
						int64_t f[6];
						const int nf = V.G.Facets(q, f);
						::dmt::U128 u;
						::dmt::S128 s;
						for (int i = 0; i < nf; i++) {
							if (f[i] == fr.cell) continue;
							u = ::dmt::Add(u, u_[(size_t)f[i]]);
							s = ::dmt::Add(s, ::dmt::MulSign(s_[(size_t)f[i]], -ie * V.G.Incidence(q, f[i])));
						}
						u_[c] = u;
						s_[c] = s;
						state_[c] = 2;
						continue;
					}
					if (epoch_[c] == cur_) continue;      // done or grey
					epoch_[c] = cur_;
					if (fr.cell == sigma) {
						u_[c] = ::dmt::U128(1);
						s_[c] = ::dmt::S128::FromInt(1);
						state_[c] = 2;
						continue;
					}
					if (!V.PairedUp(fr.cell)) {
						u_[c] = ::dmt::U128();
						s_[c] = ::dmt::S128();
						state_[c] = 2;
						continue;
					}
					state_[c] = 1;
					stack.push_back(Frame{ fr.cell, true });
					const int64_t q = V.Partner(fr.cell);
					int64_t f[6];
					const int nf = V.G.Facets(q, f);
					for (int i = 0; i < nf; i++) {
						if (f[i] == fr.cell) continue;
						const size_t k = (size_t)f[i];
						if (epoch_[k] != cur_) stack.push_back(Frame{ f[i], false });
						else if (state_[k] == 1) throw ClosedVPath(f[i]);
					}
				}
			}
		};

	} // namespace dmt
} // namespace GInt

#endif
