#ifndef GI_DMT_MORSE_COMPLEX_H
#define GI_DMT_MORSE_COMPLEX_H

// The Morse complex of a discrete gradient (Forman 1998, Sec. 8): critical
// cells as generators and the boundary
//     <d~ tau, sigma> = sum over facets e of tau of <d tau, e> *
//                        (sum of multiplicities of gradient paths e -> sigma)
// from critical (p+1)-cells tau to critical p-cells sigma, with the unsigned
// path count |Gamma| and the largest per-facet count alongside.
//
//   MorseBoundary      one matrix, rows tau in tuple order, each row's
//                      entries in first-touch order;
//   CheckMorseComplex  d~ o d~ = 0 from the signed coefficients;
//   CheckDD            d o d = 0 of the cubical complex itself;
//   BettiModP          Betti numbers of the Morse complex over Z/p.

#include "gi_dmt_gradient_view.h"
#include "gi_dmt_vpath_counter.h"
#include "dmt_core/dmt_count128.h"

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace GInt {
	namespace dmt {

		struct MorseMatrix {
			int p = 0;                          // sigma dimension; tau has p+1
			std::vector<int64_t> taus;          // tuple order
			std::vector<int64_t> rowptr;        // taus.size() + 1
			std::vector<MorseEntry> entries;
		};

		inline std::vector<int64_t> CriticalTupleOrder(const GradientView3D& V, int d) {
			std::vector<int64_t> c = V.Critical(d);
			std::sort(c.begin(), c.end(), [&](int64_t a, int64_t b) { return V.G.TupleLess(a, b); });
			return c;
		}

		// d~ from critical (p+1)-cells to critical p-cells. With levels = true the
		// memo is filled level-parallel first (identical rows).
		inline MorseMatrix MorseBoundary(const GradientView3D& V, int p, bool levels = false, int num_threads = 0) {
			MorseMatrix M;
			M.p = p;
			CounterConfig cfg;
			cfg.sinks = V.Critical(p);
			VPathCounter C(V, cfg);
			M.taus = CriticalTupleOrder(V, p + 1);
			if (levels) {
				std::vector<int64_t> starts;
				for (int64_t t : M.taus) {
					int64_t f[6];
					const int nf = V.G.Facets(t, f);
					for (int i = 0; i < nf; i++) starts.push_back(f[i]);
				}
				C.PrecomputeLevels(starts, num_threads);
			}
			std::vector<MorseEntry> row;
			M.rowptr.push_back(0);
			for (int64_t t : M.taus) {
				MorseRow(C, t, row);
				M.entries.insert(M.entries.end(), row.begin(), row.end());
				M.rowptr.push_back((int64_t)M.entries.size());
			}
			return M;
		}

		struct MorseComplexCheck {
			std::vector<MorseMatrix> B;     // B[p] maps critical (p+1)-cells to critical p-cells
			struct Bad { int p; int64_t tau; int64_t sigma2; ::dmt::S128 value; };
			std::vector<Bad> bad;           // nonzero entries of d~ o d~
			bool exact = true;              // false if a product overflowed
		};

		// d~_p o d~_{p+1} = 0 for p = 1 .. dmax-1 (as morse_signs.check_morse_complex).
		inline MorseComplexCheck CheckMorseComplex(const GradientView3D& V, int dmax = 3, int num_threads = 0) {
			MorseComplexCheck R;
			for (int p = 0; p < dmax; p++) R.B.push_back(MorseBoundary(V, p, false, num_threads));
			for (int p = 1; p < dmax; p++) {
				const MorseMatrix& hi = R.B[(size_t)p];      // (p+1) -> p
				const MorseMatrix& lo = R.B[(size_t)p - 1];  // p -> p-1
				std::unordered_map<int64_t, size_t> row_of;
				for (size_t i = 0; i < lo.taus.size(); i++) row_of[lo.taus[i]] = i;
				for (size_t i = 0; i < hi.taus.size(); i++) {
					std::vector<int64_t> order;
					std::unordered_map<int64_t, ::dmt::S128> acc;
					for (int64_t k = hi.rowptr[i]; k < hi.rowptr[i + 1]; k++) {
						const MorseEntry& e = hi.entries[(size_t)k];
						auto it = row_of.find(e.sigma);
						if (it == row_of.end()) continue;
						const size_t j = it->second;
						for (int64_t m = lo.rowptr[j]; m < lo.rowptr[j + 1]; m++) {
							const MorseEntry& e2 = lo.entries[(size_t)m];
							::dmt::S128 prod = ::dmt::Mul(e.coeff, e2.coeff);
							if (prod.overflow) R.exact = false;
							auto a = acc.find(e2.sigma);
							if (a == acc.end()) { acc[e2.sigma] = prod; order.push_back(e2.sigma); }
							else a->second = ::dmt::Add(a->second, prod);
						}
					}
					for (int64_t s2 : order) {
						const ::dmt::S128& v = acc[s2];
						if (v.overflow) R.exact = false;
						if (!v.IsZero()) R.bad.push_back(MorseComplexCheck::Bad{ p, hi.taus[i], s2, v });
					}
				}
			}
			return R;
		}

		// Nonzero entries of d o d over the cubical complex of a grid.
		inline int64_t CheckDD(const ::dmt::Cubical3D& G) {
			int64_t bad = 0;
			for (int64_t c = 0; c < G.NumCells(); c++) {
				if (G.Dim(c) < 2) continue;
				std::unordered_map<int64_t, int> acc;
				int64_t f[6];
				const int nf = G.Facets(c, f);
				for (int i = 0; i < nf; i++) {
					const int i1 = G.Incidence(c, f[i]);
					int64_t g[6];
					const int ng = G.Facets(f[i], g);
					for (int j = 0; j < ng; j++) acc[g[j]] += i1 * G.Incidence(f[i], g[j]);
				}
				for (auto& kv : acc)
					if (kv.second) bad++;
			}
			return bad;
		}

		inline uint32_t ModP(const ::dmt::S128& v, uint32_t p) {
			const ::dmt::U128 m = ::dmt::Abs(v);
			const uint64_t two64 = (0xFFFFFFFFFFFFFFFFull % p + 1) % p;   // 2^64 mod p
			uint64_t r = ((m.hi % p) * two64 + (m.lo % p)) % p;
			if (v.Negative() && r) r = p - r;
			return (uint32_t)r;
		}

		inline uint64_t PowMod(uint64_t a, uint64_t e, uint64_t p) {
			uint64_t r = 1 % p;
			a %= p;
			while (e) {
				if (e & 1) r = r * a % p;
				a = a * a % p;
				e >>= 1;
			}
			return r;
		}

		// Rank over Z/p of a sparse matrix given as columns of (row, value).
		inline int64_t RankModP(int64_t nrows, std::vector<std::vector<std::pair<int64_t, uint32_t>>> cols, uint32_t p) {
			// Column elimination with pivot rows; columns stored as dense maps.
			std::vector<std::unordered_map<int64_t, uint64_t>> pivot_col;  // by pivot row
			std::unordered_map<int64_t, size_t> pivot_of_row;
			int64_t rank = 0;
			for (auto& col : cols) {
				std::unordered_map<int64_t, uint64_t> v;
				for (auto& e : col)
					if (e.second % p) v[e.first] = (v[e.first] + e.second) % p;
				while (true) {
					// lowest row with a nonzero entry
					int64_t low = -1;
					for (auto& kv : v)
						if (kv.second && (low < 0 || kv.first < low)) low = kv.first;
					if (low < 0) break;
					auto it = pivot_of_row.find(low);
					if (it == pivot_of_row.end()) {
						pivot_of_row[low] = pivot_col.size();
						pivot_col.push_back(v);
						rank++;
						break;
					}
					const auto& pc = pivot_col[it->second];
					const uint64_t f = v[low] * PowMod(pc.at(low), p - 2, p) % p;
					for (auto& kv : pc) {
						uint64_t& x = v[kv.first];
						x = (x + p - f * kv.second % p) % p;
					}
					for (auto i = v.begin(); i != v.end();) {
						if (i->second == 0) i = v.erase(i);
						else ++i;
					}
				}
			}
			(void)nrows;
			return rank;
		}

		// Betti numbers b0..b3 of the Morse complex over Z/p, from its own
		// boundary matrices (exact modular arithmetic).
		inline std::vector<int64_t> BettiModP(const GradientView3D& V, uint32_t prime = 2147483647u, int num_threads = 0) {
			std::vector<int64_t> n(4, 0), rank(5, 0);
			for (int d = 0; d < 4; d++) n[(size_t)d] = (int64_t)V.Critical(d).size();
			for (int p = 1; p <= 3; p++) {
				const MorseMatrix M = MorseBoundary(V, p - 1, false, num_threads);
				std::unordered_map<int64_t, int64_t> row_index;
				const std::vector<int64_t> rows = CriticalTupleOrder(V, p - 1);
				for (size_t i = 0; i < rows.size(); i++) row_index[rows[i]] = (int64_t)i;
				std::vector<std::vector<std::pair<int64_t, uint32_t>>> cols(M.taus.size());
				for (size_t j = 0; j < M.taus.size(); j++)
					for (int64_t k = M.rowptr[j]; k < M.rowptr[j + 1]; k++) {
						const MorseEntry& e = M.entries[(size_t)k];
						const uint32_t v = ModP(e.coeff, prime);
						if (v) cols[j].push_back({ row_index[e.sigma], v });
					}
				rank[(size_t)p] = RankModP((int64_t)rows.size(), cols, prime);
			}
			std::vector<int64_t> b(4);
			for (int d = 0; d < 4; d++) b[(size_t)d] = n[(size_t)d] - rank[(size_t)d] - rank[(size_t)d + 1];
			return b;
		}

	} // namespace dmt
} // namespace GInt

#endif
