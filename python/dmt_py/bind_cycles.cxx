// Critical-set-preserving changes of a gradient (alternating-cycle flips) and
// the greedy searches built on them. Every mutating call changes the gradient
// in place and bumps its version.

#include "bind_common.h"
#include "bind_gradient.h"
#include "bigint.h"

#include "gi_dmt_acyclic.h"
#include "gi_dmt_alternating_cycles.h"
#include "gi_dmt_cycle_hunts.h"

#include <memory>
#include <random>

namespace dmt_py {

	namespace {

		using GInt::dmt::Move;

		py::tuple MoveTuple(const Cubical3D& G, const Move& m) {
			return py::make_tuple(TupleFromCell(G, m.b), TupleFromCell(G, m.v), TupleFromCell(G, m.a), TupleFromCell(G, m.q));
		}

		py::list MoveList(const Cubical3D& G, const std::vector<Move>& v) {
			py::list l;
			for (const Move& m : v) l.append(MoveTuple(G, m));
			return l;
		}

		Move MoveFrom(const Cubical3D& G, py::handle t) {
			py::sequence s = py::reinterpret_borrow<py::sequence>(t);
			if (py::len(s) != 4) throw py::value_error("a move is (b, v, a, q)");
			return Move{ CellFromTuple(G, s[0]), CellFromTuple(G, s[1]), CellFromTuple(G, s[2]), CellFromTuple(G, s[3]) };
		}

		std::vector<int64_t> Cells(const Cubical3D& G, py::handle cells) {
			std::vector<int64_t> v;
			for (auto c : cells) v.push_back(CellFromTuple(G, c));
			return v;
		}

		py::set CellSet(const Cubical3D& G, const std::unordered_set<int64_t>& s) {
			py::set out;
			for (int64_t c : s) out.add(TupleFromCell(G, c));
			return out;
		}

		py::tuple ScoreTuple(const GInt::dmt::PairScore& s) { return py::make_tuple(ToPy(s.paths), ToPy(s.coeff)); }

		GInt::dmt::Zone ZoneFrom(const Cubical3D& G, py::handle zone) {
			std::unordered_set<int64_t> s;
			for (auto c : zone) s.insert(CellFromTuple(G, c));
			GInt::dmt::Zone z;
			z.Assign(s);
			return z;
		}

		struct PyRng { std::mt19937_64 r; };

	} // namespace

	void BindCycles(py::module_& m) {
		py::class_<PyRng>(m, "Rng", "A std::mt19937_64 for find_cycle.")
			.def(py::init([](uint64_t seed) { return PyRng{ std::mt19937_64(seed) }; }), py::arg("seed") = 0);

		m.def("moves", [](std::shared_ptr<PyGradient> g) { return MoveList(g->G, GInt::dmt::Moves(g->View())); },
			py::arg("gradient"), "All 4-cycle moves (b, v, a, q): b paired down to v, a paired up to q, b < q, v < a.");
		m.def("local_moves", [](std::shared_ptr<PyGradient> g, py::handle zone) {
			auto z = Cells(g->G, zone);
			std::sort(z.begin(), z.end());
			return MoveList(g->G, GInt::dmt::LocalMoves(g->View(), z));
			}, py::arg("gradient"), py::arg("zone"));
		m.def("interior_moves", [](std::shared_ptr<PyGradient> g, int64_t margin, py::object zone) {
			if (zone.is_none()) return MoveList(g->G, GInt::dmt::InteriorMoves(g->View(), margin, nullptr));
			auto z = Cells(g->G, zone);
			std::sort(z.begin(), z.end());
			return MoveList(g->G, GInt::dmt::InteriorMoves(g->View(), margin, &z));
			}, py::arg("gradient"), py::arg("margin"), py::arg("zone") = py::none());
		m.def("apply_move", [](std::shared_ptr<PyGradient> g, py::handle mv) {
			auto V = g->View();
			GInt::dmt::ApplyMove(V, MoveFrom(g->G, mv));
			g->version++;
			}, py::arg("gradient"), py::arg("move"));
		m.def("undo_move", [](std::shared_ptr<PyGradient> g, py::handle mv) {
			auto V = g->View();
			GInt::dmt::UndoMove(V, MoveFrom(g->G, mv));
			g->version++;
			}, py::arg("gradient"), py::arg("move"));
		m.def("creates_cycle_at", [](std::shared_ptr<PyGradient> g, py::handle cells) {
			auto v = Cells(g->G, cells);
			return GInt::dmt::CreatesCycleAt(g->View(), v.data(), (int)v.size());
			}, py::arg("gradient"), py::arg("cells"),
			"After a flip that re-matched `cells`, does a closed V-path exist? (sound, local)");

		m.def("find_cycle", [](std::shared_ptr<PyGradient> g, py::handle zone, PyRng& rng, int kmin, int kmax, int tries) -> py::object {
			const auto z = ZoneFrom(g->G, zone);
			const auto xs = GInt::dmt::FindCycle(g->View(), z, rng.r, kmin, kmax, tries);
			if (xs.empty()) return py::none();
			py::list l;
			for (int64_t c : xs) l.append(TupleFromCell(g->G, c));
			return l;
			}, py::arg("gradient"), py::arg("zone"), py::arg("rng"), py::arg("kmin") = 2, py::arg("kmax") = 4,
			py::arg("tries") = 60, "One random alternating cycle [x_1..x_k], k in [kmin, kmax], or None.");
		m.def("cycle_ok", [](std::shared_ptr<PyGradient> g, py::handle xs) {
			return GInt::dmt::CycleOk(g->View(), Cells(g->G, xs));
			}, py::arg("gradient"), py::arg("xs"));
		m.def("apply_cycle", [](std::shared_ptr<PyGradient> g, py::handle xs) {
			auto V = g->View();
			const auto tok = GInt::dmt::ApplyCycle(V, Cells(g->G, xs));
			g->version++;
			py::list a, b;
			for (int64_t c : tok.xs) a.append(TupleFromCell(g->G, c));
			for (int64_t c : tok.ms) b.append(TupleFromCell(g->G, c));
			return py::make_tuple(a, b);
			}, py::arg("gradient"), py::arg("xs"), "Flip the cycle; returns the undo token (xs, ms).");
		m.def("undo_cycle", [](std::shared_ptr<PyGradient> g, py::tuple tok) {
			GInt::dmt::CycleToken t;
			t.xs = Cells(g->G, tok[0]);
			t.ms = Cells(g->G, tok[1]);
			auto V = g->View();
			GInt::dmt::UndoCycle(V, t);
			g->version++;
			}, py::arg("gradient"), py::arg("token"));

		m.def("pathcells", [](std::shared_ptr<PyGradient> g, py::tuple tau, py::tuple sigma) {
			return CellSet(g->G, GInt::dmt::PathCellsPair(g->View(), CellFromTuple(g->G, tau), CellFromTuple(g->G, sigma)));
			}, py::arg("gradient"), py::arg("tau"), py::arg("sigma"),
			"Cells on some gradient path from a facet of tau to sigma, plus tau and sigma.");
		m.def("neighbourhood", [](std::vector<int64_t> d, py::handle seed, int hops) {
			const Cubical3D G = Cubical3D::Make(d.at(0), d.at(1), d.at(2));
			std::unordered_set<int64_t> s;
			for (auto c : seed) s.insert(CellFromTuple(G, c));
			return CellSet(G, GInt::dmt::Neighbourhood(G, s, hops));
			}, py::arg("dims"), py::arg("seed"), py::arg("hops"));

		m.def("strangled", [](std::shared_ptr<PyGradient> g, bool require_unit, int p) {
			py::list l;
			for (const auto& c : GInt::dmt::Strangled(g->View(), require_unit, p))
				l.append(py::make_tuple(TupleFromCell(g->G, c.tau), TupleFromCell(g->G, c.sigma), ToPy(c.u), ToPy(c.coeff)));
			return l;
			}, py::arg("gradient"), py::arg("require_unit") = true, py::arg("p") = 1,
			"(tau, sigma, |Gamma|, coeff) for critical connections with more than one path\n"
			"(and coefficient +-1 when require_unit).");
		m.def("global_cost", [](std::shared_ptr<PyGradient> g, bool require_unit) {
			const auto c = GInt::dmt::GlobalCost(g->View(), require_unit);
			return py::make_tuple(c.count, ToPy(c.excess));
			}, py::arg("gradient"), py::arg("require_unit") = true);
		m.def("minima_rows", [](std::shared_ptr<PyGradient> g) {
			py::dict out;
			const auto rows = GInt::dmt::MinimaRows(g->View());
			const auto taus = GInt::dmt::CriticalTupleOrder(g->View(), 1);
			for (int64_t t : taus) {
				py::dict d;
				auto it = rows.find(t);
				if (it != rows.end())
					for (auto& e : it->second) d[TupleFromCell(g->G, e.first)] = ToPy(e.second);
				out[TupleFromCell(g->G, t)] = d;
			}
			return out;
			}, py::arg("gradient"), "1-saddle -> {minimum: coefficient}.");
		m.def("handle_slides", [](std::shared_ptr<PyGradient> before, std::shared_ptr<PyGradient> after) {
			const auto h = GInt::dmt::HandleSlides(before->View(), after->View());
			return py::make_tuple(h.saddles, h.moved, h.sign_changed);
			}, py::arg("before"), py::arg("after"),
			"(1-saddles, how many reach a different set of minima, how many the same set with\n"
			"a changed coefficient).");

		auto opts = [](int kmin, int kmax, int hops, int64_t rounds, int64_t plateau, py::object margin, uint64_t seed, bool verify) {
			GInt::dmt::HuntOptions o;
			o.kmin = kmin; o.kmax = kmax; o.hops = hops; o.rounds = rounds; o.plateau = plateau;
			o.margin = margin.is_none() ? -1 : margin.cast<int64_t>();
			o.seed = seed; o.verify = verify;
			return o;
		};

		m.def("hunt", [opts](std::shared_ptr<PyGradient> g, py::tuple tau, py::tuple sigma, int64_t rounds, uint64_t seed, bool verify) {
			const auto o = opts(2, 4, 3, rounds, 0, py::none(), seed, verify);
			const int64_t t = CellFromTuple(g->G, tau), s = CellFromTuple(g->G, sigma);
			GInt::dmt::HuntResult R;
			{
				py::gil_scoped_release release;
				auto V = g->View();
				R = GInt::dmt::Hunt(V, t, s, o);
			}
			g->version++;
			py::dict d;
			d["best"] = ScoreTuple(R.best);
			d["moves_available"] = R.moves_available;
			d["tried"] = R.tried;
			d["accepted"] = R.accepted;
			d["critical_set_preserved"] = R.critical_set_preserved;
			return d;
			}, py::arg("gradient"), py::arg("tau"), py::arg("sigma"), py::arg("rounds") = 4000, py::arg("seed") = 0,
			py::arg("verify") = false, "Descend on one pair's |Gamma| over 4-cycle moves anywhere (in place).");

		m.def("hunt_local", [opts](std::shared_ptr<PyGradient> g, py::tuple tau, py::tuple sigma, int hops, int64_t rounds,
			uint64_t seed, int64_t plateau, bool verify) {
				const auto o = opts(2, 4, hops, rounds, plateau, py::none(), seed, verify);
				const int64_t t = CellFromTuple(g->G, tau), s = CellFromTuple(g->G, sigma);
				GInt::dmt::HuntResult R;
				{
					py::gil_scoped_release release;
					auto V = g->View();
					R = GInt::dmt::HuntLocal(V, t, s, o);
				}
				g->version++;
				return py::make_tuple(ScoreTuple(R.best), R.nflip, R.critical_set_preserved);
			}, py::arg("gradient"), py::arg("tau"), py::arg("sigma"), py::arg("hops") = 3, py::arg("rounds") = 3000,
				py::arg("seed") = 0, py::arg("plateau") = 600, py::arg("verify") = false,
				"Descend on one pair's |Gamma| over 4-cycle moves near its path cells (in place).\n"
				"Returns (best, flips kept, critical set preserved).");

		m.def("hunt_cycles", [opts](std::shared_ptr<PyGradient> g, py::tuple tau, py::tuple sigma, int kmin, int kmax, int hops,
			int64_t rounds, uint64_t seed, int64_t plateau, py::object margin, bool verify) {
				const auto o = opts(kmin, kmax, hops, rounds, plateau, margin, seed, verify);
				const int64_t t = CellFromTuple(g->G, tau), s = CellFromTuple(g->G, sigma);
				GInt::dmt::HuntResult R;
				{
					py::gil_scoped_release release;
					auto V = g->View();
					R = GInt::dmt::HuntCycles(V, t, s, o);
				}
				g->version++;
				py::dict hist;
				for (auto& kv : R.hist) hist[py::int_(kv.first)] = kv.second;
				return py::make_tuple(ScoreTuple(R.best), R.nflip, hist, R.critical_set_preserved);
			}, py::arg("gradient"), py::arg("tau"), py::arg("sigma"), py::arg("kmin") = 2, py::arg("kmax") = 4,
				py::arg("hops") = 3, py::arg("rounds") = 3000, py::arg("seed") = 0, py::arg("plateau") = 600,
				py::arg("margin") = py::none(), py::arg("verify") = false,
				"Descend on one pair's |Gamma| over alternating cycles of k in [kmin, kmax] near its path\n"
				"cells, optionally only cells >= margin from the boundary (in place).\n"
				"Returns (best, flips kept, {k: flips kept}, critical set preserved).");

		m.def("hunt_all", [opts](std::shared_ptr<PyGradient> g, int64_t rounds, uint64_t seed, bool require_unit, bool verify) {
			const auto o = opts(2, 4, 3, rounds, 0, py::none(), seed, verify);
			GInt::dmt::HuntResult R;
			{
				py::gil_scoped_release release;
				auto V = g->View();
				R = GInt::dmt::HuntAll(V, o, require_unit);
			}
			g->version++;
			return py::make_tuple(py::make_tuple(R.best_cost.count, ToPy(R.best_cost.excess)), R.critical_set_preserved);
			}, py::arg("gradient"), py::arg("rounds") = 6000, py::arg("seed") = 0, py::arg("require_unit") = false,
			py::arg("verify") = false,
			"Descend on global_cost over 4-cycle moves anywhere (in place). Returns (best cost,\n"
			"critical set preserved).");
	}

} // namespace dmt_py
