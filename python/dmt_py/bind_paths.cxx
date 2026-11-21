// Gradient-path counting and the signed Morse boundary.
//
// Python shapes mirror the reference implementation: a Morse boundary is a
// dict {tau: {sigma: (|Gamma|, coefficient, best per-facet count)}} with taus
// in tuple order and each row in first-touch order; paths_from returns
// {sink: (paths, signed sum)}. Counts are exact Python ints.

#include "bind_common.h"
#include "bind_gradient.h"
#include "bigint.h"

#include "dmt/gi_dmt_free_paths.h"
#include "dmt/gi_dmt_morse_complex.h"
#include "dmt/gi_dmt_vpath_counter.h"

#include <memory>

namespace dmt_py {

	namespace {

		py::dict RowToDict(const Cubical3D& G, const std::vector<GInt::dmt::MorseEntry>& row) {
			py::dict d;
			for (const auto& e : row)
				d[TupleFromCell(G, e.sigma)] = py::make_tuple(ToPy(e.u), ToPy(e.coeff), ToPy(e.best));
			return d;
		}

		py::dict MatrixToDict(const Cubical3D& G, const GInt::dmt::MorseMatrix& M) {
			py::dict out;
			for (size_t i = 0; i < M.taus.size(); i++) {
				py::dict d;
				for (int64_t k = M.rowptr[i]; k < M.rowptr[i + 1]; k++) {
					const auto& e = M.entries[(size_t)k];
					d[TupleFromCell(G, e.sigma)] = py::make_tuple(ToPy(e.u), ToPy(e.coeff), ToPy(e.best));
				}
				out[TupleFromCell(G, M.taus[i])] = d;
			}
			return out;
		}

		std::vector<int64_t> CellList(const Cubical3D& G, py::handle cells) {
			std::vector<int64_t> v;
			for (auto c : cells) v.push_back(CellFromTuple(G, c));
			return v;
		}

		py::dict FacetRowsToDict(const Cubical3D& G, const std::vector<GInt::dmt::FacetRow>& fp) {
			py::dict out;
			for (const auto& fr : fp) {
				py::dict d;
				for (const auto& e : fr.row) d[TupleFromCell(G, e.sink)] = ToPy(e.u);
				out[TupleFromCell(G, fr.facet)] = d;
			}
			return out;
		}

		py::list CellsToList(const Cubical3D& G, const std::vector<int64_t>& v) {
			py::list l;
			for (int64_t c : v) l.append(TupleFromCell(G, c));
			return l;
		}

		struct PyPathCounter {
			std::shared_ptr<PyGradient> g;
			uint64_t version;
			std::unique_ptr<GInt::dmt::VPathCounter> C;

			void Check() const {
				if (g->version != version)
					throw std::runtime_error("the gradient changed after this PathCounter was built");
			}
		};

	} // namespace

	void BindPaths(py::module_& m) {
		// ClosedVPath -> msceer_dmt.ClosedVPathError (a RuntimeError);
		// PathCountOverflow is a std::overflow_error -> OverflowError.
		py::register_exception<GInt::dmt::ClosedVPath>(m, "ClosedVPathError", PyExc_RuntimeError);

		py::class_<PyPathCounter>(m, "PathCounter",
			"Memoised gradient-path counts from p-cells to a fixed set of sinks. The memo\n"
			"stays valid for this gradient (until it changes), these sinks and this blocked set.")
			.def(py::init([](std::shared_ptr<PyGradient> g, py::object p, py::object sinks, py::object blocked,
				bool signed_counts) {
			GInt::dmt::CounterConfig cfg;
			cfg.signed_counts = signed_counts;
			if (!sinks.is_none()) cfg.sinks = CellList(g->G, sinks);
			else if (!p.is_none()) cfg.sinks = g->View().Critical(p.cast<int>());
			else throw py::value_error("give p (sinks = critical p-cells) or sinks");
			if (!blocked.is_none()) cfg.blocked = CellList(g->G, blocked);
			PyPathCounter pc;
			pc.g = g;
			pc.version = g->version;
			pc.C.reset(new GInt::dmt::VPathCounter(g->View(), cfg));
			return pc;
				}), py::arg("gradient"), py::arg("p") = py::none(), py::arg("sinks") = py::none(),
				py::arg("blocked") = py::none(), py::arg("signed") = true)
			.def("paths_from", [](PyPathCounter& pc, py::tuple cell) {
			pc.Check();
			const auto& row = pc.C->PathsFrom(CellFromTuple(pc.g->G, cell));
			py::dict d;
			for (const auto& e : row) d[TupleFromCell(pc.g->G, e.sink)] = py::make_tuple(ToPy(e.u), ToPy(e.s));
			return d;
				}, py::arg("cell"), "{sink: (paths, signed sum)} over gradient paths from this cell.")
			.def("morse_row", [](PyPathCounter& pc, py::tuple tau) {
			pc.Check();
			std::vector<GInt::dmt::MorseEntry> row;
			GInt::dmt::MorseRow(*pc.C, CellFromTuple(pc.g->G, tau), row);
			return RowToDict(pc.g->G, row);
				}, py::arg("tau"), "{sigma: (|Gamma|, <d~ tau, sigma>, best per-facet count)}.")
			.def("precompute", [](PyPathCounter& pc, py::handle starts, int num_threads) {
			pc.Check();
			auto v = CellList(pc.g->G, starts);
			py::gil_scoped_release release;
			pc.C->PrecomputeLevels(v, num_threads);
				}, py::arg("starts"), py::arg("num_threads") = 0,
				"Fill the memo for everything reachable from `starts`, level-parallel.")
			.def("facet_sink_paths", [](PyPathCounter& pc, py::tuple src) {
			pc.Check();
			return FacetRowsToDict(pc.g->G, GInt::dmt::FacetSinkPaths(*pc.C, CellFromTuple(pc.g->G, src)));
				}, py::arg("src"),
				"{facet: {sink: paths}} from the 2-cell src (its own paired facet and blocked facets\n"
				"omitted); counts are unsigned. Build the counter with signed=False for speed.")
			.def_property_readonly("memo_bytes", [](const PyPathCounter& pc) { return pc.C->MemoBytes(); });

		// ---- free-endpoint probes and multiplicities ----

		m.def("facet_paths", [](std::shared_ptr<PyGradient> g, py::tuple s2) {
			return FacetRowsToDict(g->G, GInt::dmt::FacetPaths(g->View(), CellFromTuple(g->G, s2)));
			}, py::arg("gradient"), py::arg("s2"),
			"{facet e of s2: {critical 1-cell: paths from s2 through e}}.");

		m.def("mult_at", [](std::shared_ptr<PyGradient> g, py::tuple s2) -> py::object {
			::dmt::U128 best;
			try {
				if (GInt::dmt::MultAt(g->View(), CellFromTuple(g->G, s2), best) < 0) return py::none();
			}
			catch (const GInt::dmt::ClosedVPath&) {
				return py::none();
			}
			return ToPy(best);
			}, py::arg("gradient"), py::arg("s2"),
			"Largest per-facet path count at a critical 2-cell; None if s2 is not critical or a\n"
			"closed V-path is met.");

		m.def("mult_anywhere", [](std::shared_ptr<PyGradient> g, bool interior_only) {
			GInt::dmt::Witness w;
			::dmt::U128 best;
			{
				py::gil_scoped_release release;
				best = GInt::dmt::MultAnywhere(g->View(), interior_only, w);
			}
			py::object where = py::none();
			if (w.q >= 0)
				where = py::make_tuple(TupleFromCell(g->G, w.q), TupleFromCell(g->G, w.facet), TupleFromCell(g->G, w.sink));
			return py::make_tuple(ToPy(best), where);
			}, py::arg("gradient"), py::arg("interior_only") = true,
			"(best, (quad, facet, sink) or None): the largest per-facet count from any 2-cell to a\n"
			"critical 1-cell; the first maximum in tuple order.");

		m.def("multiplicities", [](std::shared_ptr<PyGradient> g, bool interior_only, uint64_t minmult, bool want_hist) {
			std::vector<std::pair<::dmt::U128, int64_t>> hist;
			std::vector<GInt::dmt::Multiplicity> v;
			{
				py::gil_scoped_release release;
				v = GInt::dmt::Multiplicities(g->View(), interior_only, minmult, want_hist ? &hist : nullptr);
			}
			py::list out;
			for (const auto& t : v)
				out.append(py::make_tuple(ToPy(t.n), TupleFromCell(g->G, t.s2), TupleFromCell(g->G, t.facet),
					TupleFromCell(g->G, t.s1)));
			if (!want_hist) return py::object(out);
			py::dict h;
			for (const auto& kv : hist) h[ToPy(kv.first)] = kv.second;
			return py::object(py::make_tuple(out, h));
			}, py::arg("gradient"), py::arg("interior_only") = true, py::arg("minmult") = 2,
			py::arg("want_hist") = false,
			"[(n, s2, facet, s1)] for every critical 2-cell's per-facet counts n >= minmult, sorted\n"
			"by (-n, s2, facet, s1); with want_hist, also {n: how many triples} over all n >= 1.");

		m.def("path_cells", [](std::shared_ptr<PyGradient> g, py::tuple s2, py::tuple facet, py::tuple s1) {
			const auto v = GInt::dmt::PathCells(g->View(), CellFromTuple(g->G, s2), CellFromTuple(g->G, facet),
				CellFromTuple(g->G, s1));
			py::set s;
			for (int64_t c : v) s.add(TupleFromCell(g->G, c));
			return s;
			}, py::arg("gradient"), py::arg("s2"), py::arg("facet"), py::arg("s1"),
			"Cells on the descending paths from s2 through facet to s1.");

		m.def("face_to_face", [](std::shared_ptr<PyGradient> g, bool interior, int layer) {
			GInt::dmt::FaceToFaceResult R;
			{
				py::gil_scoped_release release;
				R = GInt::dmt::FaceToFace(g->View(), interior, layer);
			}
			return py::make_tuple(ToPy(R.best_mult), ToPy(R.best_total), R.n_sources);
			}, py::arg("gradient"), py::arg("interior") = false, py::arg("layer") = 1,
			"(best per-facet count, best total, number of sources) from the interior quads under\n"
			"the top face down to the bottom face; with interior, to the interior 1-cells `layer`\n"
			"steps in, never entering a boundary cell.");

		m.def("census_best", [](std::shared_ptr<PyGradient> g) {
			std::pair<::dmt::U128, ::dmt::U128> r;
			{
				py::gil_scoped_release release;
				r = GInt::dmt::CensusBest(g->View());
			}
			return py::make_tuple(ToPy(r.first), ToPy(r.second));
			}, py::arg("gradient"),
			"(paths, largest single sink) of the busiest interior source quad under the top face,\n"
			"against the edges at the middle vertex of the bottom layer.");

		m.def("top_sources", [](std::vector<int64_t> d) {
			const Cubical3D G = Cubical3D::Make(d.at(0), d.at(1), d.at(2));
			return CellsToList(G, GInt::dmt::TopSources(G));
			}, py::arg("dims"));
		m.def("bottom_middle_sinks", [](std::vector<int64_t> d) {
			const Cubical3D G = Cubical3D::Make(d.at(0), d.at(1), d.at(2));
			return CellsToList(G, GInt::dmt::BottomMiddleSinks(G));
			}, py::arg("dims"));
		m.def("face_sinks", [](std::vector<int64_t> d, int axis, bool at_top) {
			const Cubical3D G = Cubical3D::Make(d.at(0), d.at(1), d.at(2));
			return CellsToList(G, GInt::dmt::FaceSinks(G, axis, at_top));
			}, py::arg("dims"), py::arg("axis") = 2, py::arg("at_top") = false);
		m.def("interior_sinks", [](std::vector<int64_t> d, int layer, int axis) {
			const Cubical3D G = Cubical3D::Make(d.at(0), d.at(1), d.at(2));
			return CellsToList(G, GInt::dmt::InteriorSinks(G, layer, axis));
			}, py::arg("dims"), py::arg("layer") = 1, py::arg("axis") = 2);

		m.def("morse_boundary", [](std::shared_ptr<PyGradient> g, int p, bool levels, int num_threads) {
			GInt::dmt::MorseMatrix M;
			{
				py::gil_scoped_release release;
				M = GInt::dmt::MorseBoundary(g->View(), p, levels, num_threads);
			}
			return MatrixToDict(g->G, M);
			}, py::arg("gradient"), py::arg("p"), py::arg("levels") = false, py::arg("num_threads") = 0,
			"d~ from critical (p+1)-cells to critical p-cells:\n"
			"{tau: {sigma: (|Gamma|, coefficient, best per-facet count)}}.");

		m.def("score", [](std::shared_ptr<PyGradient> g, py::tuple tau, py::tuple sigma) {
			if (!g->scorer) g->scorer = std::make_shared<GInt::dmt::PairScorer>(g->G);
			const GInt::dmt::PairScore s = g->scorer->Score(g->View(), CellFromTuple(g->G, tau), CellFromTuple(g->G, sigma));
			return py::make_tuple(ToPy(s.paths), ToPy(s.coeff));
			}, py::arg("gradient"), py::arg("tau"), py::arg("sigma"),
			"(|Gamma|, coefficient) of the connection tau -> sigma (sigma critical).");

		m.def("check_morse_complex", [](std::shared_ptr<PyGradient> g, int dmax, int num_threads) {
			GInt::dmt::MorseComplexCheck R;
			{
				py::gil_scoped_release release;
				R = GInt::dmt::CheckMorseComplex(g->View(), dmax, num_threads);
			}
			py::dict B;
			for (size_t p = 0; p < R.B.size(); p++) B[py::int_(p + 1)] = MatrixToDict(g->G, R.B[p]);
			py::list bad;
			for (const auto& b : R.bad)
				bad.append(py::make_tuple(b.p, TupleFromCell(g->G, b.tau), TupleFromCell(g->G, b.sigma2), ToPy(b.value)));
			return py::make_tuple(B, bad, R.exact);
			}, py::arg("gradient"), py::arg("dmax") = 3, py::arg("num_threads") = 0,
			"(B, bad, exact): B[p+1] = morse_boundary(p) for p < dmax, bad = nonzero entries\n"
			"(p, tau, sigma2, value) of d~ o d~, exact = no product overflowed.");

		m.def("check_dd", [](std::vector<int64_t> dims) {
			if (dims.size() != 3) throw py::value_error("dims must be (nx, ny, nz)");
			return GInt::dmt::CheckDD(Cubical3D::Make(dims[0], dims[1], dims[2]));
			}, py::arg("dims"), "Nonzero entries of d o d over the cubical complex (0 = consistent signs).");

		m.def("betti_modp", [](std::shared_ptr<PyGradient> g, uint32_t prime, int num_threads) {
			py::gil_scoped_release release;
			return GInt::dmt::BettiModP(g->View(), prime, num_threads);
			}, py::arg("gradient"), py::arg("prime") = 2147483647u, py::arg("num_threads") = 0,
			"Betti numbers b0..b3 of the Morse complex over Z/prime (exact).");
	}

} // namespace dmt_py
