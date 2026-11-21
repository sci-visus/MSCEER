// Gradient3D: a discrete gradient on the cubical complex of a 3D vertex grid,
// owned by C++ (one GradBitfield byte per cell), plus the lower-star rules
// that compute one and a tuple-keyed mapping view for Python code written
// against dict-shaped pairings.

#include "bind_common.h"
#include "bind_gradient.h"

#include "gi_dmt_acyclic.h"
#include "gi_dmt_dmf_gradient.h"
#include "gi_dmt_lowerstar_gradient.h"

#ifdef DMT_PY_GPU
#include "dgrad_gpu_api.h"
#endif

#include <memory>
#include <string>

namespace dmt_py {

	namespace {

		// device is "cpu" or "gpu". The GPU drivers (gpu_dgrad, built with
		// GPU_DGRAD_ENABLED) produce byte-identical results; when a GPU call
		// fails (no device, out of memory) the CPU driver runs instead and
		// info["device"] says which one ran.
		bool WantGpu(const std::string& device) {
			if (device == "cpu") return false;
			if (device != "gpu") throw py::value_error("device must be 'cpu' or 'gpu'");
#ifdef DMT_PY_GPU
			return true;
#else
			throw py::value_error("msceer_dmt was built without GPU support (GPU_DGRAD_ENABLED=OFF)");
#endif
		}

		std::shared_ptr<PyGradient> RunLowerStar(py::array R, GInt::dmt::LowerStarRule rule, bool journal,
			int num_threads, const std::string& device, GInt::dmt::LowerStarStats& st, std::string& used) {
			const bool gpu = WantGpu(device);
			if (gpu && GInt::dmt::RuleIsAligned(rule))
				throw py::value_error("aligned=True has no GPU kernel; use device='cpu'");
			GInt::dmt::RankField3D f = DistinctRanks(R);
			auto out = std::make_shared<PyGradient>();
			out->G = Cubical3D::Make(f.X, f.Y, f.Z);
			const size_t nc = (size_t)out->G.NumCells();
			out->bytes.assign(nc, 0);
			if (journal) {
				out->step.assign(nc, 0);
				out->part.assign(nc, 0);
			}
			used = "cpu";
			{
				py::gil_scoped_release release;
				bool done = false;
#ifdef DMT_PY_GPU
				if (gpu) {
					GInt::gpu::DmtLowerStarStats gs;
					done = GInt::gpu::ComputeLowerStarGradient3DGPU(f.ranks.data(), f.X, f.Y, f.Z,
						rule == GInt::dmt::LowerStarRule::PLS ? GInt::gpu::DmtLowerStarRule::PLS : GInt::gpu::DmtLowerStarRule::PLS_BR,
						out->bytes.data(), journal ? out->step.data() : nullptr, journal ? out->part.data() : nullptr,
						&gs, nullptr);
					if (done) {
						for (int d = 0; d < 4; d++) st.critical[d] = gs.critical[d];
						st.orphans = gs.orphans;
						st.unassigned = gs.unassigned;
						used = "gpu";
					}
				}
#else
				(void)gpu;
#endif
				if (!done)
					GInt::dmt::ComputeLowerStarGradient3D(f.ranks.data(), f.X, f.Y, f.Z, rule, out->bytes.data(),
						journal ? out->step.data() : nullptr, journal ? out->part.data() : nullptr, &st, num_threads);
			}
			if (st.unassigned != 0) throw std::runtime_error("lower-star gradient left cells unassigned");
			return out;
		}

		py::dict StatsDict(const GInt::dmt::LowerStarStats& st, const std::string& used) {
			py::dict d;
			d["critical"] = py::make_tuple(st.critical[0], st.critical[1], st.critical[2], st.critical[3]);
			d["orphans"] = st.orphans;
			d["device"] = used;
			return d;
		}

		// Tuple-keyed mapping over the gradient's pairs, mirroring a Python
		// dict that holds each pair in both directions. Critical cells are
		// absent keys. Assignment writes one side only (as dict assignment
		// does); deletion makes the cell critical.
		struct PairMap {
			std::shared_ptr<PyGradient> g;

			py::object get(py::handle key, py::object dflt) const {
				const int64_t id = CellFromTuple(g->G, key);
				const int64_t p = g->View().Partner(id);
				if (p < 0) return dflt;
				return TupleFromCell(g->G, p);
			}
		};

	} // namespace

	void BindGradient(py::module_& m) {
		py::class_<PyGradient, std::shared_ptr<PyGradient>>(m, "Gradient3D",
			"A discrete gradient on the cubical complex of a 3D vertex grid.\n"
			"Storage is one GradBitfield byte per cell (bit 0 assigned, bits 2-4 direction\n"
			"0..5 = +x,-x,+y,-y,+z,-z, 7 = critical), x fastest.")
			.def_static("pls", [](py::array R, bool journal, int num_threads, const std::string& device, bool aligned) {
			GInt::dmt::LowerStarStats st;
			std::string used;
			auto g = RunLowerStar(R, aligned ? GInt::dmt::LowerStarRule::PLS_ALIGNED : GInt::dmt::LowerStarRule::PLS,
				journal, num_threads, device, st, used);
			return py::make_tuple(g, StatsDict(st, used));
				}, py::arg("R"), py::arg("journal") = false, py::arg("num_threads") = 0, py::arg("device") = "cpu",
				py::arg("aligned") = false,
				"ProcessLowerStars (Robins, Wood, Sheppard 2011) on R[x, y, z], whose values must\n"
				"be pairwise distinct. Returns (gradient, info). With journal=True the gradient\n"
				"keeps, per cell, the index of the decision that classified it within its lower\n"
				"star (step()) and the part of the star (part()). device='gpu' runs the CUDA kernel\n"
				"(identical bytes; the CPU runs instead if the GPU call fails); info['device']\n"
				"names the device that ran.\n"
				"aligned=True is the aligned-first variant: in each lower star, right after the\n"
				"vertex is paired with its steepest edge, every quad containing that edge is paired\n"
				"with its other edge of the star before the queue order takes over, so those edges\n"
				"point the way the vertex does. CPU only.")
			.def_static("pls_br", [](py::array R, bool journal, int num_threads, const std::string& device, bool aligned) {
			GInt::dmt::LowerStarStats st;
			std::string used;
			auto g = RunLowerStar(R, aligned ? GInt::dmt::LowerStarRule::PLS_BR_ALIGNED : GInt::dmt::LowerStarRule::PLS_BR,
				journal, num_threads, device, st, used);
			return py::make_tuple(g, StatsDict(st, used));
				}, py::arg("R"), py::arg("journal") = false, py::arg("num_threads") = 0, py::arg("device") = "cpu",
				py::arg("aligned") = false,
				"Boundary-restricted ProcessLowerStars: each lower star is split into its cells on\n"
				"the domain boundary (part 0) and its interior cells (part 1), collapsed separately,\n"
				"seeding the vertex only in the part containing it. Returns (gradient, info).\n"
				"device and aligned as for pls(); aligned acts in the part that holds the vertex.")
			.def_static("dmf", [](py::array R, bool second_pass, bool want_order, bool sequential, int num_threads,
				const std::string& device) {
			const bool gpu = WantGpu(device);
			if (gpu && sequential) throw py::value_error("sequential=True is the CPU reference loop; use device='cpu'");
			GInt::dmt::RankField3D f = DistinctRanks(R);
			auto g = std::make_shared<PyGradient>();
			g->G = Cubical3D::Make(f.X, f.Y, f.Z);
			const size_t nc = (size_t)g->G.NumCells();
			g->bytes.assign(nc, 0);
			g->stage.assign(nc, 0);
			if (want_order) g->order.assign(nc, 0);
			GInt::dmt::DMFOptions opt;
			opt.second_pass = second_pass;
			opt.sequential_reference = sequential;
			opt.num_threads = num_threads;
			GInt::dmt::DMFStats st;
			std::string used = "cpu";
			{
				py::gil_scoped_release release;
				bool done = false;
#ifdef DMT_PY_GPU
				if (gpu) {
					GInt::gpu::DmtDMFStats gs;
					done = GInt::gpu::ComputeDMFGradient3DGPU(f.ranks.data(), f.X, f.Y, f.Z, second_pass, g->bytes.data(),
						g->stage.data(), want_order ? g->order.data() : nullptr, &gs, nullptr);
					if (done) {
						st.guard1 = gs.guard1; st.guard2 = gs.guard2; st.alg1 = gs.alg1; st.alg2 = gs.alg2;
						for (int k = 0; k < 4; k++) st.critical[k] = gs.critical[k];
						used = "gpu";
					}
				}
#else
				(void)gpu;
#endif
				if (!done)
					GInt::dmt::ComputeDMFGradient3D(f.ranks.data(), f.X, f.Y, f.Z, opt, g->bytes.data(), g->stage.data(),
						want_order ? g->order.data() : nullptr, &st);
			}
			py::dict d;
			d["device"] = used;
			d["guard1"] = st.guard1;
			d["guard2"] = st.guard2;
			d["alg1"] = st.alg1;
			d["alg2"] = st.alg2;
			d["critical"] = py::make_tuple(st.critical[0], st.critical[1], st.critical[2], st.critical[3]);
			return py::make_tuple(g, d);
				}, py::arg("R"), py::arg("second_pass") = true, py::arg("want_order") = false,
				py::arg("sequential") = false, py::arg("num_threads") = 0, py::arg("device") = "cpu",
				"Shivashankar-Natarajan ASSIGNGRADIENT (+ ASSIGNGRADIENT2 if second_pass) on R[x, y, z],\n"
				"values pairwise distinct. Returns (gradient, info) with the guard and claim counts.\n"
				"stage() is 1/2 for cells paired by Algorithm 1/2 and 0 for critical cells;\n"
				"with want_order=True, order() is each cell's position in the (dimension, key)\n"
				"order. sequential=True runs the literal ordered loop instead of the per-dimension\n"
				"parallel passes (same result; CPU only). device as for pls().")
			.def_static("from_pair_dict", [](std::vector<int64_t> dims, py::dict pairs) {
			if (dims.size() != 3) throw py::value_error("dims must be (nx, ny, nz)");
			auto g = std::make_shared<PyGradient>();
			g->G = Cubical3D::Make(dims[0], dims[1], dims[2]);
			g->bytes.assign((size_t)g->G.NumCells(), ::dmt::MakeGradByte(::dmt::kPairCritical));
			auto V = g->View();
			for (auto kv : pairs) {
				const int64_t a = CellFromTuple(g->G, kv.first), b = CellFromTuple(g->G, kv.second);
				if (!V.SetPartnerOneSided(a, b)) throw py::value_error("paired cells must be Hasse-adjacent");
			}
			if (GInt::dmt::CheckConsistency(V) != 0)
				throw py::value_error("pairs do not form a matching (every pair must appear in both directions)");
			return g;
				}, py::arg("dims"), py::arg("pairs"),
				"Gradient from a dict holding each pair in both directions; every other cell\n"
				"is critical.")
			.def_property_readonly("dims", [](const PyGradient& g) { return py::make_tuple(g.G.n[0], g.G.n[1], g.G.n[2]); })
			.def_property_readonly("cell_dims", [](const PyGradient& g) { return py::make_tuple(g.G.c[0], g.G.c[1], g.G.c[2]); })
			.def_property_readonly("version", [](const PyGradient& g) { return g.version; },
				"Incremented by every change to the pairing.")
			.def("bytes", [](py::object self) {
			PyGradient& g = self.cast<PyGradient&>();
			return CellArrayView<uint8_t>(g.G, g.bytes.data(), self);
				}, "Zero-copy (Cx, Cy, Cz) view of the gradient bytes, indexed [cx, cy, cz].")
			.def("step", [](py::object self) -> py::object {
			PyGradient& g = self.cast<PyGradient&>();
			if (g.step.empty()) return py::none();
			return CellArrayView<uint8_t>(g.G, g.step.data(), self);
				}, "Per-cell decision index within its lower star (journal=True only); 255 = orphan.")
			.def("part", [](py::object self) -> py::object {
			PyGradient& g = self.cast<PyGradient&>();
			if (g.part.empty()) return py::none();
			return CellArrayView<uint8_t>(g.G, g.part.data(), self);
				}, "Per-cell part of its lower star (journal=True only): 0 boundary/whole, 1 interior.")
			.def("stage", [](py::object self) -> py::object {
			PyGradient& g = self.cast<PyGradient&>();
			if (g.stage.empty()) return py::none();
			return CellArrayView<uint8_t>(g.G, g.stage.data(), self);
				}, "DMF only: 1/2 for cells paired by Algorithm 1/2, 0 for critical cells.")
			.def("order", [](py::object self) -> py::object {
			PyGradient& g = self.cast<PyGradient&>();
			if (g.order.empty()) return py::none();
			return CellArrayView<int64_t>(g.G, g.order.data(), self);
				}, "DMF with want_order=True: position of each cell in the (dimension, key) order.")
			.def("partner", [](const PyGradient& g) {
			py::array_t<int64_t> out({ (py::ssize_t)g.G.c[0], (py::ssize_t)g.G.c[1], (py::ssize_t)g.G.c[2] });
			auto w = out.mutable_unchecked<3>();
			const auto V = g.View();
			for (int64_t z = 0; z < g.G.c[2]; z++)
				for (int64_t y = 0; y < g.G.c[1]; y++)
					for (int64_t x = 0; x < g.G.c[0]; x++)
						w(x, y, z) = V.Partner(g.G.Id(x, y, z));
			return out;
				}, "(Cx, Cy, Cz) partner ids, -1 for critical cells.")
			.def("critical", [](const PyGradient& g, int dim) {
			auto ids = g.View().Critical(dim);
			SortTupleOrder(g.G, ids);
			return CoordsArray(g.G, ids);
				}, py::arg("dim") = -1, "(N, 3) critical cells (of one dimension, or all), in tuple order.")
			.def("critical_counts", [](const PyGradient& g) {
			int64_t c[4] = { 0, 0, 0, 0 };
			const auto V = g.View();
			for (int64_t i = 0; i < g.G.NumCells(); i++)
				if (V.IsCritical(i)) c[g.G.Dim(i)]++;
			return py::make_tuple(c[0], c[1], c[2], c[3]);
				})
			.def("critical_set", [](const PyGradient& g) {
			py::set s;
			const auto V = g.View();
			for (int64_t i = 0; i < g.G.NumCells(); i++)
				if (V.IsCritical(i)) s.add(TupleFromCell(g.G, i));
			return s;
				}, "Python set of critical cells (tuples).")
			.def("pair_dict", [](const PyGradient& g) {
			py::dict d;
			const auto V = g.View();
			for (int64_t i = 0; i < g.G.NumCells(); i++) {
				const int64_t p = V.Partner(i);
				if (p >= 0) d[TupleFromCell(g.G, i)] = TupleFromCell(g.G, p);
			}
			return d;
				}, "Python dict {cell: partner} holding each pair in both directions.")
			.def("pair_map", [](std::shared_ptr<PyGradient> g) { return PairMap{ g }; },
				"Live tuple-keyed MutableMapping over the pairing (no copy).")
			.def("is_acyclic", [](const PyGradient& g) {
			py::gil_scoped_release release;
			return GInt::dmt::IsAcyclic(g.View());
				})
			.def("inconsistent_cells", [](const PyGradient& g) { return GInt::dmt::CheckConsistency(g.View()); },
				"Number of cells whose pairing is not a valid matching (0 = valid).")
			.def("copy", [](const PyGradient& g) { return std::make_shared<PyGradient>(g); })
			.def("__eq__", [](const PyGradient& a, const PyGradient& b) {
			return a.G.n[0] == b.G.n[0] && a.G.n[1] == b.G.n[1] && a.G.n[2] == b.G.n[2] && a.bytes == b.bytes;
				});

		py::class_<PairMap>(m, "PairMap",
			"Tuple-keyed view of a Gradient3D's pairing, like a dict holding each pair in both\n"
			"directions. Critical cells are absent. m[a] = b writes a's side only; del m[a]\n"
			"makes a critical.")
			.def("__getitem__", [](const PairMap& pm, py::handle k) {
			py::object r = pm.get(k, py::none());
			if (r.is_none()) throw py::key_error(py::str(k).cast<std::string>());
			return r;
				})
			.def("get", &PairMap::get, py::arg("key"), py::arg("default") = py::none())
			.def("__contains__", [](const PairMap& pm, py::handle k) {
			py::sequence s = py::reinterpret_borrow<py::sequence>(k);
			if (py::len(s) != 3) return false;
			const int64_t x = s[0].cast<int64_t>(), y = s[1].cast<int64_t>(), z = s[2].cast<int64_t>();
			if (!pm.g->G.InDomain(x, y, z)) return false;
			return !pm.g->View().IsCritical(pm.g->G.Id(x, y, z));
				})
			.def("__setitem__", [](PairMap& pm, py::handle k, py::handle v) {
			const int64_t a = CellFromTuple(pm.g->G, k), b = CellFromTuple(pm.g->G, v);
			if (!pm.g->View().SetPartnerOneSided(a, b)) throw py::value_error("paired cells must be Hasse-adjacent");
			pm.g->version++;
				})
			.def("__delitem__", [](PairMap& pm, py::handle k) {
			const int64_t a = CellFromTuple(pm.g->G, k);
			auto V = pm.g->View();
			if (V.IsCritical(a)) throw py::key_error(py::str(k).cast<std::string>());
			V.SetCritical(a);
			pm.g->version++;
				})
			.def("__len__", [](const PairMap& pm) {
			int64_t n = 0;
			const auto V = pm.g->View();
			for (int64_t i = 0; i < pm.g->G.NumCells(); i++)
				if (!V.IsCritical(i)) n++;
			return n;
				})
			.def("__iter__", [](const PairMap& pm) {
			py::list keys;
			const auto V = pm.g->View();
			const auto& G = pm.g->G;
			for (int64_t x = 0; x < G.c[0]; x++)
				for (int64_t y = 0; y < G.c[1]; y++)
					for (int64_t z = 0; z < G.c[2]; z++)
						if (!V.IsCritical(G.Id(x, y, z))) keys.append(py::make_tuple(x, y, z));
			return keys.attr("__iter__")();
				}, "Keys in tuple order (x slowest).")
			.def_property_readonly("gradient", [](const PairMap& pm) { return pm.g; });
	}

} // namespace dmt_py
