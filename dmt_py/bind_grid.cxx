// Grid3D and ranks(): geometry of the cubical complex and tie-breaking.

#include "bind_common.h"

namespace dmt_py {

	struct PyGrid {
		Cubical3D G;
	};

	void BindGrid(py::module_& m) {
		py::class_<PyGrid>(m, "Grid3D",
			"Cubical complex of an nx x ny x nz vertex grid in doubled coordinates.\n"
			"Cell ids are MSCEER linear ids cx + Cx*(cy + Cy*cz).")
			.def(py::init([](std::vector<int64_t> dims) {
			if (dims.size() != 3) throw py::value_error("dims must be (nx, ny, nz)");
			for (auto d : dims) if (d < 1) throw py::value_error("dims must be positive");
			return PyGrid{ Cubical3D::Make(dims[0], dims[1], dims[2]) };
				}), py::arg("dims"))
			.def_property_readonly("dims", [](const PyGrid& g) { return py::make_tuple(g.G.n[0], g.G.n[1], g.G.n[2]); })
			.def_property_readonly("cell_dims", [](const PyGrid& g) { return py::make_tuple(g.G.c[0], g.G.c[1], g.G.c[2]); })
			.def_property_readonly("num_cells", [](const PyGrid& g) { return g.G.NumCells(); })
			.def("id", [](const PyGrid& g, py::tuple c) { return CellFromTuple(g.G, c); }, py::arg("cell"))
			.def("cell", [](const PyGrid& g, int64_t id) { return TupleFromCell(g.G, id); }, py::arg("id"))
			.def("ids", [](const PyGrid& g, py::array coords) {
			auto ids = IdsFromCoords(g.G, coords);
			return py::array_t<int64_t>((py::ssize_t)ids.size(), ids.data());
				}, py::arg("coords"), "Ids of an (N, 3) array of cells.")
			.def("coords", [](const PyGrid& g, py::array_t<int64_t, py::array::forcecast> ids) {
			auto u = ids.unchecked<1>();
			std::vector<int64_t> v((size_t)u.shape(0));
			for (py::ssize_t i = 0; i < u.shape(0); i++) v[(size_t)i] = u(i);
			return CoordsArray(g.G, v);
				}, py::arg("ids"), "(N, 3) cells of an array of ids.")
			.def("dim", [](const PyGrid& g, py::tuple c) { return g.G.Dim(CellFromTuple(g.G, c)); }, py::arg("cell"))
			.def("bnd", [](const PyGrid& g, py::tuple c) { return g.G.OnBoundary(CellFromTuple(g.G, c)); }, py::arg("cell"),
				"True if the cell lies on the domain boundary (some coordinate 0 or maximal).")
			.def("bdist", [](const PyGrid& g, py::tuple c) { return g.G.Bdist(CellFromTuple(g.G, c)); }, py::arg("cell"),
				"Steps from the domain boundary in doubled units; 0 = on it.")
			.def("incidence", [](const PyGrid& g, py::tuple tau, py::tuple sigma) {
			return g.G.Incidence(CellFromTuple(g.G, tau), CellFromTuple(g.G, sigma));
				}, py::arg("tau"), py::arg("sigma"), "Forman's cubical incidence <d tau, sigma>.");

		m.def("ranks", [](py::array R, std::string tie) {
			int64_t dims[3];
			std::vector<int64_t> tieord;
			const int64_t* tp = nullptr;
			if (tie == "c_order") {
				if (R.ndim() != 3) throw py::value_error("field must be a 3-dimensional array R[x, y, z]");
				for (int i = 0; i < 3; i++) dims[i] = (int64_t)R.shape(i);
				tieord = COrderPositions(dims);
				tp = tieord.data();
			}
			else if (tie != "msceer") {
				throw py::value_error("tie must be 'c_order' or 'msceer'");
			}
			GInt::dmt::RankField3D f;
			if (IsIntegerArray(R)) {
				auto v = ReadField<int64_t>(R, dims);
				f = GInt::dmt::RankField3D::FromValues(v.data(), dims[0], dims[1], dims[2], tp);
			}
			else {
				auto v = ReadField<double>(R, dims);
				f = GInt::dmt::RankField3D::FromValues(v.data(), dims[0], dims[1], dims[2], tp);
			}
			py::array_t<int64_t> out({ (py::ssize_t)dims[0], (py::ssize_t)dims[1], (py::ssize_t)dims[2] });
			auto w = out.mutable_unchecked<3>();
			for (int64_t z = 0; z < dims[2]; z++)
				for (int64_t y = 0; y < dims[1]; y++)
					for (int64_t x = 0; x < dims[0]; x++)
						w(x, y, z) = (int64_t)f.ranks[(size_t)(x + dims[0] * (y + dims[1] * z))];
			return out;
			}, py::arg("R"), py::arg("tie") = "c_order",
			"Ranks 0..N-1 of R[x, y, z], ties broken by numpy C order of R ('c_order',\n"
			"a stable argsort of R.ravel()) or by MSCEER vertex id, x fastest ('msceer').");
	}

} // namespace dmt_py
