#ifndef DMT_PY_BIND_COMMON_H
#define DMT_PY_BIND_COMMON_H

// Shared helpers of the msceer_dmt bindings.
//
// Python conventions:
//   - a field is a numpy array R with shape (nx, ny, nz), indexed R[x, y, z];
//   - a cell is a tuple (cx, cy, cz) in doubled coordinates;
//   - sorted outputs use tuple order (cx, cy, cz) lexicographic, x slowest,
//     which is the order of itertools.product over the axes.
// C++ uses MSCEER ids (x fastest) internally.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "gi_basic_types.h"
#include "gi_dmt_gradient_view.h"
#include "gi_dmt_rank_field.h"
#include "dmt_core/dmt_cubical3d.h"

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

namespace dmt_py {

	namespace py = pybind11;
	using ::dmt::Cubical3D;

	// Values of R[x, y, z] in vertex-id order (x fastest).
	template <class T>
	std::vector<T> ReadField(const py::array& R, int64_t dims[3]) {
		if (R.ndim() != 3) throw py::value_error("field must be a 3-dimensional array R[x, y, z]");
		auto a = py::array_t<T, py::array::forcecast>(R);
		auto u = a.template unchecked<3>();
		for (int i = 0; i < 3; i++) dims[i] = (int64_t)u.shape(i);
		std::vector<T> v((size_t)(dims[0] * dims[1] * dims[2]));
		for (int64_t z = 0; z < dims[2]; z++)
			for (int64_t y = 0; y < dims[1]; y++)
				for (int64_t x = 0; x < dims[0]; x++)
					v[(size_t)(x + dims[0] * (y + dims[1] * z))] = u(x, y, z);
		return v;
	}

	inline bool IsIntegerArray(const py::array& R) {
		const char k = R.dtype().kind();
		return k == 'i' || k == 'u' || k == 'b';
	}

	// Position of each vertex (x fastest id) in numpy C order of R[x, y, z],
	// i.e. (x * ny + y) * nz + z. Used as the 'c_order' tie-break.
	inline std::vector<int64_t> COrderPositions(const int64_t dims[3]) {
		std::vector<int64_t> p((size_t)(dims[0] * dims[1] * dims[2]));
		for (int64_t z = 0; z < dims[2]; z++)
			for (int64_t y = 0; y < dims[1]; y++)
				for (int64_t x = 0; x < dims[0]; x++)
					p[(size_t)(x + dims[0] * (y + dims[1] * z))] = (x * dims[1] + y) * dims[2] + z;
		return p;
	}

	inline std::string VertexStr(const int64_t dims[3], int64_t vid) {
		const int64_t x = vid % dims[0], y = (vid / dims[0]) % dims[1], z = vid / (dims[0] * dims[1]);
		return "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
	}

	// Exact ranks of a field whose values must be pairwise distinct.
	inline GInt::dmt::RankField3D DistinctRanks(const py::array& R) {
		int64_t dims[3];
		GInt::dmt::RankField3D f;
		int64_t a = -1, b = -1;
		bool dup;
		if (IsIntegerArray(R)) {
			auto v = ReadField<int64_t>(R, dims);
			dup = GInt::dmt::RankField3D::FirstDuplicate(v.data(), (int64_t)v.size(), a, b);
			if (!dup) f = GInt::dmt::RankField3D::FromValues(v.data(), dims[0], dims[1], dims[2]);
		}
		else {
			auto v = ReadField<double>(R, dims);
			dup = GInt::dmt::RankField3D::FirstDuplicate(v.data(), (int64_t)v.size(), a, b);
			if (!dup) f = GInt::dmt::RankField3D::FromValues(v.data(), dims[0], dims[1], dims[2]);
		}
		if (dup)
			throw py::value_error("field values must be pairwise distinct; vertices " + VertexStr(dims, a) +
				" and " + VertexStr(dims, b) + " are equal (use msceer_dmt.ranks(R, tie=...) to break ties)");
		if (dims[0] < 1 || dims[1] < 1 || dims[2] < 1) throw py::value_error("empty field");
		if (dims[0] * dims[1] * dims[2] > (int64_t)0xFFFFFFFFll) throw py::value_error("too many vertices for uint32 ranks");
		return f;
	}

	inline int64_t CellFromTuple(const Cubical3D& G, const py::handle& h) {
		py::sequence s = py::reinterpret_borrow<py::sequence>(h);
		if (py::len(s) != 3) throw py::value_error("a cell is a tuple (cx, cy, cz)");
		const int64_t x = s[0].cast<int64_t>(), y = s[1].cast<int64_t>(), z = s[2].cast<int64_t>();
		if (!G.InDomain(x, y, z)) throw py::key_error("cell outside the grid");
		return G.Id(x, y, z);
	}

	inline py::tuple TupleFromCell(const Cubical3D& G, int64_t id) {
		int64_t q[3];
		G.Coords(id, q);
		return py::make_tuple(q[0], q[1], q[2]);
	}

	// Tuple order (x slowest) of cell ids.
	inline bool TupleLess(const Cubical3D& G, int64_t a, int64_t b) {
		int64_t qa[3], qb[3];
		G.Coords(a, qa);
		G.Coords(b, qb);
		if (qa[0] != qb[0]) return qa[0] < qb[0];
		if (qa[1] != qb[1]) return qa[1] < qb[1];
		return qa[2] < qb[2];
	}

	inline void SortTupleOrder(const Cubical3D& G, std::vector<int64_t>& ids) {
		std::sort(ids.begin(), ids.end(), [&](int64_t a, int64_t b) { return TupleLess(G, a, b); });
	}

	inline py::array_t<int32_t> CoordsArray(const Cubical3D& G, const std::vector<int64_t>& ids) {
		py::array_t<int32_t> out({ (py::ssize_t)ids.size(), (py::ssize_t)3 });
		auto w = out.mutable_unchecked<2>();
		for (size_t i = 0; i < ids.size(); i++) {
			int64_t q[3];
			G.Coords(ids[i], q);
			for (int k = 0; k < 3; k++) w((py::ssize_t)i, k) = (int32_t)q[k];
		}
		return out;
	}

	// Ids from an (N, 3) array of doubled coordinates.
	inline std::vector<int64_t> IdsFromCoords(const Cubical3D& G, const py::array& coords) {
		auto a = py::array_t<int64_t, py::array::forcecast>(coords);
		if (a.ndim() != 2 || a.shape(1) != 3) throw py::value_error("coordinates must have shape (N, 3)");
		auto u = a.unchecked<2>();
		std::vector<int64_t> ids((size_t)u.shape(0));
		for (py::ssize_t i = 0; i < u.shape(0); i++) {
			if (!G.InDomain(u(i, 0), u(i, 1), u(i, 2))) throw py::value_error("cell outside the grid");
			ids[(size_t)i] = G.Id(u(i, 0), u(i, 1), u(i, 2));
		}
		return ids;
	}

	// A per-cell array shaped (Cx, Cy, Cz) so that a[cx, cy, cz] is the cell,
	// viewing x-fastest storage without a copy.
	template <class T>
	py::array_t<T> CellArrayView(const Cubical3D& G, T* data, py::handle base) {
		std::vector<py::ssize_t> shape = { (py::ssize_t)G.c[0], (py::ssize_t)G.c[1], (py::ssize_t)G.c[2] };
		std::vector<py::ssize_t> strides = { (py::ssize_t)sizeof(T), (py::ssize_t)(sizeof(T) * G.c[0]),
			(py::ssize_t)(sizeof(T) * G.c[0] * G.c[1]) };
		return py::array_t<T>(shape, strides, data, base);
	}

	void BindGrid(py::module_& m);
	void BindGradient(py::module_& m);
	void BindPaths(py::module_& m);

} // namespace dmt_py

#endif
