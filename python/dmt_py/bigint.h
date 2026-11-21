#ifndef DMT_PY_BIGINT_H
#define DMT_PY_BIGINT_H

// U128 / S128 path counts as exact Python ints (public C API only).

#include <pybind11/pybind11.h>

#include "dmt_core/dmt_count128.h"

namespace dmt_py {

	namespace py = pybind11;

	inline py::int_ ToPy(const ::dmt::U128& u) {
		if (u.overflow) throw py::value_error("path count overflowed 128 bits");
		if (u.hi == 0) return py::int_((unsigned long long)u.lo);
		py::object hi = py::reinterpret_steal<py::object>(PyLong_FromUnsignedLongLong(u.hi));
		py::object lo = py::reinterpret_steal<py::object>(PyLong_FromUnsignedLongLong(u.lo));
		py::object sh = py::reinterpret_steal<py::object>(PyLong_FromLong(64));
		py::object t = py::reinterpret_steal<py::object>(PyNumber_Lshift(hi.ptr(), sh.ptr()));
		if (!t) throw py::error_already_set();
		py::object r = py::reinterpret_steal<py::object>(PyNumber_Or(t.ptr(), lo.ptr()));
		if (!r) throw py::error_already_set();
		return py::reinterpret_borrow<py::int_>(r);
	}

	inline py::int_ ToPy(const ::dmt::S128& s) {
		if (s.overflow) throw py::value_error("signed path count overflowed 128 bits");
		if (s.FitsI64()) return py::int_((long long)s.ToI64());
		py::int_ m = ToPy(::dmt::Abs(s));
		if (!s.Negative()) return m;
		py::object r = py::reinterpret_steal<py::object>(PyNumber_Negative(m.ptr()));
		if (!r) throw py::error_already_set();
		return py::reinterpret_borrow<py::int_>(r);
	}

} // namespace dmt_py

#endif
