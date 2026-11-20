// msceer_dmt: Python module for the discrete Morse theory kernels.
//
// Bindings are split by area across bind_*.cxx; this file defines the module,
// its version, and wires the pieces together.

#include <pybind11/pybind11.h>

#include "bind_common.h"
#include "dmt_core/dmt_hd.h"

namespace py = pybind11;

static const char* kVersion = "0.5.0";

static bool selftest() {
	return dmt::GradPairCode(dmt::MakeGradByte(3)) == 3 &&
		dmt::GradCritical(dmt::MakeGradByte(dmt::kPairCritical));
}

PYBIND11_MODULE(msceer_dmt, m) {
	m.doc() = "Discrete Morse theory kernels from MSCEER: lower-star gradient rules, "
		"V-path counting, the signed Morse boundary and alternating-cycle search.";
	m.def("version", []() { return std::string(kVersion); },
		"Version string of the msceer_dmt bindings.");
	m.def("selftest", &selftest,
		"Round-trips the gradient byte encoding; True on success.");

	dmt_py::BindGrid(m);
	dmt_py::BindGradient(m);
	dmt_py::BindPaths(m);
	dmt_py::BindCycles(m);

	// PairMap behaves as a MutableMapping (keys, items, values, ==, ...).
	py::module_::import("collections.abc").attr("MutableMapping").attr("register")(m.attr("PairMap"));
}
