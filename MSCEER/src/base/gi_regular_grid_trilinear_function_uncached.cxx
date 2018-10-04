#include "base/gi_regular_grid_trilinear_function_uncached.h"

namespace GInt {

	const FLOATTYPE UncachedRegularGridTrilinearFunction::kRKCoefficients[5][9] = {
			{ 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			{ -1.0 / 2.0, 0, 1.0 / 2.0, 0, 0, 0, 0, 0, 0 },
			{ 1.0 / 12.0, -2.0 / 3.0, 0, 2.0 / 3.0, -1.0 / 12.0, 0, 0, 0, 0 },
			{ -1.0 / 60.0, 3.0 / 20.0, -3.0 / 4.0, 0, 3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0, 0, 0 },
			{ 1.0 / 280.0, -4.0 / 105.0, 1.0 / 5.0, -4.0 / 5.0, 0, 4.0 / 5.0, -1.0 / 5.0, 4.0 / 105.0, -1.0 / 280.0 }
	};

};