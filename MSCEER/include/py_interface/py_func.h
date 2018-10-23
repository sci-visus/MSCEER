/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef PY_FUNC_H
#define PY_FUNC_H

#include <vector>
#include <set>
#include <queue>
#include <time.h>

#include "base/gi_basic_types.h"
#include "base/gi_labeling.h"
#include "base/gi_discrete_gradient_labeling.h"
#include "base/gi_topological_regular_grid_2d.h"
#include "base/gi_regular_grid_2d.h"

#include "py_interface/py_basic_types.h"

// wrapper class for meshes
class Func {

public:
	virtual Coord SampleGrad(Coord p) { return{ 0, 0 }; }
	virtual float SampleImage(Coord p) { return 0; }
	virtual float MaxVal() { return 0; }
	virtual float MinVal() { return 0; }
};



#endif