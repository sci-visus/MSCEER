/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#include "base/gi_regular_grid_2d.h"


namespace GInt {

	const Vec2l RegularGrid2D::kNeighborOffsets4[4] = {
		Vec2l(1, 0), Vec2l(-1, 0),
		Vec2l(0, 1), Vec2l(0, -1)
	};

	const Vec2l RegularGrid2D::kNeighborOffsets8[8] = {
		Vec2l(1, 1), Vec2l(0, 1), Vec2l(-1, 1),
		Vec2l(1, 0), Vec2l(-1, 0),
		Vec2l(1, -1), Vec2l(0, -1), Vec2l(-1, -1),
	};
}