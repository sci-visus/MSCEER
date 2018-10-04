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