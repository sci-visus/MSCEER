#include "base/gi_regular_grid_3d.h"


namespace GInt {
	const Vec3l RegularGrid3D::kNeighborOffsets6[6] = {
		Vec3l(1, 0, 0), Vec3l(-1, 0, 0),
		Vec3l(0, 1, 0), Vec3l(0, -1, 0),
		Vec3l(0, 0, 1), Vec3l(0, 0, -1)
	};

	const Vec3l RegularGrid3D::kNeighborOffsets26[26] = {
		Vec3l(1, 1, 1), Vec3l(0, 1, 1), Vec3l(-1, 1, 1),
		Vec3l(1, 0, 1), Vec3l(0, 0, 1), Vec3l(-1, 0, 1),
		Vec3l(1, -1, 1), Vec3l(0, -1, 1), Vec3l(-1, -1, 1),

		Vec3l(1, 1, 0), Vec3l(0, 1, 0), Vec3l(-1, 1, 0),
		Vec3l(1, 0, 0), Vec3l(-1, 0, 0),
		Vec3l(1, -1, 0), Vec3l(0, -1, 0), Vec3l(-1, -1, 0),

		Vec3l(1, 1, -1), Vec3l(0, 1, -1), Vec3l(-1, 1, -1),
		Vec3l(1, 0, -1), Vec3l(0, 0, -1), Vec3l(-1, 0, -1),
		Vec3l(1, -1, -1), Vec3l(0, -1, -1), Vec3l(-1, -1, -1),

	};
}