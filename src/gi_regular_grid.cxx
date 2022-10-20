#include "gi_regular_grid.h"


using namespace GInt;
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
	Vec3l(1, 0, 0),                 Vec3l(-1, 0, 0),
	Vec3l(1, -1, 0), Vec3l(0, -1, 0), Vec3l(-1, -1, 0),

	Vec3l(1, 1, -1), Vec3l(0, 1, -1), Vec3l(-1, 1, -1),
	Vec3l(1, 0, -1), Vec3l(0, 0, -1), Vec3l(-1, 0, -1),
	Vec3l(1, -1, -1), Vec3l(0, -1, -1), Vec3l(-1, -1, -1),

};

const Vec2l RegularGrid2D::kNeighborOffsets4[4] = {
	Vec2l(1, 0), Vec2l(-1, 0),
	Vec2l(0, 1), Vec2l(0, -1)
};

const Vec2l RegularGrid2D::kNeighborOffsets8[8] = {
	Vec2l(1, 1), Vec2l(0, 1), Vec2l(-1, 1),
	Vec2l(1, 0), Vec2l(-1, 0),
	Vec2l(1, -1), Vec2l(0, -1), Vec2l(-1, -1),
};