/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#include "base/gi_topological_regular_grid_3x3x3.h"


namespace GInt {

	const int Explicit3x3x3SmallRegularGrid::m_adjacent_cell_offsets[27] = {
		-31, -30, -29, -26, -25, -24, -21, -20, -19, -6, -5, -4, -1, 0, 1, 4, 5, 6, 19, 20, 21, 24, 25, 26, 29, 30, 31
	};


	const int Explicit3x3x3SmallRegularGrid::VertNumFromCellidArray[125] = {
		0, -1, 1, -1, 2, -1, -1, -1, -1, -1, 3, -1, 4, -1, 5, -1, -1,
		-1, -1, -1, 6, -1, 7, -1, 8, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, 9, -1, 10, -1, 11, -1, -1, -1, -1, -1, 12, -1, 13, -1, 14,
		-1, -1, -1, -1, -1, 15, -1, 16, -1, 17, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, 18, -1, 19, -1, 20, -1, -1, -1, -1, -1, 21, -1, 22,
		-1, 23, -1, -1, -1, -1, -1, 24, -1, 25, -1, 26
	};

	const int Explicit3x3x3SmallRegularGrid::DimArray[125] = {
		0, 1, 0, 1, 0, 1, 2, 1, 2, 1, 0, 1, 0, 1, 0, 1, 2, 1, 2, 1,
		0, 1, 0, 1, 0, 1, 2, 1, 2, 1, 2, 3, 2, 3, 2, 1, 2, 1, 2, 1,
		2, 3, 2, 3, 2, 1, 2, 1, 2, 1, 0, 1, 0, 1, 0, 1, 2, 1, 2, 1,
		0, 1, 0, 1, 0, 1, 2, 1, 2, 1, 0, 1, 0, 1, 0, 1, 2, 1, 2, 1,
		2, 3, 2, 3, 2, 1, 2, 1, 2, 1, 2, 3, 2, 3, 2, 1, 2, 1, 2, 1,
		0, 1, 0, 1, 0, 1, 2, 1, 2, 1, 0, 1, 0, 1, 0, 1, 2, 1, 2, 1,
		0, 1, 0, 1, 0
	};

	const int Explicit3x3x3SmallRegularGrid::FacetsIterator::FacetsIdList[125][7] =
	{
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 2, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 4, 2, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 10, 0, 0, 0, 0, 0 },
		{ 4, 7, 5, 11, 1, 0, 0 },
		{ 2, 12, 2, 0, 0, 0, 0 },
		{ 4, 9, 7, 13, 3, 0, 0 },
		{ 2, 14, 4, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 12, 10, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 14, 12, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 20, 10, 0, 0, 0, 0 },
		{ 4, 17, 15, 21, 11, 0, 0 },
		{ 2, 22, 12, 0, 0, 0, 0 },
		{ 4, 19, 17, 23, 13, 0, 0 },
		{ 2, 24, 14, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 22, 20, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 24, 22, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 50, 0, 0, 0, 0, 0 },
		{ 4, 27, 25, 51, 1, 0, 0 },
		{ 2, 52, 2, 0, 0, 0, 0 },
		{ 4, 29, 27, 53, 3, 0, 0 },
		{ 2, 54, 4, 0, 0, 0, 0 },
		{ 4, 35, 25, 55, 5, 0, 0 },
		{ 6, 32, 30, 36, 26, 56, 6 },
		{ 4, 37, 27, 57, 7, 0, 0 },
		{ 6, 34, 32, 38, 28, 58, 8 },
		{ 4, 39, 29, 59, 9, 0, 0 },
		{ 2, 60, 10, 0, 0, 0, 0 },
		{ 4, 37, 35, 61, 11, 0, 0 },
		{ 2, 62, 12, 0, 0, 0, 0 },
		{ 4, 39, 37, 63, 13, 0, 0 },
		{ 2, 64, 14, 0, 0, 0, 0 },
		{ 4, 45, 35, 65, 15, 0, 0 },
		{ 6, 42, 40, 46, 36, 66, 16 },
		{ 4, 47, 37, 67, 17, 0, 0 },
		{ 6, 44, 42, 48, 38, 68, 18 },
		{ 4, 49, 39, 69, 19, 0, 0 },
		{ 2, 70, 20, 0, 0, 0, 0 },
		{ 4, 47, 45, 71, 21, 0, 0 },
		{ 2, 72, 22, 0, 0, 0, 0 },
		{ 4, 49, 47, 73, 23, 0, 0 },
		{ 2, 74, 24, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 52, 50, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 54, 52, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 60, 50, 0, 0, 0, 0 },
		{ 4, 57, 55, 61, 51, 0, 0 },
		{ 2, 62, 52, 0, 0, 0, 0 },
		{ 4, 59, 57, 63, 53, 0, 0 },
		{ 2, 64, 54, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 62, 60, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 64, 62, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 70, 60, 0, 0, 0, 0 },
		{ 4, 67, 65, 71, 61, 0, 0 },
		{ 2, 72, 62, 0, 0, 0, 0 },
		{ 4, 69, 67, 73, 63, 0, 0 },
		{ 2, 74, 64, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 72, 70, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 74, 72, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 100, 50, 0, 0, 0, 0 },
		{ 4, 77, 75, 101, 51, 0, 0 },
		{ 2, 102, 52, 0, 0, 0, 0 },
		{ 4, 79, 77, 103, 53, 0, 0 },
		{ 2, 104, 54, 0, 0, 0, 0 },
		{ 4, 85, 75, 105, 55, 0, 0 },
		{ 6, 82, 80, 86, 76, 106, 56 },
		{ 4, 87, 77, 107, 57, 0, 0 },
		{ 6, 84, 82, 88, 78, 108, 58 },
		{ 4, 89, 79, 109, 59, 0, 0 },
		{ 2, 110, 60, 0, 0, 0, 0 },
		{ 4, 87, 85, 111, 61, 0, 0 },
		{ 2, 112, 62, 0, 0, 0, 0 },
		{ 4, 89, 87, 113, 63, 0, 0 },
		{ 2, 114, 64, 0, 0, 0, 0 },
		{ 4, 95, 85, 115, 65, 0, 0 },
		{ 6, 92, 90, 96, 86, 116, 66 },
		{ 4, 97, 87, 117, 67, 0, 0 },
		{ 6, 94, 92, 98, 88, 118, 68 },
		{ 4, 99, 89, 119, 69, 0, 0 },
		{ 2, 120, 70, 0, 0, 0, 0 },
		{ 4, 97, 95, 121, 71, 0, 0 },
		{ 2, 122, 72, 0, 0, 0, 0 },
		{ 4, 99, 97, 123, 73, 0, 0 },
		{ 2, 124, 74, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 102, 100, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 104, 102, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 110, 100, 0, 0, 0, 0 },
		{ 4, 107, 105, 111, 101, 0, 0 },
		{ 2, 112, 102, 0, 0, 0, 0 },
		{ 4, 109, 107, 113, 103, 0, 0 },
		{ 2, 114, 104, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 112, 110, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 114, 112, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 120, 110, 0, 0, 0, 0 },
		{ 4, 117, 115, 121, 111, 0, 0 },
		{ 2, 122, 112, 0, 0, 0, 0 },
		{ 4, 119, 117, 123, 113, 0, 0 },
		{ 2, 124, 114, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 122, 120, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 124, 122, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 }
	};

	const int Explicit3x3x3SmallRegularGrid::CofacetsIterator::CoFacetsIdList[125][7] =
	{
		{ 3, 1, 5, 25, 0, 0, 0 },
		{ 2, 6, 26, 0, 0, 0, 0 },
		{ 4, 3, 1, 7, 27, 0, 0 },
		{ 2, 8, 28, 0, 0, 0, 0 },
		{ 3, 3, 9, 29, 0, 0, 0 },
		{ 2, 6, 30, 0, 0, 0, 0 },
		{ 1, 31, 0, 0, 0, 0, 0 },
		{ 3, 8, 6, 32, 0, 0, 0 },
		{ 1, 33, 0, 0, 0, 0, 0 },
		{ 2, 8, 34, 0, 0, 0, 0 },
		{ 4, 11, 15, 5, 35, 0, 0 },
		{ 3, 16, 6, 36, 0, 0, 0 },
		{ 5, 13, 11, 17, 7, 37, 0 },
		{ 3, 18, 8, 38, 0, 0, 0 },
		{ 4, 13, 19, 9, 39, 0, 0 },
		{ 2, 16, 40, 0, 0, 0, 0 },
		{ 1, 41, 0, 0, 0, 0, 0 },
		{ 3, 18, 16, 42, 0, 0, 0 },
		{ 1, 43, 0, 0, 0, 0, 0 },
		{ 2, 18, 44, 0, 0, 0, 0 },
		{ 3, 21, 15, 45, 0, 0, 0 },
		{ 2, 16, 46, 0, 0, 0, 0 },
		{ 4, 23, 21, 17, 47, 0, 0 },
		{ 2, 18, 48, 0, 0, 0, 0 },
		{ 3, 23, 19, 49, 0, 0, 0 },
		{ 2, 26, 30, 0, 0, 0, 0 },
		{ 1, 31, 0, 0, 0, 0, 0 },
		{ 3, 28, 26, 32, 0, 0, 0 },
		{ 1, 33, 0, 0, 0, 0, 0 },
		{ 2, 28, 34, 0, 0, 0, 0 },
		{ 1, 31, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 33, 31, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 1, 33, 0, 0, 0, 0, 0 },
		{ 3, 36, 40, 30, 0, 0, 0 },
		{ 2, 41, 31, 0, 0, 0, 0 },
		{ 4, 38, 36, 42, 32, 0, 0 },
		{ 2, 43, 33, 0, 0, 0, 0 },
		{ 3, 38, 44, 34, 0, 0, 0 },
		{ 1, 41, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 43, 41, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 1, 43, 0, 0, 0, 0, 0 },
		{ 2, 46, 40, 0, 0, 0, 0 },
		{ 1, 41, 0, 0, 0, 0, 0 },
		{ 3, 48, 46, 42, 0, 0, 0 },
		{ 1, 43, 0, 0, 0, 0, 0 },
		{ 2, 48, 44, 0, 0, 0, 0 },
		{ 4, 51, 55, 75, 25, 0, 0 },
		{ 3, 56, 76, 26, 0, 0, 0 },
		{ 5, 53, 51, 57, 77, 27, 0 },
		{ 3, 58, 78, 28, 0, 0, 0 },
		{ 4, 53, 59, 79, 29, 0, 0 },
		{ 3, 56, 80, 30, 0, 0, 0 },
		{ 2, 81, 31, 0, 0, 0, 0 },
		{ 4, 58, 56, 82, 32, 0, 0 },
		{ 2, 83, 33, 0, 0, 0, 0 },
		{ 3, 58, 84, 34, 0, 0, 0 },
		{ 5, 61, 65, 55, 85, 35, 0 },
		{ 4, 66, 56, 86, 36, 0, 0 },
		{ 6, 63, 61, 67, 57, 87, 37 },
		{ 4, 68, 58, 88, 38, 0, 0 },
		{ 5, 63, 69, 59, 89, 39, 0 },
		{ 3, 66, 90, 40, 0, 0, 0 },
		{ 2, 91, 41, 0, 0, 0, 0 },
		{ 4, 68, 66, 92, 42, 0, 0 },
		{ 2, 93, 43, 0, 0, 0, 0 },
		{ 3, 68, 94, 44, 0, 0, 0 },
		{ 4, 71, 65, 95, 45, 0, 0 },
		{ 3, 66, 96, 46, 0, 0, 0 },
		{ 5, 73, 71, 67, 97, 47, 0 },
		{ 3, 68, 98, 48, 0, 0, 0 },
		{ 4, 73, 69, 99, 49, 0, 0 },
		{ 2, 76, 80, 0, 0, 0, 0 },
		{ 1, 81, 0, 0, 0, 0, 0 },
		{ 3, 78, 76, 82, 0, 0, 0 },
		{ 1, 83, 0, 0, 0, 0, 0 },
		{ 2, 78, 84, 0, 0, 0, 0 },
		{ 1, 81, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 83, 81, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 1, 83, 0, 0, 0, 0, 0 },
		{ 3, 86, 90, 80, 0, 0, 0 },
		{ 2, 91, 81, 0, 0, 0, 0 },
		{ 4, 88, 86, 92, 82, 0, 0 },
		{ 2, 93, 83, 0, 0, 0, 0 },
		{ 3, 88, 94, 84, 0, 0, 0 },
		{ 1, 91, 0, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 2, 93, 91, 0, 0, 0, 0 },
		{ 0, 0, 0, 0, 0, 0, 0 },
		{ 1, 93, 0, 0, 0, 0, 0 },
		{ 2, 96, 90, 0, 0, 0, 0 },
		{ 1, 91, 0, 0, 0, 0, 0 },
		{ 3, 98, 96, 92, 0, 0, 0 },
		{ 1, 93, 0, 0, 0, 0, 0 },
		{ 2, 98, 94, 0, 0, 0, 0 },
		{ 3, 101, 105, 75, 0, 0, 0 },
		{ 2, 106, 76, 0, 0, 0, 0 },
		{ 4, 103, 101, 107, 77, 0, 0 },
		{ 2, 108, 78, 0, 0, 0, 0 },
		{ 3, 103, 109, 79, 0, 0, 0 },
		{ 2, 106, 80, 0, 0, 0, 0 },
		{ 1, 81, 0, 0, 0, 0, 0 },
		{ 3, 108, 106, 82, 0, 0, 0 },
		{ 1, 83, 0, 0, 0, 0, 0 },
		{ 2, 108, 84, 0, 0, 0, 0 },
		{ 4, 111, 115, 105, 85, 0, 0 },
		{ 3, 116, 106, 86, 0, 0, 0 },
		{ 5, 113, 111, 117, 107, 87, 0 },
		{ 3, 118, 108, 88, 0, 0, 0 },
		{ 4, 113, 119, 109, 89, 0, 0 },
		{ 2, 116, 90, 0, 0, 0, 0 },
		{ 1, 91, 0, 0, 0, 0, 0 },
		{ 3, 118, 116, 92, 0, 0, 0 },
		{ 1, 93, 0, 0, 0, 0, 0 },
		{ 2, 118, 94, 0, 0, 0, 0 },
		{ 3, 121, 115, 95, 0, 0, 0 },
		{ 2, 116, 96, 0, 0, 0, 0 },
		{ 4, 123, 121, 117, 97, 0, 0 },
		{ 2, 118, 98, 0, 0, 0, 0 },
		{ 3, 123, 119, 99, 0, 0, 0 }
	};

}