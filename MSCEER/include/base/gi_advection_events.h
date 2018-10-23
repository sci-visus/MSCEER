/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef ADVECTION_EVENTS_H
#define ADVECTION_EVENTS_H

namespace GInt {
	enum ADVECTION_EVENT {
		NONE, 
		OUT_OF_VOXEL,
		OVER_MAX_ITERATIONS, 
		LOW_GRADIENT, 
		HIT_EXTREMUM, 
		HIT_PREASSIGNED,
		OTHER
	};
}
#endif