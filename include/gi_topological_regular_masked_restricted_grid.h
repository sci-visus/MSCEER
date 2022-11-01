/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef TOPOLOGICAL_REGULAR_MASKED_RESTRICTED_GRID_H
#define TOPOLOGICAL_REGULAR_MASKED_RESTRICTED_GRID_H


#include <map>
#include <set>
#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_regular_grid.h"
#include "gi_topological_regular_grid.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_topological_regular_masked_grid.h"
#include "gi_labeling.h"

namespace GInt {


    class TopologicalRegularMaskedRestrictedGrid : public TopologicalRegularMaskedGrid, public TopologicalRegularGridRestricted {
    public:
        TopologicalRegularMaskedRestrictedGrid(RegularGrid3D* base_grid) :
            TopologicalRegularMaskedGrid(base_grid),
            TopologicalRegularGridRestricted(base_grid),
            TopologicalRegularGrid3D(base_grid) {}
    };


}


#endif
