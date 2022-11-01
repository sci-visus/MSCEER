/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef EXTRACT_GRAPH_H
#define EXTRACT_GRAPH_H

#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <queue>
#include <time.h>
#include <fstream>

#include "gi_basic_types.h"
#include "gi_timing.h"
#include "gi_topological_explicit_mesh_function.h"

#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"

#include "gi_topological_utility_functions.h"
#include "gi_numeric_integrator_expanding_region_stop.h" // not where comparer should be
#include "gi_topological_regular_masked_grid.h"
#include "gi_topological_regular_masked_restricted_grid.h"

#include "gi_max_vertex_labeling.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_bifiltration_pairing.h"
#include "gi_morse_smale_complex_basic.h"

// using a regular grid
typedef GInt::RegularGrid2D GridType;
typedef GInt::TopologicalRegularGrid2D MeshType;
typedef GInt::RegularGridBilinearFunction GridFuncType;
typedef GInt::DiscreteGradientLabeling<MeshType> GradType;

// for the topology we will use the lazy version of the topology-function, so no extra
// storage of computation
typedef GInt::LazyMaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
typedef GInt::TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> TopoFuncType;

// the Morse-Smale complex type based on these prior types
typedef GInt::MorseSmaleComplexBasic<float, MeshType, TopoFuncType, GradType> MscType;

#endif
