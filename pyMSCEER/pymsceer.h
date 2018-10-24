/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/


#ifndef PYMSCEER_H
#define PYMSCEER_H


#include "base/gi_ondemand_accurate_grad_builder.h"
#include "py_grad_builder.h"
#include "py_mesh.h"
#include "py_discrete_grad.h"

class TopologyContext {
public:
	Mesh* mesh;
	DiscreteGrad* dgrad;
};

TopologyContext* makeContext(int x, int y, int z);

DiscreteGradientBuilder* makeDgrad();

void makegrad(char* fname, int x, int y, int z);

#endif