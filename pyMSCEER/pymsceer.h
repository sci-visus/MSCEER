
#ifndef PYMSCEER_H
#define PYMSCEER_H


#include "base/gi_ondemand_accurate_grad_builder.h"
#include "py_interface\py_grad_builder.h"
#include "py_interface/py_mesh.h"
#include "py_interface/py_discrete_grad.h"

class TopologyContext {
public:
	Mesh* mesh;
	DiscreteGrad* dgrad;
};

TopologyContext* makeContext(int x, int y, int z);

DiscreteGradientBuilder* makeDgrad();

void makegrad(char* fname, int x, int y, int z);

#endif