/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/


 #include <iostream>
 #include "stdio.h"
 #include "stdlib.h"
#include "pymsceer.h"

 #include <cmath>

TopologyContext* makeContext(int x, int y, int z) {
	TopologyContext* t = new TopologyContext();
	if (z == 1) {
		t->mesh = new RegularGrid2D(x, y);
		t->dgrad = new GridDiscreteGrad2D();
	}
	else {
		t->mesh = new RegularGrid3D(x, y, z);
		t->dgrad = new GridDiscreteGrad3D();
	}
	return t;

}
DiscreteGradientBuilder* makeDgrad() {
	return new DiscreteGradientBuilder();
}

void makegrad(char* fname, int x, int y, int z) {

	

	GInt::OndemandDiscreteGradientBuilder* test = new GInt::OndemandDiscreteGradientBuilder();
	test->BuildGradient(x, y, z, fname, 1, 1, 1, 1, 0.0);


}


int main( int argc, char** argv) {
 return EXIT_SUCCESS;
 }