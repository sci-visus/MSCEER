#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif

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
#include "gi_modified_robins.h"
#include "gi_morse_smale_complex_basic.h"

typedef GInt::RegularGrid3D GridType;

using namespace GInt;

int main(int argc, char** argv) {

	// READ IN THE COMMAND LINE ARGUMENTS
	int X, Y, Z;
	std::string filename;
	if (argc < 5) { printf("Usage: X Y Z filename threshold\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
        filename = std::string(argv[4]);
        if (filename.find("pancakes.raw") == std::string::npos) 
	  filename = std::string(argv[4])+".pancakes.raw";
	float threshold;
	sscanf(argv[5], "%f", &threshold);
  
  float* mask_field = NULL;
  
  mask_field = new float[X*Y*Z];
  FILE* fin = fopen(filename.c_str(), "rb");
  fread(mask_field, sizeof(float), X*Y*Z, fin);
  fclose(fin);
  
	GridType* underlying_grid;
	// set up structures to navigate grid, and load the 3d image
  underlying_grid = new GridType(GInt::Vec3i(X, Y, Z), GInt::Vec3b(0, 0, 0));
	
  printf("loaded field function\n");
  
  VolumeConnectedComponents cc(underlying_grid);
  DenseLabeling<char> *maskvol = new DenseLabeling<char>(underlying_grid->NumElements());
  
  for(INDEX_TYPE i=0; i< underlying_grid->NumElements();i++)
    maskvol->SetLabel(i, mask_field[i] > threshold);
    
  cc.PerformUnionFind(maskvol);
  
  int* output = new int[underlying_grid->NumElements()];
 
  cc.mIDVol->ReMapIds(output);
  
  std::string outfilename = std::string(filename+".components.raw");
  printf("saving output to %s\n", outfilename.c_str());
  ofstream outFile;
  outFile.open(outfilename.c_str(), ios::out|ios::binary);
  outFile.write((char*)output, underlying_grid->NumElements()*sizeof(int));
  outFile.close();
  
  printf("dump connected comps\n");

	return 0;

}


