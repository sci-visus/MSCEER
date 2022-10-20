//#include "integrate3.hpp"
//#include "vector2.hpp"
//#include <stdio.h>
#include <vector>
#include <set>
#include <queue>
#include <time.h>
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

#define USEMAXV
//using namespace GInt;
//typedef IndexCompareLessThan Comparer;
//typedef RegularGrid3D GridType;
//typedef TopologicalRegularGrid3D MeshType;
//typedef RegularGridTrilinearFunction GridFuncType;
//#ifndef USEMAXV
//typedef TopologicalExplicitDenseMeshFunction<MeshType, float> TopoFuncType;
//#else
//typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
//typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType , GridFuncType, float> TopoFuncType;
//#endif
//typedef DiscreteGradientLabeling<MeshType> GradType;
//typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType> RobinsType;
//typedef MyRobins<MeshType, MaxVLType, GradType> RobinsTypeO;

using namespace GInt;
typedef RegularGrid3D GridType;
typedef RegularGridTrilinearFunction GridFuncType;
//typedef UncachedRegularGridTrilinearFunction GridFuncType;
typedef TopologicalRegularGrid3D MeshType;
typedef DiscreteGradientLabeling<MeshType> GradType;
//typedef UncachedMaximumVertexLabeling<GridType, GridFuncType> MaxVLType;
//typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
//typedef RegularGridMaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
typedef RegularGridMaxMinVertexLabeling3D<MeshType, GridFuncType> MaxVLType;
typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType, 5, 4> RobinsType;
typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> TopoFuncType;


int main(int argc, char** argv) {


	int X, Y, Z;

	if (argc < 4) { printf( "Usage: X Y Z \n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	
	char fname[1024];
	sprintf(fname, "Offsets_%dx%dx%d.txt", X, Y, Z);
	FILE* fout = fopen(fname, "w");

	GridType* m_grid;
	MeshType *m_tgrid;
	GradType *labeling;

	// format: [global activity name] [task] [start] [end] [dration]


	// -- start timing IO
	m_grid = new GridType(Vec3l((long long) X, (long long) Y, (long long) Z), Vec3b(0,0,0));
	m_tgrid = new MeshType(m_grid);

	fprintf( fout,"NumGridCells: %d\n", m_grid->NumElements());
	fprintf( fout,"NumMeshCells: %d\n", m_tgrid->numCells());

	int num_cells = m_tgrid->numCells();

	// print dimension array
	fprintf( fout,"dimension array:\n");
	fprintf( fout,"{ ");
	MeshType::AllCellsIterator allit(m_tgrid);
	for (allit.begin(); allit.valid(); ) {
		INDEX_TYPE id = allit.value();
		fprintf( fout,"%d", m_tgrid->dimension(id));
		allit.advance();
		if (allit.valid()) fprintf( fout,", ");
	}
	fprintf( fout," }\n");

	// print VertexNumberFromCellID array
	fprintf( fout,"VertexNumberFromCellID array:\n");
	fprintf( fout,"{ ");
	for (allit.begin(); allit.valid();) {
		INDEX_TYPE id = allit.value();
		if (m_tgrid->dimension(id) != 0) fprintf( fout,"-1");
		else
			fprintf( fout,"%d", m_tgrid->VertexNumberFromCellID(id));
		allit.advance();
		if (allit.valid()) fprintf( fout,", ");
	}
	fprintf( fout," }\n");

	// print Facetsiterator array
	fprintf( fout,"FacetsIterator2d array [%d][7]:\n");
	fprintf( fout,"{\n");
	for (allit.begin(); allit.valid();) {
		INDEX_TYPE id = allit.value();
		int pos = 1;
		INDEX_TYPE negs[7];
		MeshType::FacetsIterator fit(m_tgrid);
		for (fit.begin(id); fit.valid();fit.advance()) {
			INDEX_TYPE nid = fit.value();
			negs[pos++] = nid;
		}
		negs[0] = pos - 1;
		while (pos < 7) {
			negs[pos++] = 0;
		}
		fprintf( fout,"{ ");
		for (int i = 0; i < 7; i++){
			fprintf( fout,"%d", negs[i]);
			if (i != 6) fprintf( fout,", ");
		}
		fprintf( fout,"}");
		allit.advance();
		if (allit.valid()) fprintf( fout,",\n");
	}
	fprintf( fout,"\n}\n");

	// print cofacetsiterator array
	fprintf( fout,"CoFacetsIterator2d array [%d][7]:\n");
	fprintf( fout,"{\n");
	for (allit.begin(); allit.valid();) {
		INDEX_TYPE id = allit.value();
		int pos = 1;
		INDEX_TYPE negs[7];
		MeshType::CofacetsIterator fit(m_tgrid);
		for (fit.begin(id); fit.valid(); fit.advance()) {
			INDEX_TYPE nid = fit.value();
			negs[pos++] = nid;
		}
		negs[0] = pos - 1;
		while (pos < 7) {
			negs[pos++] = 0;
		}
		fprintf( fout,"{ ");
		for (int i = 0; i < 7; i++){
			fprintf( fout,"%d", negs[i]);
			if (i != 6) fprintf( fout,", ");
		}
		fprintf( fout,"}");
		allit.advance();
		if (allit.valid()) fprintf( fout,",\n");
	}
	fprintf( fout,"\n}\n");
	

	fprintf(fout, "Adjacent nonboundary offsets:\n");
	MeshType::AdjacentCellsIterator adjit(m_tgrid);
	fprintf(fout, "{ ");
	for (int i = 0; i < 27; i++) {
		fprintf(fout, "%d", m_tgrid->get27NeighborOffset(i));
		if (i < 26) fprintf(fout, ", ");
	}
	fprintf(fout, " }\n");


	fclose(fout);
	return 1;


}


