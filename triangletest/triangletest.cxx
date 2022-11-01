/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

//#include "integrate3.hpp"
//#include "vector2.hpp"
//#include <stdio.h>
#include <vector>
#include <set>
#include <queue>
#include <time.h>
#include "gi_timing.h"
#include "gi_topological_simplicial_complex.h"
#include "gi_max_vertex_labeling.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_discrete_gradient_labeling.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_bifiltration_pairing.h"
using namespace GInt;


typedef TopologicalSimplicialComplex2d MeshType;
typedef SimplicialComplexCoordinateFunction GridFuncType;

typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> TopoFuncType;

typedef DiscreteGradientLabeling<MeshType> GradType;
typedef MorseSmaleComplexBasic<FLOATTYPE, MeshType, TopoFuncType, GradType> MSCType;




MeshType* mesh;



void TestMesh(MeshType* m) {



	MeshType::AllCellsIterator ait(m);
	for (ait.begin(); ait.valid(); ait.advance()) {
		INDEX_TYPE cid = ait.value();

		MeshType::FacetsIterator fit(m);
		for (fit.begin(cid); fit.valid(); fit.advance()) {
			INDEX_TYPE fid = fit.value();
			bool has = false;
			MeshType::CofacetsIterator cfit(m);
			for (cfit.begin(fid); cfit.valid(); cfit.advance()) {
				INDEX_TYPE cfid = cfit.value();
				if (cfid == cid) has = true;
			}
			if (!has) 
				printf("inconsistency!!\n");
		}

		MeshType::CofacetsIterator cfit(m);
		for (cfit.begin(cid); cfit.valid(); cfit.advance()) {
			INDEX_TYPE cfid = cfit.value();
			bool has = false;

			MeshType::FacetsIterator fit(m);
			for (fit.begin(cfid); fit.valid(); fit.advance()) {
				INDEX_TYPE fid = fit.value();
				if (fid == cid) has = true;
			}
			if (!has) 
				printf("inconsistency2!!\n");
		}
	}


		MeshType::DCellsIterator tris(m, 2);
		for (tris.begin(); tris.valid(); tris.advance()) {

			INDEX_TYPE trid = tris.value();
			std::set<INDEX_TYPE> verts;

			MeshType::CellVerticesIterator viter(m);
			for (viter.begin(trid); viter.valid(); viter.advance()) {
				INDEX_TYPE vid = viter.value();

				if (m->dimension(vid) != 0) {
					printf("ERROR: vertex iterator non-vertex cells\n");
				}

				verts.insert(vid);

				BYTE_TYPE vo = m->CompressVertexOffsetToByte(trid, vid, viter);


				INDEX_TYPE nvid = m->UncompressByteToVertexOffset(trid, vo);
				if (nvid != vid) printf("ERROR: mismatch: %d -> %d -> %d\n", (int)vid, vo, (int)nvid);

			}

			if (verts.size() != 3) {
				printf("ERROR: triangle has %d verts\n", verts.size());
			}

		}
		printf("done test\n");

}

void TestMaxVlabelMesh(MeshType* m, MaxVLType* ml, GridFuncType* f) {



	MeshType::DCellsIterator tris(m, 2);
	for (tris.begin(); tris.valid(); tris.advance()) {

		INDEX_TYPE trid = tris.value();

		INDEX_TYPE trimaxv = ml->Cell2Vertex(trid);

		MeshType::CellVerticesIterator viter(m);
		for (viter.begin(trid); viter.valid(); viter.advance()) {
			INDEX_TYPE vid = viter.value();

			if (vid == trimaxv) continue;
			
			if (f->IsGreater(vid, trimaxv)) printf("ERROR: maxvl sucks!\n");

		}
	}
	printf("done test\n");

}

int main(int argc, char** argv) {


	if (argc < 2) { printf("Usage: filename\n"); return 0; }

	printf("starting testing:\n");
	mesh = new MeshType();
	mesh->LoadFromObj(argv[1]);

	TestMesh(mesh);

	printf("read mesh!\n");

	GridFuncType* cfunc1 = new GridFuncType(mesh, 0);
	MaxVLType* mvl1 = new MaxVLType(mesh, cfunc1);
	mvl1->ComputeOutput();

	printf("testing maxvl1\n");
	TestMaxVlabelMesh(mesh, mvl1, cfunc1);





	GridFuncType* cfunc2 = new GridFuncType(mesh, 1);
	MaxVLType* mvl2 = new MaxVLType(mesh, cfunc2);
	mvl2->ComputeOutput();



	BifiltrationPairing<MeshType, MaxVLType, MaxVLType>* pp = new BifiltrationPairing<MeshType, MaxVLType, MaxVLType>(mesh, mvl1, mvl2);
	pp->ComputePairing();


	return 1;


}


