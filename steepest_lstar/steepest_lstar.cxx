/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/


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
#include "gi_fast_robins_noalloc.h"

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
typedef SlidingWindowRobinsNoalloc < GridType, GridFuncType, MeshType, MaxVLType, GradType> NewRobinsType;

int main(int argc, char** argv) {

	ThreadedTimer timer(1);
	timer.StartGlobal();

	int X, Y, Z;
	int per_x, per_y, per_z;
	std::string filename;

	if (argc < 8) { printf("Usage: X Y Z filename per_x per_y per_z  [numthreads=max_on_machine]\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
	sscanf(argv[5], "%d", &per_x);
	sscanf(argv[6], "%d", &per_y);
	sscanf(argv[7], "%d", &per_z);
	int numthreads = omp_get_max_threads();
	if (argc >= 9)
		sscanf(argv[8], "%d", &numthreads);
	omp_set_num_threads(numthreads);

	printf("%d, %d, %d, %d, %d, %d, %d\n", X, Y, Z, per_x, per_y, per_z, numthreads);


	GridType* m_grid;
	//RegularGridTrilinearFunction* m_func;
	GridFuncType* m_func2;
	MeshType *m_tgrid;
	TopoFuncType* m_topofunc;
	GradType *labeling;
	TopologicalGradientUsingAlgorithms<MeshType, TopoFuncType, GradType>* topo_algs;

	// format: [global activity name] [task] [start] [end] [dration]


	// -- start timing IO
	m_grid = new GridType(Vec3l((long long) X, (long long) Y, (long long) Z), Vec3b(per_x, per_y, per_z));
	m_func2 = new GridFuncType(m_grid);
	m_func2->LoadImageFromFile(filename.c_str());
	std::chrono::steady_clock::time_point g_start_time = std::chrono::steady_clock::now();;
	std::chrono::steady_clock::time_point task_start_time;
	std::chrono::steady_clock::time_point prior_time;
	std::chrono::steady_clock::time_point now_time;

	// added by Harsh

	char timingname[2048];
	sprintf(timingname, "%s.steep.gtime.txt", argv[4]);
	FILE* gtiming = fopen(timingname, "w");
	prior_time = std::chrono::steady_clock::now();


	// START GRAD INTEGRATION --------------
	prior_time = std::chrono::steady_clock::now();
	task_start_time = std::chrono::steady_clock::now();

	//m_func2->Negate();
	printf("loaded cont function\n");
	m_tgrid = new MeshType(m_grid);
//	m_tgrid->set_restriction(norestrict);


#ifndef USEMAXV
	printf("made tgrid\n");
	m_topofunc = new TopoFuncType();
	m_topofunc->setMeshAndAllocate(m_tgrid);
	m_topofunc->copyVertexValuesFromGridFunction(m_func2);
	m_topofunc->setCellValuesMaxOfVerts();
	printf("made topofunc\n");
#else
	//printf("############  TESTING TIMING OF new MAXVL ############ \n");
	//MaxVLType_tst* test_newmaxv =
	//	new MaxVLType_tst(m_tgrid, m_func2);
	//test_newmaxv->ComputeOutput();
	//printf("############ DONE TESTING ############ \n");

	MaxVLType* maxv = new MaxVLType(m_tgrid, m_func2);
	maxv->ComputeOutput();


	//printf("############  TESTING Similarity OF new MAXVL ############ \n");
	//MeshType::AllCellsIterator ait(m_tgrid);
	//for (ait.begin(); ait.valid(); ait.advance()) {
	//	INDEX_TYPE cellid = ait.value();
	//	INDEX_TYPE v0 = maxv->Cell2HighestVertex(cellid);
	//	INDEX_TYPE v1 = test_newmaxv->Cell2HighestVertex(cellid);
	//	if (v0 != v1) {
	//		printf("ERROR %llu: %llu != %llu\n", cellid, v0, v1);
	//		Vec3l c0; m_tgrid->cellid2Coords(cellid, c0);
	//		c0.PrintInt();
	//	}
	//}

	//printf("############ DONE TESTING ############ \n");

	m_topofunc = new TopoFuncType();
	m_topofunc->setMeshAndFuncAndMaxVLabeling(m_tgrid, m_func2, maxv);
#endif




	//DenseLabeling<char>* norestrict = new DenseLabeling<char>(m_tgrid->numCells());
	//norestrict->SetAll(0);

	labeling = new GradType(m_tgrid);
	printf("created dgrad struct\n");

	//for (int i = 0; i < 5; i++){
		labeling->ClearAllGradient();
		now_time = std::chrono::steady_clock::now();
#if 1
		NewRobinsType* trobins =
			new NewRobinsType(m_grid, m_func2, m_tgrid, maxv, labeling);
		trobins->ComputePairing_sliding();
#else

		GradType* labeling2 = new GradType(m_tgrid);
		labeling2->ClearAllGradient();

		RobinsType* trobins_o =
			new RobinsType(m_tgrid, maxv, labeling2);
		trobins_o->ComputePairing();

		NewRobinsType* trobins =
			new NewRobinsType(m_grid, m_func2, m_tgrid, maxv, labeling);
		trobins->ComputePairing_sliding();
		//return 1;
		MeshType::AllCellsIterator allit(m_tgrid);
		for (allit.begin(); allit.valid(); allit.advance()) {
			INDEX_TYPE cid = allit.value();
			INDEX_TYPE p1 = labeling->getPair(cid);
			INDEX_TYPE p2 = labeling2->getPair(cid);
			if (m_tgrid->dimension(cid) == 0 &&  p1 != p2){
				printf("ERROR grads dont match %lld : %lld != %lld , boundary = %d!!\n", cid, p1, p2, m_tgrid->boundaryValue(cid)  );
				Vec3l cv, p1v, p2v;
				m_tgrid->cellid2Coords(cid, cv);
				m_tgrid->cellid2Coords(p1, p1v);
				m_tgrid->cellid2Coords(p2, p2v);
				printf("[%f :: %f] ",
					m_func2->SampleImage(m_tgrid->VertexNumberFromCellID(maxv->Cell2LowestVertex(cid))),
					m_func2->SampleImage(m_tgrid->VertexNumberFromCellID(maxv->Cell2HighestVertex(cid))));
				cv.PrintInt(); 
				printf("[%f :: %f] ",
					m_func2->SampleImage(m_tgrid->VertexNumberFromCellID(maxv->Cell2LowestVertex(p1))),
					m_func2->SampleImage(m_tgrid->VertexNumberFromCellID(maxv->Cell2HighestVertex(p1))));
				p1v.PrintInt(); 
				printf("[%f :: %f] ",
					m_func2->SampleImage(m_tgrid->VertexNumberFromCellID(maxv->Cell2LowestVertex(p2))),
					m_func2->SampleImage(m_tgrid->VertexNumberFromCellID(maxv->Cell2HighestVertex(p2))));
				p2v.PrintInt();

			}
		}

#endif
		printf("robins1 %d\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
//		delete trobins;

//	}
	printf("computed robins\n");

	now_time = std::chrono::steady_clock::now();

	fprintf(gtiming, "Total Overall %d %d %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(g_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count());
	printf( "Total Overall %d %d %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(g_start_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count());
	//return 1;

	//for (int i = 0; i < 5; i++){
	//	labeling->ClearAllGradient();
	//	now_time = std::chrono::steady_clock::now();
	//	RobinsTypeO* trobins =
	//		new RobinsTypeO(m_tgrid, maxv, norestrict, labeling);
	//	trobins->ComputePairing();
	//	printf("robins2 %d\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
	//	delete trobins;

	//}
	//printf("computed robins\n");


	fclose(gtiming);

	printf("setting dim of asc manifold for fast loading\n");
	topo_algs = new TopologicalGradientUsingAlgorithms<MeshType, TopoFuncType, GradType>(m_topofunc, m_tgrid, labeling);
	topo_algs->setAscendingManifoldDimensions();
	//topo_algs->CheckGradientConsistency();

	//DenseLabeling<unsigned char>* tt = new DenseLabeling<unsigned char>(m_tgrid->numCells());
	//MeshType::AllCellsIterator ait(m_tgrid);
	//for (ait.begin(); ait.valid(); ait.advance()) {
	//	INDEX_TYPE id = ait.value();
	//	tt->SetLabel(id, labeling->getDimAscMan(id));
	//}
	//tt->OutputToFile("DIMS.raw");
	//topo_algs->count_critical_points(4);
	printf("writing output\n");
	char gradname[2048];
	sprintf(gradname, "%s.grad", argv[4]);
	labeling->outputToFile(gradname);

	return 1;

	
	return 1;


}


