
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
#include "gi_modified_robins.h"
#include "gi_fast_robins_noalloc.h"

#define USEMAXV

using namespace GInt;
typedef RegularGrid3D GridType;
typedef SparseRegularGridTrilinearFunction GridFuncType;
typedef TopologicalRegularMaskedGrid MeshType;
typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;

typedef DiscreteGradientLabeling<MeshType, IndirectLabeling<GradBitfield> > GradType;
typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType, 5, 4> RobinsType;
typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> TopoFuncType;
typedef SlidingWindowRobinsNoalloc < GridType, GridFuncType, MeshType, MaxVLType, GradType> NewRobinsType;

#include <unordered_map>
template<typename dtype>
class SparseBlockedArray {
protected:
	class SubBlock {
		GInt::Vec3i m_dims;
		GInt::Vec3i m_ghost;
		dtype* m_data;
	};
	GInt::Vec3i m_global_dims;
	GInt::Vec3i m_block_dims;
	GInt::Vec3i m_num_blocks_axis;
	GInt::Vec3i m_ghost;

};


int main(int argc, char** argv) {

	ThreadedTimer timer(1);
	timer.StartGlobal();

	int X, Y, Z;
	int per_x, per_y, per_z;
	std::string filename_raw;
	std::string filename_mask;

	if (argc < 6) { printf("Usage: X Y Z filename_raw filename_mask\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename_raw = std::string(argv[4]);
	filename_mask = std::string(argv[5]);

	printf("starting compute....\n");
	printf("%d, %d, %d\n", X, Y, Z);


	GridType* underlying_grid;
	//RegularGridTrilinearFunction* m_func;
	GridFuncType* sparse_grid_function;
	MeshType *sparse_mesh;
	MaxVLType* sparse_max_vert_labeling;
	TopoFuncType* m_topofunc;
	GradType *labeling;
	TopologicalGradientUsingAlgorithms<MeshType, TopoFuncType, GradType>* topo_algs;

	// format: [global activity name] [task] [start] [end] [dration]


	// -- start timing IO
	underlying_grid = new GridType(Vec3l((long long) X, (long long) Y, (long long) Z), Vec3b(0, 0, 0));
	
	printf(" -- loading sparse function\n");
	sparse_grid_function = new GridFuncType(underlying_grid);
	sparse_grid_function->LoadImageFromFloatAndMaskFile(filename_raw.c_str(), filename_mask.c_str());
	printf(" -- loaded %d values\n", sparse_grid_function->GetSparseMap().GetNumLabels());



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

	//sparse_grid_function->Negate();
	
	sparse_mesh = new MeshType(underlying_grid);
	
	// first set tester to ALL cells so we can fill mask
	sparse_mesh->SetTester(new CellTesterDefaultTrue());
	auto MaskMap = sparse_grid_function->GetSparseMap();
	
	printf("Number of Input Vertices: %d\n", MaskMap.GetNumLabels());
	printf(" -- recording input vertices in mesh map\n");
	// make a cellid to mask map
	SparseLabeling<BYTE_TYPE>* mesh_mask = new SparseLabeling<BYTE_TYPE>();
	
	// faster to first add all vertices to the mask
	for (SparseLabeling<FLOATTYPE>::SparseKeyIterator it = MaskMap.begin(); it != MaskMap.end(); it++) {
		// each "key" is the vertex number of existing vertex
		auto vertex_number = *it;
		auto cell_id = sparse_mesh->CellIDFromVertexNumber(vertex_number);
		// add vertex to mask
		mesh_mask->SetLabel(cell_id, 1);
	}
	printf(" -- extending to higher dim cells mesh map\n");

	MeshType::CofacetsIterator cfit(sparse_mesh);
	MeshType::CellVerticesIterator vit(sparse_mesh);
	// add higher dimensional cells to mask
	for (SparseLabeling<FLOATTYPE>::SparseKeyIterator it = MaskMap.begin(); it != MaskMap.end(); it++) {
		// each "key" is the vertex number of existing vertex
		auto vertex_number = *it;
		auto cell_id = sparse_mesh->CellIDFromVertexNumber(vertex_number);
		// check each cofacet
		for (cfit.begin(cell_id); cfit.valid(); cfit.advance()) {
			auto cofacet_id = cfit.value();
			bool has_missing = false;
			for (vit.begin(cofacet_id); vit.valid(); vit.advance()) {
				auto vertex_id = vit.value();
				if (!mesh_mask->Has(vertex_id)) {
					has_missing = true;
					break;
				}				
			}
			if (!has_missing) {
				mesh_mask->SetLabel(cofacet_id, 1);
			}
		}
	}

	printf("Number of Interior Cells: %d\n", mesh_mask->GetNumLabels());

	// now apply the mask to the masked grid
	sparse_mesh->SetTester(new CellTesterSparseLabeled<BYTE_TYPE>(mesh_mask));

	// now create a max vl 
	printf(" -- computing max vertex labeling\n");
	sparse_max_vert_labeling = new MaxVLType(sparse_mesh, sparse_grid_function);
	sparse_max_vert_labeling->ComputeOutput();
	printf(" -- done\n");


	INDEX_TYPE num_cells = sparse_mesh->numCells();

//	sparse_mesh->set_restriction(norestrict);

//
//#ifndef USEMAXV
//	printf("made tgrid\n");
//	m_topofunc = new TopoFuncType();
//	m_topofunc->setMeshAndAllocate(sparse_mesh);
//	m_topofunc->copyVertexValuesFromGridFunction(sparse_grid_function);
//	m_topofunc->setCellValuesMaxOfVerts();
//	printf("made topofunc\n");
//#else
//	//printf("############  TESTING TIMING OF new MAXVL ############ \n");
//	//MaxVLType_tst* test_newmaxv =
//	//	new MaxVLType_tst(sparse_mesh, sparse_grid_function);
//	//test_newmaxv->ComputeOutput();
//	//printf("############ DONE TESTING ############ \n");
//
//	MaxVLType* maxv = new MaxVLType(sparse_mesh, sparse_grid_function);
//	maxv->ComputeOutput();
//
//
//	//printf("############  TESTING Similarity OF new MAXVL ############ \n");
//	//MeshType::AllCellsIterator ait(sparse_mesh);
//	//for (ait.begin(); ait.valid(); ait.advance()) {
//	//	INDEX_TYPE cellid = ait.value();
//	//	INDEX_TYPE v0 = maxv->Cell2HighestVertex(cellid);
//	//	INDEX_TYPE v1 = test_newmaxv->Cell2HighestVertex(cellid);
//	//	if (v0 != v1) {
//	//		printf("ERROR %llu: %llu != %llu\n", cellid, v0, v1);
//	//		Vec3l c0; sparse_mesh->cellid2Coords(cellid, c0);
//	//		c0.PrintInt();
//	//	}
//	//}
//
//	//printf("############ DONE TESTING ############ \n");
//
//	m_topofunc = new TopoFuncType();
//	m_topofunc->setMeshAndFuncAndMaxVLabeling(sparse_mesh, sparse_grid_function, maxv);
//#endif
//
//
//
//
//	//DenseLabeling<char>* norestrict = new DenseLabeling<char>(sparse_mesh->numCells());
//	//norestrict->SetAll(0);
//
//	labeling = new GradType(sparse_mesh);
//	printf("created dgrad struct\n");
//
//	//for (int i = 0; i < 5; i++){
//		labeling->ClearAllGradient();
//		now_time = std::chrono::steady_clock::now();
//#if 1
//		NewRobinsType* trobins =
//			new NewRobinsType(underlying_grid, sparse_grid_function, sparse_mesh, maxv, labeling);
//		trobins->ComputePairing_sliding();
//#else
//
//		GradType* labeling2 = new GradType(sparse_mesh);
//		labeling2->ClearAllGradient();
//
//		RobinsType* trobins_o =
//			new RobinsType(sparse_mesh, maxv, labeling2);
//		trobins_o->ComputePairing();
//
//		NewRobinsType* trobins =
//			new NewRobinsType(underlying_grid, sparse_grid_function, sparse_mesh, maxv, labeling);
//		trobins->ComputePairing_sliding();
//		//return 1;
//		MeshType::AllCellsIterator allit(sparse_mesh);
//		for (allit.begin(); allit.valid(); allit.advance()) {
//			INDEX_TYPE cid = allit.value();
//			INDEX_TYPE p1 = labeling->getPair(cid);
//			INDEX_TYPE p2 = labeling2->getPair(cid);
//			if (sparse_mesh->dimension(cid) == 0 &&  p1 != p2){
//				printf("ERROR grads dont match %lld : %lld != %lld , boundary = %d!!\n", cid, p1, p2, sparse_mesh->boundaryValue(cid)  );
//				Vec3l cv, p1v, p2v;
//				sparse_mesh->cellid2Coords(cid, cv);
//				sparse_mesh->cellid2Coords(p1, p1v);
//				sparse_mesh->cellid2Coords(p2, p2v);
//				printf("[%f :: %f] ",
//					sparse_grid_function->SampleImage(sparse_mesh->VertexNumberFromCellID(maxv->Cell2LowestVertex(cid))),
//					sparse_grid_function->SampleImage(sparse_mesh->VertexNumberFromCellID(maxv->Cell2HighestVertex(cid))));
//				cv.PrintInt(); 
//				printf("[%f :: %f] ",
//					sparse_grid_function->SampleImage(sparse_mesh->VertexNumberFromCellID(maxv->Cell2LowestVertex(p1))),
//					sparse_grid_function->SampleImage(sparse_mesh->VertexNumberFromCellID(maxv->Cell2HighestVertex(p1))));
//				p1v.PrintInt(); 
//				printf("[%f :: %f] ",
//					sparse_grid_function->SampleImage(sparse_mesh->VertexNumberFromCellID(maxv->Cell2LowestVertex(p2))),
//					sparse_grid_function->SampleImage(sparse_mesh->VertexNumberFromCellID(maxv->Cell2HighestVertex(p2))));
//				p2v.PrintInt();
//
//			}
//		}
//
//#endif
//		printf("robins1 %d\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
////		delete trobins;
//
////	}
//	printf("computed robins\n");
//
//	now_time = std::chrono::steady_clock::now();
//
//	fprintf(gtiming, "Total Overall %d %d %d\n",
//		std::chrono::duration_cast<std::chrono::milliseconds>(g_start_time - g_start_time).count(),
//		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
//		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count());
//	printf( "Total Overall %d %d %d\n",
//		std::chrono::duration_cast<std::chrono::milliseconds>(g_start_time - g_start_time).count(),
//		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
//		std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count());
//	//return 1;
//
//	//for (int i = 0; i < 5; i++){
//	//	labeling->ClearAllGradient();
//	//	now_time = std::chrono::steady_clock::now();
//	//	RobinsTypeO* trobins =
//	//		new RobinsTypeO(sparse_mesh, maxv, norestrict, labeling);
//	//	trobins->ComputePairing();
//	//	printf("robins2 %d\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
//	//	delete trobins;
//
//	//}
//	//printf("computed robins\n");
//
//
//	fclose(gtiming);
//
//	printf("setting dim of asc manifold for fast loading\n");
//	topo_algs = new TopologicalGradientUsingAlgorithms<MeshType, TopoFuncType, GradType>(m_topofunc, sparse_mesh, labeling);
//	topo_algs->setAscendingManifoldDimensions();
//	//topo_algs->CheckGradientConsistency();
//
//	//DenseLabeling<unsigned char>* tt = new DenseLabeling<unsigned char>(sparse_mesh->numCells());
//	//MeshType::AllCellsIterator ait(sparse_mesh);
//	//for (ait.begin(); ait.valid(); ait.advance()) {
//	//	INDEX_TYPE id = ait.value();
//	//	tt->SetLabel(id, labeling->getDimAscMan(id));
//	//}
//	//tt->OutputToFile("DIMS.raw");
//	//topo_algs->count_critical_points(4);
//	printf("writing output\n");
//	char gradname[2048];
//	sprintf(gradname, "%s.grad", argv[4]);
//	labeling->outputToFile(gradname);
//
//	return 1;
//
//	
	return 1;


}


