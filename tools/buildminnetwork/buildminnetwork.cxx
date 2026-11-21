//#include "integrate3.hpp"
//#include "vector2.hpp"
//#include <stdio.h>
#include <vector>
#include <set>
#include <queue>
#include <time.h>
#include "gi_timing.h"
#include "gi_topological_explicit_mesh_function.h"
#include "gi_topological_region_growing_simple_gradient_builder.h"
#include "gi_topological_convergent_gradient_builder.h"
#include "gi_robin_labeling.h"
#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_isolated_region_remover.h"
#include "gi_topological_utility_functions.h"
#include "gi_numeric_integrator_expanding_region_stop.h" // not where comparer should be
#include "gi_topological_regular_masked_grid.h"
#include "gi_topological_regular_masked_restricted_grid.h"
#include "gi_dataflow_objects.h"
#include "gi_dfo_morse_smale_selectors.h"
#include "gi_max_vertex_labeling.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_basic_geometry.h"

#define USEMAXV

using namespace GInt;
typedef IndexCompareLessThan Comparer;
typedef RegularGrid GridType;
typedef TopologicalRegularGrid MeshType;
typedef RegularGridTrilinearFunction GridFuncType;
#ifndef USEMAXV
typedef TopologicalExplicitDenseMeshFunction<MeshType, float> TopoFuncType;
#else
typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> TopoFuncType;
#endif
typedef DiscreteGradientLabeling<MeshType> GradType;
typedef MorseSmaleComplexBasic<FLOATTYPE, MeshType, TopoFuncType, GradType> MSCType;

int main(int argc, char** argv) {

	ThreadedTimer timer(1);
	timer.StartGlobal();

	int X, Y, Z;
	int per_x, per_y, per_z;
	std::string filename;

	if (argc < 8) { printf("Usage: X Y Z filename per_x per_y per_z  [outputdebug=0]\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
	sscanf(argv[5], "%d", &per_x);
	sscanf(argv[6], "%d", &per_y);
	sscanf(argv[7], "%d", &per_z);
	int outputdebug = 0;
	if (argc >= 9)
		sscanf(argv[8], "%d", &outputdebug);

	GridType* m_grid;
	//RegularGridTrilinearFunction* m_func;
	GridFuncType* m_func2;
	MeshType *m_tgrid;
	TopoFuncType* m_topofunc;
	GradType *labeling;
	TopologicalGradientUsingAlgorithms<MeshType, TopoFuncType, GradType>* topo_algs;

	// format: [global activity name] [task] [start] [end] [dration]


	// -- start timing IO
	m_grid = new GridType(Vec3i(X, Y, Z), Vec3b(per_x, per_y, per_z));
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
	MaxVLType* maxv = new MaxVLType(m_tgrid, m_func2);
	maxv->ComputeOutput();
	m_topofunc = new TopoFuncType();
	m_topofunc->setMeshAndFuncAndMaxVLabeling(m_tgrid, m_func2, maxv);
#endif




	DenseLabeling<char>* norestrict = new DenseLabeling<char>(m_tgrid->numCells());
	norestrict->SetAll(1);

	//	int num_cells = m_tgrid->numCells();
	//#pragma omp parallel
	//	{
	//		int num_threads = omp_get_num_threads();
	//		int thread_num = omp_get_thread_num();
	//
	//		std::vector<INDEX_TYPE> partition;
	//		ArrayIndexPartitioner::EvenChunkSplit(m_tgrid->numCells(), num_threads, partition);
	//
	//		GridType2::TopologicalRegularGrid::DCellsIterator vit(m_tgrid, 0, partition[thread_num], partition[thread_num + 1]);
	//		for (vit.begin(); vit.valid(); vit.advance()) {
	//			if (true /* THIS IS WHERE TRHESHOLD TEST GOES*/) {
	//				INDEX_TYPE vid = vit.value();
	//				norestrict->SetLabel(vid, 0);
	//
	//
	//				GridType2::TopologicalRegularGrid::AdjacentCellsIterator cviter(m_tgrid);
	//
	//				for (cviter.begin(vid); cviter.valid(); cviter.advance()) {
	//					norestrict->SetLabel(cviter.value(), 0);
	//				}
	//			}
	//		}
	//	}


		//CellTesterLaberInput* cellrestrict = new CellTesterLaberInput();
		//cellrestrict->SetLabeling(norestrict);
		//m_tgrid->SetTester(cellrestrict);
		//m_tgrid->set_restriction(norestrict);
		//printf("made norestrict\n");
	char gradname[1024];
	sprintf(gradname, "%s.grad", argv[4]);

	labeling = new GradType(m_tgrid);
	if (!labeling->load_from_file(gradname)) {
		labeling->ClearAllGradient();
		printf("created dgrad struct\n");
		RobinsLabelingAlgorithm<MeshType, TopoFuncType>* trobins =
			new RobinsLabelingAlgorithm<MeshType, TopoFuncType>(m_topofunc, m_tgrid, labeling);
		trobins->compute_output();
		printf("computed robins\n");

		now_time = std::chrono::steady_clock::now();

		fprintf(gtiming, "Total Overall %d %d %d\n",
			std::chrono::duration_cast<std::chrono::milliseconds>(g_start_time - g_start_time).count(),
			std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count(),
			std::chrono::duration_cast<std::chrono::milliseconds>(now_time - g_start_time).count());


		fclose(gtiming);

		topo_algs = new TopologicalGradientUsingAlgorithms<MeshType, TopoFuncType, GradType>(m_topofunc, m_tgrid, labeling);
		topo_algs->setAscendingManifoldDimensions();
		//topo_algs->CheckGradientConsistency();
		labeling->outputToFile(gradname);
	}

	MSCType* msc = new MSCType(labeling, m_tgrid, m_topofunc);
	msc->SetBuildArcGeometry(GInt::Vec3b(true, false, false));
	msc->ComputeFromGrad();
	msc->set_output_cancellation_records("pers.txt");

	msc->ComputeHierarchy(0.00900372717419840);
	printf("Done cancelling\n");

	//FILE* fnodes = fopen("nodesout.txt", "w");
	//MSCType::LivingNodesIterator nodes(msc);
	//for (nodes.begin(); nodes.valid(); nodes.advance()) {

	//	auto n = msc->getNode(nodes.value());
	//	
	//	fprintf(fnodes, "%d %f %f\n", n.dim, n.value, m_topofunc->cellValue(n.cellindex));
	//}
	//fclose(fnodes);
	const float pvals[6] = { 0.0001f, 0.0002f, 0.0003f, 0.0005f, 0.0007f, 0.0009f };
	const float fvals[1] = { 0.0001f };

	float fval = 0.01;
	for (auto pval : pvals) {
		msc->SetSelectPersAbs(pval);

		for (auto fval : fvals) {

			int counter = 0;
			map<pair<INDEX_TYPE, INDEX_TYPE>, INDEX_TYPE> mmmap;
			MSCType::LivingArcsIterator arcs(msc);
			for (arcs.begin(); arcs.valid(); arcs.advance()) {
				vector<INDEX_TYPE> res;
				auto t_arc = msc->getArc(arcs.value());
				auto t_node = msc->getNode(t_arc.upper);
				 

				//if (msc->getNode(msc->getArc(arcs.value()).lower).boundary != 0) continue;
				//printf("vals = %f, bdry = %d\n", t_node.value, t_node.boundary);
				if (t_node.dim == 1 && t_node.value > fval) {
					INDEX_TYPE other = -1;
					MSCType::SurroundingLivingArcsIterator ait(msc);
					for (ait.begin(t_arc.upper); ait.valid(); ait.advance()) {
						auto tt_arc = msc->getArc(ait.value());
						if (tt_arc.dim != 0) continue;
						if (tt_arc.lower == t_arc.lower) continue;
						other = tt_arc.lower;
					}

					pair<INDEX_TYPE, INDEX_TYPE> p;
					if (other < t_arc.lower) {
						p = pair<INDEX_TYPE, INDEX_TYPE>(other, t_arc.lower);
					}
					else {
						p = pair<INDEX_TYPE, INDEX_TYPE>(t_arc.lower, other);
					}
					if (mmmap.count(p) == 0) {
						mmmap[p] = t_arc.upper;
					counter++;
					}
				}

			}


			char fnameout[1024];
			sprintf(fnameout, "linesoup_%f_%d.bin", fval, counter);
			FILE* linesoup = fopen(fnameout, "wb");
			//for (arcs.begin(); arcs.valid(); arcs.advance()) {
			for(auto pp : mmmap) {
				INDEX_TYPE sad = pp.second;
				MSCType::SurroundingLivingArcsIterator ait(msc);
				for (ait.begin(sad); ait.valid(); ait.advance()) {
					INDEX_TYPE aid = ait.value();

					//					auto t_arc = msc->getArc(arcs.value());
					auto t_arc = msc->getArc(aid);


					vector<INDEX_TYPE> res;

					auto t_node = msc->getNode(t_arc.upper);
					//if (msc->getNode(msc->getArc(arcs.value()).lower).boundary != 0) continue;
					//printf("vals = %f, bdry = %d\n", t_node.value, t_node.boundary);
					if (t_node.dim == 1/* && t_node.value > fval*/) {
						msc->fillArcGeometry(aid, res);
						Line3d line;
						for (auto id : res) {
							if (m_topofunc->cellValue(id) < fval) continue;
							Vec3l coords;
							m_tgrid->cellid2Coords(id, coords);
							//coords.PrintInt();
							Vec3d dcoords = coords;
							dcoords *= 0.5;
							//dcoords.PrintFloat();
							line.AddToEnd(dcoords);
						}
						line.IrreversableSmooth(25);

						const vector<Vec3d>& ll = line.GetLine();
						for (int i = 0; i < ll.size() - 1; i++) {

							Vec3d c1f = ll[i];
							Vec3d c2f = ll[i + 1];

							Vec3d diff = c1f - c2f;
							if (diff.Mag() > 10.0) {
								printf("edge segment  %f in arc %d: \n", diff.Mag(), arcs.value());
								c1f.PrintFloat();
								c2f.PrintFloat();
							}
						}

						line.WriteLineToFileBin(linesoup);
					}
				}

			}
			fclose(linesoup);
			printf("wrote %d lines\n", counter);



			//int counter = 0;
			//MSCType::LivingArcsIterator arcs(msc);
			//for (arcs.begin(); arcs.valid(); arcs.advance()) {
			//	vector<INDEX_TYPE> res;
			//	auto t_arc = msc->getArc(arcs.value());
			//	auto t_node = msc->getNode(t_arc.lower);
			//	//if (msc->getNode(msc->getArc(arcs.value()).lower).boundary != 0) continue;
			//	//printf("vals = %f, bdry = %d\n", t_node.value, t_node.boundary);
			//	if (t_node.dim == 0 && t_node.value > fval) {
			//		counter++;
			//	}

			//}


			//char fnameout[1024];
			//sprintf(fnameout, "linesoup_%f_%d.bin", fval, counter);
			//FILE* linesoup = fopen(fnameout, "wb");
			//for (arcs.begin(); arcs.valid(); arcs.advance()) {
			//	vector<INDEX_TYPE> res;
			//	auto t_arc = msc->getArc(arcs.value());
			//	auto t_node = msc->getNode(t_arc.lower);
			//	//if (msc->getNode(msc->getArc(arcs.value()).lower).boundary != 0) continue;
			//	//printf("vals = %f, bdry = %d\n", t_node.value, t_node.boundary);
			//	if (t_node.dim == 0 && t_node.value > fval) {
			//		msc->fillArcGeometry(arcs.value(), res);
			//		Line3d line;
			//		for (auto id : res) {
			//			Vec3l coords;
			//			m_tgrid->cellid2Coords(id, coords);
			//			//coords.PrintInt();
			//			Vec3d dcoords = coords;
			//			dcoords *= 0.5;
			//			//dcoords.PrintFloat();
			//			line.AddToEnd(dcoords);
			//		}
			//		line.IrreversableSmooth(25);

			//		const vector<Vec3d>& ll = line.GetLine();
			//		for (int i = 0; i < ll.size() - 1; i++) {

			//			Vec3d c1f = ll[i];
			//			Vec3d c2f = ll[i + 1];

			//			Vec3d diff = c1f - c2f;
			//			if (diff.Mag() > 10.0) {
			//				printf("edge segment  %f in arc %d: \n", diff.Mag(), arcs.value());
			//				c1f.PrintFloat();
			//				c2f.PrintFloat();
			//			}
			//		}

			//		line.WriteLineToFileBin(linesoup);
			//	}

			//}
			//fclose(linesoup);
			//printf("wrote %d lines\n", counter);



		}
	}




	return 1;










	//MSCType* msc = new MSCType(labeling, m_tgrid, m_topofunc);
	//msc->ComputeFromGrad();
	//msc->ComputeHierarchy(10.0);

	//DFO_VaryingValuatorDouble* pers_select = new DFO_VaryingValuatorDouble();
	//pers_select->SetValue(5.0);

	//DFO_MSCHierarchyInterface<MSCType>* hierarchy = new DFO_MSCHierarchyInterface<MSCType>();
	//hierarchy->SetMsc(msc);
	//hierarchy->SetPersistenceValuatorDFO(pers_select);

	//DFO_MSCLivingNodes<MSCType>* nodes = new DFO_MSCLivingNodes<MSCType>();
	//nodes->SetMscDFO(hierarchy);

	//DFO_MSCLivingArcs<MSCType>* arcs = new DFO_MSCLivingArcs<MSCType>();
	//arcs->SetMscDFO(hierarchy);

	//DFO_MSCNodeIndexSelect<MSCType>* mins = new DFO_MSCNodeIndexSelect<MSCType>();
	//mins->AddInputDFO(nodes);
	//mins->SetMscDFO(hierarchy);
	//mins->SetPassingIndices(true, false, false, false);

	//DFO_MSCArcsLowerIndexSelect<MSCType>* arcs01 = new DFO_MSCArcsLowerIndexSelect<MSCType>();
	//arcs01->AddInputDFO(arcs);
	//arcs01->SetMscDFO(hierarchy);
	//arcs01->SetPassingIndices(true, false, false);

	//DFO_MSCIncidentArcs<MSCType>* inc_arcs = new DFO_MSCIncidentArcs<MSCType>();
	//inc_arcs->AddInputDFO(mins);
	//inc_arcs->SetMscDFO(hierarchy);
	//inc_arcs->SetPassingArcs(false, true);

	//DFO_MSCIncidentNodes<MSCType>* inc_nodes = new DFO_MSCIncidentNodes<MSCType>();
	//inc_nodes->AddInputDFO(inc_arcs);
	//inc_nodes->SetMscDFO(hierarchy);
	//inc_nodes->SetPassingNodes(true, false);

	//DFO_IdToPosition* id_mapper = new DFO_IdToPosition();
	//id_mapper->SetInput(inc_nodes);
	//id_mapper->compute_output();

	//DFO_MSCNodePositionValuator<MSCType>* pos_valuator = new DFO_MSCNodePositionValuator<MSCType>();
	//pos_valuator->SetMscDFO(hierarchy);
	//
	//DFO_IdValueArray<Vec3d>* pos_array = new DFO_IdValueArray<Vec3d>();
	//pos_array->SetValuatorInput(pos_valuator);
	//pos_array->SetMapperInput(id_mapper);
	//
	//DFO_MSCNodeFValValuator<MSCType, MSCType::ScalarType>* fval_valuator = 
	//	new DFO_MSCNodeFValValuator<MSCType, MSCType::ScalarType>();
	//fval_valuator->SetMscDFO(hierarchy);

	//DFO_IdValueArray<MSCType::ScalarType>* fval_array = new DFO_IdValueArray<MSCType::ScalarType>();
	//fval_array->SetValuatorInput(fval_valuator);
	//fval_array->SetMapperInput(id_mapper);

	//fval_array->compute_output();

	//DFO_IdToPosition* arc_id_mapper = new DFO_IdToPosition();
	//arc_id_mapper->SetInput(inc_arcs);
	////arc_id_mapper->compute_output();

	//DFO_MSCArcGeomListValuator<MSCType>* arc_geom_valuator =
	//	new DFO_MSCArcGeomListValuator<MSCType>();
	//arc_geom_valuator->SetMscDFO(hierarchy);
	//
	//DFO_IdValueArray<vector<INDEX_TYPE>>* arc_geom_array = new DFO_IdValueArray<vector<INDEX_TYPE>>();
	//arc_geom_array->SetValuatorInput(arc_geom_valuator);
	//arc_geom_array->SetMapperInput(arc_id_mapper);

	//arc_geom_array->compute_output();

	//DFO_MSCArcLowerFValValuator<MSCType, MSCType::ScalarType>* arc_lower_value_valuator =
	//	new DFO_MSCArcLowerFValValuator<MSCType, MSCType::ScalarType>();
	//arc_lower_value_valuator->SetMscDFO(hierarchy);
	//arc_lower_value_valuator->compute_output();
	//DFO_MSCArcUpperFValValuator<MSCType, MSCType::ScalarType>* arc_upper_value_valuator =
	//	new DFO_MSCArcUpperFValValuator<MSCType, MSCType::ScalarType>();
	//arc_upper_value_valuator->SetMscDFO(hierarchy);
	//arc_upper_value_valuator->compute_output();


	//const vector<INDEX_TYPE>* arcgeoms = arc_geom_array->GetArray();

	return 1;


}


