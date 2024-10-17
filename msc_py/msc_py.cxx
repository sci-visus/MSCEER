#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif
#include "msc_py.h"

using namespace GInt;

typedef GInt::MorseSmaleComplexBasic<float, Accurate2D::MeshType, Accurate2D::MeshFuncType, Accurate2D::GradType> MyMscType;

struct MSCInstance {
	Accurate2D::DiscreteGradientBuilder* dgb;
	Accurate2D::GridType* grid;
	Accurate2D::GridFuncType* gridfunc;
	Accurate2D::MeshType* mesh;
	Accurate2D::MeshFuncType* meshfunc;
	Accurate2D::GradType* grad;
	MyMscType* msc;
	int mX;
	int mY;
	float* frawdata;
	int* base_labeling_asc2;
	int* base_labeling_dsc2;
	float select_persistence;

};

std::vector<MSCInstance> g_msc_instances;



int MakeMSCInstance() {
	g_msc_instances.push_back(MSCInstance());
	auto msc_id = g_msc_instances.size() - 1;
	MSCInstance& msci = g_msc_instances[msc_id];
	msci.dgb = new Accurate2D::DiscreteGradientBuilder();
	msci.frawdata = NULL;
	msci.base_labeling_asc2 = NULL;
	msci.base_labeling_dsc2 = NULL;
	msci.select_persistence = 0;
	msci.mX = -1;
	msci.mY = -1;
	return msc_id;
}


// User supplies 3 numpy 2d arrays -- laplacian, gaussian, raw
// computes msc on first, samples all and outputs soups and values files
void ComputeMSC(int msc_id, py::array_t<float> raw) {

	MSCInstance& msci = g_msc_instances[msc_id];

	Py_BEGIN_ALLOW_THREADS
		auto ra = raw.unchecked<2>();
	// do dumb copy
	// in python y moves fastest
	auto pX = ra.shape(0);
	auto pY = ra.shape(1);
	auto mX = pY;
	auto mY = pX;
	msci.mX = mX;
	msci.mY = mY;

	//GridType* underlying_grid = new GridType({ la.shape(0), la.shape(1) }, { false, false });

	msci.frawdata = new float[pX * pY];



	//py::print("px:", pX, "py", pY);
	//py::print("mx:", mX, "mY:", mY);
	for (py::ssize_t i = 0; i < pX; i++)
		for (py::ssize_t j = 0; j < pY; j++) {
			msci.frawdata[j + i * mX] = ra(i, j);
		}

	msci.dgb->SetFloadArrayAndDims(mX, mY, msci.frawdata);
	msci.dgb->SetNeededAccuracy(1, 1);
	msci.dgb->SetParallelism(1);

	msci.dgb->ComputeDiscreteGradient();

	msci.grid = msci.dgb->GetGrid();
	msci.gridfunc = msci.dgb->GetGridFunc();
	msci.mesh = msci.dgb->GetTopoMesh();
	msci.meshfunc = msci.dgb->GetMeshFunc();
	msci.grad = msci.dgb->GetGrad();

	
	msci.msc = new MyMscType(msci.grad, msci.mesh, msci.meshfunc);
	msci.msc->SetBuildArcGeometry(Vec3b(true, true, true)); // we only need geometric realizations of 2d saddle-max arcs
	msci.msc->ComputeFromGrad();

	// get persistence to simplify to:
	float maxval = msci.gridfunc->GetMaxValue();
	float minval = msci.gridfunc->GetMinValue();

	float pers_limit = 0.2 * (maxval - minval);
	

	msci.msc->ComputeHierarchy(pers_limit);
	// get persistence to get to the number of maxima requested
	msci.msc->SetSelectPersAbs(pers_limit);

	Py_END_ALLOW_THREADS

}

void SetMSCPersistence(int msc_id, float value) {
	MSCInstance& msci = g_msc_instances[msc_id];
	msci.select_persistence = value;
	msci.msc->SetSelectPersAbs(value);
}

py::array_t<int> GetAsc2Manifolds(int msc_id) {
	MSCInstance& msci = g_msc_instances[msc_id];

	if (msci.base_labeling_asc2 == NULL) {
		msci.base_labeling_asc2 = new int[msci.grid->NumElements()];
		msci.msc->SetSelectPersAbs(-1);
	
		for (int i = 0; i < msci.grid->NumElements(); i++) msci.base_labeling_asc2[i] = -1;
		std::vector<INT_TYPE> nodes;
		std::unordered_map<INT_TYPE, int> nid2graphnodeid;
		INDEX_TYPE count_ids = 0;
		MyMscType::LivingNodesIterator nit(msci.msc);
		for (nit.begin(); nit.valid(); nit.advance()) {
			auto nid = nit.value();
			if (msci.msc->getNode(nid).dim != 0) continue;
			std::set<INDEX_TYPE> manifold;
			msci.msc->fillGeometry(nid, manifold, true);

			for (auto id : manifold) {
				if (msci.mesh->dimension(id) != 0) continue;
				msci.base_labeling_asc2[msci.mesh->VertexNumberFromCellID(id)] = nid;
			}
		}
	}
	msci.msc->SetSelectPersAbs(msci.select_persistence);
	std::unordered_map<INT_TYPE, int> remap;
	MyMscType::LivingNodesIterator nit(msci.msc);
	for (nit.begin(); nit.valid(); nit.advance()) {
		auto nid = nit.value();
		if (msci.msc->getNode(nid).dim != 0) continue;
		std::set<INT_TYPE> constituents;
		msci.msc->GatherNodes(nid, constituents, true);
		for (auto oid : constituents) {
			remap[oid] = nid;
		}
	}
	// now we have a map from the base node id to the representative
	py::print("allocating array");
	py::array_t<int> arr({ msci.mY, msci.mX });
	int* ptr = arr.mutable_data();
	//for (int i = 0; i < msci.mX; i++) {
	//	for (int j = 0; j < msci.mY; j++) {
	//		ptr[j+i*msci.mY] = remap[msci.base_labeling_asc2[msci.grid->Index2d({ i,j })]];
	//	}
	//}
	for (int i = 0; i < msci.mX * msci.mY; i++) {
		ptr[i] = remap[msci.base_labeling_asc2[i]];
	}

	return arr;
}
py::array_t<int> GetDsc2Manifolds(int msc_id) {
	MSCInstance& msci = g_msc_instances[msc_id];

	if (msci.base_labeling_dsc2 == NULL) {
		msci.base_labeling_dsc2 = new int[msci.grid->NumElements()];
		msci.msc->SetSelectPersAbs(-1);

		for (int i = 0; i < msci.grid->NumElements(); i++) msci.base_labeling_dsc2[i] = -1;
		std::vector<INT_TYPE> nodes;
		std::unordered_map<INT_TYPE, int> nid2graphnodeid;
		INDEX_TYPE count_ids = 0;
		MyMscType::LivingNodesIterator nit(msci.msc);
		for (nit.begin(); nit.valid(); nit.advance()) {
			auto nid = nit.value();
			if (msci.msc->getNode(nid).dim != 2) continue;
			std::set<INDEX_TYPE> manifold;
			msci.msc->fillGeometry(nid, manifold, false);

			for (auto id : manifold) {
				if (msci.mesh->dimension(id) != 2) continue;
				msci.base_labeling_dsc2[msci.mesh->VertexNumberFromCellID(id)] = nid;
			}
		}
	}
	msci.msc->SetSelectPersAbs(msci.select_persistence);
	std::unordered_map<INT_TYPE, int> remap;
	MyMscType::LivingNodesIterator nit(msci.msc);
	for (nit.begin(); nit.valid(); nit.advance()) {
		auto nid = nit.value();
		if (msci.msc->getNode(nid).dim != 2) continue;
		std::set<INT_TYPE> constituents;
		msci.msc->GatherNodes(nid, constituents, false);
		for (auto oid : constituents) {
			remap[oid] = nid;
		}
	}
	// now we have a map from the base node id to the representative
	py::print("allocating array");
	py::array_t<int> arr({ msci.mY, msci.mX });
	int* ptr = arr.mutable_data();
	//for (int i = 0; i < msci.mX; i++) {
	//	for (int j = 0; j < msci.mY; j++) {
	//		ptr[j+i*msci.mY] = remap[msci.base_labeling_asc2[msci.grid->Index2d({ i,j })]];
	//	}
	//}
	for (int i = 0; i < msci.mX * msci.mY; i++) {
		ptr[i] = remap[msci.base_labeling_dsc2[i]];
	}
	return arr;
}

py::dict GetCriticalPoints(int msc_id) {
	MSCInstance& msci = g_msc_instances[msc_id];

	py::dict res;
	std::set<INT_TYPE> living_node_ids;
	MyMscType::LivingNodesIterator nit(msci.msc);
	for (nit.begin(); nit.valid(); nit.advance()) {
		auto nid = nit.value();
		living_node_ids.insert(nid);
	}
	int num_living = living_node_ids.size();
	
	py::array_t<float> arr_xcoord({ num_living });
	float* ptr_xcoord = arr_xcoord.mutable_data();
	py::array_t<float> arr_ycoord({ num_living });
	float* ptr_ycoord = arr_ycoord.mutable_data();
	py::array_t<int> arr_index({ num_living });
	int* ptr_index = arr_index.mutable_data();
	py::array_t<int> arr_dim({ num_living });
	int* ptr_dim = arr_dim.mutable_data();

	py::array_t<float> arr_value({ num_living });
	float* ptr_value = arr_value.mutable_data();

	int counter = 0;
	for (auto id : living_node_ids) {
		auto& node = msci.msc->getNode(id);
		GInt::Vec2l coords;
		msci.mesh->cellid2Coords(node.cellindex, coords);
		GInt::Vec2f fcoords = coords;
		fcoords *= 0.5; // to grid indices
		ptr_xcoord[counter] = fcoords[0];
		ptr_ycoord[counter] = fcoords[1];
		ptr_index[counter] = id;
		ptr_dim[counter] = node.dim;
		ptr_value[counter] = node.value;
		counter++;
	}
	res[py::str("id")] = arr_index;
	res[py::str("x")] = arr_xcoord;
	res[py::str("y")] = arr_ycoord;
	res[py::str("dim")] = arr_dim;
	res[py::str("value")] = arr_value;
	return res;
}
PYBIND11_MODULE(msc_py, m) {
	m.def("MakeMSCInstance", &MakeMSCInstance, "Make an instance of a Morse-Smale complex container");
	m.def("ComputeMSC", &ComputeMSC, "Supply an msc id, and a 2d float numpy array, this computes discrete gradient, MSC, and hierarchy up to 20% of range");
	m.def("SetMSCPersistence", &SetMSCPersistence, "Supply an msc id, set the current persistence to value");
	m.def("GetAsc2Manifolds", &GetAsc2Manifolds, "create the 2d regions (basins) image at current persistence");
	m.def("GetDsc2Manifolds", &GetDsc2Manifolds, "create the 2d regions (mountains) image at current persistence");
	m.def("GetCriticalPoints", &GetCriticalPoints, "get a dictionary of values for living critical points");
}


