#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif
#include "msc_test.h"

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
void ComputeMSC(int msc_id, float* raw, int pX, int pY) {

	MSCInstance& msci = g_msc_instances[msc_id];

	auto mX = pY;
	auto mY = pX;
	msci.mX = mX;
	msci.mY = mY;

	//GridType* underlying_grid = new GridType({ la.shape(0), la.shape(1) }, { false, false });

	msci.frawdata = new float[pX * pY];

	for (int i = 0; i < pX; i++)
		for (int j = 0; j < pY; j++) {
			msci.frawdata[j + i * mX] = raw[i + j * pX];
		}

	msci.dgb->SetFloadArrayAndDims(mX, mY, msci.frawdata);
	msci.dgb->SetNeededAccuracy(1, 1);
	msci.dgb->SetParallelism(8);

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

}

void SetMSCPersistence(int msc_id, float value) {
	MSCInstance& msci = g_msc_instances[msc_id];
	msci.select_persistence = value;
	msci.msc->SetSelectPersAbs(value);
}

int* GetMyAsc2Manifolds(int msc_id) {
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
	
	int* arr = new int[ msci.mY * msci.mX ];
	int* ptr = arr;
	for (int i = 0; i < msci.mX; i++) {
		for (int j = 0; j < msci.mY; j++) {
			ptr[j+i*msci.mY] = remap[msci.base_labeling_asc2[msci.grid->Index2d({ i,j })]];
		}
	}
	return arr;
}

int main() {
	auto id = MakeMSCInstance();

	FILE* fin = fopen("C:\\Users\\jediati\\Desktop\\JEDIATI\\data\\prostate_he_1255x778.raw", "rb");
	float* data = new float[1255 * 778];
	fread(data, sizeof(float), 1255 * 778, fin);
	fclose(fin);

	ComputeMSC(id, data, 778, 1255);

	GetMyAsc2Manifolds(id);

	return 1;

}


