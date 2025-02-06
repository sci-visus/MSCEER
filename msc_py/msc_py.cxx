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


void LineAlg(int x1, int y1, int x2, int y2, std::vector<std::pair<int, int>>& line)
{
	int x, y, dx, dy, dx1, dy1, px, py, xe, ye, i;
	dx = x2 - x1;
	dy = y2 - y1;
	dx1 = abs(dx);
	dy1 = abs(dy);
	px = 2 * dy1 - dx1;
	py = 2 * dx1 - dy1;
	if (dy1 <= dx1)
	{
		if (dx >= 0)
		{
			x = x1;
			y = y1;
			xe = x2;
		}
		else
		{
			x = x2;
			y = y2;
			xe = x1;
		}
		line.push_back(std::pair<int, int>(x, y));
		for (i = 0; x < xe; i++)
		{
			x = x + 1;
			if (px < 0)
			{
				px = px + 2 * dy1;
			}
			else
			{
				if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0))
				{
					y = y + 1;
				}
				else
				{
					y = y - 1;
				}
				px = px + 2 * (dy1 - dx1);
			}
			line.push_back(std::pair<int, int>(x, y));
		}
	}
	else
	{
		if (dy >= 0)
		{
			x = x1;
			y = y1;
			ye = y2;
		}
		else
		{
			x = x2;
			y = y2;
			ye = y1;
		}
		line.push_back(std::pair<int, int>(x, y));
		for (i = 0; y < ye; i++)
		{
			y = y + 1;
			if (py <= 0)
			{
				py = py + 2 * dx1;
			}
			else
			{
				if ((dx < 0 && dy < 0) || (dx > 0 && dy > 0))
				{
					x = x + 1;
				}
				else
				{
					x = x - 1;
				}
				py = py + 2 * (dx1 - dy1);
			}
			line.push_back(std::pair<int, int>(x, y));
		}
	}
}

class RegionGraph {
public:
	struct RGNode {
		INT_TYPE id;  // External ID
		std::vector<INT_TYPE> edges;
	};

	struct RGEdge {
		INT_TYPE id;  // External ID
		size_t n1;    // Zero-based index of node 1 in m_nodes
		size_t n2;    // Zero-based index of node 2 in m_nodes
		std::vector<INDEX_TYPE> cells;
	};

protected:
	std::vector<RGNode> m_nodes;
	std::unordered_map<INT_TYPE, size_t> m_id_map;  // Maps external IDs to vector indices
	std::vector<RGEdge> m_edges;                    // Vector of edges

public:
	RegionGraph() {}

	void reset() {
		m_nodes.clear();
		m_id_map.clear();
		m_edges.clear();
	}

	// Add a node with an external ID and map it to the internal index
	void addNode(INT_TYPE external_id) {
		if (m_id_map.find(external_id) == m_id_map.end()) {
			RGNode newNode;
			newNode.id = external_id;
			m_nodes.push_back(newNode);
			m_id_map[external_id] = m_nodes.size() - 1;
		}
	}

	// Add an edge with external node IDs and map them to internal indices
	void addEdge(INT_TYPE edge_id, INT_TYPE external_id1, INT_TYPE external_id2) {
		if (m_id_map.find(external_id1) != m_id_map.end() && m_id_map.find(external_id2) != m_id_map.end()) {
			RGEdge newEdge;
			newEdge.id = edge_id;
			newEdge.n1 = m_id_map[external_id1];  // Zero-based index of node 1
			newEdge.n2 = m_id_map[external_id2];  // Zero-based index of node 2
			m_edges.push_back(newEdge);

			// Optionally, update the edges list of the nodes
			m_nodes[newEdge.n1].edges.push_back(edge_id);
			m_nodes[newEdge.n2].edges.push_back(edge_id);
		}
		else {
			// Handle the case where one or both nodes do not exist
			throw std::runtime_error("One or both nodes not found");
		}
	}

	// Get a reference to a node by external ID
	RGNode& getNode(INT_TYPE external_id) {
		return m_nodes[m_id_map.at(external_id)];
	}
};

struct Node {
	int id;
	std::vector<std::pair<float, float>> geometry; // List of 2D points
	std::vector<int> edges; // Connected node IDs
	Node() {}
	Node(int id_, std::vector<std::pair<float, float>> geometry_, std::vector<int> edges_)
		: id(id_), geometry(std::move(geometry_)), edges(std::move(edges_)) {}
};

struct Edge {
	int id;
	int from, to;
	std::vector<std::pair<float, float>> geometry;
	Edge() {}
	Edge(int f, int t) : from(f), to(t) {}
};



struct MSCInstance {
	Accurate2D::DiscreteGradientBuilder* dgb;
	Accurate2D::GridType* grid;
	Accurate2D::GridFuncType* gridfunc;
	Accurate2D::MeshType* mesh;
	Accurate2D::MeshFuncType* meshfunc;
	Accurate2D::GradType* grad;
	MyMscType* msc;
	GInt::Geometric2DGraph* geom_line_graph;
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
	msci.geom_line_graph = NULL;
	msci.select_persistence = 0;
	msci.mX = -1;
	msci.mY = -1;
	return msc_id;
}


// User supplies 3 numpy 2d arrays -- laplacian, gaussian, raw
// computes msc on first, samples all and outputs soups and values files
void ComputeMSC(int msc_id, py::array_t<float> raw, bool AccurateASC, bool AccurateDSC) {

	MSCInstance& msci = g_msc_instances[msc_id];

	if (msci.frawdata != NULL) delete[] msci.frawdata;
	if (msci.base_labeling_dsc2 != NULL) delete[] msci.base_labeling_dsc2;
	if (msci.base_labeling_asc2 != NULL) delete[] msci.base_labeling_asc2;
	msci.frawdata = NULL;
	msci.base_labeling_asc2 = NULL;
	msci.base_labeling_dsc2 = NULL;


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
	msci.dgb->SetNeededAccuracy(AccurateASC, AccurateDSC);
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

// compute the line graph for the current persistence threshold - the bool specifies valleys or ridges
void ComputePolylineGraph(int msc_id, bool use_valleys) {

	MSCInstance& msci = g_msc_instances[msc_id];

	if (msci.geom_line_graph != NULL) delete msci.geom_line_graph;

	Py_BEGIN_ALLOW_THREADS

	MeshCellsGraph* graph;
	if (use_valleys) {
		graph = GInt::BuildMeshCellsGraphFromMSCValleys<MyMscType, MeshType>(msci.msc, msci.mesh);
	}
	else {
		graph = GInt::BuildMeshCellsGraphFromMSCRidges<MyMscType, MeshType>(msci.msc, msci.mesh);
	}
	msci.geom_line_graph = GInt::BuildGeometricGraphFromMeshGraph<MeshType>(graph, msci.mesh, 10);
	delete graph;

	Py_END_ALLOW_THREADS

		// now attach to msci

}

std::tuple<std::vector<Node>, std::vector<Edge>> GetGraph(int msc_id) {
	MSCInstance& msci = g_msc_instances[msc_id];
	
	std::vector<Node> nodes;
	GInt::Geometric2DGraph::vertex_iterator vit(msci.geom_line_graph);
	for (vit.begin(); vit.valid(); vit.advance()) {
		auto vid = vit.value();
		const auto& v = msci.geom_line_graph->GetVertex(vid);

		Node n;
		n.id = v.vid;
		n.edges.insert(n.edges.begin(), v.edges.begin(), v.edges.end());
		n.geometry.push_back({ v.store[0], v.store[1] });
		nodes.push_back(n);
	}
	


	std::vector<Edge> edges;
	GInt::Geometric2DGraph::edge_iterator eit(msci.geom_line_graph);
	for (eit.begin(); eit.valid(); eit.advance()) {
		auto eid = eit.value();
		const auto& ge = msci.geom_line_graph->GetEdge(eid);

		Edge e;
		e.id = ge.eid;
		e.from = ge.v1;
		e.to = ge.v2;
		for (const auto& p : ge.store->GetLine()) {
			e.geometry.push_back({ p[0], p[1] });
		}
		edges.push_back(e);
	}
	return std::make_tuple(nodes, edges);
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
	py::class_<Node>(m, "Node")
		.def(py::init<int, std::vector<std::pair<float, float>>, std::vector<int>>())
		.def_readwrite("id", &Node::id)
		.def_readwrite("geometry", &Node::geometry)
		.def_readwrite("edges", &Node::edges);

	py::class_<Edge>(m, "Edge")
		.def(py::init<int, int>())
		.def_readwrite("from_", &Edge::from)
		.def_readwrite("geometry", &Edge::geometry)
		.def_readwrite("to", &Edge::to);

	m.def("GetGraph", &GetGraph, "return a tuple of list of nodes, list of edges. each is a struct with ids and geometry.");
	m.def("MakeMSCInstance", &MakeMSCInstance, "Make an instance of a Morse-Smale complex container");
	m.def("ComputeMSC", &ComputeMSC, "Supply an msc id, and a 2d float numpy array, this computes discrete gradient, MSC, and hierarchy up to 20% of range");
	m.def("SetMSCPersistence", &SetMSCPersistence, "Supply an msc id, set the current persistence to value");
	m.def("GetAsc2Manifolds", &GetAsc2Manifolds, "create the 2d regions (basins) image at current persistence");
	m.def("GetDsc2Manifolds", &GetDsc2Manifolds, "create the 2d regions (mountains) image at current persistence");
	m.def("GetCriticalPoints", &GetCriticalPoints, "get a dictionary of values for living critical points");
	m.def("ComputePolylineGraph", &ComputePolylineGraph, "compute the geometric line graph of the msc ridges or valleys");
}


