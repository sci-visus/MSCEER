#ifdef _WIN32
//#include <io.h>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "gi_basic_types.h"
#include "../msc_2d_lib/msc_2d_lib.h"

namespace py = pybind11;
using namespace GInt;

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
	std::unique_ptr<GInt::Msc2D::Msc2D> msc;
	int mX;
	int mY;
	float select_persistence;

};

std::vector<MSCInstance> g_msc_instances;



int MakeMSCInstance() {
	g_msc_instances.push_back(MSCInstance());
	auto msc_id = g_msc_instances.size() - 1;
	MSCInstance& msci = g_msc_instances[msc_id];
	msci.msc.reset(new GInt::Msc2D::Msc2D());
	msci.select_persistence = 0;
	msci.mX = -1;
	msci.mY = -1;
	return msc_id;
}


// User supplies 3 numpy 2d arrays -- laplacian, gaussian, raw
// computes msc on first, samples all and outputs soups and values files
void ComputeMSC(
	int msc_id,
	py::array_t<float> raw,
	bool AccurateASC,
	bool AccurateDSC,
	float base_persistence_abs,
	int num_threads) {

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
	std::vector<float> raw_data((size_t)pX * (size_t)pY, 0.0f);

	float maxval = -std::numeric_limits<float>::infinity();
	float minval = std::numeric_limits<float>::infinity();
	for (py::ssize_t i = 0; i < pX; i++)
		for (py::ssize_t j = 0; j < pY; j++) {
			const float v = ra(i, j);
			raw_data[(size_t)(j + i * mX)] = v;
			if (v > maxval) maxval = v;
			if (v < minval) minval = v;
		}

	GInt::Msc2D::Msc2D::ComputeOptions options;
	const int requested_threads = (num_threads < 1) ? 1 : num_threads;
	options.builderMode = (requested_threads == 1)
		? GInt::Msc2D::Msc2D::BuilderMode::Serial
		: GInt::Msc2D::Msc2D::BuilderMode::Partitioned;
	options.requestedParallelism = requested_threads;
	options.basePersistenceAbs = base_persistence_abs;
	options.cancelPersistenceAbs = -1.0f;
	options.accurateAsc = AccurateASC;
	options.accurateDsc = AccurateDSC;
	msci.msc->compute(raw_data.data(), (int)pX, (int)pY, options);

	const float pers_limit = 0.1f * (maxval - minval);
	msci.select_persistence = pers_limit;
	msci.msc->setPersistence(pers_limit);

	Py_END_ALLOW_THREADS

}

// compute the line graph for the current persistence threshold - the bool specifies valleys or ridges
void ComputePolylineGraph(int msc_id, bool use_valleys) {

	MSCInstance& msci = g_msc_instances[msc_id];

	Py_BEGIN_ALLOW_THREADS
	msci.msc->setPersistence(msci.select_persistence);
	msci.msc->computePolylineGraph(use_valleys);
	Py_END_ALLOW_THREADS

}

std::tuple<std::vector<Node>, std::vector<Edge>> GetGraph(int msc_id) {
	MSCInstance& msci = g_msc_instances[msc_id];
	std::vector<Node> nodes;
	std::vector<Edge> edges;
	const GInt::Msc2D::Graph g = msci.msc->graph();
	for (size_t i = 0; i < g.nodes.size(); ++i) {
		const GInt::Msc2D::Node& src = g.nodes[i];
		Node n;
		n.id = src.id;
		n.edges = src.edges;
		for (size_t pi = 0; pi < src.geometry.size(); ++pi) {
			n.geometry.push_back({ src.geometry[pi].x, src.geometry[pi].y });
		}
		nodes.push_back(n);
	}

	for (size_t i = 0; i < g.edges.size(); ++i) {
		const GInt::Msc2D::Edge& src = g.edges[i];
		Edge e;
		e.id = src.id;
		e.from = src.from;
		e.to = src.to;
		for (size_t pi = 0; pi < src.geometry.size(); ++pi) {
			e.geometry.push_back({ src.geometry[pi].x, src.geometry[pi].y });
		}
		edges.push_back(e);
	}
	return std::make_tuple(nodes, edges);
}

void SetMSCPersistence(int msc_id, float value) {
	MSCInstance& msci = g_msc_instances[msc_id];
	msci.select_persistence = value;
	msci.msc->setPersistence(value);
}

py::array_t<int> GetAsc2Manifolds(int msc_id) {
	MSCInstance& msci = g_msc_instances[msc_id];
	msci.msc->setPersistence(msci.select_persistence);
	const GInt::Msc2D::LabelImage img = msci.msc->ascending2Manifolds();
	py::array_t<int> arr({ img.height, img.width });
	int* ptr = arr.mutable_data();
	for (int i = 0; i < img.width * img.height; i++) {
		ptr[i] = img.labels[(size_t)i];
	}
	return arr;
}
py::array_t<int> GetDsc2Manifolds(int msc_id) {
	MSCInstance& msci = g_msc_instances[msc_id];
	msci.msc->setPersistence(msci.select_persistence);
	const GInt::Msc2D::LabelImage img = msci.msc->descending2Manifolds();
	py::array_t<int> arr({ img.height, img.width });
	int* ptr = arr.mutable_data();
	for (int i = 0; i < img.width * img.height; i++) {
		ptr[i] = img.labels[(size_t)i];
	}
	return arr;
}

py::dict GetCriticalPoints(int msc_id) {
	MSCInstance& msci = g_msc_instances[msc_id];

	py::dict res;
	const std::vector<GInt::Msc2D::CriticalPoint> points = msci.msc->criticalPoints();
	const int num_living = static_cast<int>(points.size());
	
	py::array_t<float> arr_xcoord(num_living);
	float* ptr_xcoord = arr_xcoord.mutable_data();
	py::array_t<float> arr_ycoord(num_living);
	float* ptr_ycoord = arr_ycoord.mutable_data();
	py::array_t<int> arr_index(num_living);
	int* ptr_index = arr_index.mutable_data();
	py::array_t<int> arr_dim(num_living);
	int* ptr_dim = arr_dim.mutable_data();

	py::array_t<float> arr_value(num_living);
	float* ptr_value = arr_value.mutable_data();

	for (int i = 0; i < num_living; ++i) {
		const GInt::Msc2D::CriticalPoint& cp = points[(size_t)i];
		ptr_xcoord[i] = cp.x;
		ptr_ycoord[i] = cp.y;
		ptr_index[i] = cp.id;
		ptr_dim[i] = cp.dim;
		ptr_value[i] = cp.value;
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
		.def_readwrite("id", &Edge::id) 
		.def_readwrite("from_", &Edge::from)
		.def_readwrite("geometry", &Edge::geometry)
		.def_readwrite("to", &Edge::to);

	m.def("GetGraph", &GetGraph, "return a tuple of list of nodes, list of edges. each is a struct with ids and geometry.");
	m.def("MakeMSCInstance", &MakeMSCInstance, "Make an instance of a Morse-Smale complex container");
	m.def("ComputeMSC", &ComputeMSC,
		"Supply an msc id and a 2d float numpy array; computes partitioned MSC by default. If num_threads==1, serial builder is used.",
		py::arg("msc_id"),
		py::arg("raw"),
		py::arg("AccurateASC") = true,
		py::arg("AccurateDSC") = true,
		py::arg("base_persistence_abs") = 0.0f,
		py::arg("num_threads") = 8);
	m.def("SetMSCPersistence", &SetMSCPersistence, "Supply an msc id, set the current persistence to value");
	m.def("GetAsc2Manifolds", &GetAsc2Manifolds, "create the 2d regions (basins) image at current persistence");
	m.def("GetDsc2Manifolds", &GetDsc2Manifolds, "create the 2d regions (mountains) image at current persistence");
	m.def("GetCriticalPoints", &GetCriticalPoints, "get a dictionary of values for living critical points");
	m.def("ComputePolylineGraph", &ComputePolylineGraph, "compute the geometric line graph of the msc ridges or valleys");
}


