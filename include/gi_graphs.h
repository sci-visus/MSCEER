/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef GI_GRAPHS_H
#define GI_GRAPHS_H

#include <vector>
#include "gi_basic_types.h"
#include "gi_basic_geometry.h"
#include "gi_array_index_partition.h"

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <omp.h>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <iterator>
// This file is a safe space to put experimental classes and functionality.

namespace GInt {




	// these are the methods and iterators used for static polymorphism - in the templated implementation of a graph traversal algorithm
	// any graph MUST implement ALL these methods with the same function prototype,
	// as well as have a vertex_iterator (iterate over all nodes of the graph) and neighbor_iterator (iterate over all nodes adjacent to a 
	// given node). furthermore, all function prototypes must match as well
	class GraphInterface {
	public:
		//INDEX_TYPE NumElements() { return 0; }


		class vertex_iterator {
			GraphInterface* mG;
		public:
			vertex_iterator(GraphInterface* g) : mG(g) {}
			void begin(INDEX_TYPE a) {}
			void advance() {}
			bool valid() { return false; }
			INDEX_TYPE value() const { return 0; }
		};
		class neighbor_iterator {
			GraphInterface* mG;
		public:
			neighbor_iterator(GraphInterface* g) : mG(g) {}
			void begin(INDEX_TYPE a) {}
			void advance() {}
			bool valid() { return false; }
			INDEX_TYPE value() const { return 0; }
		};
	};

	template<class VertexStore, class EdgeStore>
	class UndirectedGraph {
	public:
		struct vertex {
			vertex() {}
/*			bool operator==(const vertex& other) {
				return vid == other.vid;
			}	*/		
			INDEX_TYPE vid;
			std::vector<INDEX_TYPE> edges;
			VertexStore store;
			vertex(INDEX_TYPE id, VertexStore s) : vid(id), store(s) {}
		};

		struct edge {
			edge() {}
			//bool operator==(const edge& other) {
			//	return eid == other.eid && v1 == other.v1 && v2 == other.v2;
			//}
			edge(INDEX_TYPE e, INDEX_TYPE a, INDEX_TYPE b, EdgeStore s) :
				eid(e), v1(a), v2(b), store(s) {}
			INDEX_TYPE eid;
			INDEX_TYPE v1;
			INDEX_TYPE v2;
			EdgeStore store;
			INDEX_TYPE other_vertex(INDEX_TYPE v) {
				if (v1 == v) return v2;
				return v1;
			}
		};
	protected:

		std::unordered_map<INDEX_TYPE, vertex> mVertices;
		std::unordered_map<INDEX_TYPE, edge> mEdges;
		INDEX_TYPE m_edgecount;
	public:
		//INDEX_TYPE NumElements() { return 0; }
		UndirectedGraph() : m_edgecount(0) {}

		void AddVertex(INDEX_TYPE id, VertexStore s) {
			mVertices[id] = vertex(id, s);
		}
		int AddEdge(INDEX_TYPE vid1, INDEX_TYPE vid2, EdgeStore s) {
			edge e(m_edgecount, vid1, vid2, s );
			mEdges[m_edgecount] = e;

			if (mVertices.count(vid1) == 0 || mVertices.count(vid2) == 0) {
				printf("Warning: AddEdge -- count of verts %d = %d, count of verts %d = %d\n",
					vid1, mVertices.count(vid1),
					vid2, mVertices.count(vid2));
			}

			mVertices[vid1].edges.push_back(m_edgecount);
			mVertices[vid2].edges.push_back(m_edgecount);
			m_edgecount++;
			return m_edgecount - 1;
		}

		void CheckConnections() {
			int num_doubles = 0;
			for (auto& np : mVertices) {
				std::unordered_map<INDEX_TYPE, INDEX_TYPE> vertices;
				for (auto eid : np.second.edges) {
					auto ovid = OtherVertex(eid, np.first);
					if (vertices.count(ovid) != 0) {
						printf("double connection between %lld -- %lld\n", np.first, ovid);
						if (vertices[ovid] == eid) {
							for (auto& v : ((Line2d*)GetEdgeData(eid))->GetLine()) {
								v.PrintFloat(); printf(" ");
							}
							printf("\nother:\n");
							for (auto& v : ((Line2d*)GetEdgeData(vertices[ovid]))->GetLine()) {
								v.PrintFloat(); printf(" ");
							}
							printf("\n");
						}
						else {
							printf("eh, the eid's don't match: %d %d\n", vertices[ovid], eid);
						}
						num_doubles++;
					
					

						
						
					}
					vertices[ovid] = eid;
				}
			}
			if (num_doubles > 0) {
				printf("found %d multi arcs out of %d\n", num_doubles, mEdges.size());
			}
		}


		VertexStore GetVertexData(INDEX_TYPE vid) { return mVertices[vid].store; }
		EdgeStore GetEdgeData(INDEX_TYPE eid) { return mEdges[eid].store; }
		const edge& GetEdge(INDEX_TYPE eid) { return mEdges[eid]; }
		const vertex& GetVertex(INDEX_TYPE vid) { return mVertices[vid]; }

		INDEX_TYPE OtherVertex(INDEX_TYPE eid, INDEX_TYPE vid) {
			auto& e = mEdges[eid];
			if (e.v1 == vid) return e.v2;
			return e.v1;
		}

		int NumVertices() const { return mVertices.size(); }
		int NumEdges() const { return mEdges.size(); }

		class vertex_iterator {
		protected:
			UndirectedGraph<VertexStore, EdgeStore>* mG;
			typename std::unordered_map<INDEX_TYPE, vertex>::iterator it;
		public:
			vertex_iterator(UndirectedGraph<VertexStore, EdgeStore>* g) : mG(g) {}
			void begin() {
				it = mG->mVertices.begin();
			}
			void advance() {
				it++;
			}
			bool valid() { return it != mG->mVertices.end(); }
			INDEX_TYPE value() const { return (*it).first; }
		};
		class neighbor_iterator {
		protected:
			UndirectedGraph<VertexStore, EdgeStore>* mG;
			std::vector<INDEX_TYPE>::iterator it;
			std::vector<INDEX_TYPE>::iterator endit;
		public:
			neighbor_iterator(UndirectedGraph<VertexStore, EdgeStore>* g) : mG(g) {}
			void begin(INDEX_TYPE a) { 
				auto& e = mG->mVertices[a].edges;
				it = e.begin(); 
				endit = e.end();
			}
			void advance() { it++;  }
			bool valid(INDEX_TYPE a) { return it != endit; }
			INDEX_TYPE value() const { return *it; }
		};
		class edge_iterator {
		protected:
			UndirectedGraph<VertexStore, EdgeStore>* mG;
			typename std::unordered_map<INDEX_TYPE, edge>::iterator it;
		public:
			edge_iterator(UndirectedGraph<VertexStore, EdgeStore>* g) : mG(g) {}
			void begin() {
				it = mG->mEdges.begin();
			}
			void advance() {
				it++;
			}
			bool valid() { return it != mG->mEdges.end(); }
			INDEX_TYPE value() const { return (*it).first; }
		};
	};

	// a graph type for representing 1-skeletal subsets of the cells of a mesh
	// for example: the quads and hexes of a ridge-like structure
	typedef UndirectedGraph<INDEX_TYPE, std::vector<INDEX_TYPE>*> MeshCellsGraph;

	// now a graph to store pure geometry. vertices are 3d points, edges
	// store lines
	typedef UndirectedGraph<Vec3f, Line3d*> Geometric3DGraph;
	typedef UndirectedGraph<Vec2f, Line2d*> Geometric2DGraph;

	template<class MeshType>
	UndirectedGraph<
		typename MeshType::FCoordType, 
		Line<typename MeshType::FCoordType >* >*
        BuildGeometricGraphFromMeshGraph(MeshCellsGraph* in_graph, MeshType* mesh, int smooth=0) {
        //printf("i get called\n");
		UndirectedGraph<typename MeshType::FCoordType, Line<typename MeshType::FCoordType>* >* out_graph =
			new UndirectedGraph<typename MeshType::FCoordType, Line<typename MeshType::FCoordType>* >();
		int c1 = 0;
		MeshCellsGraph::vertex_iterator vit(in_graph);
		for (vit.begin(); vit.valid(); vit.advance()) {
			INDEX_TYPE cid = in_graph->GetVertexData(vit.value());
			typename MeshType::ICoordType coords;
			mesh->cellid2Coords(cid, coords);
			//printf("%llu -> ", cid); coords.PrintInt();
			typename MeshType::FCoordType dc = coords * 0.5f;
			//printf("%llu -> ", cid); dc.PrintFloat();
			out_graph->AddVertex(cid, dc);
		}
		//MeshCellsGraph::edge_iterator eit(in_graph);
		//for (eit.begin(); eit.valid(); eit.advance()) {
		// HUGE HACK
		for (int eid = 0; eid < in_graph->NumEdges(); eid++) {
			auto &ev = in_graph->GetEdge(eid); // get vector of indexes of this line
			ConstrainedLine<typename MeshType::FCoordType>* line = new ConstrainedLine<typename MeshType::FCoordType>();
			for (auto id : *(ev.store)) {
				typename MeshType::ICoordType coords;
				mesh->cellid2Coords(id, coords);
				typename MeshType::FCoordType dc = coords * 0.5f;
				line->AddToEnd(dc);
			}
			line->IrreversableSmooth(smooth);
			int edgenum = out_graph->AddEdge(ev.v1, ev.v2, line);
			if (edgenum != ev.eid) {
				printf("WHOA- geom graph %d != %d\n", edgenum, ev.eid);
			}
		}
		return out_graph;
	}



	// this class builds graph, wherethe faccet/cofacet operator of the underlying
	// mesh is used to connect cells that have been added. degree-2 cells are interior
	// to the graph edges, and others become vertices. 
	// additionally, it is possible to specify cells that are forced to become vertices
	template <class MeshType>
	class MeshCellsGraphBuilder {

	protected:
		MeshType* m_mesh;
		std::unordered_map<INDEX_TYPE, int> m_cell_counter; // # of existing cofacets/facents for the INDEX_TYPE
		std::vector<INDEX_TYPE> m_critical;
	public:
		typedef GInt::MeshCellsGraph GraphType;

		void AddToForcedVertices(INDEX_TYPE cellid) { m_critical.push_back(cellid); }
		void AddCellIndex(INDEX_TYPE cellid) {
			m_cell_counter[cellid] = 0; // set to zero because we are counting the number of exiting adjacent cofacet/facets - we don't know this yet
		}
		MeshCellsGraphBuilder(MeshType* mesh) : m_mesh(mesh) {}
	protected:

		void fill_path(INDEX_TYPE startid, INDEX_TYPE nextid, std::vector<INDEX_TYPE>& path) {
			//printf("fill path ccalled\n");
			path.push_back(startid);
			//printf("%lld\n", startid);
			while (m_cell_counter[nextid] == 2) {
				//printf("%lld\n", nextid);
				path.push_back(nextid);
				bool skip = false;
				typename MeshType::FacetsIterator fit(m_mesh);
				for (fit.begin(nextid); fit.valid(); fit.advance()) {
					INDEX_TYPE facetid = fit.value();
					if (m_cell_counter.count(facetid) != 0 && facetid != path[path.size() - 2]) {
						nextid = facetid;
						skip = true;
						break;
					}
				}
				if (skip) continue;
				typename MeshType::CofacetsIterator cfit(m_mesh);
				for (cfit.begin(nextid); cfit.valid(); cfit.advance()) {
					INDEX_TYPE cofacetid = cfit.value();
					if (m_cell_counter.count(cofacetid) != 0 && cofacetid != path[path.size() - 2]) {
						nextid = cofacetid;
						break;
					}
				}
			}

			//printf("%lld\n", nextid);
			path.push_back(nextid);
			//if (path.size() < 3) printf("path filled %d\n", path.size());

		}
	public:
		GraphType* ComputeStuff() {
			//printf("counting stuff\n");
			
			// try parallel:
#pragma omp parallel
			{			
				INDEX_TYPE t_num_elements = m_cell_counter.bucket_count();
				std::vector<INDEX_TYPE> t_partition;
				int num_threads = omp_get_num_threads();
				ArrayIndexPartitioner::EvenChunkSplit(t_num_elements, num_threads, t_partition);
				int thread_num = omp_get_thread_num();
				INDEX_TYPE num_cells_computed = 0;
				// parallelize over the buckets in the unordered_map
				for (auto bucket_num = t_partition[thread_num]; bucket_num < t_partition[thread_num + 1]; bucket_num++) {
					for (auto idpair_list_iterator = m_cell_counter.begin(bucket_num); idpair_list_iterator != m_cell_counter.end(bucket_num); idpair_list_iterator++) {
						auto& idpair = *idpair_list_iterator;
						num_cells_computed++;
						const INDEX_TYPE& cellid = idpair.first;
						typename MeshType::FacetsIterator fit(m_mesh);
						for (fit.begin(cellid); fit.valid(); fit.advance()) {
							INDEX_TYPE facetid = fit.value();
							idpair.second += m_cell_counter.count(facetid);
						}
						typename MeshType::CofacetsIterator cfit(m_mesh);
						for (cfit.begin(cellid); cfit.valid(); cfit.advance()) {
							INDEX_TYPE cofacetid = cfit.value();
							idpair.second += m_cell_counter.count(cofacetid);
						}

					}
				}
#pragma omp critical
				{
					printf("thread %d did %llu cells\n", thread_num, num_cells_computed);
				}
			}
			printf("counted incience\n");
			// count all incidence between cells

			//for (auto& idpair : m_cell_counter) {
			//	const INDEX_TYPE& cellid = idpair.first;
			//	typename MeshType::FacetsIterator fit(m_mesh);
			//	for (fit.begin(cellid); fit.valid(); fit.advance()) {
			//		INDEX_TYPE facetid = fit.value();
			//		idpair.second += m_cell_counter.count(facetid);
			//	}
			//	typename MeshType::CofacetsIterator cfit(m_mesh);
			//	for (cfit.begin(cellid); cfit.valid(); cfit.advance()) {
			//		INDEX_TYPE cofacetid = cfit.value();
			//		idpair.second += m_cell_counter.count(cofacetid);
			//	}
			//}
			for (auto cid : m_critical) { m_cell_counter[cid] = 1; }
			GraphType* m_res;
			m_res = new GraphType();
			// run a flood fill 
			// -- change the counter to -1 for a degree-2 cell that has been traversed
			// -- change the counter to -2 for a non-degree-2 cell that has been made critical
			printf("now adding vertices\n");
			int c1 = 0;
			for (auto& idpair : m_cell_counter) {
				const INDEX_TYPE& cellid = idpair.first;

				if (idpair.second < 0) {
					// this has been added already, so skip
					continue;
				}

				if (idpair.second == 2) {
					// this will be added later, so we can skip
					continue;
				}

				// so this is a NODE in the graph, and we start a traversal

				m_res->AddVertex(cellid, cellid); c1++;
				idpair.second = -2;

			}
			printf("now starting sweep...\n");
			// now all node have been added to graph
			// for each node, build star of lines starting at that node. mark each line so the same arc is not added twice
			GraphType::vertex_iterator vit(m_res);
			int c2 = 0;
			for (vit.begin(); vit.valid(); vit.advance()) {
				INDEX_TYPE cellid = vit.value();

				// start traversals;


				// now do UF on degree two stuff
				typename MeshType::FacetsIterator fit(m_mesh);
				for (fit.begin(cellid); fit.valid(); fit.advance()) {
					INDEX_TYPE facetid = fit.value();
					// skip cells that are not part of the structure or --have been seen--
					if (m_cell_counter.count(facetid) == 0 /*|| m_cell_counter[facetid] == -1*/) continue;

					// build a line
					std::vector<INDEX_TYPE> line;
					fill_path(cellid, facetid, line);
					//for (int i = 1; i < line.size() - 1; i++) m_cell_counter[line[i]] = -1; // mark the line so it does not get added twice
					//		
					// only add in down index direction
					if (line.front() < line.back()) {
						// add to graph
						std::vector<INDEX_TYPE>* toadd = new std::vector<INDEX_TYPE>(line.begin(), line.end());
						//printf("adding edge: %d\n", toadd->size());
						m_res->AddEdge(toadd->front(), toadd->back(), toadd); c2++;
					}
				}
				typename MeshType::CofacetsIterator cfit(m_mesh);
				for (cfit.begin(cellid); cfit.valid(); cfit.advance()) {
					INDEX_TYPE cofacetid = cfit.value();
					// skip cells that are not part of the structure or have been seen
					if (m_cell_counter.count(cofacetid) == 0 || m_cell_counter[cofacetid] == -1) continue;

					// build a line
					std::vector<INDEX_TYPE> line;
					fill_path(cellid, cofacetid, line);
					//for (int i = 1; i < line.size() - 1; i++) m_cell_counter[line[i]] = -1; // mark the line so it does not get added twice
																							// add to graph
					if (line.front() < line.back()) {
						std::vector<INDEX_TYPE>* toadd = new std::vector<INDEX_TYPE>(line.begin(), line.end());
						//printf("adding edge: %d\n", toadd->size());
						m_res->AddEdge(toadd->front(), toadd->back(), toadd); c2++;
					}
				}
			}
			printf("added %d verts and %d edges\n", c1, c2);
			return m_res;
		}

	};
	//class GraphMergeTree {

	//protected:
	//	UndirectedGraph<int, int>* mGraph;

	//public:

	//	struct IdComparer {
	//		IdComparer();
	//		bool operator()(INDEX_TYPE a, INDEX_TYPE b) {
	//			return false;// graph->Before(a, b);
	//		}
	//	};
	//	GraphMergeTree() {}

	//	bool Before(INDEX_TYPE id1, INDEX_TYPE id2) { return id1 < id2;  }
	//	float GetValue(INDEX_TYPE id) { return 0; }

	//	void ComputeMergeTree() {

	//		//priority_queue<INDEX_TYPE, vector<INDEX_TYPE>, Before> sorted_order;
	//		//map<INDEX_TYPE, INDEX_TYPE> myUF;

	//		//UndirectedGraph<int, int>::vertex_iterator vit(mGraph);
	//		//for (vit.begin(); vit.valid(); vit.advance()) {
	//		//	INDEX_TYPE id = vit.value();
	//		//	sorted_order.push(pair<float, INDEX_TYPE>(GetValue(id), id));
	//		//}

	//		//while (!sorted_order.empty()) {
	//		//	INDEX_TYPE cid = sorted_order.top().second;
	//		//}

	//	}



	//};


	class GraphForMLLabeling
	{
	protected:
		std::string m_basename;

	public:

		typedef int MLGGeomIdx;
		typedef int MLGNodeIdx;

		struct MLGEdge {
			MLGGeomIdx edge_geom_idx;		// the GEOMETRY index - id in to the geom array
			MLGNodeIdx node_1_idx;		// the NODE index - the id into the nodes array
			MLGNodeIdx node_2_idx;		// the NODE index - the id into the nodes array
			bool operator!=(const MLGEdge& other) const {
				return edge_geom_idx != other.edge_geom_idx ||
					node_1_idx != other.node_1_idx ||
					node_2_idx != other.node_2_idx;
			}
		};
		struct GeomFeat {
			int geom_feat_id;		// should be the index in to the geom array which this element appears
			int dim;		// the implicit dimension of the set of points
			std::vector<GInt::Vec2f> points;

			bool operator!=(const GeomFeat& other) const {
				if (geom_feat_id != other.geom_feat_id ||
					dim != other.dim ||
					points.size() != other.points.size()) return true;
				//for (int i = 0; i < points.size(); i++) {
				//	if (points[i] != other.points[i]) {
				//		points[i].PrintFloat();
				//		other.points[i].PrintFloat();
				//		return true;
				//	}
				//}
				return false;
			}
		};

		std::vector<MLGGeomIdx> MLG_nodes; // a list of the indices in MLG_geometric_features that are "nodes"
		std::vector<MLGEdge> MLG_edges; // a list of triples, index in MLG_geometric_features of "edge" geometry, indices in MLG_nodes arrays of two endpoints
		std::vector<GeomFeat> MLG_geometric_features; // the actual geometric features
		std::unordered_map<MLGGeomIdx, MLGNodeIdx> MLG_geomidx_2_nodeidx; // if a geometric feature is a node, this gives its index in the node list

		void PrintComparison(GraphForMLLabeling* other) {
			printf("COMPARING GRAPHS -----> \n");
			if (other->MLG_nodes.size() != this->MLG_nodes.size()) {
				printf(" --> ERROR: nodes sizes don't match: %d -- %d\n", this->MLG_nodes.size(), other->MLG_nodes.size());
			}
			else {
				printf(" ++ graph nodes sizes match, comparing nodes\n");
				for (int i = 0; i < MLG_nodes.size(); i++) {
					if (MLG_nodes[i] != other->MLG_nodes[i]) {
						printf("  --> error node[%d] does not match %d -- %d\n", i, MLG_nodes[i], other->MLG_nodes[i]);
					}
				}
				printf(" ++ done\n");
			}
			if (other->MLG_edges.size() != this->MLG_edges.size()) {
				printf(" --> ERROR: edges sizes don't match: %d -- %d\n", this->MLG_edges.size(), other->MLG_edges.size());
			}
			else {
				printf(" ++ graph edges sizes match, comparing edges\n");
				for (int i = 0; i < MLG_edges.size(); i++) {
					if (MLG_edges[i] != other->MLG_edges[i]) {
						printf("  --> error edge[%d] does not match (%d, %d, %d) -- (%d, %d, %d)\n", i, 
							MLG_edges[i].edge_geom_idx, MLG_edges[i].node_1_idx, MLG_edges[i].node_2_idx,
							other->MLG_edges[i].edge_geom_idx, other->MLG_edges[i].node_1_idx, other->MLG_edges[i].node_2_idx);
					}
				}
				printf(" ++ done\n");
			}
			if (other->MLG_geometric_features.size() != this->MLG_geometric_features.size()) {
				printf(" --> ERROR: geoms sizes don't match: %d -- %d\n", this->MLG_geometric_features.size(), other->MLG_geometric_features.size());
			}
			else {
				printf(" ++ graph geoms sizes match, comparing geoms\n");
				for (int i = 0; i < MLG_geometric_features.size(); i++) {
					if (MLG_geometric_features[i] != other->MLG_geometric_features[i]) {
						printf("  --> error geom[%d] does not match (%d, %d, %d) -- (%d, %d, %d)\n", i,
							MLG_geometric_features[i].geom_feat_id, MLG_geometric_features[i].dim, MLG_geometric_features[i].points.size(),
							other->MLG_geometric_features[i].geom_feat_id, other->MLG_geometric_features[i].dim, other->MLG_geometric_features[i].points.size());
					}
				}
				printf(" ++ done\n");
			}
			if (other->MLG_geomidx_2_nodeidx.size() != this->MLG_geomidx_2_nodeidx.size()) {
				printf(" --> ERROR: geomsmaps sizes don't match: %d -- %d\n", this->MLG_geomidx_2_nodeidx.size(), other->MLG_geomidx_2_nodeidx.size());
			}
			else {
				printf(" ++ graph geomsmaps sizes match, comparing geoms\n");
				for (const auto& p : this->MLG_geomidx_2_nodeidx) {
					if (other->MLG_geomidx_2_nodeidx.count(p.first) == 0) {
						printf("  --> error other map  does not have key,value %d, %d\n", p.first, p.second);
					} else if (other->MLG_geomidx_2_nodeidx[p.first] != p.second) {
						printf("  --> error other[%d] does not match %d -- %d\n", p.first, p.second, other->MLG_geomidx_2_nodeidx[p.first]);
					}
				}
				printf(" ++ done\n");
			}
			printf("done\n");



		}

		int CreateGeomFeat(int gfid, int dim, const std::vector<GInt::Vec2f>& points) {
			int num = MLG_geometric_features.size();
			MLG_geometric_features.push_back(GeomFeat());
			auto& g = MLG_geometric_features.back();
			g.geom_feat_id = gfid;
			g.dim = dim;
			g.points = points;
			return num;
		}

		int CreateNode(int gfid) {
			int num = MLG_nodes.size();
			MLG_geomidx_2_nodeidx[gfid] = num;
			MLG_nodes.push_back(gfid);
			return num;
		}

		int CreateEdge(int egid, int nid1, int nid2) {
			int num = MLG_edges.size();
			MLG_edges.push_back({ egid, nid1, nid2 });
			return num;
		}



		void CreateFrom2DPolylineGraph(GInt::Geometric2DGraph* plg) {
			MLG_nodes.clear();
			MLG_edges.clear();
			MLG_geometric_features.clear();
			MLG_geomidx_2_nodeidx.clear();

			std::unordered_map<int, MLGGeomIdx> polylineid_2_geomidx;
			// each polyline will become a node

			GInt::Geometric2DGraph::edge_iterator eit(plg);
			for (eit.begin(); eit.valid(); eit.advance()) {
				auto* line = plg->GetEdgeData(eit.value());
				//line->GetLine(); // vector of Vec2f write heere!!!
				// iterate through verices and edges
				auto& polyline_as_edge = plg->GetEdge(eit.value()); //what arcs broken into

				auto geom_feature_idx = (MLGGeomIdx) MLG_geometric_features.size();
				polylineid_2_geomidx[polyline_as_edge.eid] = geom_feature_idx;
				// now make a geometry piece for this node
				MLG_geometric_features.push_back(GeomFeat());
				auto& g = MLG_geometric_features.back();
				g.geom_feat_id = geom_feature_idx;
				g.dim = 1;
				g.points = polyline_as_edge.store->GetLine();

				// add this gid to the set of nodes
				MLG_geomidx_2_nodeidx[geom_feature_idx] = MLG_nodes.size();
				MLG_nodes.push_back(geom_feature_idx);
			}
			
			std::unordered_map<MLGNodeIdx, MLGGeomIdx> plgvertexid_2_geomidx;

			// first make the 0-dimensional geometry records for each old vertex
			GInt::Geometric2DGraph::vertex_iterator vit(plg);
			for (vit.begin(); vit.valid(); vit.advance()) {
				//const GInt::Vec2f& node_position = geometric_graph->GetVertexData(vit.value());
				auto plg_vertex_id = vit.value();

				auto geom_feature_idx = MLG_geometric_features.size();
				plgvertexid_2_geomidx[plg_vertex_id] = geom_feature_idx;
				// add node geom
				MLG_geometric_features.push_back(GeomFeat());
				auto& g = MLG_geometric_features.back();
				g.geom_feat_id = geom_feature_idx;
				g.dim = 0;
				g.points.push_back(plg->GetVertexData(plg_vertex_id));
			}

			// each old vertex adds in a bunch of arcs 
			for (vit.begin(); vit.valid(); vit.advance()) {
				//const GInt::Vec2f& node_position = geometric_graph->GetVertexData(vit.value());
				auto& plg_vertex = plg->GetVertex(vit.value());
				auto plg_vertex_id = vit.value();

				auto mlg_geom_idx = plgvertexid_2_geomidx[plg_vertex_id];
				if (plg_vertex.edges.size() == 1) {
					auto& plg_polyline_0 = plg->GetEdge(plg_vertex.edges[0]);
					auto mlg_polyline_0_as_node_idx = MLG_geomidx_2_nodeidx[polylineid_2_geomidx[plg_polyline_0.eid]];
					MLG_edges.push_back({ mlg_geom_idx, mlg_polyline_0_as_node_idx, -1 });
				}
				else {
					for (int i = 0; i < plg_vertex.edges.size(); i++) {
						for (int j = i + 1; j < plg_vertex.edges.size(); j++) {
							auto& plg_polyline_i = plg->GetEdge(plg_vertex.edges[i]);
							auto& plg_polyline_j = plg->GetEdge(plg_vertex.edges[j]);
							auto mlg_polyline_i_as_node_idx = MLG_geomidx_2_nodeidx[polylineid_2_geomidx[plg_polyline_i.eid]];
							auto mlg_polyline_j_as_node_idx = MLG_geomidx_2_nodeidx[polylineid_2_geomidx[plg_polyline_j.eid]];
							MLG_edges.push_back({ mlg_geom_idx, mlg_polyline_i_as_node_idx, mlg_polyline_j_as_node_idx });
						}
					}
				}
			}
			printf("CreateFrom2DPolylineGraph created %d nodes, %d arcs, %d geoms!\n", MLG_nodes.size(), MLG_edges.size(), MLG_geometric_features.size());
		}


		void Write(std::string basename) {
			m_basename = basename;
			std::ofstream nodes_file(basename + ".mlg_nodes.txt");
			for (auto geom_idx : MLG_nodes) {
				nodes_file << geom_idx << "\n";
			}
			nodes_file.close();

			std::ofstream edges_file(basename + ".mlg_edges.txt");
			for (auto& mlg_edge : MLG_edges) {
				edges_file << mlg_edge.edge_geom_idx << " " << mlg_edge.node_1_idx << " " << mlg_edge.node_2_idx << "\n";
			}
			edges_file.close();

			std::ofstream geom_file(basename + ".mlg_geom.txt");
			for (auto& g : MLG_geometric_features) {
				geom_file << g.geom_feat_id << " " << g.dim;
				for (auto& v : g.points) {
					geom_file << " " << v[0] << " " << v[1];
				}
				geom_file << "\n";
			}
			geom_file.close();
		}

		void Read(std::string basename) {
			m_basename = basename;
			MLG_nodes.clear();
			MLG_edges.clear();
			MLG_geometric_features.clear();

			// where to place the read line
			std::string line;
			std::ifstream geom_file(basename + ".mlg_geom.txt");
			while (std::getline(geom_file, line)) {
				std::stringstream s_stream(line); //create string stream from the string
				// fill vector of strings
				std::vector<std::string> results{ std::istream_iterator<std::string>{s_stream}, std::istream_iterator<std::string>{} };

				MLG_geometric_features.push_back(GeomFeat());
				auto& g = MLG_geometric_features.back();
				g.geom_feat_id = std::stoi(results[0]);
				g.dim = std::stoi(results[1]);
				for (int i = 2; i < results.size(); i += 2) {
					g.points.push_back({ std::stof(results[i]), std::stof(results[i + 1]) });
				}
			}
			geom_file.close();

			std::ifstream nodes_file(basename + ".mlg_nodes.txt");
			while (std::getline(nodes_file, line)) {
				MLGGeomIdx geomidx = std::stoi(line);
				MLG_geomidx_2_nodeidx[geomidx] = MLG_nodes.size();
				MLG_nodes.push_back(geomidx);
			}
			nodes_file.close();

			std::ifstream edges_file(basename + ".mlg_edges.txt");
			while (std::getline(edges_file, line)) {
				
				std::stringstream s_stream(line); //create string stream from the string
				// fill vector of strings
				std::vector<std::string> results{ std::istream_iterator<std::string>{s_stream}, std::istream_iterator<std::string>{} };
				MLG_edges.push_back({ std::stoi(results[0]), std::stoi(results[1]),std::stoi(results[2]) });
			}
			edges_file.close();
			printf("MLG Read: %d geoms, %d nodes, %d edges\n", MLG_geometric_features.size(), MLG_nodes.size(), MLG_edges.size());
		}

		GInt::Geometric2DGraph* Create2DPolylineGraph() {
			GInt::Geometric2DGraph* graph = new GInt::Geometric2DGraph();

			int numv = 0;
			int numpl = 0;

			std::unordered_map<int, int> feat_id_2_pos;

			for (auto& g : MLG_geometric_features) {
				if (g.dim == 0) {
					// add vertex
					graph->AddVertex(numv, g.points[0]);
					feat_id_2_pos[g.geom_feat_id] = numv;
					numv++;
				}
				else {
					numpl++;
				}
			}
			printf("# of stuff dim0 = %d, dimother = %d, size=%d\n", numv, numpl, graph->NumVertices());

			std::vector<std::pair<MLGGeomIdx, MLGGeomIdx>> polyline_nodeidx_2_vertex_geomidx;
			
			// iterate over the nodes (which are polylines)
			for (int i = 0; i < MLG_nodes.size(); i++) {
				polyline_nodeidx_2_vertex_geomidx.push_back({ -1, -1 });
			}
			// for each edge, mark the polyline endpoints with the vertex
			for (auto& mlg_edge : MLG_edges) {
				// add polyline
				auto polyline_node_id1 = mlg_edge.node_1_idx;
				auto polyline_node_id2 = mlg_edge.node_2_idx;

				auto vertex_geom_id = mlg_edge.edge_geom_idx;

				auto& plpair1 = polyline_nodeidx_2_vertex_geomidx[polyline_node_id1];
				if (plpair1.first == -1) {
					plpair1.first = vertex_geom_id;
				}
				else if (plpair1.first == vertex_geom_id) {
					// do nothing
				}
				else if (plpair1.second == -1) {
					plpair1.second = vertex_geom_id;
				}
				else if (plpair1.second == vertex_geom_id) {
					// do nothing
				}
				else {
					printf("WHujuuuuuuoooooopsiedaisy1!! plid1=%d, p0=%d p1=%d, a.eid=%d\n", polyline_node_id1, plpair1.first, plpair1.second, mlg_edge.edge_geom_idx);
				}
				if (polyline_node_id2 != -1) {
					auto& plpair2 = polyline_nodeidx_2_vertex_geomidx[polyline_node_id2];
					if (plpair2.first == -1) {
						plpair2.first = vertex_geom_id;
					}
					else if (plpair2.first == vertex_geom_id) {
						// do nothing
					}
					else if (plpair2.second == -1) {
						plpair2.second = vertex_geom_id;
					}
					else if (plpair2.second == vertex_geom_id) {
						// do nothing
					}
					else {
						printf("WHujuuuuuuoooooopsiedaisy2!! plid2=%d, p0=%d p1=%d, a.eid=%d\n", polyline_node_id2, plpair2.first, plpair2.second, mlg_edge.edge_geom_idx);
					}
				}
			}

			for (int i = 0; i < MLG_nodes.size(); i++) {
				auto& p = polyline_nodeidx_2_vertex_geomidx[i];
				if (p.first == -1) {
					printf("#### error: polyline %d p.first == -1\n", i);
				}				
				if (p.second == -1) {
					printf("#### error: polyline %d p.second == -1\n", i);
				}
				if (MLG_geometric_features[p.first].dim != 0) {
					printf("#### error: polyline %d p.first dim == %d\n", i, MLG_geometric_features[p.first].dim);
				}
				if (MLG_geometric_features[p.second].dim != 0) {
					printf("#### error: polyline %d p.second dim == %d\n", i, MLG_geometric_features[p.second].dim);
				}
			}

			// note we will add polylines in the order of the nodes, 
			// so node index i will be GetEdge(i) in the polylinegraph
			for (int i = 0; i < polyline_nodeidx_2_vertex_geomidx.size(); i++) {
				auto polyline_id = MLG_nodes[i];
				auto v1_id = polyline_nodeidx_2_vertex_geomidx[i].first;
				auto v2_id = polyline_nodeidx_2_vertex_geomidx[i].second;

				if (v1_id < MLG_nodes.size() || v2_id < MLG_nodes.size()) {
					printf("Warning: Create2DPolylineGraph -- v1id = %d, v2id= %d\n", v1_id, v2_id);
				}


				GInt::Line2d* points = new GInt::Line2d();
				for (auto& v : MLG_geometric_features[MLG_nodes[polyline_id]].points) points->AddToEnd(v);
				graph->AddEdge(feat_id_2_pos[ MLG_geometric_features[v1_id].geom_feat_id],
					feat_id_2_pos[ MLG_geometric_features[v2_id].geom_feat_id],
					points);
			}

			printf("Done: made %d vertices, %d polylines\n", graph->NumVertices(), graph->NumEdges());

			return graph;
		}

	};

}
#endif
