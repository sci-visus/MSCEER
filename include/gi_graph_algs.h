#ifndef GI_GRAPH_ALGS_H
#define GI_GRAPH_ALGS_H

#include <vector>
#include "gi_basic_types.h"
#include "gi_basic_geometry.h"
#include "gi_array_index_partition.h"
#include "gi_graphs.h"

#include <map>
#include <set>
#include <queue>
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
	//typedef float DType;/*typename DType,*/ 
	//typedef GInt::UndirectedGraph<int, int> GraphType; /*<class GraphType>*/
	// path is vector of vertex_id, from_edge_id, from target,downedge, to source,-1 
	template <typename DType, class GraphType>
	DType ShortestPath(GraphType* graph, INDEX_TYPE source, INDEX_TYPE target,
		const std::unordered_map<INDEX_TYPE, DType>& edge_weights, std::vector<std::pair<INDEX_TYPE, INDEX_TYPE>>& path) {
		// do priority queue based dijkstra
		std::unordered_map<INDEX_TYPE, std::pair<DType, INDEX_TYPE>> visited;
		// map distances to {from_edge_id, node_id}
			// Define the custom comparator function
		typedef  std::pair<DType, std::pair<INDEX_TYPE, INDEX_TYPE>> queueitem;
		struct Compare {
			bool operator()(const queueitem& a, const queueitem& b) const {
				return a > b; // Compare in reverse order for lowest value first
			}
		};
		std::priority_queue<queueitem, Compare> queue;
		queue.push({ 0, {-1, source} });
		while (!queue.empty()) {
			auto top = queue.top(); queue.pop();
			INDEX_TYPE current = top.second.second;
			DType dist = top.first;
			if (visited.count(current) != 0) {
				// this has been seen with a lower priority, so do nothing
				continue;
			}

			visited[current] = { dist, top.second.first }; // add dist, from edge to visited map
			
			if (current == target) {
				// build the path
				INDEX_TYPE tmpid = target;
				while (tmpid != source) {
					path.push_back({ tmpid, visited[tmpid].second });
				}
				path.push_back({ source, -1 });
				return dist;
			}

			// add neighbors
			for (auto neid : graph->GetVertex(current).edges) {
				auto& edge = graph->GetEdge(neid);
				INDEX_TYPE next_v = graph->OtherVertex(neid, current);
				if (visited.count(next_v) != 0) continue; // already seen, so ignore
				queue.push({ edge_weights[neid] + dist, {neid, next_v} }); // enqueue next
			}
		}
		return 0; // no path
	}

	// shortest path from one set of nodes to another
	template <typename DType, class GraphType>
	DType ShortestPath(GraphType* graph,
		const std::set<INDEX_TYPE>& sources,
		const std::set<INDEX_TYPE>& target,
		const std::unordered_map<INDEX_TYPE, DType>& edge_weights, std::vector<std::pair<INDEX_TYPE, INDEX_TYPE>>& path) {
		// do priority queue based dijkstra
		std::unordered_map<INDEX_TYPE, std::pair<DType, INDEX_TYPE>> visited;
		// map distances to {from_edge_id, node_id}
			// Define the custom comparator function
		typedef  std::pair<DType, std::pair<INDEX_TYPE, INDEX_TYPE>> queueitem;
		struct Compare {
			bool operator()(const queueitem& a, const queueitem& b) const {
				return a > b; // Compare in reverse order for lowest value first
			}
		};
		std::priority_queue<queueitem, Compare> queue;
		for (auto source : sources) queue.push({ 0, {-1, source} });
		while (!queue.empty()) {
			auto top = queue.top(); queue.pop();
			INDEX_TYPE current = top.second.second;
			DType dist = top.first;
			if (visited.count(current) != 0) {
				// this has been seen with a lower priority, so do nothing
				continue;
			}

			visited[current] = { dist, top.second.first }; // add dist, from edge to visited map

			if (target.count(current) != 0) {
				// build the path
				INDEX_TYPE tmpid = current;
				while (sources.count(tmpid) == 0) {
					path.push_back({ tmpid, visited[tmpid].second });
				}
				path.push_back({ tmpid, -1 });
				return dist;
			}

			// add neighbors
			for (auto neid : graph->GetVertex(current).edges) {
				auto& edge = graph->GetEdge(neid);
				INDEX_TYPE next_v = graph->OtherVertex(neid, current);
				if (visited.count(next_v) != 0) continue; // already seen, so ignore
				queue.push({ edge_weights[neid] + dist, {neid, next_v} }); // enqueue next
			}
		}
		return 0; // no path
	}


	template <typename DType>
	struct SourceItem {
		INDEX_TYPE source_id;
		INDEX_TYPE down_edge_id;
		DType distance;
	};
	// shortest path from one set of nodes to another
	// result is a map from
	template <typename DType, class GraphType>
	void NearestSources(GraphType* graph,
		const std::set<INDEX_TYPE>& sources,
		const std::unordered_map<INDEX_TYPE, DType>& edge_weights, 
		std::unordered_map<INDEX_TYPE, SourceItem<DType>>& res) {
		// do priority queue based dijkstra
		std::unordered_map<INDEX_TYPE, std::pair<DType, INDEX_TYPE>> visited;
		// map distances to {from_edge_id, node_id}
			// Define the custom comparator function
		//typedef  std::pair<DType, std::pair<INDEX_TYPE, INDEX_TYPE>> queueitem;
		//struct queueitem {
		//	DType dist;
		//	INDEX_TYPE vert_id;
		//	INDEX_TYPE from_edge_id;
		//	queueitem() : dist(0), vert_id(0), from_edge_id(0) {}
		//	queueitem(DType d, INDEX_TYPE v, INDEX_TYPE e) : dist(d), vert_id(v), from_edge_id(e) {}
		//	bool operator()(const queueitem& a, const queueitem& b) const {
		//		return a.dist > b.dist; // Compare in reverse order for lowest value first
		//	}
		//	queueitem& operator=(const queueitem& a) { 
		//		this->dist = a.dist; 
		//		this->from_edge_id = a.from_edge_id; 
		//		this->vert_id = a.vert_id; 
		//		return *this;
		//	}
		//};
		typedef  std::pair<DType, std::pair<INDEX_TYPE, INDEX_TYPE>> queueitem;
		struct Compare {
			bool operator()(const queueitem& a, const queueitem& b) const {
				return a.first > b.first; // Compare in reverse order for lowest value first
			}
		};
		std::priority_queue<queueitem, std::vector<queueitem>, Compare> queue;
		for (auto source : sources) queue.push({ (DType) 0, {source, -1ll} });
		while (!queue.empty()) {
			queueitem top = queue.top(); queue.pop();
			INDEX_TYPE current = top.second.first;
			DType dist = top.first;
			if (visited.count(current) != 0) {
				// this has been seen with a lower priority, so do nothing
				continue;
			}

			visited[current] = { dist, top.second.second }; // add dist, from edge to visited map
			// add neighbors
			for (INDEX_TYPE neid : graph->GetVertex(current).edges) {
				auto& edge = graph->GetEdge(neid);
				INDEX_TYPE next_v = graph->OtherVertex(neid, current);
				if (visited.count(next_v) != 0) continue; // already seen, so ignore
				queue.push({ (DType)((edge_weights.find(neid))->second + dist),  {(INDEX_TYPE)next_v, (INDEX_TYPE)neid} }); // enqueue next
			}
		}
		// now visited has downward map. 
		for (auto source : sources) res[source] = { source, -1, 0 };
		for (const auto& vp : visited) {
			if (res.count(vp.first) != 0) continue;
			std::vector<INDEX_TYPE> path;
			INDEX_TYPE current = vp.first;
			while (res.count(current) == 0) {
				path.push_back(current);
				current = graph->OtherVertex(visited[current].second, current); // trace down path
			}
			// we are at labeled branch
			const auto& source = res[current];
			for (auto id : path) {
				res[id] = { source.source_id, visited[id].second, visited[id].first };
			}
		}
	}

}
#endif
