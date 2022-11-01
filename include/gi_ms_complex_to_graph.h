/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef GI_MSC_2_GRAPHS_H
#define GI_MSC_2_GRAPHS_H

#include <vector>
#include "gi_basic_types.h"
#include "gi_basic_geometry.h"
#include "gi_graphs.h"

#include <map>
#include <unordered_map>

// This file is a safe space to put experimental classes and functionality.

namespace GInt {

	template<class MscType, class MeshType> 
	MeshCellsGraph* BuildMeshCellsGraphFromMSCRidges(MscType* msc, MeshType* topological_grid, bool ignore_forced = false) {

		MeshCellsGraphBuilder<MeshType> rg(topological_grid);
		typename MscType::LivingArcsIterator ait(msc);
		for (ait.begin(); ait.valid(); ait.advance()) {
			INT_TYPE aid = ait.value();
			GInt::arc<float>& a = msc->getArc(aid);
			GInt::node<float>& saddle = msc->getNode(a.lower);
			if (saddle.dim != topological_grid->maxDim() - 1) continue;
			vector<INDEX_TYPE> v;
			msc->fillArcGeometry(aid, v);
			for (auto id : v) rg.AddCellIndex(id);
			if (!ignore_forced) {
				rg.AddToForcedVertices(saddle.cellindex);
				rg.AddToForcedVertices(msc->getNode(a.upper).cellindex);
			}
			else {
				// get 2 nodes
				INT_TYPE sad_max_arcs[2];
				auto count_arcs = msc->GetLiving2SaddleExtremaArcs(a.lower, sad_max_arcs);
				if (count_arcs < 2) {
					rg.AddToForcedVertices(saddle.cellindex);
				}
				else {
					if (msc->arcUpperNode(sad_max_arcs[0]) == msc->arcUpperNode(sad_max_arcs[1])) {
						rg.AddToForcedVertices(saddle.cellindex);
					}
				}
			}
		}
		return rg.ComputeStuff();
	}

	template<class MscType, class MeshType>
	MeshCellsGraph* BuildMeshCellsGraphFromMSCValleys(MscType* msc, MeshType* topological_grid, bool ignore_forced = false) {

		MeshCellsGraphBuilder<MeshType> rg(topological_grid);
		typename MscType::LivingArcsIterator ait(msc);
		for (ait.begin(); ait.valid(); ait.advance()) {
			INT_TYPE aid = ait.value();
			GInt::arc<float>& a = msc->getArc(aid);
			GInt::node<float>& minimum = msc->getNode(a.lower);
			if (minimum.dim != 0) continue;
			vector<INDEX_TYPE> v;
			msc->fillArcGeometry(aid, v);
			for (auto id : v) rg.AddCellIndex(id);
			if (!ignore_forced) {
				rg.AddToForcedVertices(minimum.cellindex);
				rg.AddToForcedVertices(msc->getNode(a.upper).cellindex);
			}
			else {
				// get 2 nodes
				INT_TYPE sad_max_arcs[2];
				auto count_arcs = msc->GetLiving1SaddleExtremaArcs(a.upper, sad_max_arcs);
				if (count_arcs < 2) {
					rg.AddToForcedVertices(msc->getNode(a.upper).cellindex);
				}
				else {
					if (msc->arcLowerNode(sad_max_arcs[0]) == msc->arcLowerNode(sad_max_arcs[1])) {
						rg.AddToForcedVertices(msc->getNode(a.upper).cellindex);
					}
				}
			}
		}
		return rg.ComputeStuff();
	}



	template<class MscType, class MeshType, class MeshFuncType>
	MeshCellsGraph* BuildMeshCellsGraphFromMSCRidgesWithCuts(MscType* msc, MeshType* topological_grid, MeshFuncType* func, std::vector<float>& cut_points) {

		MeshCellsGraphBuilder<MeshType> rg(topological_grid);
		typename MscType::LivingArcsIterator ait(msc);
		for (ait.begin(); ait.valid(); ait.advance()) {
			INT_TYPE aid = ait.value();
			GInt::arc<float>& a = msc->getArc(aid);
			GInt::node<float>& saddle = msc->getNode(a.lower);
			GInt::node<float>& maximum = msc->getNode(a.upper);
			if (saddle.dim != topological_grid->maxDim() - 1) continue;
			vector<INDEX_TYPE> v;
			msc->fillArcGeometry(aid, v);
			std::vector<float> values;
			values.reserve(v.size());
			for (auto id : v) {
				rg.AddCellIndex(id);
				values.push_back(func->cellValue(id));
			}
			
			rg.AddToForcedVertices(saddle.cellindex);
			rg.AddToForcedVertices(msc->getNode(a.upper).cellindex);
			
			if (v.size() < 1) {
				printf("BuildMeshCellsGraphFromMSCRidgesWithCuts should never get here, %d\n", v.size());
				continue;
			}

			for (int i = 1; i < v.size(); i++) {
				for (auto val : cut_points) {
					// if the sign is different
					if ((values[i - 1] -val) * (values[i] - val) < 0) {

						rg.AddToForcedVertices(v[i-1]);
						//printf("adding %d, val = %f--%f\n", v[i], values[i - 1], values[i]);
					}
				}
			}



		}
		return rg.ComputeStuff();
	}

	template<class MscType, class MeshType, class MeshFuncType>
	MeshCellsGraph* BuildMeshCellsGraphFromMSCValleysWithCuts(MscType* msc, MeshType* topological_grid, MeshFuncType* func, std::vector<float>& cut_points) {

		MeshCellsGraphBuilder<MeshType> rg(topological_grid);
		typename MscType::LivingArcsIterator ait(msc);
		for (ait.begin(); ait.valid(); ait.advance()) {
			INT_TYPE aid = ait.value();
			GInt::arc<float>& a = msc->getArc(aid);
			GInt::node<float>& minimum = msc->getNode(a.lower);
			if (minimum.dim != 0) continue;
			vector<INDEX_TYPE> v;
			msc->fillArcGeometry(aid, v);
			std::vector<float> values;
			values.reserve(v.size());
			for (auto id : v) {
				rg.AddCellIndex(id);
				values.push_back(func->cellValue(id));
			}

			rg.AddToForcedVertices(minimum.cellindex);
			rg.AddToForcedVertices(msc->getNode(a.upper).cellindex);
			
			if (v.size() < 1) {
				printf("BuildMeshCellsGraphFromMSCValleysWithCuts should never get here, %d\n", v.size());
				continue;
			}

			for (int i = 1; i < v.size(); i++) {
				for (auto val : cut_points) {
					if ((values[i - 1] - val) * (values[i] - val) < 0) {
						rg.AddToForcedVertices(v[i]);
					}
				}
			}
		}
		return rg.ComputeStuff();
	}


}
#endif
