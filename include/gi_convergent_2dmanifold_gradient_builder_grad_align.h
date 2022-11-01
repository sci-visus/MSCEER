/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef GI_CONVERGENT_2DMANIFOLD_GRAD_BUILDER_GRAD_ALIGN
#define GI_CONVERGENT_2DMANIFOLD_GRAD_BUILDER_GRAD_ALIGN


#include "gi_convergent_2dmanifold_gradient_builder.h"
// This file is a safe space to put experimental classes and functionality.
#include "gi_max_vertex_labeling.h"

namespace GInt {

	template <class DependencyGraphType, class MeshType, class TopoFuncType, class GridFuncType, class MaxVLabelingType>
	class Convergent2ManifoldGradientBuilderGradAlign : public TopologicalGraphExpandingRegionsBTiming<SmallMemberDist, DependencyGraphType> {
	public:


		Vec3d GridCoords(INDEX_TYPE id) {
			Vec3l ccords;
			mGrid->cellid2Coords(id, ccords);
			Vec3d dcoords = ccords;
			dcoords *= 0.5;
			return dcoords;
		}


		Vec3d IdGrad(Vec3d coords) {
			return mBaseFunc->InterpolatedGrad(coords);
		}

		MeshType* mGrid;
		TopoFuncType* mFunc;
		DenseLabeling<char>* mLabel;
		MaxVLabelingType* mMaxLabel;
		//DiscreteGradientLabeling* mGrad;
		GridFuncType* mBaseFunc;

		Convergent2ManifoldGradientBuilderGradAlign(DependencyGraphType* graph,
			MeshType* grid,
			TopoFuncType* func,
			MaxVLabelingType* max_label,
			GridFuncType* base_func,
			// NEED REGULAR MESH FUNC
			// NEED MAXV LABELING
			//DiscreteGradientLabeling* grad)
			DenseLabeling<char>* label)
			: mGrid(grid), mFunc(func), mMaxLabel(max_label), mBaseFunc(base_func),  mLabel(label), TopologicalGraphExpandingRegionsBTiming(graph)
		{
		}

		//TopologicalExplicitDenseMeshFunction* GetFunction() { return m_func; }
		virtual void mdCombine(float scale, const SmallMemberDist& a, SmallMemberDist& res) {
			for (int i = 0; i < a.pairs.size(); i++) {
				bool has = false;
				for (int j = 0; j < res.pairs.size(); j++) {
					if (a.pairs[i].id == res.pairs[j].id) {
						has = true;
						res.pairs[j].prob += a.pairs[i].prob * scale;
						break;
					}
				}
				if (!has) {
					SmallPair p;
					p.id = a.pairs[i].id;
					p.prob = scale * a.pairs[i].prob;
					// TOTAL HACK!!!!!
					if (p.prob > 0.1);
						res.pairs.push_back(p);
				}
			}
		}
		
		SmallMemberDist& getPrior(INDEX_TYPE id, RegionState* r) {
			if (r->prior_data.count(id) == 0) {
				if (m_prior_data.count(id) == 0) {
					//printf("whoathere got here\n");
					return r->prior_data[id];
				}
				return m_prior_data[id];
			}
			return r->prior_data[id];
		}
		

		void CombineFacetMDs(int count, INDEX_TYPE* edges, SmallMemberDist& result, RegionState* region, Vec3d& coord) {
			if (count == 1) {
				result = getPrior(edges[0], region);
				
				coord = GridCoords(edges[0]);
				return;
			}
			coord = Vec3d(0, 0, 0); 
			float oneoversize = 1.0f / count;
			//printf("preres: ");
			//for (auto p : result.pairs) {
			//	printf("(%llu, %f) ", p.id, p.prob);
			//}
			//printf("\n");
			for (int i = 0; i < count; i++) {
			//printf("c %f: ", oneoversize);
			//for (auto p : getPrior(edges[i], region).pairs) {
			//	printf("(%llu, %f) ", p.id, p.prob);
			//}
			//printf("\n");
				mdCombine(oneoversize, getPrior(edges[i], region), result);
				coord += (GridCoords(edges[i]) * oneoversize);
			}		
			//printf("return "); coord.PrintFloat();
			//TOtal hack
				//if (result.pairs.size() > 5) {
				//	std::sort(result.pairs.begin(), result.pairs.end());
				//	result.pairs.erase(result.pairs.begin() + 5, result.pairs.end());
				//}

		}

		bool IsHIGHESTEdgeInQuadAndGather(INDEX_TYPE quad_id, INDEX_TYPE edge_id, INDEX_TYPE* edges, int& count) {
			count = 0;
			INDEX_TYPE target_vertex = mMaxLabel->Cell2HighestVertex(edge_id);
			if (target_vertex != mMaxLabel->Cell2HighestVertex(quad_id)) return false;
			//printf("ever here\n");
			/*typename */MeshType::FacetsIterator t_neighboring_primal_edge_iterator(mGrid);
			for (t_neighboring_primal_edge_iterator.begin(quad_id);
				t_neighboring_primal_edge_iterator.valid();
				t_neighboring_primal_edge_iterator.advance()){
				INDEX_TYPE t_edge_id = t_neighboring_primal_edge_iterator.value();
				if (t_edge_id == edge_id || this->mGrid->restrictionLabel(t_edge_id) == 0) continue;
				//if (target_vertex != mMaxLabel->Cell2HighestVertex(t_edge_id)) continue; // No! because the edge already has same highest vertex as quad!
				//if (mFunc->lessThan(target_vertex, mMaxLabel->Cell2HighestVertex(t_edge_id))) continue
				//printf("everhere2\n");
				if (this->AIsGreaterThanB( t_edge_id, edge_id)) return false; // might use different SOS - this might work though...
				//printf("everhere3\n"); 

				edges[count++] = t_edge_id;
			}
			return true;

		}



		virtual float mdDist(const SmallMemberDist& a, const SmallMemberDist& b) {
			float result = 0.0f;

			for (int i = 0; i < a.pairs.size(); i++) {
				bool has = false;
				for (int j = 0; j < b.pairs.size(); j++) {
					if (a.pairs[i].id == b.pairs[j].id) has = true;
				}
				if (!has) printf("a[%d] does not exist in b, should never get here\n", i);
			}

			for (int j = 0; j < b.pairs.size(); j++) {
				bool has = false;
				for (int i = 0; i < a.pairs.size(); i++) {
					if (a.pairs[i].id == b.pairs[j].id) {
						has = true;
						//float tmp = a.pairs[i].prob * b.pairs[j].prob;
						//if (tmp > result) result = tmp;
						double tmp = ((double) a.pairs[i].prob) - ((double) b.pairs[j].prob);
						result += tmp*tmp;
						//break;
					}
				}
				if (!has) {
					result += b.pairs[j].id *  b.pairs[j].id;
				}
			}
			return sqrt( result);
		}
		virtual float mdMax(const SmallMemberDist& a, const SmallMemberDist& b) {
			float result = 0.0f;
			for (int i = 0; i < a.pairs.size(); i++) {
				for (int j = 0; j < b.pairs.size(); j++) {
					if (a.pairs[i].id == b.pairs[j].id) {
						float tmp = a.pairs[i].prob * b.pairs[j].prob;
						if (tmp > result) result = tmp;
						//result += a.pairs[i].prob * b.pairs[j].prob;
						//break;
					}
				}
			}
			return result;
		}

		virtual float mdDot3(const SmallMemberDist& a, const SmallMemberDist& b) {
			double result = 0.0f;
			for (int i = 0; i < a.pairs.size(); i++) {
				for (int j = 0; j < b.pairs.size(); j++) {
					if (a.pairs[i].id == b.pairs[j].id) {
						//float tmp = a.pairs[i].prob * b.pairs[j].prob;
						//if (tmp > result) result = tmp;
						result += ((double) a.pairs[i].prob) * ((double) b.pairs[j].prob);
						//break;
					}
				}
			}
			return result;
		}
		void MakeCritical(SmallMemberDist& md, INDEX_TYPE id) {
			//printf("MAKING CRITICAL\n");
			SmallPair p; p.id = id; p.prob = 1.0f;
			md.pairs.push_back(p);
			MarkCritical(id);
		}

		//void GlobalMark(INDEX_TYPE id){
		//	mLabel->SetLabel(id, mLabel->GetLabel(id) + 1);
		//}
		void MarkCritical(INDEX_TYPE id) {
			mLabel->SetLabel(id, mLabel->GetLabel(id) + 1);
			//mGrad->setCritical(id, true);
			//mGrad->setAssigned(id, true);
			
		}

		void MarkPair(INDEX_TYPE edgeid, INDEX_TYPE toid) {
			mLabel->SetLabel(edgeid, mLabel->GetLabel(edgeid) + 1);
			mLabel->SetLabel(toid, mLabel->GetLabel(toid) + 1);
			//mGrad->setAssigned(edgeid, true);
			//mGrad->setAssigned(toid, true);
			//mGrad->setPair(edgeid, toid);
			//mGrad->setPair(toid, edgeid);
		}

		virtual SmallMemberDist DoWork(INDEX_TYPE id, RegionState* region, int asdf) {
			//printf("%f\n", mFunc->cellValue(id));
			//std::vector<INDEX_TYPE> negs;
			//int mval = 0;
			//printf("dowork2!\n");
			//FixedMemberDist res;
			INDEX_TYPE other_edges[4];
			int count_quads = 0;
			//FixedMemberDists quads_md[4];
			FLOATTYPE quad_diffs[4];
			INDEX_TYPE quad_ids[4];
			vector<SmallMemberDist> tstore;
			SmallMemberDist mydist;
			vector<Vec3d> coords;
			//FLOATTYPE my_value = mFunc->cellValue(id);
			// FIRST GATHER NEIGHBORING DISTRIBUTIONS ---------------------------
			//// FOR EACH QUAD, decide if its in the lower link of the edge
			//// if so, combine the MDs ofall lower neighboring edges
			// store the "drop" in the quad_diffs
			/*typename*/ MeshType::CofacetsIterator t_adjacent_primal_quad_iterator(mGrid);
			/*typename*/ MeshType::FacetsIterator t_neighboring_primal_edge_iterator(mGrid);
			for (t_adjacent_primal_quad_iterator.begin(id);
				t_adjacent_primal_quad_iterator.valid();
				t_adjacent_primal_quad_iterator.advance()){
				INDEX_TYPE t_quad_id = t_adjacent_primal_quad_iterator.value();

				int count_others;
				//float lowest_val;
				if (IsHIGHESTEdgeInQuadAndGather(t_quad_id, id, other_edges, count_others) && count_others > 0) {
					quad_ids[tstore.size()] = t_quad_id;
					tstore.push_back(SmallMemberDist());
					SmallMemberDist& md = tstore[tstore.size()-1];	
					md.pairs.clear();

//#pragma omp critical
//					{

					Vec3d coord;
					CombineFacetMDs(count_others, other_edges, md, region, coord);
					coords.push_back(coord);
					//	printf("combine: ");
					//		for (auto p : md.pairs) {
					//		printf("(%llu, %f) ", p.id, p.prob);
					//		}
					//		printf("\n");
					//}
				}
			}
			//printf("%f: ", my_value); for (int iii = 0; iii < tstore.size(); iii++) printf("%f ", quad_diffs[iii]); printf("\n");
			if (tstore.size() == 0) {
				// make critical
				MakeCritical(mydist, id);
				//
				//GlobalMark(id);
				return mydist;
				// RETURN HERE
			}
			else if (tstore.size() == 1) {
				// there is only one downwards direction, so pick it
				mydist = tstore[0];
				//GlobalMark(id);
				//GlobalMark(quad_ids[0]);
				MarkPair(id, quad_ids[0]);
				return mydist;
				// RETURN here
			}
			


			// NOW COMPUTE MY MEMBER DISTRIBUTION ------------------------------
			// we have options, so pick the best one!

			Vec3d cell_coords = GridCoords(id);
			Vec3d cell_grad = IdGrad(cell_coords);

			//printf("have %d options from "); cell_coords.PrintFloat();
			//printf("grad = "); cell_grad.PrintFloat();
			
			FLOATTYPE totalval = 0;
			for (int i = 0; i < tstore.size(); i++) {
				Vec3d othercoord = coords[i];
				float align_val = (cell_coords - othercoord).Dot(cell_grad);
				//printf("%d, align = %f, ", i, align_val); othercoord.PrintFloat();
				//if (mMesh->dimension(cellid) > 0) printf("%f ", align_val);
				if (align_val < 0) align_val = 0;
				totalval += align_val;
				quad_diffs[i] = align_val;
			}

				// now we need to rescale the probabilities for combination
				if (totalval <= 0) {
					// just give them even weight
					for (int i = 0; i < tstore.size(); i++) {
						quad_diffs[i] = quad_diffs[i] / tstore.size();
					}
				}
				else {
					for (int i = 0; i < tstore.size(); i++) {
						quad_diffs[i] = quad_diffs[i] / totalval;
					}
				}

				for (int i = 0; i < tstore.size(); i++) {
					mdCombine(quad_diffs[i], tstore[i], mydist);
				}



			// NOW MYDIST has the member distribution we need! ------------------
				int maxloc = 0;
				std::vector<float> tress; // the similarity of candidate[i]
				float maxval = mdDot3(tstore[0], mydist);
				//float maxval = mdMax(tstore[0], mydist);
				//float maxval = mdDist(tstore[0], mydist);
				tress.push_back(maxval);
				//printf("%d=%.2f\n", candidates[0], minval);
				for (int i = 1; i < tstore.size(); i++) {

					float otherval = mdDot3(tstore[i], mydist);
					//float otherval = mdMax(tstore[i], mydist);
					//float otherval = mdDist(tstore[i], mydist);
					//printf("%d=%.2f\n", candidates[i], otherval);
					tress.push_back(otherval);
					if (otherval > maxval ||
					//if (otherval < maxval ||
							(otherval == maxval && quad_diffs[i] > quad_diffs[maxloc])) {
						maxval = otherval;
						maxloc = i;
					}
				}

//#pragma omp critical 
//				{
//					for (auto p : mydist.pairs) {
//						printf("(%llu, %f) ", p.id, p.prob);
//					}
//					printf("=\n");
//					for (int qi = 0; qi < tstore.size(); qi++) {
//						printf("%f * ", quad_diffs[qi]);
//
//						for (auto p : tstore[qi].pairs) {
//							printf("(%llu, %f) ", p.id, p.prob);
//						}
//
//						if (qi == maxloc) printf(" <== ");
//						printf(" -> %f\n", tress[qi]);
//					}
//					printf("\n");
//				}

				// now maxloc has most similar
				// MAKE the pair
				// edge id to quad id
				// now have to pair the rest in the lower star!
				MarkPair(id, quad_ids[maxloc]);
				//GlobalMark(id);
				//GlobalMark(quad_ids[maxloc]);
				//printf("maxloc = %d\n", maxloc);
				return mydist;

		}

		virtual SmallMemberDist DoInitialWork(INDEX_TYPE id, RegionState* region) {
			SmallMemberDist res;
			if (HasNoLowerNeighbors(id)) {
				MakeCritical(res, id);
				printf("CRITICAL FROM HERE\n");
				//
				//GlobalMark(id);
			}
			else {
				res = DoWork(id, region, 0);
			}
			return res;// region->prior_data[id] = 0;
		}

	};

}
#endif
