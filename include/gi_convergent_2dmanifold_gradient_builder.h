/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef GI_EXPERIMENTAL4_H
#define GI_EXPERIMENTAL4_H

#include <vector>
#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_labeling.h"
#include "gi_regular_grid.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_topological_regular_grid.h"
#include "gi_discrete_gradient_labeling.h"
#include "gi_topological_utility_functions.h"
#include "gi_strictly_numeric_integrator.h"
#include "gi_numeric_integrator_region_stop.h"
#include "gi_numeric_integrator_expanding_region_stop.h"
#include "gi_timing.h"
#include "gi_labeling_to_bounary_labeling.h"
#include "gi_topological_explicit_mesh_function.h"
#include "gi_topological_region_growing_simple_gradient_builder.h"
#include "gi_topological_convergent_gradient_builder.h"
#include "gi_robin_labeling.h"
#include "gi_adaptive_in_quad_euler_advector.h"
#include "gi_numeric_integrator_2d_restricted_expanding_region.h"
#include "gi_topological_2d_restricted_expanding_regions.h"
#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_isolated_region_remover.h"
#include "gi_topological_utility_functions.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_morse_smale_complex_restricted.h"
#include "gi_kdtree.h"
#include "gi_msc_selectors.h"
#include "gi_numeric_streamline_integrator.h"
#include "gi_experimental.h"

//#include "concurrentqueue.h"
#include <thread>
#include <mutex>

#include "gi_experimental.h"
#include "gi_experimental3b.h"
#include "gi_experimental5.h"
// This file is a safe space to put experimental classes and functionality.

namespace GInt {


	//struct IdFPair {
	//	INDEX_TYPE id;
	//	float prob;
	//};
	//struct FixedMemberDist {
	//	IdFPair dists[5];
	//	BYTE_TYPE count;
	//	//BYTE_TYPE lowest;
	//};
	struct SmallPair {
		INDEX_TYPE id;
		float prob;
		bool operator<(const SmallPair& rhs) const { return prob > rhs.prob;  }
	};
	struct SmallMemberDist {
		std::vector<SmallPair> pairs;
	};
	template <class DependencyGraphType>
	class Convergent2ManifoldGradientBuilder : public TopologicalGraphExpandingRegionsBTiming<SmallMemberDist, DependencyGraphType> {
	public:

		TopologicalRegularGridRestricted* mGrid;
		TopologicalExplicitDenseMeshFunction<TopologicalRegularGridRestricted, float>* mFunc;
		DenseLabeling<char>* mLabel;
		//DiscreteGradientLabeling* mGrad;

		Convergent2ManifoldGradientBuilder(DependencyGraphType* graph,
			TopologicalRegularGridRestricted* grid,
			TopologicalExplicitDenseMeshFunction<TopologicalRegularGridRestricted, float>* func,
			//DiscreteGradientLabeling* grad)
			DenseLabeling<char>* label)
			: mGrid(grid), mFunc(func), /*mGrad(grad),*/  mLabel(label), TopologicalGraphExpandingRegionsBTiming(graph)
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
					printf("whoathere got here\n");
					return r->prior_data[id];
				}
				return m_prior_data[id];
			}
			return r->prior_data[id];
		}
		

		void CombineFacetMDs(int count, INDEX_TYPE* edges, SmallMemberDist& result, RegionState* region) {
			if (count == 1) {
				result = getPrior(edges[0], region);
				return;
			}
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
			}				
			//TOtal hack
				if (result.pairs.size() > 5) {
					std::sort(result.pairs.begin(), result.pairs.end());
					result.pairs.erase(result.pairs.begin() + 5, result.pairs.end());
				}

		}

		bool IsHIGHESTEdgeInQuadAndGather(INDEX_TYPE quad_id, INDEX_TYPE edge_id, INDEX_TYPE* edges, float& lowest_val, int& count) {
			count = 0;
			/*typename */TopologicalRegularGridRestricted::FacetsIterator t_neighboring_primal_edge_iterator(mGrid);
			for (t_neighboring_primal_edge_iterator.begin(quad_id);
				t_neighboring_primal_edge_iterator.valid();
				t_neighboring_primal_edge_iterator.advance()){
				INDEX_TYPE t_edge_id = t_neighboring_primal_edge_iterator.value();
				if (t_edge_id == edge_id || this->mGrid->restrictionLabel(t_edge_id) == 0) continue;
				if (this->AIsGreaterThanB( t_edge_id, edge_id)) return false;
				
				// find the lowest edge value
				if (count == 0) {
					lowest_val = mFunc->cellValue(t_edge_id);
				}
				else {
					FLOATTYPE fnew = mFunc->cellValue(t_edge_id);
					if (lowest_val > fnew) lowest_val = fnew;
				}
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

			FLOATTYPE my_value = mFunc->cellValue(id);
			// FIRST GATHER NEIGHBORING DISTRIBUTIONS ---------------------------
			//// FOR EACH QUAD, decide if its in the lower link of the edge
			//// if so, combine the MDs ofall lower neighboring edges
			// store the "drop" in the quad_diffs
			/*typename*/ TopologicalRegularGridRestricted::CofacetsIterator t_adjacent_primal_quad_iterator(mGrid);
			/*typename*/ TopologicalRegularGridRestricted::FacetsIterator t_neighboring_primal_edge_iterator(mGrid);
			for (t_adjacent_primal_quad_iterator.begin(id);
				t_adjacent_primal_quad_iterator.valid();
				t_adjacent_primal_quad_iterator.advance()){
				INDEX_TYPE t_quad_id = t_adjacent_primal_quad_iterator.value();

				int count_others;
				float lowest_val;
				if (IsHIGHESTEdgeInQuadAndGather(t_quad_id, id, other_edges, lowest_val, count_others) && count_others > 0) {
					quad_diffs[tstore.size()] = my_value - lowest_val;
					quad_ids[tstore.size()] = t_quad_id;
					tstore.push_back(SmallMemberDist());
					SmallMemberDist& md = tstore[tstore.size()-1];	
					md.pairs.clear();

//#pragma omp critical
//					{

					CombineFacetMDs(count_others, other_edges, md, region);

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
				FLOATTYPE totalval = 0;
				for (int i = 0; i < tstore.size(); i++) {
					totalval += quad_diffs[i];
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
			MakeCritical(res, id);
			//
			//GlobalMark(id);
			return res;// region->prior_data[id] = 0;
		}

	};

}
#endif
