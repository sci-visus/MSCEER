/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef MSC_SELECTORS_H
#define MSC_SELECTORS_H

#include "gi_strictly_numeric_integrator.h"
#include "gi_numeric_integrator_region_stop.h"
#include "gi_numeric_integrator_expanding_region_stop.h"
#include "gi_timing.h"
#include "gi_labeling_to_bounary_labeling.h"
#include "gi_topological_explicit_mesh_function.h"
#include "gi_topological_region_growing_simple_gradient_builder.h"
#include "gi_topological_convergent_gradient_builder.h"
//#include "gi_robin_labeling.h"
#include "gi_adaptive_in_quad_euler_advector.h"
#include "gi_numeric_integrator_2d_restricted_expanding_region.h"
#include "gi_topological_2d_restricted_expanding_regions.h"
#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"
#include "gi_isolated_region_remover.h"
//#include "gi_topological_utility_functions.h"
#include "gi_morse_smale_complex_restricted.h"
#include "gi_numeric_streamline_integrator.h"

namespace GInt {


	class Selector {

	protected:

		vector<Selector*> parents;
		vector<Selector*> children;

		// when this is called parents have been computed
		virtual void SelectorAction() {}
		bool output_valid;
		size_t add_child(Selector* s) {
			children.push_back(s);
			return children.size();
		}
	public:

		set<INT_TYPE> output;

		Selector() : output_valid(false) {}
		virtual void invalidate() {
			if (!output_valid) return; // to prevent infinite looping by accident
			output_valid = false;
			for (auto it = children.begin(); it != children.end(); it++) (*it)->invalidate();
		}

		size_t add_parent(Selector* s) {
			s->add_child(this);
			parents.push_back(s);
			invalidate();
			return parents.size();
		}

		void compute_output() {
			if (output_valid) return;
			output.clear();
			for (auto it = parents.begin(); it != parents.end(); it++) {
				(*it)->compute_output();
			}

			SelectorAction();

			output_valid = true;
			return;


		}
	};

	template<class MSCType>
	class MSCSelectorLivingNodes : public Selector {
	protected:
		MSCType* mMsc;
		virtual void SelectorAction() {
			//printf("MSCSelectorLivingNodes::SelectorAction() computing...\n");
			for (int i = 0; i < mMsc->numNodes(); i++) {
				if (mMsc->isNodeAlive(i)) this->output.insert(i);
			}
			//printf("MSCSelectorLivingNodes::SelectorAction produced %d indices\n", output.size());
		}
	public:
		MSCSelectorLivingNodes(MSCType* msc) : mMsc(msc) {}
	};

	template<class MSCType>
	class MSCSelectorNodeIndex : public Selector {
	protected:
		MSCType* mMsc;
		DIM_TYPE mIndex;
		virtual void SelectorAction() {
			//printf("MSCSelectorNodeIndex::SelectorAction() computing...\n");
			for (auto pit = parents.begin(); pit != parents.end(); pit++) {
				for (auto it = (*pit)->output.begin(); it != (*pit)->output.end(); it++) {
					if (mMsc->getNode(*it).dim == mIndex) this->output.insert(*it);
				}
			}
			//printf("MSCSelectorNodeIndex::SelectorAction produced %d indices\n", output.size());
		}
	public:
		MSCSelectorNodeIndex(MSCType* msc, DIM_TYPE index) : mMsc(msc), mIndex(index) {}
	};

	template<class MSCType>
	class MSCSelectorRepresentative1Saddle : public Selector {
	protected:
		MSCType* mMsc;
		DIM_TYPE mIndex;
		pair<INT_TYPE, INT_TYPE> getExtremumPair(INT_TYPE saddle) {
            typename MSCType::SurroundingLivingArcsIterator sit(mMsc);
			INT_TYPE extrema[2]; int numext = 0;
			for (sit.begin(saddle); sit.valid(); sit.advance()) {
				//printf("asdf %d\n", sit.value());
				INT_TYPE aid = sit.value();
				arc<FLOATTYPE>& a = mMsc->getArc(aid);
				if (a.upper == saddle) extrema[numext++] = a.lower;
			}
			//printf("done numext=%d\n", numext);
			if (numext == 1) {
				return pair<INT_TYPE, INT_TYPE>(extrema[0], extrema[0]);
			}
			else {
				if (extrema[0] < extrema[1]) {
					return pair<INT_TYPE, INT_TYPE>(extrema[0], extrema[1]);
				} 
				else {
					return pair<INT_TYPE, INT_TYPE>(extrema[1], extrema[0]);
				}

			}

		}
		bool lessthan(INT_TYPE a, INT_TYPE b) {
			node<FLOATTYPE>& na = mMsc->getNode(a);
			node<FLOATTYPE>& nb = mMsc->getNode(b);
			if (na.value == nb.value) return a < b;
			return na.value < nb.value;
		}

		virtual void SelectorAction() {
			//printf("MSCSelectorRepresentative1Saddle::SelectorAction() computing...\n");
			map<pair<INT_TYPE, INT_TYPE>, INT_TYPE> saddlemap;

			for (auto pit = parents.begin(); pit != parents.end(); pit++) {
				for (auto it = (*pit)->output.begin(); it != (*pit)->output.end(); it++) {
					
					INT_TYPE saddle = *it;
					//printf("doing %d\n", saddle);
					pair<INT_TYPE, INT_TYPE> p = getExtremumPair(saddle);

					if (saddlemap.count(p) == 0) {
						saddlemap[p] = saddle;
					}
					else {
						INT_TYPE othersaddle = saddlemap[p];
						if (lessthan(saddle, othersaddle)) {
							saddlemap[p] = saddle;
						}
					}
				}
			}
			//printf("b\n");
			for (auto it = saddlemap.begin(); it != saddlemap.end(); it++) {
				if ((*it).first.first != (*it).first.second)
					output.insert((*it).second);
			}
			//printf("MSCSelectorRepresentative1Saddle::SelectorAction produced %d indices\n", output.size());
		}

	public:
		MSCSelectorRepresentative1Saddle(MSCType* msc) : mMsc(msc) {}
	};

}

#endif
