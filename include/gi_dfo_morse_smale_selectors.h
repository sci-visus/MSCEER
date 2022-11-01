/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef GI_DFO_MSC_SELECTORS_H
#define GI_DFO_MSC_SELECTORS_H

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
#include "gi_dataflow_objects.h"

namespace GInt {

	template<class MSCType>
	class DFO_MSCHierarchyInterface : public DataflowObject {
	protected:
		MSCType* mMsc;
		DFO_VaryingValuatorDouble* mPersistenceValue;
		INDEX_TYPE mMscSelectPosition;
		virtual void Action() {
			if (mMsc == NULL) {
				printf("RUNTIME_ERROR::DFO_MSCHierarchyInterface::Action - mMsc is NULL during evaluation, doing nothing\n");
				return;
			}
			if (mPersistenceValue == NULL) {
				printf("RUNTIME_WARNING::DFO_MSCHierarchyInterface::Action - Persistence valuator is NULL, using -1 persistence default\n");
				mMsc->SetSelectPersAbs(-1);
			}
			else {
				mMsc->SetSelectPersAbs(mPersistenceValue->GetValue());
			}
		}	
	public:
		DFO_MSCHierarchyInterface() : mMsc(NULL), mPersistenceValue(NULL) {}
		
		virtual size_t add_parent(DataflowObject* s) {
			printf("STATIC_ERROR::DFO_MSCHierarchyInterface::add_parent - use SetMsc() and SetPersistenceValuatorDFO() instead!!\n");
			return 0;
		}

		virtual void SetMsc(MSCType* msc) { 
			mMsc = msc; 
		}
		virtual MSCType* GetMsc() { return mMsc; }
		virtual void SetPersistenceValuatorDFO(DFO_VaryingValuatorDouble* pers_val) {
			mPersistenceValue = pers_val;
			DataflowObject::add_parent(pers_val);
		}
	};

	template<class MSCType>
	class DFO_MSCNodeFilter : public IdSetDataflowObject {
	protected:
		DFO_MSCHierarchyInterface<MSCType>* mMsc;
		virtual void Action() {
			output_nodes.clear();
			if (mMsc == NULL) {
				printf("RUNTIME_ERROR::DFO_MSCNodeFilter::Action - mMsc is NULL! Use SetMscDFO() to set it up. doing nothing\n");
				return;
			}
			NodeFilterAction();
		}
		vector<DataflowObject*> input_objects; // the parents that are considered as "input"
		virtual void NodeFilterAction() = 0;
		set<INT_TYPE> output_nodes;
	public:
		DFO_MSCNodeFilter() : mMsc(NULL) {}
		virtual void SetMscDFO(DFO_MSCHierarchyInterface<MSCType>* msc) {
			mMsc = msc;
			add_parent(msc);
		}
		virtual void AddInputDFO(DataflowObject* s) {
			add_parent(s);
			input_objects.push_back(s);
		}
		virtual const set<INT_TYPE>&  GetOutput() const {
			return this->output_nodes;
		}

	};

	template<class MSCType>
	class DFO_MSCArcFilter : public IdSetDataflowObject {
	protected:
		DFO_MSCHierarchyInterface<MSCType>* mMsc;
		virtual void Action() {
			output_arcs.clear();
			if (mMsc == NULL) {
				printf("RUNTIME_ERROR::DFO_MSCNodeFilter::Action - mMsc is NULL! Use SetMscDFO() to set it up. doing nothing\n");
				return;
			}
			ArcFilterAction();
		}
		vector<DataflowObject*> input_objects; // the parents that are considered as "input"
		virtual void ArcFilterAction() = 0;
		set<INT_TYPE> output_arcs;
	public:
		DFO_MSCArcFilter() : mMsc(NULL) {}
		virtual void SetMscDFO(DFO_MSCHierarchyInterface<MSCType>* msc) {
			mMsc = msc;
			add_parent(msc);
		}
		virtual void AddInputDFO(DataflowObject* s) {
			add_parent(s);
			input_objects.push_back(s);
		}
		virtual const set<INT_TYPE>& GetOutput() const {
			return this->output_arcs;
		}

	};

	template<class MSCType>
	class DFO_MSCLivingNodes : public DFO_MSCNodeFilter<MSCType> {
	protected:
		virtual void NodeFilterAction() {
			MSCType* msc = this->mMsc->GetMsc();
			for (INT_TYPE i = 0; i < msc->numNodes(); i++) {
				if (msc->isNodeAlive(i)) this->output_nodes.insert(i);
			}
		}
	public:
		DFO_MSCLivingNodes() {}

		virtual void AddInputDFO(DataflowObject* s) {
			printf("STATIC_WARNING::DFO_MSCLivingNodesDFO::AddInputDFO - this function has no effect!!\n");
		}

	};

	template<class MSCType>
	class DFO_MSCLivingArcs : public DFO_MSCArcFilter<MSCType> {
	protected:
		virtual void ArcFilterAction() {
			MSCType* msc = this->mMsc->GetMsc();
			for (INT_TYPE i = 0; i < msc->numArcs(); i++) {
				if (msc->isArcAlive(i)) this->output_arcs.insert(i);
			}
		}
	public:
		DFO_MSCLivingArcs() {}

		virtual void AddInputDFO(DataflowObject* s) {
			printf("STATIC_WARNING::DFO_MSCLivingArcsDFO::AddInputDFO - this function has no effect!!\n");
		}

	};

	template<class MSCType>
	class DFO_MSCNodeIndexSelect : public DFO_MSCNodeFilter<MSCType> {
	protected:
		virtual void NodeFilterAction() {
			printf("DFO_MSCNodeIndexSelect::NodeFilterAction - computing\n");
			for (auto input : this->input_objects) {
				for (auto id : dynamic_cast<DFO_MSCNodeFilter<MSCType>*>(input)->GetOutput()) {
					if (mPass[mMsc->GetMsc()->getNode(id).dim]) this->output_nodes.insert(id);
				}
			}
		}

		bool mPass[4];
	public:

		void SetPassingIndices(bool minima, bool saddle_1, bool saddle_2, bool maxima) {
			invalidate();
			mPass[0] = minima;
			mPass[1] = saddle_1;
			mPass[2] = saddle_2;
			mPass[3] = maxima;
		}
	};

	template<class MSCType>
	class DFO_MSCArcsLowerIndexSelect : public DFO_MSCArcFilter<MSCType> {
	protected:
		virtual void ArcFilterAction() {
			printf("DFO_MSCArcsLowerIndexSelect::ArcFilterAction - computing\n");
			for (auto input : this->input_objects) {
				for (auto id : dynamic_cast<DFO_MSCArcFilter<MSCType>*>(input)->GetOutput()) {
					if (mPass[mMsc->GetMsc()->getArc(id).dim]) this->output_arcs.insert(id);
				}
			}
		}

		bool mPass[3];
	public:

		void SetPassingIndices(bool arcs01, bool arcs12, bool arcs23) {
			invalidate();
			mPass[0] = arcs01;
			mPass[1] = arcs12;
			mPass[2] = arcs23;
		}
	};

	template<class MSCType>
	class DFO_MSCIncidentArcs : public DFO_MSCArcFilter<MSCType> {
	protected:
		virtual void ArcFilterAction() {
			printf("DFO_MSCIncidentArcs::ArcFilterAction - computing\n");
			MSCType* msc = this->mMsc->GetMsc();

			for (auto input : this->input_objects) {
				for (auto id : dynamic_cast<DFO_MSCNodeFilter<MSCType>*>(input)->GetOutput()) {
					//if (mPass[mMsc->GetMsc()->getArc(id).dim]) this->output_arcs.insert(id);

					for (INT_TYPE aid = msc->firstIncidentLivingArc(id);
						msc->isValidArcId(aid);
						aid = msc->nextIncidentLivingArc(aid, id)) {
						if (msc->arcLowerNode(aid) == id) {
							if (mPassLower) this->output_arcs.insert(aid);
						}
						else {
							if (mPassUpper) this->output_arcs.insert(aid);
						}
						
					}

				}
			}
		}
		bool mPassLower;
		bool mPassUpper;
	public:
		DFO_MSCIncidentArcs() : mPassLower(true), mPassUpper(true) {}

		void SetPassingArcs(bool lower, bool upper) {
			mPassLower = lower;
			mPassUpper = upper;
		}
	};

	template<class MSCType>
	class DFO_MSCIncidentNodes : public DFO_MSCNodeFilter<MSCType> {
	protected:
		virtual void NodeFilterAction() {
			printf("DFO_MSCIncidentNodes::NodeFilterAction - computing\n");
			MSCType* msc = this->mMsc->GetMsc();
			for (auto input : this->input_objects) {
				for (auto id : dynamic_cast<DFO_MSCArcFilter<MSCType>*>(input)->GetOutput()) {
					if (mPassLower) this->output_nodes.insert(msc->arcLowerNode(id));
					if (mPassUpper) this->output_nodes.insert(msc->arcUpperNode(id));
				}
			}
		}
		bool mPassLower;
		bool mPassUpper;
	public:
		DFO_MSCIncidentNodes() : mPassLower(true), mPassUpper(true) {}
		void SetPassingNodes(bool lower, bool upper) {
			mPassLower = lower;
			mPassUpper = upper;
		}
	};

	// every valuator of elements of the morse-smale complex has to have the ms complex
	// as something to refer to. all valuators will be subclasses of this class, which simply
	// provides storage and set method for the complex
	template<class MSCType, typename DType>
	class DFO_MSCIdValuator : public DFO_IdValuator < DType > {
	protected:
		DFO_MSCHierarchyInterface<MSCType>* mMsc;
	public:
		virtual void SetMscDFO(DFO_MSCHierarchyInterface<MSCType>* msc) {
			mMsc = msc;
			add_parent(msc);
		}
	};

	template<class MSCType>
	class DFO_MSCNodePositionValuator : public DFO_MSCIdValuator<MSCType, Vec3d> {
	public:
		virtual Vec3d GetValue(INT_TYPE id) const {
			INDEX_TYPE ind = this->mMsc->GetMsc()->getNode(id).cellindex;
			Vec3l resl;
			this->mMsc->GetMsc()->GetMesh()->cellid2Coords(ind, resl);
			return resl;
		}
	};

	// access to function values of the msc, Dtype must match MSCType's scalar type
	template<class MSCType, typename DType>
	class DFO_MSCNodeFValValuator : public DFO_MSCIdValuator<MSCType, DType> {
	public:
		virtual DType GetValue(INT_TYPE id) const {
			return this->mMsc->GetMsc()->getNode(id).value;
		}
	};

	template<class MSCType, typename DType>
	class DFO_MSCArcLowerFValValuator : public DFO_MSCIdValuator<MSCType, DType> {
	public:
		virtual DType GetValue(INT_TYPE id) const {
			 MSCType* msc = this->mMsc->GetMsc();
			INT_TYPE lower_node_id = msc->getArc(id).lower;
			return msc->getNode(lower_node_id).value;
		}
	};
	template<class MSCType, typename DType>
	class DFO_MSCArcUpperFValValuator : public DFO_MSCIdValuator<MSCType, DType> {
	public:
		virtual DType GetValue(INT_TYPE id) const {
			 MSCType* msc = this->mMsc->GetMsc();
			INT_TYPE lower_node_id = msc->getArc(id).upper;
			return msc->getNode(lower_node_id).value;
		}
	};
	template<class MSCType, typename DType>
	class DFO_MSCArcPersistenceFValValuator : public DFO_MSCIdValuator<MSCType, DType> {
	public:
		virtual DType GetValue(INT_TYPE id) const {
			MSCType* msc = this->mMsc->GetMsc();
			return msc->getArc(id).persistence;
		}
	};

	template<class MSCType>
	class DFO_MSCArcLowerDimValuator : public DFO_MSCIdValuator<MSCType, DIM_TYPE> {
	public:
		virtual DIM_TYPE GetValue(INT_TYPE id) const {
			MSCType* msc = this->mMsc->GetMsc();
			return msc->getArc(id).dim;
		}
	};

	template<class MSCType>
	class DFO_MSCArcBoundaryValuator : public DFO_MSCIdValuator<MSCType, BOUNDARY_TYPE> {
	public:
		virtual BOUNDARY_TYPE GetValue(INT_TYPE id) const {
			MSCType* msc = this->mMsc->GetMsc();
			return msc->getArc(id).boundary;
		}
	};
	template<class MSCType>
	class DFO_MSCNodeCellIndexValuator : public DFO_MSCIdValuator<MSCType, INDEX_TYPE> {
	public:
		virtual INDEX_TYPE GetValue(INT_TYPE id) const {
			return this->mMsc->GetMsc()->getNode(id).cellindex;
		}
	};
	template<class MSCType>
	class DFO_MSCNodeDimValuator : public DFO_MSCIdValuator<MSCType, DIM_TYPE> {
	public:
		virtual DIM_TYPE GetValue(INT_TYPE id) const {
			return this->mMsc->GetMsc()->getNode(id).dim;
		}
	};
	template<class MSCType>
	class DFO_MSCNodeBoundaryValuator : public DFO_MSCIdValuator<MSCType, BOUNDARY_TYPE> {
	public:
		virtual BOUNDARY_TYPE GetValue(INT_TYPE id) const {
			return this->mMsc->GetMsc()->getNode(id).boundary;
		}
	};
	template<class MSCType>
	class DFO_MSCArcGeomListValuator : public DFO_MSCIdValuator<MSCType, vector<INDEX_TYPE> > {
	public:
		virtual  vector<INDEX_TYPE> GetValue(INT_TYPE id) const {
			vector<INDEX_TYPE> v;
			this->mMsc->GetMsc()->fillArcGeometry(id, v);
			return v;
		}
	};
	template<class MSCType>
	class DFO_MSCNodeAscManValuator : public DFO_MSCIdValuator<MSCType, set<INDEX_TYPE> > {
	public:
		virtual  set<INDEX_TYPE> GetValue(INT_TYPE id) const {
			set<INDEX_TYPE> v;
			this->mMsc->GetMsc()->fillGeometry(id, v, true);
			return v;
		}
	};
	template<class MSCType>
	class DFO_MSCNodeDscManValuator : public DFO_MSCIdValuator<MSCType, set<INDEX_TYPE> > {
	public:
		virtual  set<INDEX_TYPE> GetValue(INT_TYPE id) const {
			set<INDEX_TYPE> v;
			this->mMsc->GetMsc()->fillGeometry(id, v, false);
			return v;
		}
	};
};

#endif
