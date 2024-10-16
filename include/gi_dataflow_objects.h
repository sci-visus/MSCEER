#ifndef GI_DATA_FLOW_OBJECTS_H
#define GI_DATA_FLOW_OBJECTS_H

#include <type_traits>

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
#include "gi_regular_grid_trilinear_function.h"


namespace GInt {

	// a base class for objects having the property that they are part of a directed acyclic graph 
	// for computing, where each object has inputs, outputs, and does some action to transform input 
	// data to output data. Furthermore, each object can be either in a state of valid, or invalid,
	// in order to propagate changes, e.g. in parameters, downstream. Each subclass implements Action
	// which ultimately will do the work of the subclass, and is called by compute_output().
	// -- each dataflow object owns its output - children will only access by const reference
	// -- a child can assume that when its Action() call is called, all its parents are valid
	class DataflowObject {
	protected:
		vector<DataflowObject*> parents; // any data that this object may depend on
		vector<DataflowObject*> children;  // any objects that may depend on this object
		// the actual work of a dataflow object happens in teh action call.
		// when this is called parents have been computed
		virtual void Action() {};
		bool output_valid;		
		size_t add_child(DataflowObject* s) {
			children.push_back(s);
			invalidate();
			return children.size();
		}
	public:
		virtual size_t add_parent(DataflowObject* s) {
			s->add_child(this);
			parents.push_back(s);
			return parents.size();
		}	
		DataflowObject() : output_valid(false) {}
		// invalidate the result of this, and all downstream objects
		virtual void invalidate() {
			if (!output_valid) return; // to prevent infinite looping by accident
			output_valid = false;
			for (auto it = children.begin(); it != children.end(); it++) (*it)->invalidate();
		}
		// force computation of this, and all objects it depends on
		void compute_output() {
			if (output_valid) return;
			for (auto it = parents.begin(); it != parents.end(); it++) {
				(*it)->compute_output();
			}
			Action(); // do local work for class
			output_valid = true;
			return;
		}
	};

	// a special class for dataflow objects, where the result of the action 
	// populates a set of integers - typically a set of indices. 
	class IdSetDataflowObject : public DataflowObject {
	public:
		virtual const set<INT_TYPE>& GetOutput() const = 0;
	};
	


	// utility dataflow object which takes an index set (which can be sparse in the index space)
	// and creates a map to a dense 0...size index set. the id-to-position map itself is available through the 
	// getoutputmap call, and position-to-id mapping through the getoutputlist call.
	class DFO_IdToPosition : public DataflowObject {
	protected: 

		IdSetDataflowObject* mInput;

		virtual void Action() {
			mPosMap.clear();
			INT_TYPE pos = 0;
			for (auto id : mInput->GetOutput()) {
				mIdList.push_back(id);
				mPosMap[id] = pos++;
			}
		}
		unordered_map<INT_TYPE, INT_TYPE> mPosMap;
		vector<INT_TYPE> mIdList;
	public:
		void SetInput(IdSetDataflowObject* input) {
			if (parents.size() > 0) {
				printf("STATIC_WARNING::DFO_IdToPosition::SetInput - only uses one input, ignoring call\n");
				return;
			}
			this->add_parent(input);
			mInput = input;
		}

		virtual const unordered_map<INT_TYPE, INT_TYPE>& GetOutputMap() const {
			return mPosMap;
		}
		virtual const vector<INT_TYPE>& GetOutputList() const {
			return mIdList;
		}
	};


	template<typename DType>
	class DFO_Valuator : public DataflowObject {

	public:
		virtual DType GetValue() const = 0;
	};

	template<typename DType>
	class DFO_ConstValuator : public DFO_Valuator<DType> {
	protected:
		const DType mValue;
	public:
		DFO_ConstValuator(DType value) : mValue(value){}
		virtual DType GetValue() const {
			return mValue;
		}
	};

	typedef DFO_ConstValuator<float> DFO_ConstFloat;
	typedef DFO_ConstValuator<double> DFO_ConstDouble;
	typedef DFO_ConstValuator<int> DFO_ConstInt;
	typedef DFO_ConstValuator<unsigned int> DFO_ConstUInt;
	typedef DFO_ConstValuator<INDEX_TYPE> DFO_ConstIndexType;
	typedef DFO_ConstValuator<BOUNDARY_TYPE> DFO_ConstBoundaryType;
	typedef DFO_ConstValuator<DIM_TYPE> DFO_ConstDimType;


	template<typename DType>
	class DFO_VaryingValuator : public DFO_Valuator<DType> {
	protected:
		DType m_value;
	public:
		virtual DType GetValue() const { return m_value; }
		virtual void SetValue(DType value) {
			invalidate();
			m_value = value;
		}
	};

	typedef DFO_VaryingValuator<double> DFO_VaryingValuatorDouble;
	typedef DFO_VaryingValuator<float> DFO_VaryingValuatorFloat;
	typedef DFO_VaryingValuator<int> DFO_VaryingValuatorInt;

	template<typename DType1, typename DType2>
	class DFO_UnaryValuator : public DataflowObject {
	public:
		virtual DType2 GetValue(DType1 input) const = 0;
	};
	
	template<typename DType>
	class DFO_IdValuator : public DFO_UnaryValuator<INT_TYPE, DType> {
	public:
		virtual DType GetValue(INT_TYPE id) const = 0;
	};

	
	


	template<class ValueType>
	class DFO_IdValueArray : public DataflowObject {
	protected:
		ValueType* mValues;

		void DestroyInternals( size_t idsize) {}
		//void DestroyInternals(ValueType* values, size_t idsize) {}

		//inline typename std::enable_if<std::is_pointer<ValueType>::value, void>::type DestroyInterals(size_t idsize)
		//{
		//	for (size_t i = 0; i < idsize; i++) delete mValues[i];
		//}

		virtual void Action() {
			//ADD TESTS HERE
			const vector<INT_TYPE>& ids = mIdMapper->GetOutputList();
			size_t idsize = ids.size();

			if (mValues != NULL) {
				DestroyInternals(idsize);
				delete[] mValues;
				mValues = new ValueType[idsize];
			}

#pragma omp parallel for
			for (int i = 0; i < idsize; i++) {
				mValues[i] = mValuator->GetValue(ids[i]);
			}
		}
		DFO_IdToPosition* mIdMapper;
		DFO_IdValuator<ValueType>* mValuator;
	public:
		DFO_IdValueArray() : mValues(NULL), mIdMapper(NULL), mValuator(NULL) {}

		virtual void SetMapperInput(DFO_IdToPosition* id_mapper) {
			mIdMapper = id_mapper;
			this->add_parent(mIdMapper);
			invalidate();
		}
		virtual void SetValuatorInput(DFO_IdValuator<ValueType>* id_valuator) {
			mValuator = id_valuator;
			this->add_parent(mValuator);
			invalidate();
		}

		size_t GetPosition(INT_TYPE id) const { return mIdMapper->GetOutputMap()[id]; }
		const ValueType* GetArray() const { return mValues; }

	};

	template<class TopoGridType>
	class DFO_TopoGrid : public DataflowObject {
		TopoGridType* mGrid;
	public:
		virtual void SetTopoGrid(TopoGridType* grid) {
			mGrid = grid;
			invalidate();
		}
		const TopoGridType* GetGrid() const {
			return mGrid;
		}
	};


	template<class TopoGridType, typename DType>
	class DFO_TopoIdConverter : public DFO_UnaryValuator<INDEX_TYPE, DType>{
	protected:
		DFO_TopoGrid<TopoGridType>* mTopoGridDFO;
	public:
		virtual void SetTopoGridDFO(DFO_TopoGrid<TopoGridType>* grid) {
			mTopoGridDFO = grid;
		}
		virtual DType GetValue(INDEX_TYPE id) const = 0;
	};

	template<class TopoGridType>
	class DFO_TopoIdToGridIdValuator : public DFO_TopoIdConverter<TopoGridType, INDEX_TYPE> {
	public:
		virtual INDEX_TYPE GetValue(INDEX_TYPE id) const {
			return mTopoGridDFO->GetGrid()->VertexNumberFromCellID(id);
		}
	};

	template<class TopoGridType>
	class DFO_TopoIdToTopoCoordsValuator : public DFO_TopoIdConverter<TopoGridType, Vec3d> {
	public:
		virtual Vec3d GetValue(INDEX_TYPE id) const {
			Vec3d res;
			mTopoGridDFO->GetGrid()->centroid(id, &res);
			return res;
		}
	};

	template<class TopoGridType>
	class DFO_TopoIdBoundaryValuator : public DFO_TopoIdConverter<TopoGridType, BOUNDARY_TYPE> {
	public:
		virtual BOUNDARY_TYPE GetValue(INDEX_TYPE id) const {
			return mTopoGridDFO->GetGrid()->boundaryValue(id);
		}
	};
	template<class TopoGridType>
	class DFO_TopoIdDimValuator : public DFO_TopoIdConverter<TopoGridType, BOUNDARY_TYPE> {
	public:
		virtual BOUNDARY_TYPE GetValue(INDEX_TYPE id) const {
			return mTopoGridDFO->GetGrid()->dimension(id);
		}
	};
	typedef DFO_IdValuator<bool> DFO_IdTest;



	template<typename DType>
	class DFO_FunctionValueSampler : public DFO_UnaryValuator<INDEX_TYPE, DType> {
	};


	template<typename DType>
	class DFO_RegularGridFunctionSampler : public DFO_FunctionValueSampler<DType> {
	protected:
		RegularGridTrilinearFunction* mFunc;
	public:
		virtual void SetInputFunction(RegularGridTrilinearFunction* func) {
			mFunc = func;
			invalidate();
		}
		virtual DType GetValue(INDEX_TYPE pos) const {
			return mFunc->SampleImage(pos);
		}

	};


	template<typename DType>
	class DFO_RegularGridTopoFunctionSampler : public DFO_FunctionValueSampler<DType> {
	protected:
		TopologicalExplicitDenseMeshFunction* mFunc;
	public:
		virtual void SetInputFunction(TopologicalExplicitDenseMeshFunction* func) {
			mFunc = func;
			invalidate();
		}
		virtual DType GetValue(INDEX_TYPE pos) {
			return mFunc->cellValue(pos);
		}

	};
	//template<class GridType>
	template<typename DType>
	class DFO_FunctionInterpolatingSampler : public DFO_UnaryValuator<Vec3d, DType> {
	};

	template<typename DType>
	class DFO_RegularGridTrilinearFunctionSampler : public DFO_FunctionInterpolatingSampler<DType> {
	protected:
		RegularGridTrilinearFunction* mFunc;
	public:
		virtual void SetInputFunction(RegularGridTrilinearFunction* func) {
			mFunc = func;
			invalidate();
		}
		virtual DType GetValue(Vec3d pos) {
			return mFunc->SampleImage(pos);
		}

	};

	template<typename Vec3d>
	class DFO_RegularGridTrilinearGradSampler : public DFO_FunctionInterpolatingSampler<Vec3d> {
	protected:
		RegularGridTrilinearFunction* mFunc;
	public:
		virtual void SetInputFunction(RegularGridTrilinearFunction* func) {
			mFunc = func;
			invalidate();
		}
		virtual Vec3d GetValue(Vec3d pos) {
			return mFunc->SampleGrad(pos);
		}
	};
	//template<class MSCType>
	//class MSCSelectorLivingNodes : public Selector {
	//protected:
	//	MSCType* mMsc;
	//	virtual void SelectorAction() {
	//		//printf("MSCSelectorLivingNodes::SelectorAction() computing...\n");
	//		for (int i = 0; i < mMsc->numNodes(); i++) {
	//			if (mMsc->isNodeAlive(i)) this->output.insert(i);
	//		}
	//		//printf("MSCSelectorLivingNodes::SelectorAction produced %d indices\n", output.size());
	//	}
	//public:
	//	MSCSelectorLivingNodes(MSCType* msc) : mMsc(msc) {}
	//};

	//template<class MSCType>
	//class MSCSelectorNodeIndex : public Selector {
	//protected:
	//	MSCType* mMsc;
	//	DIM_TYPE mIndex;
	//	virtual void SelectorAction() {
	//		//printf("MSCSelectorNodeIndex::SelectorAction() computing...\n");
	//		for (auto pit = parents.begin(); pit != parents.end(); pit++) {
	//			for (auto it = (*pit)->output.begin(); it != (*pit)->output.end(); it++) {
	//				if (mMsc->getNode(*it).dim == mIndex) this->output.insert(*it);
	//			}
	//		}
	//		//printf("MSCSelectorNodeIndex::SelectorAction produced %d indices\n", output.size());
	//	}
	//public:
	//	MSCSelectorNodeIndex(MSCType* msc, DIM_TYPE index) : mMsc(msc), mIndex(index) {}
	//};

	//template<class MSCType>
	//class MSCSelectorRepresentative1Saddle : public Selector {
	//protected:
	//	MSCType* mMsc;
	//	DIM_TYPE mIndex;
	//	pair<INT_TYPE, INT_TYPE> getExtremumPair(INT_TYPE saddle) {
 //           typename MSCType::SurroundingArcsIterator sit(mMsc);
	//		INT_TYPE extrema[2]; int numext = 0;
	//		for (sit.begin(saddle); sit.valid(); sit.advance()) {
	//			//printf("asdf %d\n", sit.value());
	//			INT_TYPE aid = sit.value();
	//			arc<FLOATTYPE>& a = mMsc->getArc(aid);
	//			if (a.upper == saddle) extrema[numext++] = a.lower;
	//		}
	//		//printf("done numext=%d\n", numext);
	//		if (numext == 1) {
	//			return pair<INT_TYPE, INT_TYPE>(extrema[0], extrema[0]);
	//		}
	//		else {
	//			if (extrema[0] < extrema[1]) {
	//				return pair<INT_TYPE, INT_TYPE>(extrema[0], extrema[1]);
	//			} 
	//			else {
	//				return pair<INT_TYPE, INT_TYPE>(extrema[1], extrema[0]);
	//			}

	//		}

	//	}
	//	bool lessthan(INT_TYPE a, INT_TYPE b) {
	//		node<FLOATTYPE>& na = mMsc->getNode(a);
	//		node<FLOATTYPE>& nb = mMsc->getNode(b);
	//		if (na.value == nb.value) return a < b;
	//		return na.value < nb.value;
	//	}

	//	virtual void SelectorAction() {
	//		//printf("MSCSelectorRepresentative1Saddle::SelectorAction() computing...\n");
	//		map<pair<INT_TYPE, INT_TYPE>, INT_TYPE> saddlemap;

	//		for (auto pit = parents.begin(); pit != parents.end(); pit++) {
	//			for (auto it = (*pit)->output.begin(); it != (*pit)->output.end(); it++) {
	//				
	//				INT_TYPE saddle = *it;
	//				//printf("doing %d\n", saddle);
	//				pair<INT_TYPE, INT_TYPE> p = getExtremumPair(saddle);

	//				if (saddlemap.count(p) == 0) {
	//					saddlemap[p] = saddle;
	//				}
	//				else {
	//					INT_TYPE othersaddle = saddlemap[p];
	//					if (lessthan(saddle, othersaddle)) {
	//						saddlemap[p] = saddle;
	//					}
	//				}
	//			}
	//		}
	//		//printf("b\n");
	//		for (auto it = saddlemap.begin(); it != saddlemap.end(); it++) {
	//			if ((*it).first.first != (*it).first.second)
	//				output.insert((*it).second);
	//		}
	//		//printf("MSCSelectorRepresentative1Saddle::SelectorAction produced %d indices\n", output.size());
	//	}

	//public:
	//	MSCSelectorRepresentative1Saddle(MSCType* msc) : mMsc(msc) {}
	//};

}

#endif
