#ifndef MC_LIGHT_GEOM_MSC
#define MC_LIGHT_GEOM_MSC

#include <stdio.h>
#include <vector>
#include <map>
#include <queue>
#include <set>


#ifdef __APPLE__
#include <unordered_map>
#else
#include <unordered_map>
#endif

#include "gi_basic_types.h"
#include "gi_topological_explicit_mesh_function.h"
#include "gi_discrete_gradient_labeling.h"
#include "gi_topological_regular_grid_restricted.h"


#define MAX_CELLID_VALUE 999999999
//#define WRITEREMAP( x ) (x==-1?-1:x+1)
//#define datap unsigned long long
//#define datap unsigned long long
#define UINT_INFTY -1
#define INT_INFTY 4294967295 /2
#define NULLID 4294967295/2

using namespace std;

// a class that computes an MSC (in parallel) with NO arc geometry, and a hierarchy (in serial) storing the merging of manifolds
namespace GInt {

	//typedef DiscreteGradientLabeling GRAD_TYPE;
	//typedef TopologicalRegularGrid3D MESH_TYPE;
	//typedef TopologicalExplicitDenseMeshFunction FUNC_TYPE;
	//typedef float SCALAR_TYPE;

	struct msbitfield
	{
		//// rehash this
		unsigned char dim : 3;
		unsigned char boundary : 1;
		unsigned char f1 : 1;
		unsigned char f2 : 1;
		unsigned char f3 : 2;
	};


	// nodes start out with -1 merged_manifolds - means that you use discrete grad to fill in the merged_manifold geometry
	// as cancellations go, we add to global merged_manifold arrays
	//
	struct merged_manifold {
		INT_TYPE merged[2]; // the merged_manifold ids
		INT_TYPE basenode; // the node id
		INT_TYPE mergetime;
	};


	// original arcs get their geometry explicitly, while new arcs from merging get different ones
	struct arc_base_geometry {
		vector<INDEX_TYPE> geometry;
	};
	struct arc_merged_geometry {
		INT_TYPE fields[3];
	};

	template<typename SCALAR_TYPE>
	struct node
	{
		INDEX_TYPE cellindex;
		INT_TYPE firstarc;
		INT_TYPE destroyed; // the time this node is cancelled.
		INT_TYPE amanifoldid; // set to -1 if this is the base
		INT_TYPE dmanifoldid; // set to -1 if this is the base
		unsigned short numarcs;
		unsigned short numlower;
		SCALAR_TYPE value;
		DIM_TYPE dim;
		BOUNDARY_TYPE boundary;
	};


	template<typename SCALAR_TYPE>
	struct arc
	{
		INT_TYPE lower; // node
		INT_TYPE lower_next; //arc - INT_INFTY is the null arc!
		INT_TYPE upper; // node
		INT_TYPE upper_next; //arc
		INT_TYPE created;
		INT_TYPE destroyed;
		INT_TYPE geom; // if created == 0, this is an original arc, and its geometry will be found in base_geom list of MSC, else in the merged_geom list
		SCALAR_TYPE persistence;
		DIM_TYPE dim; // the .dim stores the dim of the lower endpoint
		BOUNDARY_TYPE boundary;
	};
 
	template<typename SCALAR_TYPE, class MESH_TYPE, class FUNC_TYPE, class GRAD_TYPE>
	class MorseSmaleComplexBasic {
	public:


		typedef SCALAR_TYPE ScalarType;
		typedef MESH_TYPE MeshType;
		typedef FUNC_TYPE FuncType;
		typedef GRAD_TYPE GradType;

	protected:
		struct sortedEdge {
			SCALAR_TYPE persistence;
			int countweight;
			INT_TYPE ep;

			bool operator()(const sortedEdge& _Left, const sortedEdge& _Right) {
				if (_Left.persistence < _Right.persistence) return false;
				if (_Left.persistence > _Right.persistence) return true;
				if (_Left.countweight < _Right.countweight) return false;
				if (_Left.countweight > _Right.countweight) return true;
				return _Left.ep > _Right.ep;
			}
		};

		priority_queue<sortedEdge, vector< sortedEdge>, sortedEdge > edges_to_cancel;
	public:
		struct cancellation_record {
			int index;
			SCALAR_TYPE persistence;
			SCALAR_TYPE lval;
			SCALAR_TYPE uval;
			SCALAR_TYPE persPerc;
			INT_TYPE arcp;
			int boundary;
		};
	protected:

		vector<cancellation_record> mCRecords;
	public:
		const vector<cancellation_record>& GetCancellationRecords() {
			return mCRecords;
		}
	protected:

		bool output_cancellation_records;

		const char* mCRecordFName;
		float m_temp_perc_pers;
		//INDEX_TYPE select_persistence;
		INDEX_TYPE num_destroyed;
		GRAD_TYPE* mGrad;
		MESH_TYPE* mMesh;
		FUNC_TYPE* mFunc;
		Vec3b mBuildArcGeometry;
		Vec3b mBuildArcAtAll;

		vector<node<SCALAR_TYPE> > nodes;
		vector<arc<SCALAR_TYPE> > arcs;
		vector<arc_base_geometry> arc_base_geoms;
		vector<arc_merged_geometry> arc_merge_geoms;
		vector<merged_manifold> mans;

#ifdef WIN32
		std::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#else
#ifdef __APPLE__
		//__gnu_cxx::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
		std::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#else
		std::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#endif
#endif

		// this is NOT thread safe - uses references to items in a vector that could be reallocated
		inline void connectArc(INT_TYPE arcID, arc<SCALAR_TYPE>& a, node<SCALAR_TYPE>& lower, node<SCALAR_TYPE>& upper) {
			a.lower_next = lower.firstarc;
			lower.firstarc = arcID;
			a.upper_next = upper.firstarc;
			upper.firstarc = arcID;

			lower.numarcs++;
			upper.numarcs++;
			upper.numlower++;

		}

		void connectArc(INT_TYPE arcID) {

			arc<SCALAR_TYPE>& a = arcs[arcID];
			node<SCALAR_TYPE>& lower = nodes[a.lower];
			node<SCALAR_TYPE>& upper = nodes[a.upper];
			connectArc(arcID, a, lower, upper);
		}
		const merged_manifold& getManifold(INT_TYPE m) const {
			return mans[m];
		}

	public:

		INT_TYPE numArcs() const {
			return (INT_TYPE)arcs.size();
		}
		INT_TYPE numNodes() const  {
			return (INT_TYPE)nodes.size();
		}
		const node<SCALAR_TYPE>& getNode(INT_TYPE e) const {
			return nodes[e];
		}
		node<SCALAR_TYPE>& getNode(INT_TYPE e) {
			return nodes[e];
		}

		const arc<SCALAR_TYPE>& getArc(INT_TYPE e) const {
			return arcs[e];
		}
		arc<SCALAR_TYPE>& getArc(INT_TYPE e) {
			return arcs[e];
		}
		MESH_TYPE* const GetMesh() const {
			return mMesh;
		}

		void set_output_cancellation_records(const char* fname) {
			output_cancellation_records = true;
			mCRecordFName = fname;
		}


		MorseSmaleComplexBasic(GRAD_TYPE* grad,
			MESH_TYPE* mesh,
			FUNC_TYPE* func) :
			mGrad(grad), mMesh(mesh), mFunc(func),
			mBuildArcAtAll(Vec3b(true, true, true)), mBuildArcGeometry(Vec3b(true, true, true))
		{
			num_destroyed = 0;
			select_persistence = 0;
			output_cancellation_records = false;
			for (int i = 0; i < 4; i++) {
				countIDs[i] = 0;
				countBoundaryIds[i] = 0;
			}
		}

		// must be calleb before computefromgrad
		void SetBuildArcs(Vec3b v) { mBuildArcAtAll = v; }
		void SetBuildArcGeometry(Vec3b v) { mBuildArcGeometry = v; }

	protected:

		INT_TYPE createManifold(INT_TYPE baseId, INT_TYPE baseManId, INT_TYPE mergeManId, INT_TYPE mergetime) {
			INT_TYPE manid = mans.size();
			merged_manifold m;
			m.basenode = baseId;
			m.merged[0] = baseManId;
			m.merged[1] = mergeManId;
			m.mergetime = mergetime;
			mans.push_back(m);
			return manid;
		}

		INT_TYPE createNode(INDEX_TYPE cellID) {
			return createNode(cellID, mFunc->cellValue(cellID));
		}
		INT_TYPE createNode(INDEX_TYPE cellID, SCALAR_TYPE value) {
				node<SCALAR_TYPE> tn;
			INT_TYPE listID = nodes.size();
			//if (nodes.size() == nodes.capacity()) printf("about to expand capacity - nodes createnode!\b");
			nodes.push_back(tn);
			nodemap[cellID] = listID;

			//node<SCALAR_TYPE> &n = nodes->operator [](listID);
			node<SCALAR_TYPE> &n = nodes[listID];
			n.cellindex = cellID;

			n.destroyed = INT_INFTY;
			n.firstarc = NULLID;
			n.boundary = mMesh->boundaryValue(cellID);
			//printf("n.boundary = %d\n", n.boundary);
			//if (n.boundary) {
			//	INDEX_TYPE coords[3];
			//	sgg->getCoords(coords, cellID);
			//	printf("(%d, %d, %d) = %d\n",
			//	   (int) coords[0], (int) coords[1], (int) coords[2], n.boundary);
			//}
			n.dim = mMesh->dimension(cellID);
			this->countIDs[n.dim]++;
			if (n.boundary) this->countBoundaryIds[n.dim]++;
			n.numarcs = 0;
			n.numlower = 0;
			n.value = value;

			n.amanifoldid = createManifold(listID, -1, -1, 0);
			n.dmanifoldid = createManifold(listID, -1, -1, 0);

			return listID;
		}

		//// create an arc connecting lower to upper <- cell id's from gradient
		//INT_TYPE createArc(arc<SCALAR_TYPE> &a){
		//	INT_TYPE listID = (INT_TYPE)arcs.push_back(a);
		//	//printf("%f - %f\n", (float) a.persistence, arcs[listID].persistence);
		//	return listID;
		//}

		INT_TYPE createArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID, std::vector<INDEX_TYPE>& geometry) {

			arc<SCALAR_TYPE> at;
			INT_TYPE listID = arcs.size();
			arcs.push_back(at);
			INT_TYPE geomID = this->arc_base_geoms.size();
			arc_base_geoms.push_back(arc_base_geometry());

			arc<SCALAR_TYPE>& a = arcs[listID];
			arc_base_geometry& ge = arc_base_geoms[geomID];
			ge.geometry.insert(ge.geometry.begin(), geometry.begin(), geometry.end());



			a.created = 0;
			a.destroyed = INT_INFTY;
			a.lower = nodemap[lowerCellID];
			a.upper = nodemap[upperCellID];
			a.geom = geomID;
			if (nodes[a.upper].value < nodes[a.lower].value)
				printf("ERROR: upper (%f, %d) - lower (%f, %d)\n",
					float(nodes[a.upper].value), nodes[a.upper].dim,
					float(nodes[a.lower].value), nodes[a.lower].dim);

			a.persistence = nodes[a.upper].value - nodes[a.lower].value;
			a.dim = nodes[a.lower].dim;
			connectArc(listID);

			return listID;
		}

		INT_TYPE createArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID) {
			std::vector<INDEX_TYPE> a;
			return createArc(lowerCellID, upperCellID, a);
		}
		virtual void InsertArcIntoSimplification(arc<SCALAR_TYPE>& a, INT_TYPE id) {

			if (a.persistence <= gPersThreshold) {
				sortedEdge se;

				se.persistence = a.persistence;
				se.countweight = edgeCountWeight(a);
				se.ep = id;
				edges_to_cancel.push(se);
			}
			//if (arcs.size() == arcs.capacity()) {


		}


		INT_TYPE createArc(INT_TYPE luap, INT_TYPE ma, INT_TYPE ulap, INT_TYPE ctime) {
			arc<SCALAR_TYPE>& lua = arcs[luap];
			arc<SCALAR_TYPE>& ula = arcs[ulap];
			INT_TYPE na_id = arcs.size();
			INT_TYPE na_geom_id = arc_merge_geoms.size();
			arc_merge_geoms.push_back(arc_merged_geometry());
			arc_merged_geometry& na_geom = arc_merge_geoms[na_geom_id];
			na_geom.fields[0] = luap;
			na_geom.fields[1] = ma;
			na_geom.fields[2] = ulap;

			//this->arcs.push_back(ma); // copy setting from old
			arc<SCALAR_TYPE> na;// = this->arcs[na_id];

			na.created = ctime;
			na.destroyed = INT_INFTY;
			na.lower = ula.lower;
			na.upper = lua.upper;
			na.geom = na_geom_id;

			node<SCALAR_TYPE>& nup = this->nodes[na.upper];
			node<SCALAR_TYPE>& nlo = this->nodes[na.lower];

			this->connectArc(na_id, na, nlo, nup);
			na.persistence = nup.value - nlo.value;
			na.boundary = nlo.boundary + nup.boundary;
			na.dim = nlo.dim;
			if (na.persistence < 0 && nlo.boundary == nup.boundary) printf("creatinga inversion %f, %d-%d, b%d-b%d\n",
				(float)na.persistence, nlo.dim, nup.dim, nlo.boundary, nup.boundary);

			// now insert into global sort if it has a chance of being cancelled
			InsertArcIntoSimplification(na, na_id);
			//	printf("about to expand capacity2 - arcs createnode!\n");
			//}
			this->arcs.push_back(na); // this could destroy references so put at end;
			return na_id;

		}
		//void rec_tdcr_no_geom(const INDEX_TYPE& cellid, DIM_TYPE& temp_dim, const INDEX_TYPE start) {
		//	INDEX_TYPE current = cellid;
		//	MESH_TYPE::FacetsIterator facets(mMesh);
		//	for (facets.begin(current); facets.valid(); facets.advance()) {
		//		INDEX_TYPE temp_id = facets.value();
		//		if (mGrad->getCritical(temp_id)) {
		//			//printf("adding arc: %llu\n", temp_id);
		//			createArc(temp_id, start);

		//		}
		//		else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
		//			INDEX_TYPE pair = mGrad->getPair(temp_id);
		//			if (pair != cellid && mMesh->dimension(pair) == mMesh->dimension(cellid)) {

		//				//result.push_back(pair);
		//				rec_tdcr_no_geom(pair, temp_dim, start);

		//			}
		//		}

		//	}
		//}
		//void trace_down_cells_restricted_no_geom(const INDEX_TYPE& cellid) {

		//	DIM_TYPE temp_dim = mGrad->getDimAscMan(cellid) + 1;
		//	rec_tdcr_no_geom(cellid, temp_dim, cellid);
		//}

        void rec_tdcr(const INDEX_TYPE& cellid, DIM_TYPE& temp_dim, const INDEX_TYPE start, vector<INDEX_TYPE>& geom) {
			INDEX_TYPE current = cellid;

            geom.push_back(current);
            typename MESH_TYPE::FacetsIterator facets(mMesh);
            for (facets.begin(current); facets.valid(); facets.advance()) {
				INDEX_TYPE temp_id = facets.value();
				geom.push_back(temp_id);
				if (mGrad->getCritical(temp_id)) {
#pragma omp critical
                    {
						createArc(temp_id, start, geom);
                    }

				}
				else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
					INDEX_TYPE pair = mGrad->getPair(temp_id);
					if (pair != cellid && mMesh->dimension(pair) == mMesh->dimension(cellid)) {
                        rec_tdcr(pair, temp_dim, start, geom);

					}
				}
				geom.pop_back();

            }
			geom.pop_back();
		}

        void rec_tdcr(const INDEX_TYPE& cellid, DIM_TYPE& temp_dim, const INDEX_TYPE start) {
            INDEX_TYPE current = cellid;

            typename MESH_TYPE::FacetsIterator facets(mMesh);
            for (facets.begin(current); facets.valid(); facets.advance()) {
                INDEX_TYPE temp_id = facets.value();
				if (mGrad->getCritical(temp_id)) {
#pragma omp critical
                    {
						createArc(temp_id, start);
                    }

				}
				else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
					INDEX_TYPE pair = mGrad->getPair(temp_id);
					if (pair != cellid && mMesh->dimension(pair) == mMesh->dimension(cellid)) {

                        rec_tdcr(pair, temp_dim, start);

					}
				}

            }
		}

        void trace_down_cells_restricted(INDEX_TYPE cellid) {
            vector<INDEX_TYPE> geom;
			DIM_TYPE temp_dim = mGrad->getDimAscMan(cellid) + 1;
            rec_tdcr(cellid, temp_dim, cellid, geom);
		}
        void trace_down_cells_restricted_nogeom(INDEX_TYPE cellid) {
			DIM_TYPE temp_dim = mGrad->getDimAscMan(cellid) + 1;
            rec_tdcr(cellid, temp_dim, cellid);
		}

        void AddArcs(INT_TYPE nodeId) {
			node<SCALAR_TYPE>& n = getNode(nodeId);
			if (n.dim == 0) return;
			if (mBuildArcAtAll[n.dim - 1]) {
				if (mBuildArcGeometry[n.dim - 1]) {
                    trace_down_cells_restricted(nodes[nodeId].cellindex);
				}
				else {
                    trace_down_cells_restricted_nogeom(nodes[nodeId].cellindex);
				}
			}
		}

	public:
		void ComputeFromGrad(bool restricted = false) {
			typename MESH_TYPE::AllCellsIterator t_cells(mMesh);
			printf("   -- Finding critical points...");
			fflush(stdout);
			int count = 0;
			for (t_cells.begin(); t_cells.valid(); t_cells.advance()) {
				INDEX_TYPE t_id = t_cells.value();
				if (mGrad->getCritical(t_id)) {
					createNode(t_id);
					count++;
				}
			}
			printf(" done! found %d critical points\n", count);
			printf("   -- Finding arcs...");
            fflush(stdout);
			// then add arcs
            INT_TYPE numnodes = nodes.size();
            INT_TYPE j;
#pragma omp parallel for// private(numnodes)
            for (j = 0; j < numnodes; j++) {
                AddArcs(j);
			}
            printf(" done! found %d arcs\n", this->arcs.size());
		}



	protected:

		//////////// NOW CANCELLATION STUFF

		int countUpperArcs(INT_TYPE id) const {
			return this->nodes[id].numarcs - this->nodes[id].numlower;
		}
		int countLowerArcs(INT_TYPE id) const {
			return this->nodes[id].numlower;
		}

		int edgeCountWeight(arc<SCALAR_TYPE> &a) const {
			INT_TYPE lower = a.lower;
			INT_TYPE upper = a.upper;
			int nl = countUpperArcs(lower);
			int nu = countLowerArcs(upper);

			////first option: return the number created
			return (nl - 1) * (nu - 1);

			//2nd option: return the change in the number of arcs
			return (nl - 1) * (nu - 1) - (nl + nu - 1);

			//return 1;
		}

		inline INT_TYPE current_pers() const { return select_persistence; }
		inline bool isAlive(const arc<SCALAR_TYPE>&a, INT_TYPE place) const {
			return a.created <= place && a.destroyed > place;
		}
	public:
		bool isNodeAlive(INT_TYPE n) const {
			return isAlive(this->nodes[n], select_persistence);
		}
		bool isArcAlive(INT_TYPE a) const {
			bool res = isAlive(this->arcs[a], select_persistence);
			if (res && !isNodeAlive(this->arcs[a].lower)) printf("ERROR LOWER NODE NOT ALIVE\n");
			if (res && !isNodeAlive(this->arcs[a].upper)) printf("ERROR UPPER NODE NOT ALIVE\n");
			return res;
		}

	protected:
		inline bool isAlive(const node<SCALAR_TYPE>&n, INT_TYPE place) const {
			//printf("n.destroyed = %d, place = %d\n", n.destroyed, place);
			return  n.destroyed > place;
		}
		inline bool isAlive(INT_TYPE n, INT_TYPE place) const {
			return isAlive(this->nodes[n], place);
		}


	public:
		SCALAR_TYPE GetSelectPersAbs() {
			return select_persistence;
		}
		void SetSelectPersAbs(SCALAR_TYPE value) {
			//printf("mcLightGeomMSC::SuperLightMSC::SetSelectPersAbs -> %f\n", value);
			int offset = 1;// cancel_num_to_pers.size() - 1;
			for (int i = 0; i < cancel_num_to_pers.size(); i += offset) {
				if (cancel_num_to_pers[i] > value) {
					select_persistence = i;

					return;
				}
				//while (i + offset > cancel_num_to_pers.size() - 1) offset = offset / 2;
				//while (offset > 1 && cancel_num_to_pers[i + offset] > value) offset = offset / 2;
			}
			select_persistence = cancel_num_to_pers.size() /*- 1*/;
		}

	protected:
		int cancel(INT_TYPE a) {

			int createcounter = 0;
			int deletecounter = 0;

			//int tmp1 = (int) edges_to_cancel.size();

			arc<SCALAR_TYPE>* ap = &(this->arcs[a]);
			node<SCALAR_TYPE>& lower = this->nodes[ap->lower];
			node<SCALAR_TYPE>& upper = this->nodes[ap->upper];

			//if (output_cancellation_records) {
			cancellation_record cr;
			cr.index = lower.dim;
			cr.lval = lower.value;
			cr.uval = upper.value;
			cr.persistence = ap->persistence;
			//SCALAR_TYPE diff = this->sgg->sgd->maxval - this->sgg->sgd->minval;
			//cr.persPerc = 100.0f * ap->persistence / diff;
			cr.arcp = a;
			cr.boundary = lower.boundary + upper.boundary;
			this->mCRecords.push_back(cr);
			//}



			if (lower.destroyed != INT_INFTY) printf("Error: MorseSmaleComplexBasic::cancel(%d) - lower.destroyed != INT_INFTY\n", a);
			if (upper.destroyed != INT_INFTY) printf("Error: MorseSmaleComplexBasic::cancel(%d) - upper.destroyed != INT_INFTY\n", a);
			//printf("\n");
			//printNode(ap->lower);
			//printNode(ap->upper);

			int initialguess = (lower.numarcs - lower.numlower - 1) * (upper.numlower - 1);
			int init2 = lower.numarcs + upper.numarcs - 1;

			if (ap->persistence > max_pers_so_far) max_pers_so_far = ap->persistence;



			// for each upwards arc connected to the lower node,
			// for each downwards arc connected to the upper node,
			// create a new connection from the other endpoints
			INT_TYPE la = lower.firstarc;
			while (la != NULLID) {
				arc<SCALAR_TYPE>& lap = this->arcs[la];

				// skip the arc itself
				if (la == a) {
					la = lap.lower_next; // we're guarantee this is the next
					continue;
				}

				//test if they are the same kind
				if (ap->lower != lap.lower) {
					la = lap.upper_next;
					continue;
				}
				if (!isAlive(lap, num_cancelled)) {
					la = lap.lower_next;
					continue;
				}

				INT_TYPE ua = upper.firstarc;
				while (ua != NULLID) {
					arc<SCALAR_TYPE>& uap = this->arcs[ua];

					// skip the arc itself
					if (ua == a) {
						ua = uap.upper_next; // we're guarantee this is the next
						continue;
					}

					//test if they are the same kind
					if (ap->upper != uap.upper) {
						ua = uap.lower_next;
						continue;
					}
					if (!isAlive(uap, num_cancelled)) {
						ua = uap.upper_next;
						continue;
					}

					// create the arc here!
					INT_TYPE newarc = createArc(la, a, ua, num_cancelled + 1);
					//printf("ha1\n");
					ap = &(this->arcs[a]);
					//printf("ha\n");
					createcounter++;

					ua = this->arcs[ua].upper_next;
					//printf("asdf\n");
				}

				la = this->arcs[la].lower_next;

			}

			//printf("hahan\n");
			 // go through and set the delete time on the arcs that are to be removed,
			 // and update the arc counters at the other endpoints

			 // following comments are for 0-1 cancellations
			 // first look in neighborhood of minimum = lower
			la = lower.firstarc;							// pick first arc attached to min
			while (la != NULLID) {							// while there are more arcs in neighborhood, keep looking
				arc<SCALAR_TYPE>& lap = this->arcs[la];			// get the reference to the arc attached to the min

				// skip the arc itself
				if (la == a) {
					la = lap.lower_next; // we're guarantee this is the next
					continue;
				}

				//test if they are the same kind
				if (ap->lower == lap.lower) {				// we do not want arcs that are not 0-1, so do this check
					if (isAlive(lap, num_cancelled)) {		// make sure we are only working with living arcs
						node<SCALAR_TYPE>& n = this->nodes[lap.upper]; // now "n" is the 1-saddle attached to minimum

						n.numarcs--;						// if we remove this arc, we need to decrement the number of living arcs attached to it
						n.numlower--;						// this cancellation will also remove a "lower" arc of the 1-saddle

						// create merged_manifold for merging
						// extend the descending merged_manifold of the 1-saddle by merging it with the descending merged_manifold of the cancelled 1-saddle
						INT_TYPE nmanid =
							createManifold(lap.upper, n.dmanifoldid, upper.dmanifoldid, num_cancelled + 1);
						n.dmanifoldid = nmanid;				// change the merged_manifold reference of the 1-saddle to reference this new decending merged_manifold

						lap.destroyed = num_cancelled + 1;		// we remove the old arc to the deleted min, so mark it destroyed with the num_cancelld
						deletecounter++;					// keep track of deleted arc cound
					}
					la = lap.lower_next;					// continue on the next arc around the minimum
				}
				else if (ap->lower == lap.upper) {          // if this is in fact an i-1 - i arc, we just remove the arcs around it - no merging happens
					if (isAlive(lap, num_cancelled)) {
						node<SCALAR_TYPE>& n = this->nodes[lap.lower];
						n.numarcs--;
						lap.destroyed = num_cancelled + 1;
						deletecounter++;
					}
					la = lap.upper_next;
				}
				else {
					printf("ERROR SHOULD NEVER GET HERE\n");
				}



			}
			INT_TYPE ua = upper.firstarc;
			while (ua != NULLID) {
				arc<SCALAR_TYPE>& uap = this->arcs[ua];

				// skip the arc itself
				if (ua == a) {
					ua = uap.upper_next; // we're guarantee this is the next
					continue;
				}
				//test if they are the same kind
				if (ap->upper == uap.upper) {
					if (isAlive(uap, num_cancelled)) {

						node<SCALAR_TYPE>& n = this->nodes[uap.lower]; // pick the other

						n.numarcs--;

						// create merged_manifold for merging
						INT_TYPE nmanid = createManifold(uap.lower, n.amanifoldid, lower.amanifoldid, num_cancelled + 1);
						n.amanifoldid = nmanid;

						uap.destroyed = num_cancelled + 1;
						deletecounter++;
					}
					ua = uap.upper_next;
				}
				else if (ap->upper == uap.lower) {
					if (isAlive(uap, num_cancelled)) {
						node<SCALAR_TYPE>& n = this->nodes[uap.upper];
						n.numarcs--;
						n.numlower--;
						uap.destroyed = num_cancelled + 1;
						deletecounter++;
					}
					ua = uap.lower_next;
				}
				else {
					printf("ERROR SHOULD NEVER GET HERE\n");
				}


			}
			num_cancelled++;


			ap->destroyed = num_cancelled;
			lower.destroyed = num_cancelled;
			upper.destroyed = num_cancelled;
			deletecounter++;

			//int tmp2 = (int) edges_to_cancel.size();

			//printf("initialguess: %d actual:%d del:%d td:%d tmp1:%d tmp2:%d \n", initialguess, createcounter, deletecounter,init2
			//	,tmp1, tmp2);
			//printf("gothere\n");
			return 1;
		}



		INT_TYPE nextArc(const arc<SCALAR_TYPE>& ap, INT_TYPE n) const {
			if (ap.lower == n) return ap.lower_next;
			if (ap.upper == n) return ap.upper_next;
			printf("ERROR BARF POOP\n");
			return NULLID;
		}
		int countMultiplicity(arc<SCALAR_TYPE>& ap, INT_TYPE ctime) const {
			INT_TYPE nu = ap.upper;
			INT_TYPE nl = ap.lower;
			INT_TYPE a = this->nodes[nu].firstarc;
			int counter = 0;
			while (a != NULLID) {
				const arc<SCALAR_TYPE>& nap = this->arcs[a];
				if (isAlive(nap, ctime) && nap.lower == nl && nap.upper == nu) {
					counter++;
				}
				a = nextArc(nap, nu);
			}
			return counter;
		}

		virtual bool isValid(INT_TYPE a, arc<SCALAR_TYPE>& ap) const {
			// test for boundary
			if (this->nodes[ap.lower].boundary !=
				this->nodes[ap.upper].boundary) return false;

			// 2 endpoints must be connected by exactly one arc
			if (countMultiplicity(ap, num_cancelled) != 1) return false;
			// test for inversions?

			return true;


		}

		// THIS WILL HAVE TO WORK LATER, but for now only need 0 and maxdim manifolds
		//void recCounLeaftManifolds(INT_TYPE mId, map< INT_TYPE, int >& counter) {
		//	merged_manifold &m = this->getManifold(mId);
		//	if (m->merge[0] != NULL) {
		//		recCounLeaftManifolds(m->merge[0], counter);
		//		recCounLeaftManifolds(m->merge[1], counter);
		//	}
		//	else {
		//		if (counter.count(m) == 0) {
		//			counter[m] = 1;
		//		}
		//		else {
		//			counter[m]++;
		//		}
		//	}
		//}

		void printmanifold(INT_TYPE man) const {
			merged_manifold& m = getManifold(man);
			printf("man=%d, man.base=%d, man.mergetime=%d, man.merge[0]=%d, man.merge[1]=%d\n",
				man, m.basenode, m.mergetime, m.merged[0], m.merged[1]);
		}

		INT_TYPE getActiveMan(INT_TYPE man) const {
			//printmanifold(man);
			while (getManifold(man).mergetime > this->current_pers()) {
				man = getManifold(man).merged[0];
				//printf("  --");  printmanifold(man);
			}
			//printf("ret %d\n", man);
			return man;
		}



	protected:
		void recGatherNodes(INT_TYPE mid, set<INT_TYPE>& res)  const {
			vector<INT_TYPE> expand;
			expand.push_back(mid);
			while (expand.size() > 0) {

				INT_TYPE curr = expand.back(); expand.pop_back();

				const merged_manifold& m = getManifold(curr);
				if (m.merged[0] != -1) {
					expand.push_back(m.merged[0]);
					expand.push_back(m.merged[1]);
				}
				res.insert(m.basenode);


			}


		}

		public:

		void rec_man_trace_up(INDEX_TYPE cellid, set<INDEX_TYPE>& res) const {
			res.insert(cellid);
			INDEX_TYPE current = cellid;
			DIM_TYPE cdim = mMesh->dimension(cellid);
			typename MESH_TYPE::CofacetsIterator cofacets(mMesh);
			for (cofacets.begin(current); cofacets.valid(); cofacets.advance()) {
				INDEX_TYPE temp_id = cofacets.value();
				if (mGrad->getCritical(temp_id)) continue;

				INDEX_TYPE temp_pair = mGrad->getPair(temp_id);

				if (temp_pair == cellid) continue;

				if (mMesh->dimension(temp_pair) != cdim) continue;

				rec_man_trace_up(temp_pair, res);
			}
		}
		void rec_man_trace_down(INDEX_TYPE cellid, set<INDEX_TYPE>& res) const {
			res.insert(cellid);
			INDEX_TYPE current = cellid;
			DIM_TYPE cdim = mMesh->dimension(cellid);
			typename MESH_TYPE::FacetsIterator facets(mMesh);
			for (facets.begin(current); facets.valid(); facets.advance()) {
				INDEX_TYPE temp_id = facets.value();
				if (mGrad->getCritical(temp_id)) continue;

				INDEX_TYPE temp_pair = mGrad->getPair(temp_id);

				if (temp_pair == cellid) continue;

				if (mMesh->dimension(temp_pair) != cdim) continue;

				rec_man_trace_down(temp_pair, res);
			}
		}

		protected:

		void fillUnsimplifiedGeometry(INT_TYPE nodeID, set<INDEX_TYPE>& res, bool ascending) const {

			node<SCALAR_TYPE>& n = this->getNode(nodeID);

			if (ascending) {
				rec_man_trace_up(n.cellindex, res);
			}
			else {
				rec_man_trace_down(n.cellindex, res);
			}

		}

	public:


		void GatherNodes(INT_TYPE nodeID, set<INT_TYPE>& res, bool ascending) const {

			if (!isAlive(nodeID, this->select_persistence)) {
				return;
			}
			const node<SCALAR_TYPE>& n = this->getNode(nodeID);
			INT_TYPE man;
			if (ascending) {
				man = n.amanifoldid;
			}
			else {
				man = n.dmanifoldid;
			}
			INT_TYPE manID = getActiveMan(man);
			recGatherNodes(manID, res);

		}

		void fillGeometry(INT_TYPE nodeID, set<INDEX_TYPE>& res, bool ascending) const {

			if (!isNodeAlive(nodeID)) return;

			set<INT_TYPE> nodeset;
			//printf("gothere1\n");
			GatherNodes(nodeID, nodeset, ascending);
			//printf("gothere2\n");

			for (set<INT_TYPE>::iterator it = nodeset.begin(); it != nodeset.end(); it++) {
				const node<SCALAR_TYPE>& n = this->getNode(*it);

				if (ascending) {
					rec_man_trace_up(n.cellindex, res);
				}
				else {
					rec_man_trace_down(n.cellindex, res);
				}

			}
			//printf("gothere3\n");



		}


		void fillBaseGeometry(INT_TYPE nodeID, set<INDEX_TYPE>& res, bool ascending) const {
			const node<SCALAR_TYPE>& n = this->getNode(nodeID);
			if (ascending) {
				rec_man_trace_up(n.cellindex, res);
			}
			else {
				rec_man_trace_down(n.cellindex, res);
			}
		}

	protected:

		bool get_next_to_cancel(INT_TYPE& a) {
			///printf("getnext to cancel called\n");
			while (!edges_to_cancel.empty()) {
				sortedEdge se = edges_to_cancel.top();
				edges_to_cancel.pop();
				a = se.ep;
				arc<SCALAR_TYPE>& ap = this->arcs[a];
				///printf("%u ", a);
				// is it alive in the current context
				if (!isAlive(ap, num_cancelled)) {
					//printf("adsf1\n");
					//printf("skip1->");
					///printArc(a);
					continue;
				}

				//test if it's a valid cancellation
				if (!isValid(a, ap)) {
					//printf("adsf2\n");
					continue;
				}

				int newcountweight = edgeCountWeight(ap);
				if (newcountweight > se.countweight) {
					//printf("adsf3\n");
					se.countweight = newcountweight;
					edges_to_cancel.push(se);
					continue;
				}

				if (newcountweight > 1500) {
					se.persistence += 1;
					if (se.persistence <= gPersThreshold)
						edges_to_cancel.push(se);
					continue;
				}

				//printf("getnext to cancel returned true\n");

				return true;


			}
			//printf("getnext to cancel returned false\n");
			return false;
		}

		vector<SCALAR_TYPE> cancel_num_to_pers;
		INT_TYPE select_persistence; // the persistence to select in hierarchy
		INT_TYPE num_cancelled;

		SCALAR_TYPE max_pers_so_far;
		SCALAR_TYPE gPersThreshold;

		int countIDs[4];
		int countBoundaryIds[4];
	public:

    void SetPersistanceByNumOfMaxima(int num_maxima){
      int maxcount = this->countIDs[3];
      
      for (int i = 0; i < this->mCRecords.size(); i++) {
        cancellation_record& cr = this->mCRecords[i];
      
        if(maxcount < num_maxima){
          select_persistence = i-1;
          return;
        }
        
        if (cr.index == 2) maxcount--;

      }
    }
    
		virtual void ComputeHierarchy(SCALAR_TYPE pers_limit) {
			cancel_num_to_pers.clear();
			printf(" -- Performing cancellation to %f...\n", pers_limit);
			gPersThreshold = pers_limit;
			max_pers_so_far = 0;
			num_cancelled = 0;
			// insert every arc to cancel list
			printf("  -- Adding arcs to sorter...");
			INT_TYPE mysize = (INT_TYPE) this->arcs.size();
			for (INT_TYPE i = 0; i < mysize; i++) {
				// test if it passes the hierarchy test
				//if (!hierarchy_test->slice_position_list(i, this)) continue;
				arc<SCALAR_TYPE> &a = this->arcs[i];

				if (a.persistence > gPersThreshold) continue;

				sortedEdge se;

				se.persistence = a.persistence;
				se.countweight = edgeCountWeight(a);
				se.ep = i;
				edges_to_cancel.push(se);
			}

			printf("Done!\n");
			printf("  -- Cancelling:");

			INT_TYPE a;
			//float maxv = 0;

			while (get_next_to_cancel(a) && this->arcs[a].persistence <= gPersThreshold) {
				//if (arcs[a].persistence > maxv) maxv = arcs[a].persistence;
				if (num_cancelled % 1000 == 0) {
					printf("\r  -- Cancelling: %u val=%f", num_cancelled, (float)max_pers_so_far);
				}
				//printf("\n\ncancelling: %d", num_cancelled);
				//printArc(a);

				cancel(a);
				cancel_num_to_pers.push_back(max_pers_so_far);

				//printf("Done\n");
			}
			printf("\r  -- Cancelling finished: %u val=%f\n", num_cancelled, (float)max_pers_so_far);


			//_evercomputed = true;
			if (output_cancellation_records) {
				int mincount = this->countIDs[0];
				int sad1count = this->countIDs[1];
				int sad2count = this->countIDs[2];
				int maxcount = this->countIDs[3];

				int bcount[4];
				for (int i = 0; i < 4; i++) bcount[i] = this->countBoundaryIds[i];

				FILE* fout = fopen(this->mCRecordFName, "w");
				//fprintf(fout, "pers_so_far actual_pers lval uval perc_pers lindex #min #1s #2s #max is_bdry\n", max_pers_so_far, cr.persistence,

				float max_pers_so_far = this->mCRecords[0].persistence;
				for (int i = 0; i < this->mCRecords.size(); i++) {
					cancellation_record& cr = this->mCRecords[i];
					max_pers_so_far = fmax(max_pers_so_far, cr.persistence);
					fprintf(fout, "%f %f %f %f %f %d %d %d %d %d %d %d %d %d %d\n", max_pers_so_far, cr.persistence, 
						cr.lval, cr.uval, cr.persPerc, cr.index, mincount, sad1count, sad2count, maxcount,
						 cr.boundary, bcount[0], bcount[1], bcount[2], bcount[3]);
					if (cr.index == 0) {
						mincount--; sad1count--;
					}
					else if (cr.index == 1) {
						sad1count--; sad2count--;
					}
					else if (cr.index == 2) {
						maxcount--; sad2count--;
					}
					if (cr.boundary) {
						bcount[cr.index]--; 
						bcount[cr.index + 1]--;
					}
				}
				fclose(fout);
			}
			// now fill in the lists?
		}
		void Output1SaddleRecord(const char* fname) const  {

			FILE* fsadrec = fopen(fname, "w");

			for (int i = 0; i < nodes.size(); i++) {
				node<SCALAR_TYPE>& n = nodes[i];
				if (n.dim > 1) continue;
				float persval =
					(n.destroyed < cancel_num_to_pers.size() + 1 ?
						cancel_num_to_pers[n.destroyed - 1] :
						cancel_num_to_pers[cancel_num_to_pers.size() - 1]);
				fprintf(fsadrec, "%d %d %llu %f %d %f\n", n.dim, n.boundary, n.cellindex, n.value, n.destroyed, persval);


			}
			fclose(fsadrec);


		}
	protected:
		void _fillArcGeometry(INT_TYPE aid, vector<INDEX_TYPE>& v, bool direction) const  {
			const arc<SCALAR_TYPE>& a = this->arcs[aid];
			//printf("gothere1\n");
			if (a.created == 0) {
				// this is base
				const arc_base_geometry& base = this->arc_base_geoms[a.geom];
				if (direction) {
					for (int i = 0; i < base.geometry.size(); i++) {
						if (v.size() > 0 && v[v.size() - 1] == base.geometry[i]) {
							if (v.size() > 1 && i < base.geometry.size() -1 && v[v.size()-2] == base.geometry[i+1]) v.pop_back();
						}
						else {
							v.push_back(base.geometry[i]);
						}
					}
				}
				else {
					for (int i = base.geometry.size() - 1; i >= 0; i--) {
						if (v.size() > 0 && v[v.size() - 1] == base.geometry[i]) {
							if (v.size() > 1 && i >0 && v[v.size() - 2] == base.geometry[i - 1])v.pop_back();
						}
						else {
							v.push_back(base.geometry[i]);
						}
					}
				}
				return;
			}
			//printf("gothere2\n");
			// recurse on children
			const arc_merged_geometry& m = arc_merge_geoms[a.geom];
			if (direction) {
				_fillArcGeometry(m.fields[0], v, true);
				_fillArcGeometry(m.fields[1], v, false);
				_fillArcGeometry(m.fields[2], v, true);
			}
			else {
				_fillArcGeometry(m.fields[2], v, false);
				_fillArcGeometry(m.fields[1], v, true);
				_fillArcGeometry(m.fields[0], v, false);
			}


		}
	public:
		INT_TYPE arcLowerNode(INT_TYPE aid) const {
			return getArc(aid).lower;
		}
		INT_TYPE arcUpperNode(INT_TYPE aid) const {
			return getArc(aid).upper;
		}
		INT_TYPE nextIncidentLivingArc(INT_TYPE aid, INT_TYPE nid) const {
			//printf("%d-\n", aid);
			INT_TYPE naid = nextArc(getArc(aid), nid);
			//printf("  -%d\n", naid);
			while (naid != NULLID) {
				if (isArcAlive(naid)) return naid;
				naid = nextArc(getArc(naid), nid);
				//printf("  -%d\n", naid);
			}
			return naid;
		}
		bool isValidArcId(INT_TYPE aid) const {
			return aid != NULLID;
		}

		// return count, fill array of saddle extrema arcs
		int GetLiving2SaddleExtremaArcs(INT_TYPE sadid, INT_TYPE* extrema_arcs) {
			auto first_arc = firstIncidentLivingArc(sadid);
			while (first_arc != NULLID && arcLowerNode(first_arc) != sadid) {
				first_arc = nextIncidentLivingArc(first_arc, sadid);
			}
			if (first_arc == NULLID) return 0;
			extrema_arcs[0] = first_arc;
			auto second_arc = nextIncidentLivingArc(first_arc, sadid);
			while (second_arc != NULLID && arcLowerNode(second_arc) != sadid) {
				second_arc = nextIncidentLivingArc(second_arc, sadid);
			}
			if (second_arc == NULLID) return 1;
			extrema_arcs[1] = second_arc;
			return 2;
			
		}
		// return count, fill array of saddle extrema arcs
		int GetLiving1SaddleExtremaArcs(INT_TYPE sadid, INT_TYPE* extrema_arcs) {
			auto first_arc = firstIncidentLivingArc(sadid);
			while (first_arc != NULLID && arcUpperNode(first_arc) != sadid) {
				first_arc = nextIncidentLivingArc(first_arc, sadid);
			}
			if (first_arc == NULLID) return 0;
			extrema_arcs[0] = first_arc;
			auto second_arc = nextIncidentLivingArc(first_arc, sadid);
			while (second_arc != NULLID && arcUpperNode(second_arc) != sadid) {
				second_arc = nextIncidentLivingArc(second_arc, sadid);
			}
			if (second_arc == NULLID) return 1;
			extrema_arcs[1] = second_arc;
			return 2;

		}
		INT_TYPE firstIncidentLivingArc(INT_TYPE nid) const {
			const node<SCALAR_TYPE>& n = getNode(nid);
			INT_TYPE aid = n.firstarc;
			if (!isArcAlive(aid)) return nextIncidentLivingArc(aid, nid);
			return aid;
		}

		void fillArcGeometry(INT_TYPE aid, vector<INDEX_TYPE>& v) const {
			//v.push_back(nodes[arcs[aid].upper].cellindex);
			v.clear();
			_fillArcGeometry(aid, v, true);
			//v.push_back(nodes[arcs[aid].lower].cellindex);
		}

		class SurroundingArcsIterator {
		protected:
			MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>* mMSC;
			INT_TYPE mNID;
			INT_TYPE currarc;
			INT_TYPE next_arc(INT_TYPE arcid) {
				return mMSC->nextArc(mMSC->getArc(arcid), mNID);
			}

		public:
			SurroundingArcsIterator(MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>* msc) :
				mMSC(msc) {}

			void begin(INT_TYPE nid) {
				//printf("begin %d...\n");
				mNID = nid;
				node<SCALAR_TYPE>& n = mMSC->getNode(mNID);
				currarc = n.firstarc;
			}
			bool valid() {
				//printf("valid? %d %d\n", currarc, NULLID);
				return currarc != NULLID;
			}
			INT_TYPE value() {
				return currarc;
			}
			void advance() {
				//printf("Advance? %d\n", currarc);
				if (currarc == NULLID) return;
				//printf("...\n");
				currarc = next_arc(currarc);
				//printf("advance = %d\n", currarc);
			}
		};

		class SurroundingLivingArcsIterator : public SurroundingArcsIterator {
		protected:
			bool advance_until_alive() {
				this->currarc = this->next_arc(this->currarc);

				if (this->currarc == NULLID) return false;
				while (!this->mMSC->isArcAlive(this->currarc)) {
					this->currarc = this->next_arc(this->currarc);
					if (this->currarc == NULLID) return false;
				}
				return true;
			}
		public:
			SurroundingLivingArcsIterator(MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>* msc) :
				SurroundingArcsIterator(msc) {}

			void begin(INT_TYPE nid) {
				//printf("begin %d...\n");
				SurroundingArcsIterator::begin(nid);
				if (!this->valid())
					printf("ERROR: node has no arcs!\n");
				if (!this->mMSC->isArcAlive(this->currarc)) advance_until_alive();
				//printf("begin return %d\n", currarc);
			}

			void advance() {
				//printf("Advance? %d\n", currarc);
				if (this->currarc == NULLID) return;
				//printf("...\n");
				advance_until_alive();
				//printf("advance = %d\n", currarc);
			}
		};

		class NodesIterator {
		protected:
			MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>* mMSC;
			INT_TYPE currid;

		public:
			NodesIterator(MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>* msc) :
				mMSC(msc) {}

			void begin() {
				currid = 0;
			}
			bool valid() {
				return currid < mMSC->nodes.size();
			}
			INT_TYPE value() {
				return currid;
			}
			void advance() {
				currid++;
			}
		};
		class LivingNodesIterator : public NodesIterator {
		protected:
 
			void advance_until_alive() {
				this->currid++;
				while (this->valid() && !this->mMSC->isNodeAlive(this->currid)) {
					this->currid++;
				}
			}
		public:
			LivingNodesIterator(MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>* msc) :
				NodesIterator(msc) {}

			void begin() {
				NodesIterator::begin();
				if (this->valid() && !this->mMSC->isNodeAlive(this->currid)) advance_until_alive();
			}
			void advance() {
				advance_until_alive();
			}
		};

		// replicates functionality just uses different end criteria
		class ArcsIterator {
		protected:
			MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>* mMSC;
			INT_TYPE currid;

		public:
			ArcsIterator(MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>* msc) :
				mMSC(msc) {}

			void begin() {
				currid = 0;
			}
			bool valid() {
				return currid < mMSC->arcs.size();
			}
			INT_TYPE value() {
				return currid;
			}
			void advance() {
				currid++;
			}
		};
		class LivingArcsIterator : public ArcsIterator {
		protected:
			void advance_until_alive() {
				this->currid++;
				while (this->valid() && !this->mMSC->isArcAlive(this->currid)) {
					this->currid++;
				}
			}
		public:
			LivingArcsIterator(MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>* msc) :
				ArcsIterator(msc) {}

			void begin() {
				ArcsIterator::begin();
				if ((this->valid() && !this->mMSC->isArcAlive(this->currid))) advance_until_alive();
			}
			void advance() {
				advance_until_alive();
			}
		};


		void print_complex_info(bool verbose) {
			LivingNodesIterator nit(this);
			int numnodes = 0;
			int dims[4]; for (int i = 0; i < 4; i++) dims[i] = 0;
			for (nit.begin(); nit.valid(); nit.advance()) {
				numnodes++;
				INT_TYPE nid = nit.value();
				node<SCALAR_TYPE>& n = this->getNode(nid);
				if (verbose) {
					printf("node: %lld, %d, %f\n", n.cellindex, n.dim, n.value);
				}
				dims[n.dim]++;
			}
			int numarcs = 0;
			int dims2[4]; for (int i = 0; i < 4; i++) dims2[i] = 0;
			for (int i = 0; i < 4; i++) dims2[i] = 0;
			LivingArcsIterator ait(this);
			for (ait.begin(); ait.valid(); ait.advance()) {
				INT_TYPE aid = ait.value();
				numarcs++;
				arc<SCALAR_TYPE>& a = this->getArc(aid);
				dims2[a.dim]++;
				if (verbose) {
					printf("arc: %lld--%lld %d\n", this->getNode(a.lower).cellindex, this->getNode(a.upper).cellindex, a.dim);
				}
			}
			printf("msc: %d nodes [%d, %d, %d, %d] ", numnodes, dims[0], dims[1], dims[2], dims[3]);

			printf("%d arcs [%d %d %d]\n", numarcs, dims2[0], dims2[1], dims2[2]);
		};


		struct write_node {
			INDEX_TYPE cellindex;
			INT_TYPE amanifoldid; // set to -1 if this is the base
			INT_TYPE dmanifoldid; // set to -1 if this is the base
			SCALAR_TYPE value; // KEEP this?

		};
		struct write_arc {
			INT_TYPE lower; // node
			INT_TYPE upper; // node
			INT_TYPE geom_size; // if created == 0, this is an original arc, and its geometry will be found in base_geom list of MSC, else in the merged_geom list
		};

		struct write_state {
			vector<write_node> nodes;
			vector<write_arc> arcs;
			vector<INDEX_TYPE> geom;

			void write_to_file(const char* name) {
			
				FILE* fout = fopen(name, "wb");

				int numnodes = nodes.size();
				int numarcs = arcs.size();
				int numgeom = geom.size();
				fwrite(&numnodes, sizeof(int), 1, fout);
				fwrite(&numarcs, sizeof(int), 1, fout);
				fwrite(&numgeom, sizeof(int), 1, fout);
				//printf("%d %d %d\n", numnodes, numarcs, numgeom);

				fwrite(&nodes[0], sizeof(write_node), numnodes, fout);
				fwrite(&arcs[0], sizeof(write_arc), numarcs, fout);
				fwrite(&geom[0], sizeof(INDEX_TYPE), numgeom, fout);
				fclose(fout);

			}

			void read_from_file(const char* name) {

				FILE* fin = fopen(name, "rb");

				int numnodes;
				fread(&numnodes, sizeof(int), 1, fin);
				nodes.resize(numnodes);
				int numarcs;
				fread(&numarcs, sizeof(int), 1, fin);
				arcs.resize(numarcs);
				int numgeom;
				fread(&numgeom, sizeof(int), 1, fin);
				geom.resize(numgeom);
				//printf("%d %d %d\n", numnodes, numarcs, numgeom);

				fread(&nodes[0], sizeof(write_node), numnodes, fin);
				//for (auto wn : nodes) printf("%lld\n", wn.cellindex);
				fread(&arcs[0], sizeof(write_arc), numarcs, fin);
				fread(&geom[0], sizeof(INDEX_TYPE), numgeom, fin);
				fclose(fin);

			}
		};

		void WriteComplex1Skeleton(const char* filename) {

			map<INT_TYPE, INT_TYPE> living_nodes;
		
			write_state ws;

			LivingNodesIterator nit(this);
			int pos = 0;
			for (nit.begin(); nit.valid(); nit.advance()) {
				INT_TYPE nid = nit.value();
				living_nodes[nid] = pos;
				write_node wn;
				node<SCALAR_TYPE>& n = this->getNode(nid);
				wn.cellindex = n.cellindex;
				wn.value = n.value;
				ws.nodes.push_back(wn);
				pos++;
			}

			LivingArcsIterator ait(this);
			for (ait.begin(); ait.valid(); ait.advance()) {
				INT_TYPE aid = ait.value();
				write_arc wa;
				arc<SCALAR_TYPE>& a = this->getArc(aid);
				wa.lower = living_nodes[a.lower];
				wa.upper = living_nodes[a.upper];
				vector<INDEX_TYPE> geom;
				this->fillArcGeometry(aid, geom);
				wa.geom_size = geom.size();
				ws.geom.insert(ws.geom.end(), geom.begin(), geom.end());
				ws.arcs.push_back(wa);
			}

			ws.write_to_file(filename);

		}

		void LoadComplex1Skeleton(const char* filename) {
			
			write_state ws;
			ws.read_from_file(filename);
			printf("read %d nodes, %d arcs, %d geom\n", ws.nodes.size(), ws.arcs.size(), ws.geom.size());
			for (auto wn : ws.nodes) {
				this->createNode(wn.cellindex, wn.value);
			}
			auto geom_it = ws.geom.begin();
			for (auto wa : ws.arcs) {
				int size = wa.geom_size;
				vector<INDEX_TYPE> geom(wa.geom_size);
				geom.insert(geom.begin(), geom_it, geom_it + wa.geom_size);
				geom_it += wa.geom_size;
				this->createArc(ws.nodes[wa.lower].cellindex, ws.nodes[wa.upper].cellindex, geom);
			}
		}


	};

	  

	namespace LightMSC {

		struct msbitfield
		{
			//// rehash this
			unsigned char dim : 3;
			unsigned char boundary : 1;
			unsigned char f1 : 1;
			unsigned char f2 : 1;
			unsigned char f3 : 2;
		};

		template<typename DTYPE>
		struct node
		{
			INDEX_TYPE cellindex;
			INT_TYPE firstarc;
			INT_TYPE destroyed; // the time this node is cancelled.
			unsigned short numarcs;
			unsigned short numlower;
			DTYPE value;
			msbitfield flags;
		};


		template<typename DTYPE>
		struct arc
		{
			INT_TYPE lower; // node
			INT_TYPE lower_next; //arc - INT_INFTY is the null arc!
			INT_TYPE upper; // node
			INT_TYPE upper_next; //arc
			INT_TYPE created;
			INT_TYPE destroyed;
			int multiplicity;
			DTYPE persistence;
			msbitfield flags; // the .dim stores the dim of the lower endpoint
		};

		template<typename DTYPE, class MESH_TYPE, class FUNC_TYPE, class GRAD_TYPE>
		class SuperLightMSC {
		protected:
			struct sortedEdge {
				DTYPE persistence;
				int countweight;
				INT_TYPE ep;

				bool operator()(const sortedEdge& _Left, const sortedEdge& _Right) {
					if (_Left.persistence < _Right.persistence) return false;
					if (_Left.persistence > _Right.persistence) return true;
					if (_Left.countweight < _Right.countweight) return false;
					if (_Left.countweight > _Right.countweight) return true;
					return _Left.ep > _Right.ep;
				}
			};

			priority_queue<sortedEdge, vector< sortedEdge>, sortedEdge > edges_to_cancel;

			struct cancellation_record {
				int index;
				DTYPE persistence;
				DTYPE lval;
				DTYPE uval;
				DTYPE persPerc;
				INT_TYPE arcp;
				int boundary;
			};

			vector<cancellation_record> mCRecords;
			bool output_cancellation_records;

			const char* mCRecordFName;
			float m_temp_perc_pers;
			//INDEX_TYPE select_persistence;
			INDEX_TYPE num_destroyed;
			GRAD_TYPE* mGrad;
			MESH_TYPE* mMesh;
			FUNC_TYPE* mFunc;

			vector<node<DTYPE> > nodes;
			vector<arc<DTYPE> > arcs;
#ifdef WIN32
			unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#elif defined(__APPLE__)
			unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#else
			unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
			//__gnu_cxx::hash_map<INDEX_TYPE, INT_TYPE> nodemap;
#endif
			inline void connectArc(INT_TYPE arcID, arc<DTYPE>& a, node<DTYPE>& lower, node<DTYPE>& upper) {
				a.lower_next = lower.firstarc;
				lower.firstarc = arcID;
				a.upper_next = upper.firstarc;
				upper.firstarc = arcID;

				lower.numarcs++;
				upper.numarcs++;
				upper.numlower++;

			}

			void connectArc(INT_TYPE arcID) {

				arc<DTYPE>& a = arcs[arcID];
				node<DTYPE>& lower = nodes[a.lower];
				node<DTYPE>& upper = nodes[a.upper];
				connectArc(arcID, a, lower, upper);
			}


		public:

			INT_TYPE numArcs() {
				return (INT_TYPE)arcs.size();
			}
			INT_TYPE numNodes() {
				return (INT_TYPE)nodes.size();
			}
			node<DTYPE>& getNode(INT_TYPE e) {
				return nodes[e];
			}

			arc<DTYPE>& getArc(INT_TYPE e) {
				return arcs[e];
			}


			void set_output_cancellation_records(const char* fname) {
				output_cancellation_records = true;
				mCRecordFName = fname;
			}


			SuperLightMSC(GRAD_TYPE* grad,
				MESH_TYPE* mesh,
				FUNC_TYPE* func) :
				mGrad(grad), mMesh(mesh), mFunc(func)
			{
				num_destroyed = 0;
				select_persistence = 0;
				output_cancellation_records = false;
			}


			//INT_TYPE createNode(node<DTYPE> &n) {
			//	INT_TYPE listID = (INT_TYPE)nodes.push_back(n);
			//	nodemap[n.cellindex] = listID;
			//	//this->countIDs[n.flags.dim]++;
			//	//if (!((bool)n.flags.boundary)) this->countInteriorIDs[n.flags.dim]++;
			//	return listID;
			//}

			INT_TYPE createNode(INDEX_TYPE cellID) {
				node<DTYPE> tn;
				INT_TYPE listID = nodes.size();
				//if (nodes.size() == nodes.capacity()) printf("about to expand capacity - nodes createnode!\b");
				nodes.push_back(tn);
				nodemap[cellID] = listID;

				//node<DTYPE> &n = nodes->operator [](listID);
				node<DTYPE> &n = nodes[listID];
				n.cellindex = cellID;

				n.destroyed = INT_INFTY;
				n.firstarc = NULLID;
				n.flags.boundary = mMesh->boundaryValue(cellID);
				//if (n.flags.boundary) {
				//	INDEX_TYPE coords[3];
				//	sgg->getCoords(coords, cellID);
				//	printf("(%d, %d, %d) = %d\n", 
				//	   (int) coords[0], (int) coords[1], (int) coords[2], n.flags.boundary);
				//}
				n.flags.dim = mMesh->dimension(cellID);
				this->countIDs[n.flags.dim]++;
				n.numarcs = 0;
				n.numlower = 0;
				n.value = mFunc->cellValue(cellID);

				return listID;
			}

			//// create an arc connecting lower to upper <- cell id's from gradient
			//INT_TYPE createArc(arc<DTYPE> &a){
			//	INT_TYPE listID = (INT_TYPE)arcs.push_back(a);
			//	//printf("%f - %f\n", (float) a.persistence, arcs[listID].persistence); 
			//	return listID;
			//}

			INT_TYPE createArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID, int count) {

				arc<DTYPE> at;
				INT_TYPE listID = arcs.size();
				//if (arcs.size() == arcs.capacity()) printf("about to expand capacity - arcs createnode!\b");

				arcs.push_back(at);

				arc<DTYPE>& a = arcs[listID];

				//if (count != -1) printf("creating arc with count %d\n", count);
				a.multiplicity = count;
				a.created = 0;
				a.destroyed = INT_INFTY;
				a.lower = nodemap[lowerCellID];
				a.upper = nodemap[upperCellID];
				//node<DTYPE>& ln = nodes[a.lower];
				//node<DTYPE>& un = nodes[a.upper];
				//sgg->printCell(lowerCellID); printf("->"); sgg->printCell(upperCellID); printf("\n");
				//if (ln.flags.dim != un.flags.dim - 1) printf("ERROR IN CREATE ARC\n");
				if (nodes[a.upper].value < nodes[a.lower].value)
					printf("ERROR: upper (%f, %d) - lower (%f, %d)\n",
						float(nodes[a.upper].value), nodes[a.upper].flags.dim,
						float(nodes[a.lower].value), nodes[a.lower].flags.dim);

				a.persistence = nodes[a.upper].value - nodes[a.lower].value;
				//printf("n=%f\n", (float) a.persistence);
				a.flags.dim = nodes[a.lower].flags.dim;
				// now connect the crap!
				connectArc(listID);

				//printf("creating arc: \n\t");
				//sgg->printCell(lowerCellID);
				//for (int i = 0; i < g->size(); i++) {
				//	printf("\n");
				//	sgg->printCell(g->operator [](i));
				//}


				//printf("\n\t");
				//sgg->printCell(upperCellID);
				//printf("\n\t%d\n", g->size());
				/*		for (int i=0; i < g->size(); i++) printf("%d ", g->operator [](i));
				printf("\n")*/;

				//printf("Creating arc: %d-%d, %d-%d\n", lowerCellID, upperCellID, nodes[a.lower].flags.dim, nodes[a.upper].flags.dim);

				return listID;
			}

			INT_TYPE createArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID) {
				return createArc(lowerCellID, upperCellID, -1);
			}

			INT_TYPE createArc(INT_TYPE luap, INT_TYPE ulap, INT_TYPE ctime) {
				arc<DTYPE>& lua = arcs[luap];
				arc<DTYPE>& ula = arcs[ulap];
				INT_TYPE na_id = arcs.size();
				//this->arcs.push_back(ma); // copy setting from old

				arc<DTYPE> na;// = this->arcs[na_id];
				na.created = ctime;
				na.destroyed = INT_INFTY;
				na.lower = ula.lower;
				na.upper = lua.upper;
				na.multiplicity = -1;

				node<DTYPE>& nup = this->nodes[na.upper];
				node<DTYPE>& nlo = this->nodes[na.lower];

				this->connectArc(na_id, na, nlo, nup);
				na.persistence = nup.value - nlo.value;
				na.flags.boundary = nlo.flags.boundary || nup.flags.boundary;
				na.flags.dim = nlo.flags.dim;
				if (na.persistence < 0 && nlo.flags.boundary == nup.flags.boundary) printf("creatinga inversion %f, %d-%d, b%d-b%d\n",
					(float)na.persistence, nlo.flags.dim, nup.flags.dim, nlo.flags.boundary, nup.flags.boundary);
#ifdef HACK_REST_CANCEL	
				if (na.persistence <= gPersThreshold - 1.0) {
#else 
				if (na.persistence <= gPersThreshold) {
#endif
					sortedEdge se;
#ifdef HACK_REST_CANCEL
					ulong64 nid1 = this->nodes[na.lower].cellindex;
					ulong64 nid2 = this->nodes[na.upper].cellindex;
					bool same = this->sgg->bboundary[nid1] == this->sgg->bboundary[nid2];
					if (same) {
						na.persistence = na.persistence / maxp;
					}
					else {
						na.persistence = na.persistence + 1.0;
					}
#endif
					se.persistence = na.persistence;
					se.countweight = edgeCountWeight(na);
					se.ep = na_id;
					edges_to_cancel.push(se);
				}
				//if (arcs.size() == arcs.capacity()) {
				//	printf("about to expand capacity2 - arcs createnode!\n");
				//}
				this->arcs.push_back(na);
				return na_id;

			}
			void rec_tdcr(const INDEX_TYPE& cellid, DIM_TYPE& temp_dim, const INDEX_TYPE start) {
				INDEX_TYPE current = cellid;
				typename MESH_TYPE::FacetsIterator facets(mMesh);
				for (facets.begin(current); facets.valid(); facets.advance()) {
					INDEX_TYPE temp_id = facets.value();
					if (mGrad->getCritical(temp_id)) {
						//printf("adding arc: %llu\n", temp_id);
						createArc(temp_id, start);

					}
					else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
						INDEX_TYPE pair = mGrad->getPair(temp_id);
						if (pair != cellid && mMesh->dimension(pair) == mMesh->dimension(cellid)) {

							//result.push_back(pair);
							rec_tdcr(pair, temp_dim, start);

						}
					}

				}
			}


			void rec_tdcr_build_graph(const INDEX_TYPE& cellid, DIM_TYPE& temp_dim, map<INDEX_TYPE, pair<int, int> >& localmap) {
				INDEX_TYPE current = cellid;
				typename MESH_TYPE::FacetsIterator facets(mMesh);
				for (facets.begin(current); facets.valid(); facets.advance()) {
					INDEX_TYPE temp_id = facets.value();
					if (mGrad->getCritical(temp_id)) {
						if (localmap.count(temp_id) == 0) {
							localmap[temp_id] = pair<int, int>(1, 0);
						}
						else {
							localmap[temp_id].first += 1;
						}

					}
					else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
						INDEX_TYPE p = mGrad->getPair(temp_id);
						if (p != cellid &&  mMesh->dimension(p) == mMesh->dimension(cellid)) {
							if (localmap.count(p) == 0) {
								localmap[p] = pair<int, int>(1, 0);
								rec_tdcr_build_graph(p, temp_dim, localmap);
							}
							else {
								localmap[temp_id].first += 1;
							}
						}
					}

				}
			}

			void rec_tdcr_count_arcs(const INDEX_TYPE cellid, DIM_TYPE temp_dim, const INDEX_TYPE start,
				map<INDEX_TYPE, pair<int, int> >& localmap) {
				INDEX_TYPE current = cellid;
				typename MESH_TYPE::FacetsIterator facets(mMesh);
				for (facets.begin(current); facets.valid(); facets.advance()) {
					INDEX_TYPE temp_id = facets.value();
					if (mGrad->getCritical(temp_id)) {
						// if we are the last remaining cofacet to recurse on this, then do something
						if (localmap[temp_id].first == 1) {
							int count = localmap[temp_id].second + localmap[cellid].second;
							createArc(temp_id, start, count);
						}
						else {
							localmap[temp_id].second += localmap[cellid].second;
							localmap[temp_id].first--;
						}

					}
					else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
						INDEX_TYPE p = mGrad->getPair(temp_id);
						if (p != cellid &&  mMesh->dimension(p) == mMesh->dimension(cellid)) {
							if (localmap[p].first == 1) {
								localmap[p].second += localmap[cellid].second;
								rec_tdcr_count_arcs(p, temp_dim, start, localmap);
							}
							else {
								localmap[temp_id].second += localmap[cellid].second;
								localmap[temp_id].first--;
							}
						}
					}

				}
			}
			void trace_down_1_2(const INDEX_TYPE cellid) {

				// start from a 2-cell
				DIM_TYPE temp_dim = mGrad->getDimAscMan(cellid) + 1;
				map<INDEX_TYPE, pair<int, int> > localmap;
				localmap[cellid] = pair<int, int>(0, 1);
				rec_tdcr_build_graph(cellid, temp_dim, localmap);
				// now have # incoming edges seen in local map
				rec_tdcr_count_arcs(cellid, temp_dim, cellid, localmap);

			}
			void trace_down_cells_restricted(const INDEX_TYPE& cellid) {

				DIM_TYPE temp_dim = mGrad->getDimAscMan(cellid) + 1;
				rec_tdcr(cellid, temp_dim, cellid);
			}

			void ComputeFromGrad(bool restricted = false) {
				typename MESH_TYPE::AllCellsIterator t_cells(mMesh);
				printf("Finding critical points...\n");
				int count = 0;
				for (t_cells.begin(); t_cells.valid(); t_cells.advance()) {
					INDEX_TYPE t_id = t_cells.value();
					if (mGrad->getCritical(t_id)) {
						createNode(t_id);
						count++;
					}
				}
				printf("found %d critical points\n", count);

				// then add arcs
				if (!restricted) {
					printf("beginning arc tracing\n");
					for (INT_TYPE i = 0; i < nodes.size(); i++) {
						INDEX_TYPE t_id = nodes[i].cellindex;
						if (nodes[i].flags.dim == 2) {
							this->trace_down_1_2(t_id);
						}
						else {
							this->trace_down_cells_restricted(t_id);
						}
						//if (i % 1000 == 0) 
						// printf("traced %d arcs from %d nodes %d\n", arcs.size(), i, nodes[i].flags.dim);
					}
					printf("found %d unique arcs!\n", this->arcs.size());
				}
				//else {
				//	for (typename map< INDEX_TYPE, node<DTYPE>* >::iterator nit = nodes.begin();
				//		nit != nodes.end(); nit++) {
				//		INDEX_TYPE t_id = (*nit).first;
				//		if ((*nit).second->index == 1)
				//			this->trace_down_cells_restricted(t_id);

				//	}
				//}
				//ValidateComplex();
			}




			//////////// NOW CANCELLATION STUFF

			int countUpperArcs(INT_TYPE id) {
				return this->nodes[id].numarcs - this->nodes[id].numlower;
			}
			int countLowerArcs(INT_TYPE id) {
				return this->nodes[id].numlower;
			}

			int edgeCountWeight(arc<DTYPE> &a) {
				INT_TYPE lower = a.lower;
				INT_TYPE upper = a.upper;
				int nl = countUpperArcs(lower);
				int nu = countLowerArcs(upper);

				////first option: return the number created
				return (nl - 1) * (nu - 1);

				//2nd option: return the change in the number of arcs
				return (nl - 1) * (nu - 1) - (nl + nu - 1);

				//return 1;
			}

			inline INT_TYPE current_pers() { return select_persistence; }
			inline bool isAlive(arc<DTYPE>&a, INT_TYPE place) {
				return a.created <= place && a.destroyed > place;
			}
			bool isNodeAlive(INT_TYPE n) {
				return isAlive(this->nodes[n], select_persistence);
			}
			bool isArcAlive(INT_TYPE a) {
				bool res = isAlive(this->arcs[a], select_persistence);
				if (res && !isNodeAlive(this->arcs[a].lower)) printf("ERROR LOWER NODE NOT ALIVE\n");
				if (res && !isNodeAlive(this->arcs[a].upper)) printf("ERROR UPPER NODE NOT ALIVE\n");
				return isAlive(this->arcs[a], select_persistence);
			}

			inline bool isAlive(node<DTYPE>&n, INT_TYPE place) {
				return  n.destroyed > place;
			}
			inline bool isAlive(INT_TYPE n, INT_TYPE place) {
				return isAlive(this->nodes[n], place);
			}





			int cancel(INT_TYPE a) {

				int createcounter = 0;
				int deletecounter = 0;

				//int tmp1 = (int) edges_to_cancel.size();

				arc<DTYPE>* ap = &(this->arcs[a]);
				node<DTYPE>& lower = this->nodes[ap->lower];
				node<DTYPE>& upper = this->nodes[ap->upper];

				//if (output_cancellation_records) {
				cancellation_record cr;
				cr.index = lower.flags.dim;
				cr.lval = lower.value;
				cr.uval = upper.value;
				cr.persistence = ap->persistence;
				//DTYPE diff = this->sgg->sgd->maxval - this->sgg->sgd->minval;
				//cr.persPerc = 100.0f * ap->persistence / diff;
				cr.arcp = a;
				cr.boundary = lower.flags.boundary + upper.flags.boundary;
				this->mCRecords.push_back(cr);
				//}



				if (lower.destroyed != INT_INFTY) printf("ERRORORORORRO");
				if (upper.destroyed != INT_INFTY) printf("ERRORORORadfsORRO");
				//printf("\n");
				//printNode(ap->lower);
				//printNode(ap->upper);

				int initialguess = (lower.numarcs - lower.numlower - 1) * (upper.numlower - 1);
				int init2 = lower.numarcs + upper.numarcs - 1;

				if (ap->persistence > max_pers_so_far) max_pers_so_far = ap->persistence;
				num_cancelled++;


				// for each upwards arc connected to the lower node,
				// for each downwards arc connected to the upper node,
				// create a new connection from the other endpoints
				INT_TYPE la = lower.firstarc;
				while (la != NULLID) {
					arc<DTYPE>& lap = this->arcs[la];

					// skip the arc itself
					if (la == a) {
						la = lap.lower_next; // we're guarantee this is the next
						continue;
					}

					//test if they are the same kind
					if (ap->lower != lap.lower) {
						la = lap.upper_next;
						continue;
					}
					if (!isAlive(lap, num_cancelled)) {
						la = lap.lower_next;
						continue;
					}

					INT_TYPE ua = upper.firstarc;
					while (ua != NULLID) {
						arc<DTYPE>& uap = this->arcs[ua];

						// skip the arc itself
						if (ua == a) {
							ua = uap.upper_next; // we're guarantee this is the next
							continue;
						}

						//test if they are the same kind
						if (ap->upper != uap.upper) {
							ua = uap.lower_next;
							continue;
						}
						if (!isAlive(uap, num_cancelled)) {
							ua = uap.upper_next;
							continue;
						}

						// create the arc here!
						INT_TYPE newarc = createArc(la, ua, num_cancelled);
						//printf("ha1\n");
						ap = &(this->arcs[a]);
						//printf("ha\n");
						createcounter++;

						ua = this->arcs[ua].upper_next;
						//printf("asdf\n");
					}

					la = this->arcs[la].lower_next;

				}
				//printf("hahan\n");
				// go through and set the delete time on the arcs that are to be removed,
				// and update the arc counters at the other endpoints
				la = lower.firstarc;
				while (la != NULLID) {
					arc<DTYPE>& lap = this->arcs[la];

					// skip the arc itself
					if (la == a) {
						la = lap.lower_next; // we're guarantee this is the next
						continue;
					}

					//test if they are the same kind ////;TEST TO SEE IF ALIVE!!
					if (ap->lower == lap.lower) {
						if (isAlive(lap, num_cancelled)) {
							node<DTYPE>& n = this->nodes[lap.upper];

							n.numarcs--;
							n.numlower--;
							lap.destroyed = num_cancelled;
							deletecounter++;
						}
						la = lap.lower_next;
					}
					else if (ap->lower == lap.upper) {
						if (isAlive(lap, num_cancelled)) {
							node<DTYPE>& n = this->nodes[lap.lower];
							n.numarcs--;
							lap.destroyed = num_cancelled;
							deletecounter++;
						}
						la = lap.upper_next;
					}
					else {
						printf("ERROR SHOULD NEVER GET HERE\n");
					}



				}
				INT_TYPE ua = upper.firstarc;
				while (ua != NULLID) {
					arc<DTYPE>& uap = this->arcs[ua];

					// skip the arc itself
					if (ua == a) {
						ua = uap.upper_next; // we're guarantee this is the next
						continue;
					}
					//test if they are the same kind
					if (ap->upper == uap.upper) {
						if (isAlive(uap, num_cancelled)) {

							node<DTYPE>& n = this->nodes[uap.lower]; // pick the other

							n.numarcs--;
							uap.destroyed = num_cancelled;
							deletecounter++;
						}
						ua = uap.upper_next;
					}
					else if (ap->upper == uap.lower) {
						if (isAlive(uap, num_cancelled)) {
							node<DTYPE>& n = this->nodes[uap.upper];
							n.numarcs--;
							n.numlower--;
							uap.destroyed = num_cancelled;
							deletecounter++;
						}
						ua = uap.lower_next;
					}
					else {
						printf("ERROR SHOULD NEVER GET HERE\n");
					}


				}

				ap->destroyed = num_cancelled;
				lower.destroyed = num_cancelled;
				upper.destroyed = num_cancelled;
				deletecounter++;

				//int tmp2 = (int) edges_to_cancel.size();

				//printf("initialguess: %d actual:%d del:%d td:%d tmp1:%d tmp2:%d \n", initialguess, createcounter, deletecounter,init2
				//	,tmp1, tmp2);
				//printf("gothere\n");
				return 1;
			}



			INT_TYPE nextArc(arc<DTYPE>& ap, INT_TYPE n) {
				if (ap.lower == n) return ap.lower_next;
				if (ap.upper == n) return ap.upper_next;
				printf("\nERROR BARF POOP\n");
				return 0;
			}
			int countMultiplicity(arc<DTYPE>& ap, INT_TYPE ctime) {
				if (ap.multiplicity != -1) return ap.multiplicity;
				INT_TYPE nu = ap.upper;
				INT_TYPE nl = ap.lower;
				INT_TYPE a = this->nodes[nu].firstarc;
				int counter = 0;
				while (a != NULLID) {
					arc<DTYPE>& nap = this->arcs[a];
					if (isAlive(nap, ctime) && nap.lower == nl && nap.upper == nu) {
						if (nap.multiplicity == -1) {
							counter++;
						}
						else {
							counter += nap.multiplicity;
						}
					}
					a = nextArc(nap, nu);
				}
				return counter;
			}

			bool isValid(INT_TYPE a, arc<DTYPE>& ap) {
				// test for boundary
				if (this->nodes[ap.lower].flags.boundary !=
					this->nodes[ap.upper].flags.boundary) return false;

				// 2 endpoints must be connected by exactly one arc
				if (countMultiplicity(ap, num_cancelled) != 1) return false;
				// test for inversions?

				return true;


			}

			bool get_next_to_cancel(INT_TYPE& a) {
				///printf("getnext to cancel called\n");
				while (!edges_to_cancel.empty()) {
					sortedEdge se = edges_to_cancel.top();
					edges_to_cancel.pop();
					a = se.ep;
					arc<DTYPE>& ap = this->arcs[a];
					///printf("%u ", a);
					// is it alive in the current context
					if (!isAlive(ap, num_cancelled)) {
						//printf("adsf1\n");
						//printf("skip1->"); 
						///printArc(a);
						continue;
					}

					//test if it's a valid cancellation
					if (!isValid(a, ap)) {
						//printf("adsf2\n");
						continue;
					}

					int newcountweight = edgeCountWeight(ap);
					if (newcountweight > se.countweight) {
						//printf("adsf3\n");
						se.countweight = newcountweight;
						edges_to_cancel.push(se);
						continue;
					}

					if (newcountweight > 1500) {
						se.persistence += 1;
						if (se.persistence <= gPersThreshold)
							edges_to_cancel.push(se);
						continue;
					}

					//printf("getnext to cancel returned true\n");

					return true;


				}
				//printf("getnext to cancel returned false\n");
				return false;
			}

			vector<DTYPE> cancel_num_to_pers;
			INT_TYPE select_persistence; // the persistence to select in hierarchy
			INT_TYPE num_cancelled;

			void SetSelectPersAbs(DTYPE value) {
				int offset = cancel_num_to_pers.size() - 1;
				for (int i = 0; i < cancel_num_to_pers.size(); i += offset) {
					if (cancel_num_to_pers[i] > value) {
						select_persistence = i;
						return;
					}
					while (i + offset > cancel_num_to_pers.size() - 1) offset = offset / 2;
					while (offset > 1 && cancel_num_to_pers[i + offset] > value) offset = offset / 2;
				}
				select_persistence = cancel_num_to_pers.size();
			}


			DTYPE max_pers_so_far;
			DTYPE gPersThreshold;
			void ComputeHierarchy(DTYPE pers_limit) {
				cancel_num_to_pers.clear();
				printf(" Computing instance\n");
				gPersThreshold = pers_limit;
				max_pers_so_far = 0;
				num_cancelled = 0;
				// insert every arc to cancel list
				printf("  - Adding arcs to sorter...");
				INT_TYPE mysize = (INT_TYPE) this->arcs.size();
				for (INT_TYPE i = 0; i < mysize; i++) {
					// test if it passes the hierarchy test
					//if (!hierarchy_test->slice_position_list(i, this)) continue;
					arc<DTYPE> &a = this->arcs[i];
#ifdef HACK_REST_CANCEL
					if (a.persistence > gPersThreshold - 1.0) continue;
#else
					if (a.persistence > gPersThreshold) continue;
#endif
					sortedEdge se;
#ifdef HACK_REST_CANCEL
					ulong64 nid1 = this->nodes[a.lower].cellindex;
					ulong64 nid2 = this->nodes[a.upper].cellindex;
					bool same = this->sgg->bboundary[nid1] == this->sgg->bboundary[nid2];
					if (same) {
						a.persistence = a.persistence / maxp;
					}
					else {
						a.persistence = a.persistence + 1.0;
					}
#endif
					se.persistence = a.persistence;
					se.countweight = edgeCountWeight(a);
					se.ep = i;
					edges_to_cancel.push(se);
				}

				printf("Done!\n  - Cancelling:");

				INT_TYPE a;
				//float maxv = 0;

				while (get_next_to_cancel(a) && this->arcs[a].persistence <= gPersThreshold) {
					//if (arcs[a].persistence > maxv) maxv = arcs[a].persistence;
					if (num_cancelled % 1000 == 0) {
						printf("\r    - Cancelling: %u val=%f", num_cancelled, (float)max_pers_so_far);
					}
					//printf("\n\ncancelling: %d", num_cancelled);
					//printArc(a);

					cancel(a);
					cancel_num_to_pers.push_back(max_pers_so_far);

					//printf("Done\n");
				}
				printf("\r    - Cancelling finished: %u val=%f\n", num_cancelled, (float)max_pers_so_far);


				//_evercomputed = true;
				//if (output_cancellation_records) {
				//	int mincount = this->countIDs[0];
				//	int maxcount = this->countIDs[3];
				//	int interiormincount = this->countInteriorIDs[0];
				//	FILE* fout = fopen(this->mCRecordFName, "w");
				//	for (int i = 0; i < this->mCRecords.size(); i++) {
				//		cancellation_record& cr = this->mCRecords[i];
				//		fprintf(fout, "%f %f %f %f %d %d %d %d %d\n", cr.persistence,
				//			cr.lval, cr.uval, cr.persPerc, cr.index, mincount, maxcount,
				//			interiormincount, cr.boundary);
				//		if (cr.index == 0) {
				//			mincount--;
				//			node<DTYPE>& n =
				//				this->getNode(this->getArc(cr.arcp).lower);
				//			if (!n.flags.boundary) {
				//				interiormincount--;
				//			}
				//		}
				//		if (cr.index == 2) maxcount--;
				//	}
				//	fclose(fout);
				//}
				// now fill in the lists?	
			}
			void Output1SaddleRecord(const char* fname) {

				FILE* fsadrec = fopen(fname, "w");

				for (int i = 0; i < nodes.size(); i++) {
					node<DTYPE>& n = nodes[i];
					if (n.flags.dim > 1) continue;
					float persval =
						(n.destroyed < cancel_num_to_pers.size() + 1 ?
							cancel_num_to_pers[n.destroyed - 1] :
							cancel_num_to_pers[cancel_num_to_pers.size() - 1]);
					fprintf(fsadrec, "%d %d %llu %f %d %f\n", n.flags.dim, n.flags.boundary, n.cellindex, n.value, n.destroyed, persval);


				}
				fclose(fsadrec);


			}
		};




	};

	  

	namespace LightGeomMSC {

		struct msbitfield
		{
			//// rehash this
			unsigned char dim : 3;
			unsigned char boundary : 1;
			unsigned char f1 : 1;
			unsigned char f2 : 1;
			unsigned char f3 : 2;
		};

		struct manifold {
			INT_TYPE merged[2]; // the manifold ids
			INT_TYPE basenode; // the node id
			INT_TYPE mergetime;
		};

		template<typename DTYPE>
		struct node
		{
			INDEX_TYPE cellindex;
			INT_TYPE firstarc;
			INT_TYPE destroyed; // the time this node is cancelled.
			INT_TYPE amanifoldid; // set to -1 if this is the base
			INT_TYPE dmanifoldid; // set to -1 if this is the base
			unsigned short numarcs;
			unsigned short numlower;
			DTYPE value;
			msbitfield flags;
		};


		template<typename DTYPE>
		struct arc
		{
			INT_TYPE lower; // node
			INT_TYPE lower_next; //arc - INT_INFTY is the null arc!
			INT_TYPE upper; // node
			INT_TYPE upper_next; //arc
			INT_TYPE created;
			INT_TYPE destroyed;
			int multiplicity;
			DTYPE persistence;
			msbitfield flags; // the .dim stores the dim of the lower endpoint
		};

		template<typename DTYPE, class MESH_TYPE, class FUNC_TYPE, class GRAD_TYPE>
		class SuperLightMSC {
		protected:
			struct sortedEdge {
				DTYPE persistence;
				int countweight;
				INT_TYPE ep;

				bool operator()(const sortedEdge& _Left, const sortedEdge& _Right) {
					if (_Left.persistence < _Right.persistence) return false;
					if (_Left.persistence > _Right.persistence) return true;
					if (_Left.countweight < _Right.countweight) return false;
					if (_Left.countweight > _Right.countweight) return true;
					return _Left.ep > _Right.ep;
				}
			};

			priority_queue<sortedEdge, vector< sortedEdge>, sortedEdge > edges_to_cancel;

			struct cancellation_record {
				int index;
				DTYPE persistence;
				DTYPE lval;
				DTYPE uval;
				DTYPE persPerc;
				INT_TYPE arcp;
				int boundary;
			};

			vector<cancellation_record> mCRecords;
			bool output_cancellation_records;

			const char* mCRecordFName;
			float m_temp_perc_pers;
			//INDEX_TYPE select_persistence;
			INDEX_TYPE num_destroyed;
			GRAD_TYPE* mGrad;
			MESH_TYPE* mMesh;
			FUNC_TYPE* mFunc;

			vector<node<DTYPE> > nodes;
			vector<arc<DTYPE> > arcs;
			vector<manifold> mans;

#ifdef WIN32
			unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#elif defined(__APPLE__)
			unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#else
	unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
			//__gnu_cxx::hash_map<INDEX_TYPE, INT_TYPE> nodemap;
#endif
			inline void connectArc(INT_TYPE arcID, arc<DTYPE>& a, node<DTYPE>& lower, node<DTYPE>& upper) {
				a.lower_next = lower.firstarc;
				lower.firstarc = arcID;
				a.upper_next = upper.firstarc;
				upper.firstarc = arcID;

				lower.numarcs++;
				upper.numarcs++;
				upper.numlower++;

			}

			void connectArc(INT_TYPE arcID) {

				arc<DTYPE>& a = arcs[arcID];
				node<DTYPE>& lower = nodes[a.lower];
				node<DTYPE>& upper = nodes[a.upper];
				connectArc(arcID, a, lower, upper);
			}


		public:
			void printMemory() {
				printf("SuperLightMSC memory size: \n");
				printf("\t %d nodes x nodesize %d = %d bytes\n", nodes.size(), sizeof(node<DTYPE>), nodes.size() * sizeof(node<DTYPE>));
				printf("\t %d arcs x arcsize %d = %d bytes\n", arcs.size(), sizeof(arc<DTYPE>), arcs.size() * sizeof(arc<DTYPE>));
				printf("\t %d manifolds x manifoldsize %d = %d bytes\n", mans.size(), sizeof(manifold), mans.size() * sizeof(manifold));
				printf("\t %d canrecs x canrec %d = %d bytes\n", mCRecords.size(), sizeof(cancellation_record), mCRecords.size() * sizeof(cancellation_record));

				size_t count = 0;
				for (unsigned i = 0; i < nodemap.bucket_count(); ++i) {
					size_t bucket_size = nodemap.bucket_size(i);
					if (bucket_size == 0) {
						count++;
					}
					else {
						count += bucket_size;
					}
				}

				printf("\t nodemap size: %d bytes\n", count);
			}
			INT_TYPE numArcs() {
				return (INT_TYPE)arcs.size();
			}
			INT_TYPE numNodes() {
				return (INT_TYPE)nodes.size();
			}
			node<DTYPE>& getNode(INT_TYPE e) {
				return nodes[e];
			}

			const manifold& getManifold(INT_TYPE m) const{
				return mans[m];
			}

			arc<DTYPE>& getArc(INT_TYPE e) {
				return arcs[e];
			}


			void set_output_cancellation_records(const char* fname) {
				output_cancellation_records = true;
				mCRecordFName = fname;
			}


			SuperLightMSC(GRAD_TYPE* grad,
				MESH_TYPE* mesh,
				FUNC_TYPE* func) :
				mGrad(grad), mMesh(mesh), mFunc(func)
			{
				num_destroyed = 0;
				select_persistence = 0;
				output_cancellation_records = false;

			}


			//INT_TYPE createNode(node<DTYPE> &n) {
			//	INT_TYPE listID = (INT_TYPE)nodes.push_back(n);
			//	nodemap[n.cellindex] = listID;
			//	//this->countIDs[n.flags.dim]++;
			//	//if (!((bool)n.flags.boundary)) this->countInteriorIDs[n.flags.dim]++;
			//	return listID;
			//}

			INT_TYPE createManifold(INT_TYPE baseId, INT_TYPE baseManId, INT_TYPE mergeManId, INT_TYPE mergetime) {
				INT_TYPE manid = mans.size();
				manifold m;
				m.basenode = baseId;
				m.merged[0] = baseManId;
				m.merged[1] = mergeManId;
				m.mergetime = mergetime;
				mans.push_back(m);
				return manid;
			}

			INT_TYPE createNode(INDEX_TYPE cellID) {
				node<DTYPE> tn;
				INT_TYPE listID = nodes.size();
				//if (nodes.size() == nodes.capacity()) printf("about to expand capacity - nodes createnode!\b");
				nodes.push_back(tn);
				nodemap[cellID] = listID;

				//node<DTYPE> &n = nodes->operator [](listID);
				node<DTYPE> &n = nodes[listID];
				n.cellindex = cellID;

				n.destroyed = INT_INFTY;
				n.firstarc = NULLID;
				n.flags.boundary = mMesh->boundaryValue(cellID);
				//if (n.flags.boundary) {
				//	INDEX_TYPE coords[3];
				//	sgg->getCoords(coords, cellID);
				//	printf("(%d, %d, %d) = %d\n", 
				//	   (int) coords[0], (int) coords[1], (int) coords[2], n.flags.boundary);
				//}
				n.flags.dim = mMesh->dimension(cellID);
				n.numarcs = 0;
				n.numlower = 0;
				n.value = mFunc->cellValue(cellID);

				n.amanifoldid = createManifold(listID, -1, -1, 0);
				n.dmanifoldid = createManifold(listID, -1, -1, 0);

				return listID;
			}

			//// create an arc connecting lower to upper <- cell id's from gradient
			//INT_TYPE createArc(arc<DTYPE> &a){
			//	INT_TYPE listID = (INT_TYPE)arcs.push_back(a);
			//	//printf("%f - %f\n", (float) a.persistence, arcs[listID].persistence); 
			//	return listID;
			//}

			INT_TYPE createArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID, int count) {

				arc<DTYPE> at;
				INT_TYPE listID = arcs.size();
				//if (arcs.size() == arcs.capacity()) printf("about to expand capacity - arcs createnode!\b");

				arcs.push_back(at);

				arc<DTYPE>& a = arcs[listID];

				//if (count != -1) printf("creating arc with count %d\n", count);
				a.multiplicity = count;
				a.created = 0;
				a.destroyed = INT_INFTY;
				a.lower = nodemap[lowerCellID];
				a.upper = nodemap[upperCellID];
				//node<DTYPE>& ln = nodes[a.lower];
				//node<DTYPE>& un = nodes[a.upper];
				//sgg->printCell(lowerCellID); printf("->"); sgg->printCell(upperCellID); printf("\n");
				//if (ln.flags.dim != un.flags.dim - 1) printf("ERROR IN CREATE ARC\n");
				if (nodes[a.upper].value < nodes[a.lower].value)
					printf("ERROR: upper (%f, %d) - lower (%f, %d)\n",
						float(nodes[a.upper].value), nodes[a.upper].flags.dim,
						float(nodes[a.lower].value), nodes[a.lower].flags.dim);

				a.persistence = nodes[a.upper].value - nodes[a.lower].value;
				//printf("n=%f\n", (float) a.persistence);
				a.flags.dim = nodes[a.lower].flags.dim;
				// now connect the crap!
				connectArc(listID);

				//printf("creating arc: \n\t");
				//sgg->printCell(lowerCellID);
				//for (int i = 0; i < g->size(); i++) {
				//	printf("\n");
				//	sgg->printCell(g->operator [](i));
				//}


				//printf("\n\t");
				//sgg->printCell(upperCellID);
				//printf("\n\t%d\n", g->size());
				/*		for (int i=0; i < g->size(); i++) printf("%d ", g->operator [](i));
				printf("\n")*/;

				//printf("Creating arc: %d-%d, %d-%d\n", lowerCellID, upperCellID, nodes[a.lower].flags.dim, nodes[a.upper].flags.dim);

				return listID;
			}

			INT_TYPE createArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID) {
				return createArc(lowerCellID, upperCellID, -1);
			}

			INT_TYPE createArc(INT_TYPE luap, INT_TYPE ulap, INT_TYPE ctime) {
				arc<DTYPE>& lua = arcs[luap];
				arc<DTYPE>& ula = arcs[ulap];
				INT_TYPE na_id = arcs.size();
				//this->arcs.push_back(ma); // copy setting from old

				arc<DTYPE> na;// = this->arcs[na_id];
				na.created = ctime;
				na.destroyed = INT_INFTY;
				na.lower = ula.lower;
				na.upper = lua.upper;
				na.multiplicity = -1;

				node<DTYPE>& nup = this->nodes[na.upper];
				node<DTYPE>& nlo = this->nodes[na.lower];

				this->connectArc(na_id, na, nlo, nup);
				na.persistence = nup.value - nlo.value;
				na.flags.boundary = nlo.flags.boundary || nup.flags.boundary;
				na.flags.dim = nlo.flags.dim;
				if (na.persistence < 0 && nlo.flags.boundary == nup.flags.boundary) printf("creatinga inversion %f, %d-%d, b%d-b%d\n",
					(float)na.persistence, nlo.flags.dim, nup.flags.dim, nlo.flags.boundary, nup.flags.boundary);
#ifdef HACK_REST_CANCEL	
				if (na.persistence <= gPersThreshold - 1.0) {
#else 
				if (na.persistence <= gPersThreshold) {
#endif
					sortedEdge se;
#ifdef HACK_REST_CANCEL
					ulong64 nid1 = this->nodes[na.lower].cellindex;
					ulong64 nid2 = this->nodes[na.upper].cellindex;
					bool same = this->sgg->bboundary[nid1] == this->sgg->bboundary[nid2];
					if (same) {
						na.persistence = na.persistence / maxp;
					}
					else {
						na.persistence = na.persistence + 1.0;
					}
#endif
					se.persistence = na.persistence;
					se.countweight = edgeCountWeight(na);
					se.ep = na_id;
					edges_to_cancel.push(se);
				}
				//if (arcs.size() == arcs.capacity()) {
				//	printf("about to expand capacity2 - arcs createnode!\n");
				//}
				this->arcs.push_back(na);
				return na_id;

				}
			void rec_tdcr(const INDEX_TYPE& cellid, DIM_TYPE& temp_dim, const INDEX_TYPE start) {
				INDEX_TYPE current = cellid;
				typename MESH_TYPE::FacetsIterator facets(mMesh);
				for (facets.begin(current); facets.valid(); facets.advance()) {
					INDEX_TYPE temp_id = facets.value();
					if (mGrad->getCritical(temp_id)) {
						//printf("adding arc: %llu\n", temp_id);
						createArc(temp_id, start);

					}
					else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
						INDEX_TYPE pair = mGrad->getPair(temp_id);
						if (pair != cellid && mMesh->dimension(pair) == mMesh->dimension(cellid)) {

							//result.push_back(pair);
							rec_tdcr(pair, temp_dim, start);

						}
					}

				}
			}


			void rec_tdcr_build_graph(const INDEX_TYPE& cellid, DIM_TYPE& temp_dim, map<INDEX_TYPE, pair<int, int> >& localmap) {
				INDEX_TYPE current = cellid;
				typename MESH_TYPE::FacetsIterator facets(mMesh);
				for (facets.begin(current); facets.valid(); facets.advance()) {
					INDEX_TYPE temp_id = facets.value();
					if (mGrad->getCritical(temp_id)) {
						if (localmap.count(temp_id) == 0) {
							localmap[temp_id] = pair<int, int>(1, 0);
						}
						else {
							localmap[temp_id].first += 1;
						}

					}
					else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
						INDEX_TYPE p = mGrad->getPair(temp_id);
						if (p != cellid &&  mMesh->dimension(p) == mMesh->dimension(cellid)) {
							if (localmap.count(p) == 0) {
								localmap[p] = pair<int, int>(1, 0);
								rec_tdcr_build_graph(p, temp_dim, localmap);
							}
							else {
								localmap[temp_id].first += 1;
							}
						}
					}

				}
			}

			void rec_tdcr_count_arcs(const INDEX_TYPE cellid, DIM_TYPE temp_dim, const INDEX_TYPE start,
				map<INDEX_TYPE, pair<int, int> >& localmap) {
				INDEX_TYPE current = cellid;
				typename MESH_TYPE::FacetsIterator facets(mMesh);
				for (facets.begin(current); facets.valid(); facets.advance()) {
					INDEX_TYPE temp_id = facets.value();
					if (mGrad->getCritical(temp_id)) {
						// if we are the last remaining cofacet to recurse on this, then do something
						if (localmap[temp_id].first == 1) {
							int count = localmap[temp_id].second + localmap[cellid].second;
							createArc(temp_id, start, count);
						}
						else {
							localmap[temp_id].second += localmap[cellid].second;
							localmap[temp_id].first--;
						}

					}
					else if (mGrad->getDimAscMan(temp_id) == temp_dim) {
						INDEX_TYPE p = mGrad->getPair(temp_id);
						if (p != cellid &&  mMesh->dimension(p) == mMesh->dimension(cellid)) {
							if (localmap[p].first == 1) {
								localmap[p].second += localmap[cellid].second;
								rec_tdcr_count_arcs(p, temp_dim, start, localmap);
							}
							else {
								localmap[temp_id].second += localmap[cellid].second;
								localmap[temp_id].first--;
							}
						}
					}

				}
			}
			void trace_down_1_2(const INDEX_TYPE cellid) {

				// start from a 2-cell
				DIM_TYPE temp_dim = mGrad->getDimAscMan(cellid) + 1;
				map<INDEX_TYPE, pair<int, int> > localmap;
				localmap[cellid] = pair<int, int>(0, 1);
				rec_tdcr_build_graph(cellid, temp_dim, localmap);
				// now have # incoming edges seen in local map
				rec_tdcr_count_arcs(cellid, temp_dim, cellid, localmap);

			}
			void trace_down_cells_restricted(const INDEX_TYPE& cellid) {

				DIM_TYPE temp_dim = mGrad->getDimAscMan(cellid) + 1;
				rec_tdcr(cellid, temp_dim, cellid);
			}

			// only works with dim 1 right now
			void rec_tdc_gather(const INDEX_TYPE& cellid, vector<INDEX_TYPE>& res) {
				INDEX_TYPE current = cellid;
				typename MESH_TYPE::FacetsIterator facets(mMesh);
				for (facets.begin(current); facets.valid(); facets.advance()) {
					INDEX_TYPE temp_id = facets.value();
					if (mGrad->getCritical(temp_id)) {
						//printf("adding arc: %llu\n", temp_id);
						res.push_back(temp_id);

					}
					else {
						INDEX_TYPE pair = mGrad->getPair(temp_id);
						if (pair != cellid/* && mMesh->dimension(pair) == mMesh->dimension(cellid)*/) { // can avoid check on only 0-1 connections
																										//result.push_back(pair);
							rec_tdc_gather(pair, res);

						}
					}

				}
			}

			//mutex AddNodesMutex;
			void AddArcs(INT_TYPE start) {
				//vector<pair<INDEX_TYPE, INDEX_TYPE>> pairs;

				map<pair<INDEX_TYPE, INDEX_TYPE>, pair<float, INDEX_TYPE> > endtoendmap;
				INT_TYPE i = start; 
					node<DTYPE>& n = nodes[i];

					if (n.flags.dim == 1) {

						vector<INDEX_TYPE> destinations;
						this->rec_tdc_gather(n.cellindex, destinations);

						if (destinations.size() != 2) printf("Error, %d dests found\n", destinations.size());


						INDEX_TYPE d1 = destinations[0];
						INDEX_TYPE d2 = destinations[1];

						//if (d1 == d2) continue; // no point adding self loops

						if (d2 < d1) {
							d1 = destinations[1];
							d2 = destinations[0];
						}

						pair<INDEX_TYPE, INDEX_TYPE> npair(d1, d2);
						if (endtoendmap.count(npair) == 0) {
							endtoendmap[npair] = pair<float, INDEX_TYPE>(n.value, n.cellindex);
						}
						else {
							if (endtoendmap[npair].first > n.value) {
								endtoendmap[npair] = pair<float, INDEX_TYPE>(n.value, n.cellindex);
							}
						}

						//
						//sgg->addDescendingArcsG(n.cellindex, destinations, false);
						//
						//for (int j = 0; j < destinations.size(); j++) {
						//	pairs.push_back(pair<INDEX_TYPE, INDEX_TYPE>(n.cellindex, destinations[j]));
						//}

					}
				

				//unique_lock<mutex> lock(this->AddNodesMutex);

				for (auto pp : endtoendmap) {
					this->createArc(pp.first.first, pp.second.second);
					this->createArc(pp.first.second, pp.second.second);
				}
				//for (int i = 0; i < pairs.size(); i++) {
				//	this->createArc(pairs[i].second, pairs[i].first);
				//}
			}


			void ComputeFromGrad(bool restricted = false) {
				typename MESH_TYPE::AllCellsIterator t_cells(mMesh);
				printf("Finding critical points...\n");
				int count = 0;
				for (t_cells.begin(); t_cells.valid(); t_cells.advance()) {
					INDEX_TYPE t_id = t_cells.value();
					// only add minima and saddles
					if (mGrad->getCritical(t_id) && mMesh->dimension(t_id) < 2) {
						createNode(t_id);
						count++;
					}
				}
				printf("found %d critical points\n", count);

				printf("Finding Edges...\n");
				// then add arcs
				//vector<thread> ThreadPool;
				
				for (int i = 0; i < this->nodes.size(); i++) {
					AddArcs(i);
				}
				printf("found %d Edges\n", this->arcs.size());

			}




			//////////// NOW CANCELLATION STUFF

			int countUpperArcs(INT_TYPE id) {
				return this->nodes[id].numarcs - this->nodes[id].numlower;
			}
			int countLowerArcs(INT_TYPE id) {
				return this->nodes[id].numlower;
			}

			int edgeCountWeight(arc<DTYPE> &a) {
				INT_TYPE lower = a.lower;
				INT_TYPE upper = a.upper;
				int nl = countUpperArcs(lower);
				int nu = countLowerArcs(upper);

				////first option: return the number created
				return (nl - 1) * (nu - 1);

				//2nd option: return the change in the number of arcs
				return (nl - 1) * (nu - 1) - (nl + nu - 1);

				//return 1;
			}

			inline INT_TYPE current_pers() { return select_persistence; }
			inline bool isAlive(arc<DTYPE>&a, INT_TYPE place) {
				return a.created <= place && a.destroyed > place;
			}
			bool isNodeAlive(INT_TYPE n) {
				return isAlive(this->nodes[n], select_persistence);
			}
			bool isArcAlive(INT_TYPE a) {
				bool res = isAlive(this->arcs[a], select_persistence);
				if (res && !isNodeAlive(this->arcs[a].lower)) printf("ERROR LOWER NODE NOT ALIVE\n");
				if (res && !isNodeAlive(this->arcs[a].upper)) printf("ERROR UPPER NODE NOT ALIVE\n");
				return isAlive(this->arcs[a], select_persistence);
			}

			inline bool isAlive(node<DTYPE>&n, INT_TYPE place) {
				return  n.destroyed > place;
			}
			inline bool isAlive(INT_TYPE n, INT_TYPE place) {
				return isAlive(this->nodes[n], place);
			}


			void SetSelectPersAbs(DTYPE value) {
				printf("mcLightGeomMSC::SuperLightMSC::SetSelectPersAbs -> %f\n", value);
				int offset = 1;// cancel_num_to_pers.size() - 1;
				for (int i = 0; i < cancel_num_to_pers.size(); i += offset) {
					if (cancel_num_to_pers[i] > value) {
						select_persistence = i;

						return;
					}
					//while (i + offset > cancel_num_to_pers.size() - 1) offset = offset / 2;
					//while (offset > 1 && cancel_num_to_pers[i + offset] > value) offset = offset / 2;
				}
				select_persistence = cancel_num_to_pers.size() - 1;
			}


			int cancel(INT_TYPE a) {

				int createcounter = 0;
				int deletecounter = 0;

				//int tmp1 = (int) edges_to_cancel.size();

				arc<DTYPE>* ap = &(this->arcs[a]);
				node<DTYPE>& lower = this->nodes[ap->lower];
				node<DTYPE>& upper = this->nodes[ap->upper];

				//if (output_cancellation_records) {
				cancellation_record cr;
				cr.index = lower.flags.dim;
				cr.lval = lower.value;
				cr.uval = upper.value;
				cr.persistence = ap->persistence;
				//DTYPE diff = this->sgg->sgd->maxval - this->sgg->sgd->minval;
				//cr.persPerc = 100.0f * ap->persistence / diff;
				cr.arcp = a;
				cr.boundary = lower.flags.boundary + upper.flags.boundary;
				this->mCRecords.push_back(cr);
				//}



				if (lower.destroyed != INT_INFTY) printf("ERRORORORORRO");
				if (upper.destroyed != INT_INFTY) printf("ERRORORORadfsORRO");
				//printf("\n");
				//printNode(ap->lower);
				//printNode(ap->upper);

				int initialguess = (lower.numarcs - lower.numlower - 1) * (upper.numlower - 1);
				int init2 = lower.numarcs + upper.numarcs - 1;

				if (ap->persistence > max_pers_so_far) max_pers_so_far = ap->persistence;
				num_cancelled++;


				// for each upwards arc connected to the lower node,
				// for each downwards arc connected to the upper node,
				// create a new connection from the other endpoints
				INT_TYPE la = lower.firstarc;
				while (la != NULLID) {
					arc<DTYPE>& lap = this->arcs[la];

					// skip the arc itself
					if (la == a) {
						la = lap.lower_next; // we're guarantee this is the next
						continue;
					}

					//test if they are the same kind
					if (ap->lower != lap.lower) {
						la = lap.upper_next;
						continue;
					}
					if (!isAlive(lap, num_cancelled)) {
						la = lap.lower_next;
						continue;
					}

					INT_TYPE ua = upper.firstarc;
					while (ua != NULLID) {
						arc<DTYPE>& uap = this->arcs[ua];

						// skip the arc itself
						if (ua == a) {
							ua = uap.upper_next; // we're guarantee this is the next
							continue;
						}

						//test if they are the same kind
						if (ap->upper != uap.upper) {
							ua = uap.lower_next;
							continue;
						}
						if (!isAlive(uap, num_cancelled)) {
							ua = uap.upper_next;
							continue;
						}

						// create the arc here!
						INT_TYPE newarc = createArc(la, ua, num_cancelled);
						//printf("ha1\n");
						ap = &(this->arcs[a]);
						//printf("ha\n");
						createcounter++;

						ua = this->arcs[ua].upper_next;
						//printf("asdf\n");
					}

					la = this->arcs[la].lower_next;

				}
				//printf("hahan\n");
				// go through and set the delete time on the arcs that are to be removed,
				// and update the arc counters at the other endpoints

				// following comments are for 0-1 cancellations
				// first look in neighborhood of minimum = lower
				la = lower.firstarc;							// pick first arc attached to min
				while (la != NULLID) {							// while there are more arcs in neighborhood, keep looking
					arc<DTYPE>& lap = this->arcs[la];			// get the reference to the arc attached to the min

																// skip the arc itself
					if (la == a) {
						la = lap.lower_next; // we're guarantee this is the next
						continue;
					}

					//test if they are the same kind 
					if (ap->lower == lap.lower) {				// we do not want arcs that are not 0-1, so do this check
						if (isAlive(lap, num_cancelled)) {		// make sure we are only working with living arcs
							node<DTYPE>& n = this->nodes[lap.upper]; // now "n" is the 1-saddle attached to minimum

							n.numarcs--;						// if we remove this arc, we need to decrement the number of living arcs attached to it
							n.numlower--;						// this cancellation will also remove a "lower" arc of the 1-saddle

																// create manifold for merging
																// extend the descending manifold of the 1-saddle by merging it with the descending manifold of the cancelled 1-saddle
							INT_TYPE nmanid =
								createManifold(lap.upper, n.dmanifoldid, upper.dmanifoldid, num_cancelled);
							n.dmanifoldid = nmanid;				// change the manifold reference of the 1-saddle to reference this new decending manifold

							lap.destroyed = num_cancelled;		// we remove the old arc to the deleted min, so mark it destroyed with the num_cancelld
							deletecounter++;					// keep track of deleted arc cound
						}
						la = lap.lower_next;					// continue on the next arc around the minimum
					}
					else if (ap->lower == lap.upper) {          // if this is in fact an i-1 - i arc, we just remove the arcs around it - no merging happens
						if (isAlive(lap, num_cancelled)) {
							node<DTYPE>& n = this->nodes[lap.lower];
							n.numarcs--;
							lap.destroyed = num_cancelled;
							deletecounter++;
						}
						la = lap.upper_next;
					}
					else {
						printf("ERROR SHOULD NEVER GET HERE\n");
					}



				}
				INT_TYPE ua = upper.firstarc;
				while (ua != NULLID) {
					arc<DTYPE>& uap = this->arcs[ua];

					// skip the arc itself
					if (ua == a) {
						ua = uap.upper_next; // we're guarantee this is the next
						continue;
					}
					//test if they are the same kind
					if (ap->upper == uap.upper) {
						if (isAlive(uap, num_cancelled)) {

							node<DTYPE>& n = this->nodes[uap.lower]; // pick the other

							n.numarcs--;

							// create manifold for merging
							INT_TYPE nmanid = createManifold(uap.lower, n.amanifoldid, lower.amanifoldid, num_cancelled);
							n.amanifoldid = nmanid;

							uap.destroyed = num_cancelled;
							deletecounter++;
						}
						ua = uap.upper_next;
					}
					else if (ap->upper == uap.lower) {
						if (isAlive(uap, num_cancelled)) {
							node<DTYPE>& n = this->nodes[uap.upper];
							n.numarcs--;
							n.numlower--;
							uap.destroyed = num_cancelled;
							deletecounter++;
						}
						ua = uap.lower_next;
					}
					else {
						printf("ERROR SHOULD NEVER GET HERE\n");
					}


				}

				ap->destroyed = num_cancelled;
				lower.destroyed = num_cancelled;
				upper.destroyed = num_cancelled;
				deletecounter++;

				//int tmp2 = (int) edges_to_cancel.size();

				//printf("initialguess: %d actual:%d del:%d td:%d tmp1:%d tmp2:%d \n", initialguess, createcounter, deletecounter,init2
				//	,tmp1, tmp2);
				//printf("gothere\n");
				return 1;
			}



			INT_TYPE nextArc(arc<DTYPE>& ap, INT_TYPE n) {
				if (ap.lower == n) return ap.lower_next;
				if (ap.upper == n) return ap.upper_next;
				printf("\nERROR BARF POOP\n");
				return 0;
			}
			int countMultiplicity(arc<DTYPE>& ap, INT_TYPE ctime) {
				if (ap.multiplicity != -1) return ap.multiplicity;
				INT_TYPE nu = ap.upper;
				INT_TYPE nl = ap.lower;
				INT_TYPE a = this->nodes[nu].firstarc;
				int counter = 0;
				while (a != NULLID) {
					arc<DTYPE>& nap = this->arcs[a];
					if (isAlive(nap, ctime) && nap.lower == nl && nap.upper == nu) {
						if (nap.multiplicity == -1) {
							counter++;
						}
						else {
							counter += nap.multiplicity;
						}
					}
					a = nextArc(nap, nu);
				}
				return counter;
			}

			bool isValid(INT_TYPE a, arc<DTYPE>& ap) {
				// test for boundary
				if (this->nodes[ap.lower].flags.boundary !=
					this->nodes[ap.upper].flags.boundary) return false;

				// 2 endpoints must be connected by exactly one arc
				if (countMultiplicity(ap, num_cancelled) != 1) return false;
				// test for inversions?

				return true;


			}

			// THIS WILL HAVE TO WORK LATER, but for now only need 0 and maxdim manifolds
			//void recCounLeaftManifolds(INT_TYPE mId, map< INT_TYPE, int >& counter) {
			//	manifold &m = this->getManifold(mId);
			//	if (m->merge[0] != NULL) {
			//		recCounLeaftManifolds(m->merge[0], counter);
			//		recCounLeaftManifolds(m->merge[1], counter);
			//	}
			//	else {
			//		if (counter.count(m) == 0) {
			//			counter[m] = 1;
			//		}
			//		else {
			//			counter[m]++;
			//		}
			//	}
			//}

			void printmanifold(INT_TYPE man) {
				manifold& m = getManifold(man);
				printf("man=%d, man.base=%d, man.mergetime=%d, man.merge[0]=%d, man.merge[1]=%d\n",
					man, m.basenode, m.mergetime, m.merged[0], m.merged[1]);
			}

			INT_TYPE getActiveMan(INT_TYPE man) {
				//printmanifold(man); 
				while (getManifold(man).mergetime > this->current_pers()) {
					man = getManifold(man).merged[0];
					//printf("  --");  printmanifold(man);
				}
				//printf("ret %d\n", man);
				return man;
			}

			void GatherNodes(INT_TYPE nodeID, set<INT_TYPE>& res) {

				if (!isAlive(nodeID, this->select_persistence)) {
					return;
				}
				node<DTYPE>& n = this->getNode(nodeID);
				INT_TYPE man;
				if (n.flags.dim < this->mMesh->maxDim() / 2) {
					man = n.amanifoldid;
				}
				else {
					man = n.dmanifoldid;
				}
				INT_TYPE manID = getActiveMan(man);

				if (n.flags.dim == 0 || n.flags.dim == this->mMesh->maxDim())
				{
					//set<INT_TYPE> basenodes;
					recGatherNodes(manID, res);
				}

			}

			void recGatherNodes(INT_TYPE mid, set<INT_TYPE>& res) {
				manifold& m = getManifold(mid);
				if (m.merged[0] != -1) {
					recGatherNodes(m.merged[0], res);
					recGatherNodes(m.merged[1], res);
				}
				//printf("gothere\n");
				res.insert(m.basenode);
			}

			void rec_man_trace_up(INDEX_TYPE& cellid, set<INDEX_TYPE>& res) {
				res.insert(cellid);
				INDEX_TYPE current = cellid;
				DIM_TYPE cdim = mMesh->dimension(cellid);
				typename MESH_TYPE::CofacetsIterator cofacets(mMesh);
				for (cofacets.begin(current); cofacets.valid(); cofacets.advance()) {
					INDEX_TYPE temp_id = cofacets.value();
					if (mGrad->getCritical(temp_id)) continue;

					INDEX_TYPE temp_pair = mGrad->getPair(temp_id);

					if (temp_pair == cellid) continue;

					if (mMesh->dimension(temp_pair) != cdim) continue;

					rec_man_trace_up(temp_pair, res);
				}
			}
			void rec_man_trace_down(INDEX_TYPE& cellid, set<INDEX_TYPE>& res) {
				res.insert(cellid);
				INDEX_TYPE current = cellid;
				DIM_TYPE cdim = mMesh->dimension(cellid);
				typename MESH_TYPE::FacetsIterator facets(mMesh);
				for (facets.begin(current); facets.valid(); facets.advance()) {
					INDEX_TYPE temp_id = facets.value();
					if (mGrad->getCritical(temp_id)) continue;

					INDEX_TYPE temp_pair = mGrad->getPair(temp_id);

					if (temp_pair == cellid) continue;

					if (mMesh->dimension(temp_pair) != cdim) continue;

					rec_man_trace_down(temp_pair, res);
				}
			}

			void fillUnsimplifiedGeometry(INT_TYPE nodeID, set<INDEX_TYPE>& res) {

				node<DTYPE>& n = this->getNode(nodeID);

				if (n.flags.dim == 0) {
					rec_man_trace_up(n.cellindex, res);
				}
				else if (n.flags.dim == this->mMesh->maxDim()) {
					rec_man_trace_down(n.cellindex, res);
				}

			}


			void fillGeometry(INT_TYPE nodeID, set<INDEX_TYPE>& res) {

				set<INT_TYPE> nodeset;
				GatherNodes(nodeID, nodeset);

				for (set<INT_TYPE>::iterator it = nodeset.begin(); it != nodeset.end(); it++) {
					node<DTYPE>& n = this->getNode(*it);

					if (n.flags.dim == 0) {
						rec_man_trace_up(n.cellindex, res);
					}
					else if (n.flags.dim == this->mMesh->maxDim()) {
						rec_man_trace_down(n.cellindex, res);
					}

				}


			}






			bool get_next_to_cancel(INT_TYPE& a) {
				///printf("getnext to cancel called\n");
				while (!edges_to_cancel.empty()) {
					sortedEdge se = edges_to_cancel.top();
					edges_to_cancel.pop();
					a = se.ep;
					arc<DTYPE>& ap = this->arcs[a];
					///printf("%u ", a);
					// is it alive in the current context
					if (!isAlive(ap, num_cancelled)) {
						//printf("adsf1\n");
						//printf("skip1->"); 
						///printArc(a);
						continue;
					}

					//test if it's a valid cancellation
					if (!isValid(a, ap)) {
						//printf("adsf2\n");
						continue;
					}

					int newcountweight = edgeCountWeight(ap);
					if (newcountweight > se.countweight) {
						//printf("adsf3\n");
						se.countweight = newcountweight;
						edges_to_cancel.push(se);
						continue;
					}

					if (newcountweight > 1500) {
						se.persistence += 1;
						if (se.persistence <= gPersThreshold)
							edges_to_cancel.push(se);
						continue;
					}

					//printf("getnext to cancel returned true\n");

					return true;


				}
				//printf("getnext to cancel returned false\n");
				return false;
			}

			vector<DTYPE> cancel_num_to_pers;
			INT_TYPE select_persistence; // the persistence to select in hierarchy
			INT_TYPE num_cancelled;

			DTYPE max_pers_so_far;
			DTYPE gPersThreshold;
			void ComputeHierarchy(DTYPE pers_limit) {
				cancel_num_to_pers.clear();
				printf(" Computing instance\n");
				gPersThreshold = pers_limit;
				max_pers_so_far = 0;
				num_cancelled = 0;
				// insert every arc to cancel list
				printf("  - Adding arcs to sorter...");
				INT_TYPE mysize = (INT_TYPE) this->arcs.size();
				for (INT_TYPE i = 0; i < mysize; i++) {
					// test if it passes the hierarchy test
					//if (!hierarchy_test->slice_position_list(i, this)) continue;
					arc<DTYPE> &a = this->arcs[i];
#ifdef HACK_REST_CANCEL
					if (a.persistence > gPersThreshold - 1.0) continue;
#else
					if (a.persistence > gPersThreshold) continue;
#endif
					sortedEdge se;
#ifdef HACK_REST_CANCEL
					ulong64 nid1 = this->nodes[a.lower].cellindex;
					ulong64 nid2 = this->nodes[a.upper].cellindex;
					bool same = this->sgg->bboundary[nid1] == this->sgg->bboundary[nid2];
					if (same) {
						a.persistence = a.persistence / maxp;
					}
					else {
						a.persistence = a.persistence + 1.0;
					}
#endif
					se.persistence = a.persistence;
					se.countweight = edgeCountWeight(a);
					se.ep = i;
					edges_to_cancel.push(se);
				}

				printf("Done!\n  - Cancelling:");

				INT_TYPE a;
				//float maxv = 0;

				while (get_next_to_cancel(a) && this->arcs[a].persistence <= gPersThreshold) {
					//if (arcs[a].persistence > maxv) maxv = arcs[a].persistence;
					if (num_cancelled % 1000 == 0) {
						printf("\r    - Cancelling: %u val=%f", num_cancelled, (float)max_pers_so_far);
					}
					//printf("\n\ncancelling: %d", num_cancelled);
					//printArc(a);

					cancel(a);
					cancel_num_to_pers.push_back(max_pers_so_far);

					//printf("Done\n");
				}
				printf("\r    - Cancelling finished: %u val=%f\n", num_cancelled, (float)max_pers_so_far);


				//_evercomputed = true;
				//if (output_cancellation_records) {
				//	int mincount = this->countIDs[0];
				//	int maxcount = this->countIDs[3];
				//	int interiormincount = this->countInteriorIDs[0];
				//	FILE* fout = fopen(this->mCRecordFName, "w");
				//	for (int i = 0; i < this->mCRecords.size(); i++) {
				//		cancellation_record& cr = this->mCRecords[i];
				//		fprintf(fout, "%f %f %f %f %d %d %d %d %d\n", cr.persistence,
				//			cr.lval, cr.uval, cr.persPerc, cr.index, mincount, maxcount,
				//			interiormincount, cr.boundary);
				//		if (cr.index == 0) {
				//			mincount--;
				//			node<DTYPE>& n =
				//				this->getNode(this->getArc(cr.arcp).lower);
				//			if (!n.flags.boundary) {
				//				interiormincount--;
				//			}
				//		}
				//		if (cr.index == 2) maxcount--;
				//	}
				//	fclose(fout);
				//}
				// now fill in the lists?	
			}
			void OutputCancelHistory(const char* fname) {

				FILE* fsadrec = fopen(fname, "w");
				for (int i = 0; i < cancel_num_to_pers.size(); i++) {
					fprintf(fsadrec, "%d %f\n", cancel_num_to_pers.size() - i, cancel_num_to_pers[i]);
				}
				fclose(fsadrec);

			}
			void Output1SaddleRecord(const char* fname) {

				FILE* fsadrec = fopen(fname, "w");

				for (int i = 0; i < nodes.size(); i++) {
					node<DTYPE>& n = nodes[i];
					if (n.flags.dim > 1) continue;
					float persval =
						(n.destroyed < cancel_num_to_pers.size() + 1 ?
							cancel_num_to_pers[n.destroyed - 1] :
							cancel_num_to_pers[cancel_num_to_pers.size() - 1]);
					fprintf(fsadrec, "%d %d %llu %f %d %f\n", n.flags.dim, n.flags.boundary, n.cellindex, n.value, n.destroyed, persval);


				}
				fclose(fsadrec);


			}
			};




		};

		template<typename SCALAR_TYPE, class MESH_TYPE>
		class NoFuncNoGradMSC {
		public:


			typedef SCALAR_TYPE ScalarType;
			typedef MESH_TYPE MeshType;

		protected:
			struct sortedEdge {
				SCALAR_TYPE persistence;
				int countweight;
				INT_TYPE ep;

				bool operator()(const sortedEdge& _Left, const sortedEdge& _Right) {
					if (_Left.persistence < _Right.persistence) return false;
					if (_Left.persistence > _Right.persistence) return true;
					if (_Left.countweight < _Right.countweight) return false;
					if (_Left.countweight > _Right.countweight) return true;
					return _Left.ep > _Right.ep;
				}
			};

			priority_queue<sortedEdge, vector< sortedEdge>, sortedEdge > edges_to_cancel;
		public:
			struct cancellation_record {
				int index;
				SCALAR_TYPE persistence;
				SCALAR_TYPE lval;
				SCALAR_TYPE uval;
				SCALAR_TYPE persPerc;
				INT_TYPE arcp;
				int boundary;
			};
		protected:

			vector<cancellation_record> mCRecords;
		public:
			const vector<cancellation_record>& GetCancellationRecords() {
				return mCRecords;
			}
		protected:

			bool output_cancellation_records;

			const char* mCRecordFName;
			float m_temp_perc_pers;
			//INDEX_TYPE select_persistence;
			INDEX_TYPE num_destroyed;
			MESH_TYPE* mMesh;
			Vec3b mBuildArcGeometry;
			Vec3b mBuildArcAtAll;

			vector<node<SCALAR_TYPE> > nodes;
			vector<arc<SCALAR_TYPE> > arcs;
			vector<arc_base_geometry> arc_base_geoms;
			vector<arc_merged_geometry> arc_merge_geoms;
			vector<merged_manifold> mans;

#ifdef WIN32
			std::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#else
#ifdef __APPLE__
			//__gnu_cxx::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
			std::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#else
			std::unordered_map<INDEX_TYPE, INT_TYPE> nodemap;
#endif
#endif

			// this is NOT thread safe - uses references to items in a vector that could be reallocated
			inline void connectArc(INT_TYPE arcID, arc<SCALAR_TYPE>& a, node<SCALAR_TYPE>& lower, node<SCALAR_TYPE>& upper) {
				a.lower_next = lower.firstarc;
				lower.firstarc = arcID;
				a.upper_next = upper.firstarc;
				upper.firstarc = arcID;

				lower.numarcs++;
				upper.numarcs++;
				upper.numlower++;

			}

			void connectArc(INT_TYPE arcID) {

				arc<SCALAR_TYPE>& a = arcs[arcID];
				node<SCALAR_TYPE>& lower = nodes[a.lower];
				node<SCALAR_TYPE>& upper = nodes[a.upper];
				connectArc(arcID, a, lower, upper);
			}
			const merged_manifold& getManifold(INT_TYPE m) const {
				return mans[m];
			}

		public:
			// add an arc to the complex, with geometry - vector ordered from upper node to lower node
			INT_TYPE AddArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID, vector<INDEX_TYPE>& geometry) {
				return createArc(nodemap[lowerCellID], nodemap[upperCellID], geometry);
			}
			INT_TYPE AddNode(INDEX_TYPE node_cell_id, SCALAR_TYPE value) {
				return createNode(node_cell_id, value);
			}

			INT_TYPE numArcs() const {
				return (INT_TYPE)arcs.size();
			}
			INT_TYPE numNodes() const {
				return (INT_TYPE)nodes.size();
			}
			const node<SCALAR_TYPE>& getNode(INT_TYPE e) const {
				return nodes[e];
			}
			node<SCALAR_TYPE>& getNode(INT_TYPE e) {
				return nodes[e];
			}

			const arc<SCALAR_TYPE>& getArc(INT_TYPE e) const {
				return arcs[e];
			}
			arc<SCALAR_TYPE>& getArc(INT_TYPE e) {
				return arcs[e];
			}
			MESH_TYPE* const GetMesh() const {
				return mMesh;
			}

			void set_output_cancellation_records(const char* fname) {
				output_cancellation_records = true;
				mCRecordFName = fname;
			}


			NoFuncNoGradMSC(
				MESH_TYPE* mesh) :
				 mMesh(mesh),
				mBuildArcAtAll(Vec3b(true, true, true)), mBuildArcGeometry(Vec3b(true, true, true))
			{
				num_destroyed = 0;
				select_persistence = 0;
				output_cancellation_records = false;
				for (int i = 0; i < 4; i++) {
					countIDs[i] = 0;
					countBoundaryIds[i] = 0;
				}
			}

			// must be calleb before computefromgrad
			void SetBuildArcs(Vec3b v) { mBuildArcAtAll = v; }
			void SetBuildArcGeometry(Vec3b v) { mBuildArcGeometry = v; }

		protected:

			INT_TYPE createManifold(INT_TYPE baseId, INT_TYPE baseManId, INT_TYPE mergeManId, INT_TYPE mergetime) {
				INT_TYPE manid = mans.size();
				merged_manifold m;
				m.basenode = baseId;
				m.merged[0] = baseManId;
				m.merged[1] = mergeManId;
				m.mergetime = mergetime;
				mans.push_back(m);
				return manid;
			}

			INT_TYPE createNode(INDEX_TYPE cellID, SCALAR_TYPE value) {
				node<SCALAR_TYPE> tn;
				INT_TYPE listID = nodes.size();
				//if (nodes.size() == nodes.capacity()) printf("about to expand capacity - nodes createnode!\b");
				nodes.push_back(tn);
				nodemap[cellID] = listID;

				//node<SCALAR_TYPE> &n = nodes->operator [](listID);
				node<SCALAR_TYPE> &n = nodes[listID];
				n.cellindex = cellID;

				n.destroyed = INT_INFTY;
				n.firstarc = NULLID;
				n.boundary = mMesh->boundaryValue(cellID);
				//printf("n.boundary = %d\n", n.boundary);
				//if (n.boundary) {
				//	INDEX_TYPE coords[3];
				//	sgg->getCoords(coords, cellID);
				//	printf("(%d, %d, %d) = %d\n",
				//	   (int) coords[0], (int) coords[1], (int) coords[2], n.boundary);
				//}
				n.dim = mMesh->dimension(cellID);
				this->countIDs[n.dim]++;
				if (n.boundary) this->countBoundaryIds[n.dim]++;
				n.numarcs = 0;
				n.numlower = 0;
				n.value = value;

				n.amanifoldid = createManifold(listID, -1, -1, 0);
				n.dmanifoldid = createManifold(listID, -1, -1, 0);

				return listID;
			}

			//// create an arc connecting lower to upper <- cell id's from gradient
			//INT_TYPE createArc(arc<SCALAR_TYPE> &a){
			//	INT_TYPE listID = (INT_TYPE)arcs.push_back(a);
			//	//printf("%f - %f\n", (float) a.persistence, arcs[listID].persistence);
			//	return listID;
			//}

			INT_TYPE createArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID, std::vector<INDEX_TYPE>& geometry) {

				arc<SCALAR_TYPE> at;
				INT_TYPE listID = arcs.size();
				arcs.push_back(at);
				INT_TYPE geomID = this->arc_base_geoms.size();
				arc_base_geoms.push_back(arc_base_geometry());

				arc<SCALAR_TYPE>& a = arcs[listID];
				arc_base_geometry& ge = arc_base_geoms[geomID];
				ge.geometry.insert(ge.geometry.begin(), geometry.begin(), geometry.end());



				a.created = 0;
				a.destroyed = INT_INFTY;
				a.lower = nodemap[lowerCellID];
				a.upper = nodemap[upperCellID];
				a.geom = geomID;
				if (nodes[a.upper].value < nodes[a.lower].value)
					printf("ERROR: upper (%f, %d) - lower (%f, %d)\n",
						float(nodes[a.upper].value), nodes[a.upper].dim,
						float(nodes[a.lower].value), nodes[a.lower].dim);

				a.persistence = nodes[a.upper].value - nodes[a.lower].value;
				a.dim = nodes[a.lower].dim;
				connectArc(listID);

				return listID;
			}

			INT_TYPE createArc(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID) {
				std::vector<INDEX_TYPE> a;
				return createArc(lowerCellID, upperCellID, a);
			}
			virtual void InsertArcIntoSimplification(arc<SCALAR_TYPE>& a, INT_TYPE id) {

				if (a.persistence <= gPersThreshold) {
					sortedEdge se;

					se.persistence = a.persistence;
					se.countweight = edgeCountWeight(a);
					se.ep = id;
					edges_to_cancel.push(se);
				}
				//if (arcs.size() == arcs.capacity()) {


			}


			INT_TYPE createArc(INT_TYPE luap, INT_TYPE ma, INT_TYPE ulap, INT_TYPE ctime) {
				arc<SCALAR_TYPE>& lua = arcs[luap];
				arc<SCALAR_TYPE>& ula = arcs[ulap];
				INT_TYPE na_id = arcs.size();
				INT_TYPE na_geom_id = arc_merge_geoms.size();
				arc_merge_geoms.push_back(arc_merged_geometry());
				arc_merged_geometry& na_geom = arc_merge_geoms[na_geom_id];
				na_geom.fields[0] = luap;
				na_geom.fields[1] = ma;
				na_geom.fields[2] = ulap;

				//this->arcs.push_back(ma); // copy setting from old
				arc<SCALAR_TYPE> na;// = this->arcs[na_id];

				na.created = ctime;
				na.destroyed = INT_INFTY;
				na.lower = ula.lower;
				na.upper = lua.upper;
				na.geom = na_geom_id;

				node<SCALAR_TYPE>& nup = this->nodes[na.upper];
				node<SCALAR_TYPE>& nlo = this->nodes[na.lower];

				this->connectArc(na_id, na, nlo, nup);
				na.persistence = nup.value - nlo.value;
				na.boundary = nlo.boundary + nup.boundary;
				na.dim = nlo.dim;
				if (na.persistence < 0 && nlo.boundary == nup.boundary) printf("creatinga inversion %f, %d-%d, b%d-b%d\n",
					(float)na.persistence, nlo.dim, nup.dim, nlo.boundary, nup.boundary);

				// now insert into global sort if it has a chance of being cancelled
				InsertArcIntoSimplification(na, na_id);
				//	printf("about to expand capacity2 - arcs createnode!\n");
				//}
				this->arcs.push_back(na); // this could destroy references so put at end;
				return na_id;

			}

				// old way is to trace "down" so the first of geom vector
			// is the upper dim critical poitn

		public:



		protected:

			//////////// NOW CANCELLATION STUFF

			int countUpperArcs(INT_TYPE id) const {
				return this->nodes[id].numarcs - this->nodes[id].numlower;
			}
			int countLowerArcs(INT_TYPE id) const {
				return this->nodes[id].numlower;
			}

			int edgeCountWeight(arc<SCALAR_TYPE> &a) const {
				INT_TYPE lower = a.lower;
				INT_TYPE upper = a.upper;
				int nl = countUpperArcs(lower);
				int nu = countLowerArcs(upper);

				////first option: return the number created
				return (nl - 1) * (nu - 1);

				//2nd option: return the change in the number of arcs
				return (nl - 1) * (nu - 1) - (nl + nu - 1);

				//return 1;
			}

			inline INT_TYPE current_pers() const { return select_persistence; }
			inline bool isAlive(const arc<SCALAR_TYPE>&a, INT_TYPE place) const {
				return a.created <= place && a.destroyed > place;
			}
		public:
			bool isNodeAlive(INT_TYPE n) const {
				return isAlive(this->nodes[n], select_persistence);
			}
			bool isArcAlive(INT_TYPE a) const {
				bool res = isAlive(this->arcs[a], select_persistence);
				if (res && !isNodeAlive(this->arcs[a].lower)) printf("ERROR LOWER NODE NOT ALIVE\n");
				if (res && !isNodeAlive(this->arcs[a].upper)) printf("ERROR UPPER NODE NOT ALIVE\n");
				return res;
			}

		protected:
			inline bool isAlive(const node<SCALAR_TYPE>&n, INT_TYPE place) const {
				//printf("n.destroyed = %d, place = %d\n", n.destroyed, place);
				return  n.destroyed > place;
			}
			inline bool isAlive(INT_TYPE n, INT_TYPE place) const {
				return isAlive(this->nodes[n], place);
			}


		public:
			SCALAR_TYPE GetSelectPersAbs() {
				return select_persistence;
			}
			void SetSelectPersAbs(SCALAR_TYPE value) {
				//printf("mcLightGeomMSC::SuperLightMSC::SetSelectPersAbs -> %f\n", value);
				int offset = 1;// cancel_num_to_pers.size() - 1;
				for (int i = 0; i < cancel_num_to_pers.size(); i += offset) {
					if (cancel_num_to_pers[i] > value) {
						select_persistence = i;

						return;
					}
					//while (i + offset > cancel_num_to_pers.size() - 1) offset = offset / 2;
					//while (offset > 1 && cancel_num_to_pers[i + offset] > value) offset = offset / 2;
				}
				select_persistence = cancel_num_to_pers.size() - 1;
			}

		protected:
			int cancel(INT_TYPE a) {

				int createcounter = 0;
				int deletecounter = 0;

				//int tmp1 = (int) edges_to_cancel.size();

				arc<SCALAR_TYPE>* ap = &(this->arcs[a]);
				node<SCALAR_TYPE>& lower = this->nodes[ap->lower];
				node<SCALAR_TYPE>& upper = this->nodes[ap->upper];

				//if (output_cancellation_records) {
				cancellation_record cr;
				cr.index = lower.dim;
				cr.lval = lower.value;
				cr.uval = upper.value;
				cr.persistence = ap->persistence;
				//SCALAR_TYPE diff = this->sgg->sgd->maxval - this->sgg->sgd->minval;
				//cr.persPerc = 100.0f * ap->persistence / diff;
				cr.arcp = a;
				cr.boundary = lower.boundary + upper.boundary;
				this->mCRecords.push_back(cr);
				//}



				if (lower.destroyed != INT_INFTY) printf("Error: NoFuncNoGradMSC::cancel(%d) - lower.destroyed != INT_INFTY\n", a);
				if (upper.destroyed != INT_INFTY) printf("Error: NoFuncNoGradMSC::cancel(%d) - upper.destroyed != INT_INFTY\n", a);
				//printf("\n");
				//printNode(ap->lower);
				//printNode(ap->upper);

				int initialguess = (lower.numarcs - lower.numlower - 1) * (upper.numlower - 1);
				int init2 = lower.numarcs + upper.numarcs - 1;

				if (ap->persistence > max_pers_so_far) max_pers_so_far = ap->persistence;



				// for each upwards arc connected to the lower node,
				// for each downwards arc connected to the upper node,
				// create a new connection from the other endpoints
				INT_TYPE la = lower.firstarc;
				while (la != NULLID) {
					arc<SCALAR_TYPE>& lap = this->arcs[la];

					// skip the arc itself
					if (la == a) {
						la = lap.lower_next; // we're guarantee this is the next
						continue;
					}

					//test if they are the same kind
					if (ap->lower != lap.lower) {
						la = lap.upper_next;
						continue;
					}
					if (!isAlive(lap, num_cancelled)) {
						la = lap.lower_next;
						continue;
					}

					INT_TYPE ua = upper.firstarc;
					while (ua != NULLID) {
						arc<SCALAR_TYPE>& uap = this->arcs[ua];

						// skip the arc itself
						if (ua == a) {
							ua = uap.upper_next; // we're guarantee this is the next
							continue;
						}

						//test if they are the same kind
						if (ap->upper != uap.upper) {
							ua = uap.lower_next;
							continue;
						}
						if (!isAlive(uap, num_cancelled)) {
							ua = uap.upper_next;
							continue;
						}

						// create the arc here!
						INT_TYPE newarc = createArc(la, a, ua, num_cancelled + 1);
						//printf("ha1\n");
						ap = &(this->arcs[a]);
						//printf("ha\n");
						createcounter++;

						ua = this->arcs[ua].upper_next;
						//printf("asdf\n");
					}

					la = this->arcs[la].lower_next;

				}

				//printf("hahan\n");
				// go through and set the delete time on the arcs that are to be removed,
				// and update the arc counters at the other endpoints

				// following comments are for 0-1 cancellations
				// first look in neighborhood of minimum = lower
				la = lower.firstarc;							// pick first arc attached to min
				while (la != NULLID) {							// while there are more arcs in neighborhood, keep looking
					arc<SCALAR_TYPE>& lap = this->arcs[la];			// get the reference to the arc attached to the min

																	// skip the arc itself
					if (la == a) {
						la = lap.lower_next; // we're guarantee this is the next
						continue;
					}

					//test if they are the same kind
					if (ap->lower == lap.lower) {				// we do not want arcs that are not 0-1, so do this check
						if (isAlive(lap, num_cancelled)) {		// make sure we are only working with living arcs
							node<SCALAR_TYPE>& n = this->nodes[lap.upper]; // now "n" is the 1-saddle attached to minimum

							n.numarcs--;						// if we remove this arc, we need to decrement the number of living arcs attached to it
							n.numlower--;						// this cancellation will also remove a "lower" arc of the 1-saddle

																// create merged_manifold for merging
																// extend the descending merged_manifold of the 1-saddle by merging it with the descending merged_manifold of the cancelled 1-saddle
							INT_TYPE nmanid =
								createManifold(lap.upper, n.dmanifoldid, upper.dmanifoldid, num_cancelled + 1);
							n.dmanifoldid = nmanid;				// change the merged_manifold reference of the 1-saddle to reference this new decending merged_manifold

							lap.destroyed = num_cancelled + 1;		// we remove the old arc to the deleted min, so mark it destroyed with the num_cancelld
							deletecounter++;					// keep track of deleted arc cound
						}
						la = lap.lower_next;					// continue on the next arc around the minimum
					}
					else if (ap->lower == lap.upper) {          // if this is in fact an i-1 - i arc, we just remove the arcs around it - no merging happens
						if (isAlive(lap, num_cancelled)) {
							node<SCALAR_TYPE>& n = this->nodes[lap.lower];
							n.numarcs--;
							lap.destroyed = num_cancelled + 1;
							deletecounter++;
						}
						la = lap.upper_next;
					}
					else {
						printf("ERROR SHOULD NEVER GET HERE\n");
					}



				}
				INT_TYPE ua = upper.firstarc;
				while (ua != NULLID) {
					arc<SCALAR_TYPE>& uap = this->arcs[ua];

					// skip the arc itself
					if (ua == a) {
						ua = uap.upper_next; // we're guarantee this is the next
						continue;
					}
					//test if they are the same kind
					if (ap->upper == uap.upper) {
						if (isAlive(uap, num_cancelled)) {

							node<SCALAR_TYPE>& n = this->nodes[uap.lower]; // pick the other

							n.numarcs--;

							// create merged_manifold for merging
							INT_TYPE nmanid = createManifold(uap.lower, n.amanifoldid, lower.amanifoldid, num_cancelled + 1);
							n.amanifoldid = nmanid;

							uap.destroyed = num_cancelled + 1;
							deletecounter++;
						}
						ua = uap.upper_next;
					}
					else if (ap->upper == uap.lower) {
						if (isAlive(uap, num_cancelled)) {
							node<SCALAR_TYPE>& n = this->nodes[uap.upper];
							n.numarcs--;
							n.numlower--;
							uap.destroyed = num_cancelled + 1;
							deletecounter++;
						}
						ua = uap.lower_next;
					}
					else {
						printf("ERROR SHOULD NEVER GET HERE\n");
					}


				}
				num_cancelled++;


				ap->destroyed = num_cancelled;
				lower.destroyed = num_cancelled;
				upper.destroyed = num_cancelled;
				deletecounter++;

				//int tmp2 = (int) edges_to_cancel.size();

				//printf("initialguess: %d actual:%d del:%d td:%d tmp1:%d tmp2:%d \n", initialguess, createcounter, deletecounter,init2
				//	,tmp1, tmp2);
				//printf("gothere\n");
				return 1;
			}



			INT_TYPE nextArc(const arc<SCALAR_TYPE>& ap, INT_TYPE n) const {
				if (ap.lower == n) return ap.lower_next;
				if (ap.upper == n) return ap.upper_next;
				printf("ERROR BARF POOP\n");
				return NULLID;
			}
			int countMultiplicity(arc<SCALAR_TYPE>& ap, INT_TYPE ctime) const {
				INT_TYPE nu = ap.upper;
				INT_TYPE nl = ap.lower;
				INT_TYPE a = this->nodes[nu].firstarc;
				int counter = 0;
				while (a != NULLID) {
					const arc<SCALAR_TYPE>& nap = this->arcs[a];
					if (isAlive(nap, ctime) && nap.lower == nl && nap.upper == nu) {
						counter++;
					}
					a = nextArc(nap, nu);
				}
				return counter;
			}

			virtual bool isValid(INT_TYPE a, arc<SCALAR_TYPE>& ap) const {
				// test for boundary
				if (this->nodes[ap.lower].boundary !=
					this->nodes[ap.upper].boundary) return false;

				// 2 endpoints must be connected by exactly one arc
				if (countMultiplicity(ap, num_cancelled) != 1) return false;
				// test for inversions?

				return true;


			}

			// THIS WILL HAVE TO WORK LATER, but for now only need 0 and maxdim manifolds
			//void recCounLeaftManifolds(INT_TYPE mId, map< INT_TYPE, int >& counter) {
			//	merged_manifold &m = this->getManifold(mId);
			//	if (m->merge[0] != NULL) {
			//		recCounLeaftManifolds(m->merge[0], counter);
			//		recCounLeaftManifolds(m->merge[1], counter);
			//	}
			//	else {
			//		if (counter.count(m) == 0) {
			//			counter[m] = 1;
			//		}
			//		else {
			//			counter[m]++;
			//		}
			//	}
			//}

			void printmanifold(INT_TYPE man) const {
				merged_manifold& m = getManifold(man);
				printf("man=%d, man.base=%d, man.mergetime=%d, man.merge[0]=%d, man.merge[1]=%d\n",
					man, m.basenode, m.mergetime, m.merged[0], m.merged[1]);
			}

			INT_TYPE getActiveMan(INT_TYPE man) const {
				//printmanifold(man);
				while (getManifold(man).mergetime > this->current_pers()) {
					man = getManifold(man).merged[0];
					//printf("  --");  printmanifold(man);
				}
				//printf("ret %d\n", man);
				return man;
			}



		protected:
			void recGatherNodes(INT_TYPE mid, set<INT_TYPE>& res)  const {
				vector<INT_TYPE> expand;
				expand.push_back(mid);
				while (expand.size() > 0) {

					INT_TYPE curr = expand.back(); expand.pop_back();

					const merged_manifold& m = getManifold(curr);
					if (m.merged[0] != -1) {
						expand.push_back(m.merged[0]);
						expand.push_back(m.merged[1]);
					}
					res.insert(m.basenode);


				}


			}



		public:


			void GatherNodes(INT_TYPE nodeID, set<INT_TYPE>& res, bool ascending) const {

				if (!isAlive(nodeID, this->select_persistence)) {
					return;
				}
				const node<SCALAR_TYPE>& n = this->getNode(nodeID);
				INT_TYPE man;
				if (ascending) {
					man = n.amanifoldid;
				}
				else {
					man = n.dmanifoldid;
				}
				INT_TYPE manID = getActiveMan(man);
				recGatherNodes(manID, res);

			}
		private:
			void fillGeometry(INT_TYPE nodeID, set<INDEX_TYPE>& res, bool ascending) const {

				if (!isNodeAlive(nodeID)) return;

				set<INT_TYPE> nodeset;
				//printf("gothere1\n");
				GatherNodes(nodeID, nodeset, ascending);
				//printf("gothere2\n");

				for (set<INT_TYPE>::iterator it = nodeset.begin(); it != nodeset.end(); it++) {
					const node<SCALAR_TYPE>& n = this->getNode(*it);

					// NEED TO ADD MANIFOLD GEOM HERE

					//if (ascending) {
					//	rec_man_trace_up(n.cellindex, res);
					//}
					//else {
					//	rec_man_trace_down(n.cellindex, res);
					//}

				}
				//printf("gothere3\n");



			}


			void fillBaseGeometry(INT_TYPE nodeID, set<INDEX_TYPE>& res, bool ascending) const {
				const node<SCALAR_TYPE>& n = this->getNode(nodeID);
				
				// NEED TO ADD MANIFOLD GEOM HERE
				//if (ascending) {
				//	rec_man_trace_up(n.cellindex, res);
				//}
				//else {
				//	rec_man_trace_down(n.cellindex, res);
				//}
			}

		protected:

			bool get_next_to_cancel(INT_TYPE& a) {
				///printf("getnext to cancel called\n");
				while (!edges_to_cancel.empty()) {
					sortedEdge se = edges_to_cancel.top();
					edges_to_cancel.pop();
					a = se.ep;
					arc<SCALAR_TYPE>& ap = this->arcs[a];
					///printf("%u ", a);
					// is it alive in the current context
					if (!isAlive(ap, num_cancelled)) {
						//printf("adsf1\n");
						//printf("skip1->");
						///printArc(a);
						continue;
					}

					//test if it's a valid cancellation
					if (!isValid(a, ap)) {
						//printf("adsf2\n");
						continue;
					}

					int newcountweight = edgeCountWeight(ap);
					if (newcountweight > se.countweight) {
						//printf("adsf3\n");
						se.countweight = newcountweight;
						edges_to_cancel.push(se);
						continue;
					}

					if (newcountweight > 1500) {
						se.persistence += 1;
						if (se.persistence <= gPersThreshold)
							edges_to_cancel.push(se);
						continue;
					}

					//printf("getnext to cancel returned true\n");

					return true;


				}
				//printf("getnext to cancel returned false\n");
				return false;
			}

			vector<SCALAR_TYPE> cancel_num_to_pers;
			INT_TYPE select_persistence; // the persistence to select in hierarchy
			INT_TYPE num_cancelled;

			SCALAR_TYPE max_pers_so_far;
			SCALAR_TYPE gPersThreshold;

			int countIDs[4];
			int countBoundaryIds[4];
		public:

			void SetPersistanceByNumOfMaxima(int num_maxima) {
				int maxcount = this->countIDs[3];

				for (int i = 0; i < this->mCRecords.size(); i++) {
					cancellation_record& cr = this->mCRecords[i];

					if (maxcount < num_maxima) {
						select_persistence = i - 1;
						return;
					}

					if (cr.index == 2) maxcount--;

				}
			}

			virtual void ComputeHierarchy(SCALAR_TYPE pers_limit) {
				cancel_num_to_pers.clear();
				printf(" -- Performing cancellation to %f...\n", pers_limit);
				gPersThreshold = pers_limit;
				max_pers_so_far = 0;
				num_cancelled = 0;
				// insert every arc to cancel list
				printf("  -- Adding arcs to sorter...");
				INT_TYPE mysize = (INT_TYPE) this->arcs.size();
				for (INT_TYPE i = 0; i < mysize; i++) {
					// test if it passes the hierarchy test
					//if (!hierarchy_test->slice_position_list(i, this)) continue;
					arc<SCALAR_TYPE> &a = this->arcs[i];

					if (a.persistence > gPersThreshold) continue;

					sortedEdge se;

					se.persistence = a.persistence;
					se.countweight = edgeCountWeight(a);
					se.ep = i;
					edges_to_cancel.push(se);
				}

				printf("Done!\n");
				printf("  -- Cancelling:");

				INT_TYPE a;
				//float maxv = 0;

				while (get_next_to_cancel(a) && this->arcs[a].persistence <= gPersThreshold) {
					//if (arcs[a].persistence > maxv) maxv = arcs[a].persistence;
					if (num_cancelled % 1000 == 0) {
						printf("\r  -- Cancelling: %u val=%f", num_cancelled, (float)max_pers_so_far);
					}
					//printf("\n\ncancelling: %d", num_cancelled);
					//printArc(a);

					cancel(a);
					cancel_num_to_pers.push_back(max_pers_so_far);

					//printf("Done\n");
				}
				printf("\r  -- Cancelling finished: %u val=%f\n", num_cancelled, (float)max_pers_so_far);


				//_evercomputed = true;
				if (output_cancellation_records) {
					int mincount = this->countIDs[0];
					int sad1count = this->countIDs[1];
					int sad2count = this->countIDs[2];
					int maxcount = this->countIDs[3];

					int bcount[4];
					for (int i = 0; i < 4; i++) bcount[i] = this->countBoundaryIds[i];

					FILE* fout = fopen(this->mCRecordFName, "w");
					//fprintf(fout, "pers_so_far actual_pers lval uval perc_pers lindex #min #1s #2s #max is_bdry\n", max_pers_so_far, cr.persistence,

					float max_pers_so_far = this->mCRecords[0].persistence;
					for (int i = 0; i < this->mCRecords.size(); i++) {
						cancellation_record& cr = this->mCRecords[i];
						max_pers_so_far = fmax(max_pers_so_far, cr.persistence);
						fprintf(fout, "%f %f %f %f %f %d %d %d %d %d %d %d %d %d %d\n", max_pers_so_far, cr.persistence,
							cr.lval, cr.uval, cr.persPerc, cr.index, mincount, sad1count, sad2count, maxcount,
							cr.boundary, bcount[0], bcount[1], bcount[2], bcount[3]);
						if (cr.index == 0) {
							mincount--; sad1count--;
						}
						else if (cr.index == 1) {
							sad1count--; sad2count--;
						}
						else if (cr.index == 2) {
							maxcount--; sad2count--;
						}
						if (cr.boundary) {
							bcount[cr.index]--;
							bcount[cr.index + 1]--;
						}
					}
					fclose(fout);
				}
				// now fill in the lists?
			}
			void Output1SaddleRecord(const char* fname) const {

				FILE* fsadrec = fopen(fname, "w");

				for (int i = 0; i < nodes.size(); i++) {
					node<SCALAR_TYPE>& n = nodes[i];
					if (n.dim > 1) continue;
					float persval =
						(n.destroyed < cancel_num_to_pers.size() + 1 ?
							cancel_num_to_pers[n.destroyed - 1] :
							cancel_num_to_pers[cancel_num_to_pers.size() - 1]);
					fprintf(fsadrec, "%d %d %llu %f %d %f\n", n.dim, n.boundary, n.cellindex, n.value, n.destroyed, persval);


				}
				fclose(fsadrec);


			}
		protected:
			void _fillArcGeometry(INT_TYPE aid, vector<INDEX_TYPE>& v, bool direction) const {
				const arc<SCALAR_TYPE>& a = this->arcs[aid];
				//printf("gothere1\n");
				if (a.created == 0) {
					// this is base
					const arc_base_geometry& base = this->arc_base_geoms[a.geom];
					if (direction) {
						for (int i = 0; i < base.geometry.size(); i++) {
							if (v.size() > 0 && v[v.size() - 1] == base.geometry[i]) {
								if (v.size() > 1 && i < base.geometry.size() - 1 && v[v.size() - 2] == base.geometry[i + 1]) v.pop_back();
							}
							else {
								v.push_back(base.geometry[i]);
							}
						}
					}
					else {
						for (int i = base.geometry.size() - 1; i >= 0; i--) {
							if (v.size() > 0 && v[v.size() - 1] == base.geometry[i]) {
								if (v.size() > 1 && i >0 && v[v.size() - 2] == base.geometry[i - 1])v.pop_back();
							}
							else {
								v.push_back(base.geometry[i]);
							}
						}
					}
					return;
				}
				//printf("gothere2\n");
				// recurse on children
				const arc_merged_geometry& m = arc_merge_geoms[a.geom];
				if (direction) {
					_fillArcGeometry(m.fields[0], v, true);
					_fillArcGeometry(m.fields[1], v, false);
					_fillArcGeometry(m.fields[2], v, true);
				}
				else {
					_fillArcGeometry(m.fields[2], v, false);
					_fillArcGeometry(m.fields[1], v, true);
					_fillArcGeometry(m.fields[0], v, false);
				}


			}
		public:
			INT_TYPE arcLowerNode(INT_TYPE aid) const {
				return getArc(aid).lower;
			}
			INT_TYPE arcUpperNode(INT_TYPE aid) const {
				return getArc(aid).upper;
			}
			INT_TYPE nextIncidentLivingArc(INT_TYPE aid, INT_TYPE nid) const {
				//printf("%d-\n", aid);
				INT_TYPE naid = nextArc(getArc(aid), nid);
				//printf("  -%d\n", naid);
				while (naid != NULLID) {
					if (isArcAlive(naid)) return naid;
					naid = nextArc(getArc(naid), nid);
					//printf("  -%d\n", naid);
				}
				return naid;
			}
			bool isValidArcId(INT_TYPE aid) const {
				return aid != NULLID;
			}

			INT_TYPE firstIncidentLivingArc(INT_TYPE nid) const {
				node<SCALAR_TYPE>& n = getNode(nid);
				INT_TYPE aid = n.firstarc;
				if (!isArcAlive(aid)) return nextIncidentLivingArc(aid, nid);
				return aid;
			}

			void fillArcGeometry(INT_TYPE aid, vector<INDEX_TYPE>& v) const {
				//v.push_back(nodes[arcs[aid].upper].cellindex);
				v.clear();
				_fillArcGeometry(aid, v, true);
				//v.push_back(nodes[arcs[aid].lower].cellindex);
			}

			class SurroundingArcsIterator {
			protected:
				NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE>* mMSC;
				INT_TYPE mNID;
				INT_TYPE currarc;
				INT_TYPE next_arc(INT_TYPE arcid) {
					return mMSC->nextArc(mMSC->getArc(arcid), mNID);
				}

			public:
				SurroundingArcsIterator(NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE>* msc) :
					mMSC(msc) {}

				void begin(INT_TYPE nid) {
					//printf("begin %d...\n");
					mNID = nid;
					node<SCALAR_TYPE>& n = mMSC->getNode(mNID);
					currarc = n.firstarc;
				}
				bool valid() {
					//printf("valid? %d %d\n", currarc, NULLID);
					return currarc != NULLID;
				}
				INT_TYPE value() {
					return currarc;
				}
				void advance() {
					//printf("Advance? %d\n", currarc);
					if (currarc == NULLID) return;
					//printf("...\n");
					currarc = next_arc(currarc);
					//printf("advance = %d\n", currarc);
				}
			};

			class SurroundingLivingArcsIterator : public SurroundingArcsIterator {
			protected:
				bool advance_until_alive() {
					this->currarc = this->next_arc(this->currarc);

					if (this->currarc == NULLID) return false;
					while (!this->mMSC->isArcAlive(this->currarc)) {
						this->currarc = this->next_arc(this->currarc);
						if (this->currarc == NULLID) return false;
					}
					return true;
				}
			public:
				SurroundingLivingArcsIterator(NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE>* msc) :
					SurroundingArcsIterator(msc) {}

				void begin(INT_TYPE nid) {
					//printf("begin %d...\n");
					SurroundingArcsIterator::begin(nid);
					if (!this->valid())
						printf("ERROR: node has no arcs!\n");
					if (!this->mMSC->isArcAlive(this->currarc)) advance_until_alive();
					//printf("begin return %d\n", currarc);
				}

				void advance() {
					//printf("Advance? %d\n", currarc);
					if (this->currarc == NULLID) return;
					//printf("...\n");
					advance_until_alive();
					//printf("advance = %d\n", currarc);
				}
			};

			class NodesIterator {
			protected:
				NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE>* mMSC;
				INT_TYPE currid;

			public:
				NodesIterator(NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE>* msc) :
					mMSC(msc) {}

				void begin() {
					currid = 0;
				}
				bool valid() {
					return currid < mMSC->nodes.size();
				}
				INT_TYPE value() {
					return currid;
				}
				void advance() {
					currid++;
				}
			};
			class LivingNodesIterator : public NodesIterator {
			protected:

				void advance_until_alive() {
					this->currid++;
					while (this->valid() && !this->mMSC->isNodeAlive(this->currid)) {
						this->currid++;
					}
				}
			public:
				LivingNodesIterator(NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE>* msc) :
					NodesIterator(msc) {}

				void begin() {
					NodesIterator::begin();
					if (this->valid() && !this->mMSC->isNodeAlive(this->currid)) advance_until_alive();
				}
				void advance() {
					advance_until_alive();
				}
			};

			// replicates functionality just uses different end criteria
			class ArcsIterator {
			protected:
				NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE>** mMSC;
				INT_TYPE currid;

			public:
				ArcsIterator(NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE>* msc) :
					mMSC(msc) {}

				void begin() {
					currid = 0;
				}
				bool valid() {
					return currid < mMSC->arcs.size();
				}
				INT_TYPE value() {
					return currid;
				}
				void advance() {
					currid++;
				}
			};
			class LivingArcsIterator : public ArcsIterator {
			protected:
				void advance_until_alive() {
					this->currid++;
					while (this->valid() && !this->mMSC->isArcAlive(this->currid)) {
						this->currid++;
					}
				}
			public:
				LivingArcsIterator(NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE>* msc) :
					ArcsIterator(msc) {}

				void begin() {
					ArcsIterator::begin();
					if ((this->valid() && !this->mMSC->isArcAlive(this->currid))) advance_until_alive();
				}
				void advance() {
					advance_until_alive();
				}
			};


			void print_complex_info(bool verbose) {
				LivingNodesIterator nit(this);
				int numnodes = 0;
				int dims[4]; for (int i = 0; i < 4; i++) dims[i] = 0;
				for (nit.begin(); nit.valid(); nit.advance()) {
					numnodes++;
					INT_TYPE nid = nit.value();
					node<SCALAR_TYPE>& n = this->getNode(nid);
					if (verbose) {
						printf("node: %lld, %d, %f\n", n.cellindex, n.dim, n.value);
					}
					dims[n.dim]++;
				}
				int numarcs = 0;
				int dims2[4]; for (int i = 0; i < 4; i++) dims2[i] = 0;
				for (int i = 0; i < 4; i++) dims2[i] = 0;
				LivingArcsIterator ait(this);
				for (ait.begin(); ait.valid(); ait.advance()) {
					INT_TYPE aid = ait.value();
					numarcs++;
					arc<SCALAR_TYPE>& a = this->getArc(aid);
					dims2[a.dim]++;
					if (verbose) {
						printf("arc: %lld--%lld %d\n", this->getNode(a.lower).cellindex, this->getNode(a.upper).cellindex, a.dim);
					}
				}
				printf("msc: %d nodes [%d, %d, %d, %d] ", numnodes, dims[0], dims[1], dims[2], dims[3]);

				printf("%d arcs [%d %d %d]\n", numarcs, dims2[0], dims2[1], dims2[2]);
			};


			struct write_node {
				INDEX_TYPE cellindex;
				INT_TYPE amanifoldid; // set to -1 if this is the base
				INT_TYPE dmanifoldid; // set to -1 if this is the base
				SCALAR_TYPE value; // KEEP this?

			};
			struct write_arc {
				INT_TYPE lower; // node
				INT_TYPE upper; // node
				INT_TYPE geom_size; // if created == 0, this is an original arc, and its geometry will be found in base_geom list of MSC, else in the merged_geom list
			};

			struct write_state {
				vector<write_node> nodes;
				vector<write_arc> arcs;
				vector<INDEX_TYPE> geom;

				void write_to_file(const char* name) {

					FILE* fout = fopen(name, "wb");

					int numnodes = nodes.size();
					int numarcs = arcs.size();
					int numgeom = geom.size();
					fwrite(&numnodes, sizeof(int), 1, fout);
					fwrite(&numarcs, sizeof(int), 1, fout);
					fwrite(&numgeom, sizeof(int), 1, fout);
					//printf("%d %d %d\n", numnodes, numarcs, numgeom);

					fwrite(&nodes[0], sizeof(write_node), numnodes, fout);
					fwrite(&arcs[0], sizeof(write_arc), numarcs, fout);
					fwrite(&geom[0], sizeof(INDEX_TYPE), numgeom, fout);
					fclose(fout);

				}

				void read_from_file(const char* name) {

					FILE* fin = fopen(name, "rb");

					int numnodes;
					fread(&numnodes, sizeof(int), 1, fin);
					nodes.resize(numnodes);
					int numarcs;
					fread(&numarcs, sizeof(int), 1, fin);
					arcs.resize(numarcs);
					int numgeom;
					fread(&numgeom, sizeof(int), 1, fin);
					geom.resize(numgeom);
					//printf("%d %d %d\n", numnodes, numarcs, numgeom);

					fread(&nodes[0], sizeof(write_node), numnodes, fin);
					//for (auto wn : nodes) printf("%lld\n", wn.cellindex);
					fread(&arcs[0], sizeof(write_arc), numarcs, fin);
					fread(&geom[0], sizeof(INDEX_TYPE), numgeom, fin);
					fclose(fin);

				}
			};

			void WriteComplex1Skeleton(const char* filename) {

				map<INT_TYPE, INT_TYPE> living_nodes;

				write_state ws;

				LivingNodesIterator nit(this);
				int pos = 0;
				for (nit.begin(); nit.valid(); nit.advance()) {
					INT_TYPE nid = nit.value();
					living_nodes[nid] = pos;
					write_node wn;
					node<SCALAR_TYPE>& n = this->getNode(nid);
					wn.cellindex = n.cellindex;
					wn.value = n.value;
					ws.nodes.push_back(wn);
					pos++;
				}

				LivingArcsIterator ait(this);
				for (ait.begin(); ait.valid(); ait.advance()) {
					INT_TYPE aid = ait.value();
					write_arc wa;
					arc<SCALAR_TYPE>& a = this->getArc(aid);
					wa.lower = living_nodes[a.lower];
					wa.upper = living_nodes[a.upper];
					vector<INDEX_TYPE> geom;
					this->fillArcGeometry(aid, geom);
					wa.geom_size = geom.size();
					ws.geom.insert(ws.geom.end(), geom.begin(), geom.end());
					ws.arcs.push_back(wa);
				}

				ws.write_to_file(filename);

			}

			void LoadComplex1Skeleton(const char* filename) {

				write_state ws;
				ws.read_from_file(filename);
				printf("read %d nodes, %d arcs, %d geom\n", ws.nodes.size(), ws.arcs.size(), ws.geom.size());
				for (auto wn : ws.nodes) {
					this->createNode(wn.cellindex, wn.value);
				}
				auto geom_it = ws.geom.begin();
				for (auto wa : ws.arcs) {
					int size = wa.geom_size;
					vector<INDEX_TYPE> geom(wa.geom_size);
					geom.insert(geom.begin(), geom_it, geom_it + wa.geom_size);
					geom_it += wa.geom_size;
					this->createArc(ws.nodes[wa.lower].cellindex, ws.nodes[wa.upper].cellindex, geom);
				}
			}


		};

		template<typename SCALAR_TYPE, class MESH_TYPE, class FUNC_TYPE, class GRAD_TYPE>
		class MorseSmaleComplexRestrict1SaddleCancel : public MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>
		{
		public:

			MorseSmaleComplexRestrict1SaddleCancel(GRAD_TYPE* grad,
				MESH_TYPE* mesh,
				FUNC_TYPE* func) :
				MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>(grad, mesh, func) {
			}


		protected:
			inline int find_mins_of_1saddle(INT_TYPE id_1saddle, INT_TYPE* min_ids, INT_TYPE ctime) const {
				INT_TYPE arcs_around_1saddle = this->nodes[id_1saddle].firstarc;
				int num_mins = 0;
				while (arcs_around_1saddle != NULLID) {
					const arc<SCALAR_TYPE>& arc_ref = this->arcs[arcs_around_1saddle];
					if (isAlive(arc_ref, ctime) && arc_ref.dim == 0) {
						min_ids[num_mins++] = arc_ref.lower;
					}
					arcs_around_1saddle = nextArc(arc_ref, id_1saddle);
				}
				return num_mins;
			}

			inline int count_other_connection_to_min(INT_TYPE min_id, INT_TYPE other_min_id, INT_TYPE ctime) const {
				INT_TYPE arcs_around_min = this->nodes[min_id].firstarc;
				int count_connections = 0;
				while (arcs_around_min != NULLID) {
					const arc<SCALAR_TYPE>& arc_ref = this->arcs[arcs_around_min];
					if (isAlive(arc_ref, ctime)) {
						INT_TYPE min_ids[2];
						int num_mins = find_mins_of_1saddle(arc_ref.upper, min_ids, ctime);
						if (num_mins != 2) continue;
						if (other_min_id == min_ids[0] || other_min_id == min_ids[1]) {
							count_connections++;
							if (count_connections >= 2) return count_connections;
						}
					}
					arcs_around_min = nextArc(arc_ref, min_id);
				}
				return count_connections;
			}

			int is_only_remaining_1saddle_between_two_mins(arc<SCALAR_TYPE>& ap, INT_TYPE ctime) const {
				if (ap.dim != 1) return false;
				INT_TYPE id_1saddle = ap.lower;
				//int counter = 0;
				int num_mins = 0;
				INT_TYPE min_ids[2];
				num_mins = find_mins_of_1saddle(id_1saddle, min_ids, ctime);
				if (num_mins != 2) return false;

				// have 2 distinct minima
				int count = count_other_connection_to_min(min_ids[0], min_ids[1], ctime);

				return count < 2;
			}

			virtual bool isValid(INT_TYPE a, arc<SCALAR_TYPE>& ap) const {
				// test for boundary
				if (this->nodes[ap.lower].boundary !=
					this->nodes[ap.upper].boundary) return false;

				// 2 endpoints must be connected by exactly one arc
				if (countMultiplicity(ap, this->num_cancelled) != 1) return false;
				// test for inversions?

				if (is_only_remaining_1saddle_between_two_mins(ap, this->num_cancelled)) return false;

				return true;


			}


		};



		template<typename SCALAR_TYPE, class MESH_TYPE, class FUNC_TYPE>
		class GradientlessMSC : public NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE> {
		public:
			typedef FUNC_TYPE FuncType;


		protected:

			FUNC_TYPE* mFunc;

		public:

			INT_TYPE AddNode(INDEX_TYPE node_cell_id) {
				return createNode(node_cell_id);
			}

			GradientlessMSC(
				MESH_TYPE* mesh,
				FUNC_TYPE* func) : NoFuncNoGradMSC<SCALAR_TYPE, MESH_TYPE>(mesh), mFunc(func)
			{}
		protected:
			INT_TYPE createNode(INDEX_TYPE cellID) {
				return createNode(cellID, mFunc->cellValue(cellID));
			}

		};
};
#endif


