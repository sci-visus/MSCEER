#ifndef BASIC_GEOMETRY_H	
#define BASIC_GEOMETRY_H

#include "gi_basic_types.h"
#include "gi_vectors.h"
#include <vector>

namespace GInt {


	template <typename FCoordType>
	class Line {
	protected:
		std::vector<FCoordType> m_simple_points;
		std::vector<FCoordType> m_points;
		float m_length;
	public:
		Line() : m_length(0) {}


		virtual FCoordType GetPointAlongLine(bool from_start, float dist) const {
			int id = (from_start ? 0 : m_points.size() - 1);

			if (from_start) {
				if (dist > Length()) return m_points[m_points.size() - 1];
				int id = 0;
				float distleft = dist;
				while (true) {
					if (id > m_points.size() - 2 || distleft <= 0) {
						printf("Something went wrong, %d, %f\n", id, distleft);
						return m_points[m_points.size() - 1];
					}
					float seg_length = (m_points[id] - m_points[id + 1]).Mag();
					if (seg_length >= distleft) {
						// interpolate to get point
						auto a = m_points[id];
						auto b = m_points[id + 1];
						return (b * (distleft / seg_length) + a * (1 - (distleft / seg_length)));
					}
					else {
						distleft -= seg_length;
						id++;
					}
				}
			} else {
				if (dist > Length()) return m_points[0];
				int id = m_points.size() - 1;
				float distleft = dist;
				while (true) {
					
					if (id < 0 || distleft <= 0) {
						printf("Something went wrong, %d, %f\n", id, distleft);
						return m_points[0];
					}
					float seg_length = (m_points[id] - m_points[id - 1]).Mag();
					if (seg_length >= distleft) {
						// interpolate to get point
						auto a = m_points[id];
						auto b = m_points[id - 1];
						return (b * (distleft / seg_length) + a * (1 - (distleft / seg_length)));
					}
					else {
						distleft -= seg_length;
						id--;
					}
				}
			}
		}



		virtual const std::vector<FCoordType>& GetLine() const { return m_points; }
		virtual const std::vector<FCoordType>& GetSimpleLine() const { return m_simple_points; }

        virtual const int line_size() const { return m_points.size(); }

		virtual void RecomputeLength() {
			m_length = 0;
			for (int i = 1; i < m_points.size(); i++) {
				m_length += (m_points[i] - m_points[i - 1]).Mag();
			}
		}
		virtual float Length() const {
			return m_length;
		}
		virtual void AddToEnd(FCoordType point) {
			m_points.push_back(point);
			if (m_points.size() > 1)
				m_length += (m_points[m_points.size() - 1] -
					m_points[m_points.size() - 2]).Mag();
		}

		virtual void IrreversableSmooth(int iterations) {
			int size = m_points.size();
			if (size < 3) return;
			FCoordType a, b, c;

			for (; iterations > 0; iterations--) {
				a = m_points[0]; // old a, b
				b = m_points[1];
				for (int j = 1; j < size - 1; j++) {
					c = m_points[j + 1];
					m_points[j] = a * 0.25 + b*0.5 + c * .25;
					a = b;
					b = c;
				}
			}
			RecomputeLength();
		}

		virtual bool ReadNextLineFromFileBin(FILE* f) {
			if (feof(f)) return false;
			int count;
			fread(&count, sizeof(int), 1, f);
			while (count > 0) {
				FCoordType v;
				fread(&v, sizeof(FCoordType), 1, f);
				m_points.push_back(v);
				count--;
			}
			RecomputeLength();
			return true;
		}
		virtual void WriteLineToFileBin(FILE* f) {
			int count = m_points.size();
			fwrite(&count, sizeof(int), 1, f);

			fwrite(m_points.data(), sizeof(FCoordType), count, f);
		}
		float triangle_error_sq(FCoordType a, FCoordType b, FCoordType c) const {
			auto bc = b - c;
			auto ac = a - c;
			ac.Normalize();
			auto proj_p = ac * (bc.Dot(ac));
			return (bc - proj_p).MagSq();
		}
		
		virtual void ComputeDecimatedLine(float threshold)  {
			m_simple_points.clear();
			auto thesh_sq = threshold * threshold;
			
			int num_points = this->m_points.size();
			// doubly linked list elements to live in a vector
			struct llp {
				int prev;  // -1 is no element
				int me;
				int next;	// -1 is no element
				int exists; // the iteration this dies, -1 is alive
			};
			std::vector<llp> ll;
			ll.reserve(num_points);
			for (int i = 0; i < num_points; i++) {
				ll.push_back({ i - 1, i, i + 1, -1 });
			}
			(*(ll.rbegin())).next = -1; // fix last element

			int iteration = 0; // preserve the order things were removed
			while (true) {
				// get the lowest cost living element
				// always skip first and last
				float mindist = 99999999.0f;;
				int lowest_element = -1;
				for (int i = 0; i != -1; i = ll[i].next) {
					if (i == 0 || i == num_points - 1) continue; // skip first and last - we dont simplify these
					// compute weight
					auto& n = ll[i];
					float ndist = triangle_error_sq(this->m_points[n.prev],
						this->m_points[n.me], this->m_points[n.next]);

					if (ndist < mindist) {
						mindist = ndist;
						lowest_element = n.me;
					}
				}
				//printf("stick iter %d, min dist %f, min element %d\n", iteration, mindist, lowest_element);
				// now we have checked for the lowest cost node removal

				// if it existed AND lower than our threshold, 
				// remove the node and continue
				// -- note we are using squared values to avoid square root
				if (lowest_element != -1 && mindist < thesh_sq) {
					auto& n = ll[lowest_element];
					auto& prev = ll[n.prev];
					auto& next = ll[n.next];
					n.exists = iteration++;
					prev.next = n.next;	// reconnect linked list to remove n
					next.prev = n.prev;
				}
				// otherwise we are done simplifying so break out of loop
				else {
					break;
				}
			}

			// now gather the points that are still alive
			int ni = 0;
			while (ni != -1) {
				m_simple_points.push_back(this->m_points[ll[ni].me]);
				ni = ll[ni].next;
			}
		}

	};

	typedef Line<Vec3f> Line3d;
	typedef Line<Vec2f> Line2d;

	template <typename FCoordType>
	class ConstrainedLine : public Line<FCoordType> {
	protected:
		std::vector<FCoordType> m_original_points;
	public:

		virtual void AddToEnd(FCoordType point) {
			this->m_points.push_back(point);
			m_original_points.push_back(point);
		}

		virtual void IrreversableSmooth(int iterations) {
			int size = this->m_points.size();
			if (size < 3) return;
			FCoordType a, b, c;

			for (; iterations > 0; iterations--) {
				a = this->m_points[0]; // old a, b
				b = this->m_points[1];
				for (int j = 1; j < size - 1; j++) {
					c = this->m_points[j + 1];
					this->m_points[j] = a * 0.25 + b*0.5 + c * .25;
					a = b;
					b = c;
				}
				for (int j = 1; j < size - 1; j++) {
					FCoordType diff = this->m_points[j] - m_original_points[j];
					for (int k = 0; k < diff.Dim(); k++) {
						if (diff[k] < -0.5) diff[k] = -0.5;
						if (diff[k] > 0.5) diff[k] = 0.5;
					}
					this->m_points[j] = m_original_points[j] + diff;
				}
			}
		}

		//virtual bool ReadNextLineFromFileBin(FILE* f) {
		//	if (feof(f)) return false;
		//	int count;
		//	fread(&count, sizeof(int), 1, f);
		//	while (count > 0) {
		//		Vec3f v;
		//		fread(&v, sizeof(Vec3f), 1, f);
		//		m_points.push_back(v);
		//		count--;
		//	}
		//	return true;
		//}
		//virtual void WriteLineToFileBin(FILE* f) {
		//	int count = m_points.size();
		//	fwrite(&count, sizeof(int), 1, f);

		//	fwrite(m_points.data(), sizeof(Vec3f), count, f);
		//}

	};

	// this is what you get if you take a normal, an "up" direction (assert normal not parallel to up)
	// and put a patch 2d grid on it \_\_\
	//                                \_\_\
	// follows right-hand-rule: thumb = norm, index = up, middle = left
	// the grid has m_x_samples, m_y_samples, 
	// 
	// 0,0 of the grid occurs at geometry m_origin + m_up*m_x_dist + m_left*m_y_dist
	// m_x_sampels, 0 = m_origin + m_up*m_x_dist - m_left*m_y_dist
	// 0, m_y_samples = m_origin -  m_up*m_x_dist + m_left*m_y_dist
	class PlaneGrid2D {
	protected:
		Vec3f m_origin;
		Vec3f m_norm;
		Vec3f m_up;
		Vec3f m_left;
		float m_last_x;
		float m_last_y;
		float m_x_dist;
		float m_y_dist;
		void initialize(Vec3f origin, Vec3f norm, Vec3f up,
			int numx, int numy, float radius_x, float radius_y) {
			m_origin = origin;
			m_norm = norm;
			m_up = up;
			m_last_x = numx - 1;
			m_last_y = numy - 1;
			m_x_dist = radius_x;
			m_y_dist = radius_y;

			// now compute some
			m_norm.Normalize();
			m_up.Normalize();
			auto dot = m_norm.Dot(m_up);
			if (dot > 0.95) {
				printf("warning planeGrid2D norm and up dot = %f\n", dot);
			}
			// project up onto norm
			auto proj = m_norm * dot;
			// make up orthogonal to norm
			m_up -= proj;
			m_up.Normalize();
			// now they are definitely orthogonal

			m_left = GInt::Vec3f::CrossF(m_norm, m_up);

			// now scale
			m_up *= m_y_dist;
			m_left *= m_x_dist;

			// now put origin at the corner
			// and flip up and left
			m_origin = m_origin + m_up + m_left;
			m_up *= -2.f;
			m_left *= -2.f;
		}
	public:
		PlaneGrid2D(Vec3f origin, Vec3f norm, Vec3f up, 
			int numx, int numy, float radius_x, float radius_y) {
			initialize(origin, norm, up, numx, numy, radius_x, radius_y);
		}

		// here we don't have an "up" vector - so we look at the smallest component
		// of the normal vector, 
		PlaneGrid2D(Vec3f origin, Vec3f norm,
			int numx, int numy, float radius_x, float radius_y) {
			int smallest_axis = 0;
			for (int i = 1; i < 3; i++) {
				if (norm[i] < norm[smallest_axis]) smallest_axis = i;
			}
			Vec3f up(0, 0, 0);
			up[smallest_axis] = 1;
			initialize(origin, norm, up, numx, numy, radius_x, radius_y);
		}

		Vec3f Get3DPosition(int i, int j) {
			return m_origin +
				m_up * ((float)j / m_last_y) +
				m_left * ((float)i / m_last_x);
		}
	};


	// start wiht set of "edges"
	// make set of quads
	// each vertex has list of neighbors

	//class my_surf {

	//	struct myV {
	//		Vec3f xyz;
	//		INDEX_TYPE id;
	//		int start_neg;
	//		int end_neg;
	//	};

	//	std::vector<INDEX_TYPE> mEdgeIds;
	//	std::vector<int> neighbors;
	//	std::vector<myV> mVertices;
	//	void createSurface() {

	//	}

	//public: 
	//	void add_edge_as_quad(INDEX_TYPE id) {
	//		mEdgeIds.push_back(id);
	//	}

	//};



	//template< int NUM_SIDES>
	//class abstract_surface {
	//protected:
	//	// my half-edge data structure
	//	// vertex - pointer to first half-edge		-- 4 bytes
	//	// half-edge:	tail, head vertices			-- 8 bytes
	//	//				next half-edge around edge  -- 4 bytes
	//	//				n/p half-edge around face	-- 8 bytes
	//	//				face id						-- 4 bytes
	//	// face - pointer to first half-edge		-- 4 bytes


	//	// list of vertices
	//	// map to position kept elsewhere
	//	struct abstract_vertex {
	//		INT_TYPE first_edge;
	//	};

	//	struct abstract_edge {
	//		INT_TYPE vertices[2];
	//		INT_TYPE next_edge[2];

	//		INT_TYPE first_face;
	//	};

	//	struct abstract_face {
	//		INT_TYPE vertices[NUM_SIDES];
	//		INT_TYPE next_face_vertex[NUM_SIDES];
	//	};


	//	abstract_vertex& get_vertex(INT_TYPE vid);
	//	abstract_edge& get_edge(INT_TYPE eid);
	//	abstract_face& get_face(INT_TYPE fid);

	//	INT_TYPE get_next_edge_id(INT_TYPE vid, INT_TYPE curr_eid) const {
	//		const abstract_edge& e = get_edge(curr_eid);
	//		if (e.vertices[0] == vid) return e.next_edge[0];
	//		return e.next_edge[1];
	//	}


	//	INT_TYPE fill_neighbor_edges(INT_TYPE vid, INT_TYPE* edge_list) const {
	//		int pos = 0;
	//		INT_TYPE curr_eid = get_vertex(vid).first_edge;
	//		while (curr_eid != NULLID) {
	//			edge_list[pos++] = curr_eid;
	//			curr_eid = get_next_edge_id(vid, curr_eid);
	//		}
	//		return pos;
	//	}




	//	// vertex, has pointer to first edge

	//	// list of edges
	//	// edge has next vertex for start, end, next face

	//	// face has ?? vertices-> fast rendering, normal?
	//





	//	struct geombitfield {
	//		unsigned char f1 : 1;
	//		unsigned char f2 : 1;
	//		unsigned char f3 : 1;
	//		unsigned char f4 : 5;
	//	};

	//	struct g_vertex {
	//		float coords[3];
	//		float norms[3];
	//		INT_TYPE id;
	//		geombitfield flags;
	//		set<INT_TYPE> adjacent;
	//		INDEX_TYPE insertID;

	//	};


	//	struct g_quad {
	//		geombitfield flags;
	//		INT_TYPE v1;
	//		INT_TYPE v2;
	//		INT_TYPE v3;
	//		INT_TYPE v4;
	//		float centroid[3];
	//	};

	//	INT_TYPE make_key(INT_TYPE v1, INT_TYPE v2) {
	//		INT_TYPE a = (INT_TYPE) (v1 < v2 ? v1 : v2);
	//		INT_TYPE b = (INT_TYPE) (v1 > v2 ? v2 : v1);
	//		INT_TYPE res = (INT_TYPE)a;
	//		res = (res << 16) + (INT_TYPE)b;
	//		return res;
	//	}

	//	INT_TYPE UID;
	//public:

	//	vector< g_quad > quads;
	//	vector< g_vertex<dtype> > vertices;
	//	vector< g_edge > edges;

	//	std::unordered_map< INDEX_TYPE, INT_TYPE > vertex_map;
	//	std::unordered_map< INT_TYPE, INT_TYPE > edge_map;
	//	surface() {
	//		UID = INT_INFTY / 2;
	//	}



	//	// returns the id of the vertex added, or the id of the vertex
	//	// if it alredy exists
	//	INT_TYPE add_vertex(float* coords, INDEX_TYPE gid, dtype val) {
	//		//printf("b");
	//		if (vertex_map.count(gid) > 0) {
	//			//printf("b");
	//			return vertex_map[gid];
	//		}
	//		g_vertex<dtype> v;
	//		for (int i = 0; i < 3; i++) v.coords[i] = coords[i];
	//		v.flags.f4 = 0;
	//		v.value = val;
	//		v.insertID = gid;
	//		INT_TYPE position = (INT_TYPE)vertices.size();
	//		v.id = position;
	//		vertices.push_back(v);
	//		vertex_map[gid] = position;
	//		//printf("b");
	//		return position;
	//	}

	//	INT_TYPE add_edge(INT_TYPE v1, INT_TYPE v2, INDEX_TYPE keyval) {
	//		g_edge e;
	//		e.v1 = v1; e.v2 = v2; e.flags.f4 = 1;
	//		INT_TYPE pos = (INT_TYPE)edges.size();
	//		edges.push_back(e);
	//		edge_map[keyval] = pos;
	//		return pos;
	//	}

	//	void add_quad_edge(INT_TYPE v1, INT_TYPE v2) {
	//		INT_TYPE keyval = make_key(v1, v2);
	//		INT_TYPE edgepos;
	//		if (edge_map.count(keyval) > 0) {
	//			edgepos = edge_map[keyval];
	//			g_edge& e = edges[edgepos];
	//			e.flags.f4++;
	//		}
	//		else {
	//			edgepos = add_edge(v1, v2, keyval);
	//		}
	//		vertices[v1].adjacent.insert(edgepos);
	//		vertices[v2].adjacent.insert(edgepos);
	//	}

	//	INT_TYPE add_quad(INT_TYPE v1, INT_TYPE v2, INT_TYPE v3, INT_TYPE v4) {
	//		g_quad q;
	//		q.flags.f4 = 0;
	//		q.v1 = v1; q.v2 = v2; q.v3 = v3; q.v4 = v4;
	//		INT_TYPE pos = (INT_TYPE)quads.size();
	//		quads.push_back(q);
	//		// now connect them
	//		// and increment the number of quads touching that vertex
	//		vertices[v1].flags.f4++;
	//		add_quad_edge(v1, v2);

	//		vertices[v2].flags.f4++;
	//		add_quad_edge(v2, v3);
	//		//vertices[v2].adjacent.insert(pos);
	//		////vertices[v2].adjacent.insert(v3);

	//		vertices[v3].flags.f4++;
	//		add_quad_edge(v3, v4);

	//		vertices[v4].flags.f4++;
	//		add_quad_edge(v4, v1);

	//		return pos;
	//	}

	//	void merge_surface(surface<dtype>* s) {

	//		for (int i = 0; i < s->vertices.size(); i++) {
	//			g_vertex<dtype>& v = s->vertices[i];
	//			this->add_vertex(v.coords, v.insertID, v.value);
	//		}
	//		for (int i = 0; i < s->quads.size(); i++) {
	//			g_quad& q = s->quads[i];
	//			INT_TYPE v1 = add_vertex(NULL, s->vertices[q.v1].insertID, 0);
	//			INT_TYPE v2 = add_vertex(NULL, s->vertices[q.v2].insertID, 0);
	//			INT_TYPE v3 = add_vertex(NULL, s->vertices[q.v3].insertID, 0);
	//			INT_TYPE v4 = add_vertex(NULL, s->vertices[q.v4].insertID, 0);
	//			this->add_quad(v1, v2, v3, v4);
	//		}

	//	}

	//	void dump_surf(FILE* fs) {
	//		int iheader[3];

	//		iheader[0] = vertices.size();
	//		iheader[1] = quads.size();
	//		iheader[2] = 0;
	//		fwrite(iheader, sizeof(int), fs);

	//		float vert[5];
	//		for (INT_TYPE i = 0; i < vertices.size(); i++) {
	//			g_vertex<dtype>& v = vertices[i];
	//			fwrite(v.coords, sizeof(float), fs);
	//			vert[0] = ((float)v.value);
	//			fwrite(vert, sizeof(float), 1, fs);
	//		}

	//		INT_TYPE quadv[4];

	//		for (INT_TYPE i = 0; i < quads.size(); i++) {
	//			g_quad &q = quads[i];
	//			quadv[0] = q.v1;
	//			quadv[1] = q.v2;
	//			quadv[2] = q.v3;
	//			quadv[3] = q.v4;
	//			fwrite(quadv, sizeof(INT_TYPE), 4, fs);
	//		}
	//	}



	//	void dumpInObj(char* filename) {
	//		FILE* obj = fopen(filename, "w");
	//		fprintf(obj, "# my object\n");

	//		for (INT_TYPE i = 0; i < vertices.size(); i++) {
	//			g_vertex< dtype > &v = vertices[i];
	//			fprintf(obj, "v %.4f %.4f %.4f\n", v.coords[0], v.coords[1], v.coords[2]);
	//		}

	//		////// USING TRIANGLES
	//		for (INT_TYPE i = 0; i < quads.size(); i++) {
	//			g_quad &q = quads[i];
	//			fprintf(obj, "f %d %d %d\n", q.v1 + 1, q.v2 + 1, q.v3 + 1);
	//			fprintf(obj, "f %d %d %d\n", q.v3 + 1, q.v4 + 1, q.v1 + 1);
	//			//fprintf(obj, "f %d %d %d\n", q.v3 +1, q.v4 +1, q.v1 +1);
	//		}
	//		fclose(obj);
	//	}

	//	void dumpInObjWN(char* filename) {
	//		FILE* obj = fopen(filename, "w");
	//		fprintf(obj, "# my object\n");

	//		for (INT_TYPE i = 0; i < vertices.size(); i++) {
	//			g_vertex< dtype > &v = vertices[i];
	//			fprintf(obj, "v %.4f %.4f %.4f\n", v.coords[0], v.coords[1], v.coords[2]);
	//		}
	//		for (INT_TYPE i = 0; i < vertices.size(); i++) {
	//			g_vertex< dtype > &v = vertices[i];
	//			fprintf(obj, "vn %.4f %.4f %.4f\n", v.norms[0], v.norms[1], v.norms[2]);
	//		}
	//		for (INT_TYPE i = 0; i < quads.size(); i++) {
	//			g_quad &q = quads[i];
	//			fprintf(obj, "f %d//%d %d//%d %d//%d %d//%d\n", q.v1 + 1, q.v1 + 1,
	//				q.v2 + 1, q.v2 + 1, q.v3 + 1, q.v3 + 1, q.v4 + 1, q.v4 + 1);
	//			//fprintf(obj, "f %d %d %d\n", q.v3 +1, q.v4 +1, q.v1 +1);
	//		}
	//		fclose(obj);
	//	}


	//	void compute_normals() {
	//		////// compute normal on each quad

	//		////// average normals of quads
	//		for (INT_TYPE i = 0; i < quads.size(); i++) {
	//			g_quad &q = quads[i];
	//			g_vertex<dtype>& v1 = vertices[q.v1];
	//			g_vertex<dtype>& v2 = vertices[q.v2];
	//			g_vertex<dtype>& v3 = vertices[q.v3];
	//			g_vertex<dtype>& v4 = vertices[q.v4];
	//			Vec3d vm1;
	//			vm1[0] = v1.coords[0] * 0.5f - v2.coords[0] * 0.5f;
	//			vm1[1] = v1.coords[1] * 0.5f - v2.coords[1] * 0.5f;
	//			vm1[2] = v1.coords[2] * 0.5f - v2.coords[2] * 0.5f;
	//			Vec3d vm2;
	//			vm2[0] = v1.coords[0] * 0.5f - v3.coords[0] * 0.5f;
	//			vm2[1] = v1.coords[1] * 0.5f - v3.coords[1] * 0.5f;
	//			vm2[2] = v1.coords[2] * 0.5f - v3.coords[2] * 0.5f;
	//			Vec3d res = vm1.Cross(vm2);
	//			res.Normalize();
	//			res.Abs();

	//		}
	//	}

	//	void smooth(int iterations) {

	//		for (INT_TYPE i = 0; i < vertices.size(); i++) {
	//			g_vertex< dtype > &v = vertices[i];
	//			int minnum = 11110;
	//			set<INT_TYPE>::iterator it = v.adjacent.begin();
	//			while (it != v.adjacent.end()) {
	//				g_edge &e = edges[*it];
	//				if (e.flags.f4 < minnum) minnum = e.flags.f4;
	//				it++;
	//			}
	//			v.flags.f4 = minnum;
	//		}
	//		while (iterations > 0) {
	//			iterations--;
	//			// compute centroids
	//			for (INT_TYPE i = 0; i < edges.size(); i++) {
	//				g_edge &e = edges[i];
	//				g_vertex<dtype>& v1 = vertices[e.v1];
	//				for (int j = 0; j<3; j++) e.centroid[j] = v1.coords[j];

	//				g_vertex<dtype>& v2 = vertices[e.v2];
	//				for (int j = 0; j<3; j++) e.centroid[j] += v2.coords[j];

	//				for (int j = 0; j<3; j++) e.centroid[j] *= 0.5f;
	//			}
	//			// now move vertices!
	//			for (INT_TYPE i = 0; i < vertices.size(); i++) {
	//				g_vertex< dtype > &v = vertices[i];

	//				//printf(" v=%d %d \n", v.flags.f4, v.adjacent.size());
	//				//if (v.flags.f4 <= 1) continue;
	//				//clear the position
	//				for (int j = 0; j<3; j++) v.coords[j] = 0.0f;
	//				float numadded = 0.0f;
	//				set<INT_TYPE>::iterator it = v.adjacent.begin();
	//				while (it != v.adjacent.end()) {
	//					g_edge &e = edges[*it];
	//					if (v.flags.f4 >= e.flags.f4) {
	//						for (int j = 0; j<3; j++) v.coords[j] += e.centroid[j];
	//						numadded += 1.0f;
	//					}
	//					it++;
	//				}
	//				float scale = 1.0f / numadded;
	//				for (int j = 0; j<3; j++) v.coords[j] *= scale;
	//			}
	//		}


	//	}
	//};
};

#endif
