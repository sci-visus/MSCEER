#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif
#include "extractmsc.h"
#include <time.h>

using namespace GInt;




namespace GeomSurf {

	struct g_vertex {
		INDEX_TYPE insertID;
		INT_TYPE id;
		GInt::Vec3f coords;
		std::set<INT_TYPE> adjacent;
	};


	struct g_quad {
		INT_TYPE v1;
		INT_TYPE v2;
		INT_TYPE v3;
		INT_TYPE v4;
		//float centroid[3];
	};


	struct g_edge {
		INT_TYPE v1;
		INT_TYPE v2;
		BYTE_TYPE num_quads;
		//float centroid[3];
		g_edge() {
			num_quads = 0;
		}
	};

	class surface {
	protected:
		INDEX_TYPE make_key(INT_TYPE v1, INT_TYPE v2) {
			if (v1 < v2) {
				return (((INDEX_TYPE)v1) << 32) | ((INDEX_TYPE)v2);
			}
			else {
				return (((INDEX_TYPE)v2) << 32) | ((INDEX_TYPE)v1);
			}
		}

		void add_quad_edge(INT_TYPE v1, INT_TYPE v2) {
			INDEX_TYPE edgekey = make_key(v1, v2);
			if (edge_map.count(edgekey) == 0) {
				edge_map[edgekey] = edges.size();
				g_edge e; e.v1 = v1; e.v2 = v2;
				edges.push_back(e);
			}
			INT_TYPE edgepos = edge_map[edgekey];
			edges[edgepos].num_quads++;
			vertices[v1].adjacent.insert(edgepos);
			vertices[v2].adjacent.insert(edgepos);
		}

		std::vector< g_vertex > vertices;
		std::vector< g_edge > edges;
		std::vector< g_quad > quads;

		std::unordered_map< INDEX_TYPE, INT_TYPE > vertex_map;
		std::unordered_map< INDEX_TYPE, INT_TYPE > edge_map;

		INT_TYPE other_vertex(INT_TYPE vid, INT_TYPE eid) {
			if (vid == edges[eid].v1) return edges[eid].v2;
			return edges[eid].v1;
		}
	public:

		surface() {}

		// returns the id of the vertex added, or the id of the vertex
		// if it alredy exists
		virtual INT_TYPE add_vertex(GInt::Vec3f coords, INDEX_TYPE gid) {
			if (vertex_map.count(gid) > 0) {
				return vertex_map[gid];
			}
			g_vertex v;
			v.coords = coords;
			v.insertID = gid;
			INT_TYPE position = (INT_TYPE)vertices.size();
			v.id = position;
			vertices.push_back(v);
			vertex_map[gid] = position;
			return position;
		}

		virtual int num_vertices() {
			return vertices.size();
		}
		virtual int num_quads() {
			return quads.size();
		}

		// ordered around quad
		virtual INT_TYPE add_quad(INT_TYPE v1, INT_TYPE v2, INT_TYPE v3, INT_TYPE v4) {
			g_quad q;
			q.v1 = v1; q.v2 = v2; q.v3 = v3; q.v4 = v4;
			INT_TYPE pos = (INT_TYPE)quads.size();
			quads.push_back(q);
			// now connect them
			// and increment the number of quads touching that vertex
			//vertices[v1].flags.f4++;
			add_quad_edge(v1, v2);

			//vertices[v2].flags.f4++;
			add_quad_edge(v2, v3);
			//vertices[v2].adjacent.insert(pos);
			////vertices[v2].adjacent.insert(v3);

			//vertices[v3].flags.f4++;
			add_quad_edge(v3, v4);

			//vertices[v4].flags.f4++;
			add_quad_edge(v4, v1);

			return pos;
		}

		virtual GInt::Vec3f vertex_position(INT_TYPE id) {
			printf("USING NON_SMOOTH POSITION\n");
			return vertices[id].coords;
		}

		virtual GInt::Vec3f quad_normal(INT_TYPE id) {
			GInt::Vec3f v2 = vertex_position(quads[id].v2);
			GInt::Vec3f v1 = vertex_position(quads[id].v1);
			GInt::Vec3f v2_v1 = v2 - v1;
			GInt::Vec3f v3 = vertex_position(quads[id].v3);
			GInt::Vec3f v2_v3 = v2 - v3;
			GInt::Vec3f norm1 = GInt::Vec3f::CrossF(v2_v1, v2_v3);
			GInt::Vec3f v4_v3 = vertex_position(quads[id].v4) - vertex_position(quads[id].v3);
			GInt::Vec3f v4_v1 = vertex_position(quads[id].v4) - vertex_position(quads[id].v1);
			GInt::Vec3f norm2 = GInt::Vec3f::CrossF(v4_v3, v4_v1);
			GInt::Vec3f crosssum = (norm1 + norm2) * 0.5;
			crosssum.Normalize();
			return crosssum;
		}

		// expects verts to be INT_TYPE[4]
		virtual void quad_vertices(INT_TYPE quadid, INT_TYPE* verts) {
			verts[0] = quads[quadid].v1;
			verts[1] = quads[quadid].v2;
			verts[2] = quads[quadid].v3;
			verts[3] = quads[quadid].v4;
		}

		// number of neighbors, expects to fill INT_TYPE[6]
		virtual INT_TYPE vertex_neighbors(INT_TYPE vid, INT_TYPE* verts) {
			INT_TYPE count = 0;
			for (auto eid : vertices[vid].adjacent) {
				verts[count++] = other_vertex(vid, eid);
			}
			return count;
		}

		virtual int edge_incidence(INT_TYPE v1, INT_TYPE v2) {
			return edges[edge_map[make_key(v1, v2)]].num_quads;
		}

		virtual BYTE_TYPE is_edge_boundary(INT_TYPE eid) {
			return edges[eid].num_quads == 1;
		}

	};

	class smooth_surface : public surface {
	protected:
		std::vector<GInt::Vec3f> new_vert_pos;
	public:

		virtual INT_TYPE add_vertex(GInt::Vec3f coords, INDEX_TYPE gid) {
			INT_TYPE pos = surface::add_vertex(coords, gid);
			if (pos == new_vert_pos.size()) new_vert_pos.push_back(coords);
			return pos;
		}
		virtual GInt::Vec3f vertex_position(INT_TYPE id) {
			return new_vert_pos[id];
		}

		virtual void smooth(int num) {
			//printf("%d sizes %d, %d, %d\n", num, quads.size(), new_vert_pos.size(), vertices.size());
			if (quads.size() == 0) return;
			for (int i = 0; i < vertices.size(); i++) new_vert_pos[i] = vertices[i].coords;
			std::vector<GInt::Vec3f> tmp_pos;
			tmp_pos.resize(vertices.size());
			for (int iteration = num; iteration > 0; iteration--) {

				for (int i = 0; i < vertices.size(); i++) {
					INT_TYPE negs[29];
					int n = vertex_neighbors(i, negs);
					//printf("--doing %d, %d negs:\n",i, n);
					//new_vert_pos[i].PrintFloat();
					//if (is_vertex_boundary(i)) {
					//	int cn = 0;
					//	INT_TYPE nnegs[8];
					//	for (int j = 0; j < n; j++) {
					//		if (is_vertex_boundary(negs[j])) nnegs[cn++]=negs[j];
					//	}
					//	float oneoversize = 1.0f / (1.0f + cn);
					//	GInt::Vec3f newpos = new_vert_pos[i] * oneoversize;
					//	for (int j = 0; j < cn; j++) {
					//		newpos += vertices[nnegs[j]].coords * oneoversize;
					//	}
					//	tmp_pos[i] = newpos;
					//}
					//else {
					float oneoversize = 1.0f / (1.0f + n);
					GInt::Vec3f newpos = new_vert_pos[i] * oneoversize;
					for (int j = 0; j < n; j++) {
						newpos += vertices[negs[j]].coords * oneoversize;
					}
					tmp_pos[i] = newpos;
					//}
					//tmp_pos[i].PrintFloat();
					//printf("\n");
				}
				//new_vert_pos.insert(new_vert_pos.begin(), tmp_pos.begin(), tmp_pos.end());
				new_vert_pos.swap(tmp_pos);
			}
		}

	};





}


template<class MeshType>
class EdgeSurface {
protected:

	MeshType* m_mesh;
public:
	SurfaceStatistics m_stat;

	set<EdgeSurface*> m_negs;
	EdgeSurface(MeshType* mesh) : m_mesh(mesh) {
		m_surf = new GeomSurf::smooth_surface();
	}


	GeomSurf::smooth_surface* m_surf;
	vector<INDEX_TYPE> m_edges;


	void AddEdge(INDEX_TYPE id) {

		if (m_mesh->boundaryValue(id) != 0) return;

		m_edges.push_back(id);
		// now add to m_surf
		GInt::Vec3l coords, cd_0, cd_1, cd_2, cd_3;
		m_mesh->cellid2Coords(id, coords);


		if (coords[0] % 2 == 1) {
			cd_0 = coords + GInt::Vec3l(0, -1, -1);
			cd_1 = coords + GInt::Vec3l(0, -1, 1);
			cd_2 = coords + GInt::Vec3l(0, 1, 1);
			cd_3 = coords + GInt::Vec3l(0, 1, -1);
		}
		else if (coords[1] % 2 == 1) {
			cd_0 = coords + GInt::Vec3l(-1, 0, -1);
			cd_1 = coords + GInt::Vec3l(-1, 0, 1);
			cd_2 = coords + GInt::Vec3l(1, 0, 1);
			cd_3 = coords + GInt::Vec3l(1, 0, -1);
		}
		else {
			cd_0 = coords + GInt::Vec3l(-1, -1, 0);
			cd_1 = coords + GInt::Vec3l(-1, 1, 0);
			cd_2 = coords + GInt::Vec3l(1, 1, 0);
			cd_3 = coords + GInt::Vec3l(1, -1, 0);
		}

		INT_TYPE v1 = m_surf->add_vertex(GInt::Vec3f(cd_0) * 0.5f, m_mesh->coords2Cellid(cd_0));
		INT_TYPE v2 = m_surf->add_vertex(GInt::Vec3f(cd_1) * 0.5f, m_mesh->coords2Cellid(cd_1));
		INT_TYPE v3 = m_surf->add_vertex(GInt::Vec3f(cd_2) * 0.5f, m_mesh->coords2Cellid(cd_2));
		INT_TYPE v4 = m_surf->add_vertex(GInt::Vec3f(cd_3) * 0.5f, m_mesh->coords2Cellid(cd_3));

		m_surf->add_quad(v1, v2, v3, v4);
	}

	void ComputeStatistics(int numsmooth, GInt::RegularGridTrilinearFunction* gridfunc) {
		if (m_surf->num_quads() == 0) return;
		//printf("about to smooth\n");
		m_surf->smooth(numsmooth);
		//printf("smoothed\n");

		int numverts = m_surf->num_vertices();// m_edges.size();

		m_stat.m_ave_value = 0;

		m_stat.num_verts = numverts;
		m_stat.num_quads = m_surf->num_quads();
		//m_points.reserve(m_edges.size());
		m_stat.m_directions = GInt::Vec3f(0, 0, 0);
		m_stat.m_centroid = GInt::Vec3f(0, 0, 0);
		float oneoverv = 1.0 / numverts;
		m_stat.m_orientation = GInt::Vec3f(0, 0, 0);
		// compute centroid
		m_stat.m_max_value = m_stat.m_min_value = gridfunc->InterpolatedValue(m_surf->vertex_position(0));
		for (int i = 0; i < numverts; i++) {
			m_stat.m_centroid += m_surf->vertex_position(i) * oneoverv;
			float gfv = gridfunc->InterpolatedValue(m_surf->vertex_position(i));
			m_stat.m_ave_value += oneoverv * gfv;
			m_stat.m_max_value = (m_stat.m_max_value < gfv ? gfv : m_stat.m_max_value);
			m_stat.m_min_value = (m_stat.m_min_value > gfv ? gfv : m_stat.m_min_value);
		}


		int numquads = m_surf->num_quads();
		float oneoverq = 1.0f / numquads;
		// now have centroids and points
		for (int i = 0; i < numquads; i++) {
			GInt::Vec3f norm = m_surf->quad_normal(i);
			if (norm.Dot(m_stat.m_orientation) < 0) norm *= -1.0;
			//	if (norm[1] < 0.0) norm *= -1.0;
			m_stat.m_orientation += (norm * oneoverq);
			m_stat.m_directions[1] += norm[1] * oneoverq;
		}
		//m_directions.PrintFloat();
		// EIGEN STUFF HERE
		Eigen::MatrixXf M = Eigen::MatrixXf(numverts, 3);
		for (int i = 0; i < numverts; i++) {
			GInt::Vec3f v = m_surf->vertex_position(i);
			for (int j = 0; j < 3; j++) {
				M(i, j) = v[j];
			}
		}
		Eigen::VectorXf mean = M.colwise().mean();
		//cout << "hand computed: " << m_centroid[0] << ", " << m_centroid[1] << ", " << m_centroid[2] << endl;
		//cout << mean.transpose() << endl;
		//auto mt = mean.transpose();
		Eigen::MatrixXf C = M.rowwise() - mean.transpose();
		//cout << "numedges " << numedges << " " << C.rows() << "x" << C.cols() << endl << C << endl << endl;
		auto cov = C.adjoint() * C;
		////cout << "size: " << cov.rows() << "x" << cov.cols() << endl;
		////mean = mean.transpose();
		////Eigen::MatrixXf C = M.rowwise() - mean;
		////Eigen::MatrixXf cov = C.adjoint() * C;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig(cov);

		////eig.compute(cov);
		//auto v = eig.eigenvectors();
		m_stat.m_eigenvals = eig.eigenvalues().transpose();
		m_stat.m_eigenvectors = eig.eigenvectors();
		//// columns are eigenvectors;
		////cout << "eigenV:" << endl << eig.eigenvalues() << endl << v << endl << m_eigenvectors.col(1);
	}
};

// do the vertex->edge -> other_vertex step in tracing a discrete gradient path
INDEX_TYPE one_step_down_vertex(GradType* mGrad, INDEX_TYPE cellid) {
	if (mGrad->getCritical(cellid)) {
		return cellid;
	}
	// total hack - only works for grids (non-blocked, dense), doubles offset between 
	// vertex-edge pair to get vertex-edge-vertex (discrete gradient arrow, other facet)
	INDEX_TYPE pair = mGrad->getPair(cellid);
	INDEX_TYPE next_v = cellid + 2ll * ((INDEX_TYPE)pair - cellid);
	return next_v;
}

std::unordered_map<INDEX_TYPE, int> crit_cells_2_nodes;
// for vertex ids start thru end, compute ascending manifold and set the corresponding ids in the index volume
void SetAscendingManifolds(DenseLabeling<int>* mIds, MeshType* mesh, GradType* grad, MscType* msc, GridFuncType* data, INDEX_TYPE start, INDEX_TYPE end) {

	vector<INDEX_TYPE> v_path;
	v_path.reserve(1024);
	for (INDEX_TYPE vid = start; vid < end; vid++) {
		if (mIds->GetLabel(vid) != -1) continue; // this has been seen

		// uncomment this block if you want to cut off the manifolds at a funciton value
		// this is used in geobodies, quantum chemistry, imaging to remove "background" 
		/*
			if (data->SampleImage(vid) > 0.9) {
			mIds->SetLabel(vid, -2); // mark a detected background
			continue;
		}
		*/

		// start a V-path trace 
		INDEX_TYPE v_cid = mesh->CellIDFromVertexNumber(vid);
		INDEX_TYPE current = v_cid;
		int dest_id = -1;
		while (true) {

			INDEX_TYPE vert_num = mesh->VertexNumberFromCellID(current);
			v_path.push_back(vert_num);
			
			// early termination:
			// we have seen a path before? just end here then
			if (mIds->GetLabel(vert_num) != -1) {
				dest_id = mIds->GetLabel(vert_num);
				break;
			}
			// did we hit a critical cell? then find out which node and end v-path trace
			else if (grad->getCritical(current)) {
				dest_id = crit_cells_2_nodes[current];
				break;
			}
			current = one_step_down_vertex(grad, current);
		}
		assert(dest_id != -1);
		// now paint back the path to the index volume
		for (auto vid : v_path) mIds->SetLabel(vid, dest_id);
		v_path.clear();
	}
}

void FillAscendingIndexVolumeByVPaths(DenseLabeling<int>* mIds, MeshType* mesh, GradType* grad, MscType* msc, GridFuncType* data) {

#pragma omp parallel
	{
		int numthreads = omp_get_num_threads();
		int threadnum = omp_get_thread_num();
		auto num_cells = data->GetGrid()->NumElements();
		auto startid = threadnum * num_cells / numthreads;
		auto endid = (threadnum + 1) * num_cells / numthreads;
		if (threadnum == numthreads - 1) endid = num_cells;
		SetAscendingManifolds(mIds, mesh, grad, msc, data, startid, endid);
	}
}

int main(int argc, char** argv) {

	// READ IN THE COMMAND LINE ARGUMENTS
	int X, Y, Z;
	std::string filename;
	if (argc < 6) { printf("Usage: X Y Z filename persistence WRITE?\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
	float persistence;
	sscanf(argv[5], "%f", &persistence);
  

	int writeorread = 0;
	sscanf(argv[6], "%d", &writeorread);
  std::string mask_filename;
  
  int data_dims[3]={X,Y,Z};
  printf("extractmsc persistence is %f\n", persistence);
	



	GridType* underlying_grid;
	GridFuncType* grid_function;
	MeshType *topological_grid;
	TopoFuncType* topological_grid_function;
	GradType *discrete_gradient;
	if (writeorread) {
		// set up structures to navigate grid, and load the 3d image
		underlying_grid = new GridType(GInt::Vec3i(X, Y, Z), GInt::Vec3b(0, 0, 0));
		grid_function = new GridFuncType(underlying_grid);
		grid_function->LoadImageFromFile(filename.c_str());

		printf("loaded cont function\n");

		// now set up an indexing scheme to use in a topological interpretation of the 
		// regular grid
		topological_grid = new MeshType(underlying_grid);

		// we will use a lazy-evaluation max vertex mesh function - i.e. the value of a 
		// cell in the topological grid is the maximum value of its vertices in the input
		// image
		MaxVLType* maxv = new MaxVLType(topological_grid, grid_function);
		maxv->ComputeOutput();
		topological_grid_function = new TopoFuncType();
		topological_grid_function->setMeshAndFuncAndMaxVLabeling(topological_grid, grid_function, maxv);

		// read the discrete gradient from disk
		discrete_gradient = new GradType(topological_grid);
		printf("created dgrad struct\n");
		char gradname[2048];
		sprintf(gradname, "%s.grad", argv[4]);
		discrete_gradient->load_from_file(gradname);

		// now compute the Morse-Smale complex
		auto now_time = std::chrono::steady_clock::now();
		MscType* msc = new MscType(discrete_gradient, topological_grid, topological_grid_function);
		msc->SetBuildArcGeometry(Vec3b(false, false, false)); // we only need geometric realizations of 2saddle-max arcs
		msc->SetBuildArcs(Vec3b(true, false, true)); // dont need saddle-saddle
		msc->ComputeFromGrad();
		printf("### build msc from dgrad in:  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());

		
		now_time = std::chrono::steady_clock::now();
		std::vector<int> base_mins;
		std::vector<int> base_maxs;

		printf("Gathering base mins/maxs\n");
		MscType::LivingNodesIterator criticals(msc);
		for (criticals.begin(); criticals.valid(); criticals.advance()) {
			auto n = criticals.value(); // node id
			int dim = msc->getNode(n).dim;
			if (dim == 0) base_mins.push_back(n);
			if (dim == 3) base_maxs.push_back(n);
			crit_cells_2_nodes[msc->getNode(n).cellindex] = n;
		}
		DenseLabeling<int>* basins = new DenseLabeling<int>(underlying_grid->NumElements());

		printf("painting basins\n");
		auto num_mins = base_mins.size();
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < num_mins; i++) {
			auto n = base_mins[i];
			set<INDEX_TYPE> manifold;
			msc->fillGeometry(n, manifold, true);
			for (auto id : manifold) {
				// just paint the vertices
				if (topological_grid->dimension(id) == 0) {
					auto vid = topological_grid->VertexNumberFromCellID(id);
					basins->SetLabel(vid, n);
				}
			}
		}
		printf("### painted base basins using msc ascending manifolds in:  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());


		now_time = std::chrono::steady_clock::now();
		basins->SetAll(-1);
		FillAscendingIndexVolumeByVPaths(basins, topological_grid, discrete_gradient, msc, grid_function);
		printf("### painted base basins using path paintig in:  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());

		now_time = std::chrono::steady_clock::now();

		printf("painting mountains\n");
		DenseLabeling<int>* mounts = new DenseLabeling<int>(underlying_grid->NumElements());
		mounts->SetAll(-1);

		auto num_maxs = base_maxs.size();
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < num_maxs; i++) {
			auto n = base_maxs[i];
			set<INDEX_TYPE> manifold;
			msc->fillGeometry(n, manifold, false); // get descending manifold
			for (auto id : manifold) {
				// just paint the vertices
				if (topological_grid->dimension(id) == 3) {
					auto vid = topological_grid->VertexNumberFromCellID(id);
					mounts->SetLabel(vid, n);
				}
			}
		}
		printf("### painted mountains basins using msc ascending manifolds in:  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());

		
		now_time = std::chrono::steady_clock::now();

		printf("now building hierarchy\n");
		// want to switch persistence? first build a hierarchy with max persistence we ever want to browse to...
		
		// simplify to persistence and dump record to file
		char persname[2048];
		sprintf(persname, "%s.pers", argv[4]);
		msc->set_output_cancellation_records(persname);
		msc->ComputeHierarchy(persistence);
		printf("### built msc hierarchy in:  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());


		// now "browse" -- set the selection persistene, then repaint basins
		printf("starting interactive browse"); // could here try different persistence values
		// just arbitrarily show 4 different persistences
		for (float per = 0.0; per < persistence; per += persistence*.25) {
			DenseLabeling<int>* simp_basins = new DenseLabeling<int>(underlying_grid->NumElements());
			DenseLabeling<int>* simp_mounts = new DenseLabeling<int>(underlying_grid->NumElements());
			now_time = std::chrono::steady_clock::now();
			printf("\n\n");
			// START TIMINING FOR SWITCHING PERSISTENCE
			msc->SetSelectPersAbs(per);
			printf("make critical point -> representative map\n");

			// first get a node->node maps
			unordered_map<int, int> node_to_rep;
			MscType::LivingNodesIterator criticals(msc);
			for (criticals.begin(); criticals.valid(); criticals.advance()) {
				auto n = criticals.value(); // node id
				std::set<int> result;
				if (msc->getNode(n).dim == 0) {
					msc->GatherNodes(n, result, true);
					for (auto nid : result) {
						node_to_rep[nid] = n;
					}
				}
				else if (msc->getNode(n).dim == 3) {
					msc->GatherNodes(n, result, false);
					for (auto nid : result) {
						node_to_rep[nid] = n;
					}
				}
			}

			

			node_to_rep[-1] = -1; // add this in to handle out-of-range on mountains
			auto num = underlying_grid->NumElements();
			printf("repainting from base ids\n");
#pragma omp parallel for 
			for (long long i = 0; i < num; i++) {
				simp_basins->SetLabel(i, node_to_rep.at(basins->GetLabel(i)));
				simp_mounts->SetLabel(i, node_to_rep.at(mounts->GetLabel(i)));
			}
			printf("### re-painted basins AND mountains after pers change in:  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			now_time = std::chrono::steady_clock::now();
			// TIMING TO HERE IS WHAT IS NEEDED FOR BASIN/MOUNTAIN ID MAP
			vector<INDEX_TYPE> basin_seps;
#pragma omp parallel
			{
				int numthreads = omp_get_num_threads();
				int threadnum = omp_get_thread_num();
				auto num_cells = topological_grid->numCells();
				auto startid = threadnum * num_cells / numthreads;
				auto endid = (threadnum + 1) * num_cells / numthreads;
				if (threadnum == numthreads - 1) endid = num_cells;
				vector<INDEX_TYPE> local_seps;
				MeshType::FacetsIterator facets(topological_grid);
				// iterate over edges of mesh
				MeshType::DCellsIterator cells(topological_grid, 1, startid, endid);
				for (cells.begin(); cells.valid(); cells.advance()) {
					auto edgeid = cells.value();
					facets.begin(edgeid);// move facets iterator to first vertex
					auto v1 = facets.value();
					facets.advance();
					auto v2 = facets.value();
					if (simp_basins->GetLabel(topological_grid->VertexNumberFromCellID(v1)) !=
						simp_basins->GetLabel(topological_grid->VertexNumberFromCellID(v2)))
						local_seps.push_back(edgeid);

				}
#pragma omp critical
				{
					basin_seps.insert(basin_seps.end(), local_seps.begin(), local_seps.end());
				}
			}
			printf("found %llu basin separating edges\n", basin_seps.size());
			printf("### re-computed basins sep after pers change in:  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			now_time = std::chrono::steady_clock::now();
			vector<INDEX_TYPE> mountain_seps;
#pragma omp parallel
			{
				int numthreads = omp_get_num_threads();
				int threadnum = omp_get_thread_num();
				auto num_cells = topological_grid->numCells();
				auto startid = threadnum * num_cells / numthreads;
				auto endid = (threadnum + 1) * num_cells / numthreads;
				if (threadnum == numthreads - 1) endid = num_cells;
				vector<INDEX_TYPE> local_seps;
				MeshType::CofacetsIterator cofacets(topological_grid);
				// iterate over edges of mesh
				MeshType::DCellsIterator cells(topological_grid, 2, startid, endid);
				for (cells.begin(); cells.valid(); cells.advance()) {
					auto quadid = cells.value();
					cofacets.begin(quadid);// move quads cofacets iterator to first hex
					auto v1 = cofacets.value();
					cofacets.advance();
					if (!cofacets.valid()) continue; // we are on boundary of data
					auto v2 = cofacets.value();
					if (simp_mounts->GetLabel(topological_grid->VertexNumberFromCellID(v1)) !=
						simp_mounts->GetLabel(topological_grid->VertexNumberFromCellID(v2)))
						local_seps.push_back(quadid);

				}
#pragma omp critical
				{
					mountain_seps.insert(mountain_seps.end(), local_seps.begin(), local_seps.end());
				}
			}
			printf("found %llu mountain separating quads\n", mountain_seps.size());
		
			printf("### re-painted mountains sep after pers change in:  %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - now_time).count());
			


			printf("dumping to file\n");

			char fname[1024];
			sprintf(fname, "TEST_BASINS_%.4f.raw", per);
			simp_basins->OutputToIntFile(fname);
			sprintf(fname, "TEST_SMOUNTS_%.4f.raw", per);
			simp_mounts->OutputToIntFile(fname);
		}

		basins->OutputToIntFile("TESTBASINS.raw");
		mounts->OutputToIntFile("TESTMOUNTS.raw");		

		
		//msc->WriteComplex1Skeleton("testcomplex.bin");
		//msc->print_complex_info(true);
		//
		////sanity check arcs
		//unordered_set<INT_TYPE> living_arcs;
		//MscType::LivingArcsIterator liv_arcs_it(msc);
		//for (liv_arcs_it.begin(); liv_arcs_it.valid(); liv_arcs_it.advance()) {
		//	living_arcs.insert(liv_arcs_it.value());
		//}

		//char pairsname[1024];
		//sprintf(pairsname, "%s.pairs.txt", argv[4]);
		//FILE* f_pairs = fopen(pairsname, "w");
		//for (auto arcid : living_arcs) {
		//	auto& a = msc->getArc(arcid);
		//	if (a.dim != 2) continue;
		//	if (a.boundary) continue;
		//	float low_val = msc->getNode(a.lower).value;
		//	float high_val = msc->getNode(a.upper).value;
		//	fprintf(f_pairs, "%f %f\n", low_val, high_val);
		//}
		//fclose(f_pairs);

	}
	else {
		// set up structures to navigate grid, and load the 3d image
		underlying_grid = new GridType(GInt::Vec3i(X, Y, Z), GInt::Vec3b(0, 0, 0));

		// now set up an indexing scheme to use in a topological interpretation of the 
		// regular grid
		topological_grid = new MeshType(underlying_grid);


		// now compute the Morse-Smale complex
		MscType* msc = new MscType(NULL, topological_grid, NULL);
		msc->LoadComplex1Skeleton("testcomplex.bin");
		msc->print_complex_info(true);
		msc->ComputeHierarchy(persistence);

		MscType::LivingArcsIterator ait(msc);
		for (ait.begin(); ait.valid(); ait.advance()) {
			INT_TYPE aid = ait.value();
			vector<INDEX_TYPE> geom;
			msc->fillArcGeometry(aid, geom);
		}
	}



	return 0;

}


