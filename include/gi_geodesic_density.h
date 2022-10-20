#include <stdio.h>
#include <omp.h>
#include <queue>
#include <vector>
#include <string>
#include <utility>
#include <iostream>

#include "gi_vectors.h"

namespace GInt {
	class GeodesicDensityComputer {

	protected:
		// members for initialized during the init() call
		std::vector<GInt::Vec3i> m_box_neighbors;
		std::vector<int> m_box_neighbor_offsets;
		std::vector<float> m_distance_to_neighbor;

		int m_box_extent; // the length (integer) of each side of the box
		int m_box_middle; // the integer offset to the middle of a side
		int m_box_num_verts; // number of verts in a sliding box
		int m_max_steps;		// the integer maximum distance we can walk in a box
		float m_stopping_dist; // the floating point distance we actually want to integrate
		float m_lowest_definitely_mat_val;
		float m_lowest_passable_value;
		float m_oneover_max_minus_min;
		int m_box_start_id;					// the local ids for the middle will all be the same;

		// for setting up the bfs in the box, we need distance, id pairs
		typedef std::pair<float, int> distpair;
		struct compare {
			bool operator()(const distpair& lhs, const distpair& rhs) {
				return lhs > rhs;
			}
		};

		int to_box_index(const GInt::Vec3i& v) {
			return v[0] + v[1] * m_box_extent + v[2] * m_box_extent* m_box_extent;
		}

		enum { BLANK, ENQUEUED, SEEN };


		// return a value between 0 and 1
		// This scales the distance UP -- the retuned value is multiplied by the distance to get the "time"
		// it takes to go that distance 
		// so the "slowest" we want to go is (m_max_steps - 1) -- in the case where we are at the lowest passable value
		// at the "fastest", we have slowdown factor of 0 -- is this right? YES - we will add 1 to it.
		inline float slowdown_factor(float startvalue, float endvalue) {

			// do scaling
			float ave = 0.5f * (startvalue + endvalue); // first get the material value at the current spot
			if (ave > m_lowest_definitely_mat_val) return 0;
			//ave = (ave > m_lowest_definitely_mat_val ? m_lowest_definitely_mat_val : ave);
			return (m_stopping_dist - 1.0) * (1.0f - (ave - m_lowest_passable_value) * m_oneover_max_minus_min);

			//return 1.0;
		}
		// return value dependins on distance to voxel
		inline double accumulate_function(double distance) {
			return 1.0;
		}

		std::vector<std::pair<float, int>> m_all_dist_to_vol;
		GInt::Vec3i m_xyz;
		float* m_global_data;
	public:
		// the box xyz
		GeodesicDensityComputer(GInt::Vec3i xyz, float* global_data) : m_xyz(xyz), m_global_data(global_data) {
			m_is_initialized = false;
		}

		// now provide the functions to do  BFS integral
		bool m_is_initialized;
		void initialize(float max_distance, float lowest_definitely_mat_val, float lowest_passable_value) {

			m_stopping_dist = max_distance;
			m_lowest_definitely_mat_val = lowest_definitely_mat_val;
			m_lowest_passable_value = lowest_passable_value;
			m_max_steps = (int)m_stopping_dist + 1;
			m_oneover_max_minus_min = 1.0f / (m_lowest_definitely_mat_val - m_lowest_passable_value);

			printf("using m_max_steps %d size %d %d %d, mv=[%f--%f]\n", m_max_steps, m_xyz[0], m_xyz[1], m_xyz[1], m_lowest_definitely_mat_val, m_lowest_passable_value);


			long long global_size = ((long long)m_xyz[0]) * m_xyz[1] * m_xyz[2];
			//float* valfield = new float[X*Y*Z];
			//FILE* fin = fopen(argv[4], "rb");
			//fread(valfield, sizeof(float), X*Y*Z, fin);
			//fclose(fin);

			//printf("read data\n");
			//float* newdat = new float[X*Y*Z];
			//memset(newdat, 0, sizeof(float)*X*Y*Z);
			//printf("made data structures\n");

			// set all neighbor info
			m_box_extent = (m_max_steps * 2) + 1;
			m_box_middle = m_max_steps;
			m_box_num_verts = m_box_extent*m_box_extent*m_box_extent;
			for (int i = -1; i <= 1; i++) {
				for (int j = -1; j <= 1; j++) {
					for (int k = -1; k <= 1; k++) {
						if (i == 0 && j == 0 && k == 0) continue;
						GInt::Vec3i v(i, j, k);
						int noffset = to_box_index(v);
						m_box_neighbors.push_back(v);
						m_box_neighbor_offsets.push_back(noffset);
						m_distance_to_neighbor.push_back(((GInt::Vec3d) v).Mag());
					}
				}
			}
			m_box_start_id = to_box_index(GInt::Vec3i(m_box_middle, m_box_middle, m_box_middle)); // initialize m_box_start_id with address of local middle

			printf("set up neighborhoods\n");
			// now compute volumes an unrestricted sphere would ahve
			if (m_is_initialized) {

			}

			FillDistToVol(m_all_dist_to_vol);
		}


		void initialize(const std::vector<float>& distances, float lowest_definitely_mat_val, float lowest_passable_value) {
			initialize(distances[distances.size() - 1], lowest_definitely_mat_val, lowest_passable_value); // set things up wiht max sized box
		}

	protected:
		// fill the vector with the distance to volume pairs for a perfect ball - no slowdowns
		void FillDistToVol(std::vector<std::pair<float, int>>& dist_to_vol) {
			unsigned char* localmarked = new unsigned char[m_box_num_verts];
			memset(localmarked, BLANK, m_box_num_verts);
			int startid = m_box_start_id;
			int current_dist_index = 0;

			double sum = 0;
			std::vector<distpair> avec;
			avec.reserve(m_box_num_verts * 6);
			std::priority_queue<distpair, std::vector<distpair>, compare> mypq(compare(), avec);

			mypq.push(distpair(0.0f, startid));
			localmarked[startid] = ENQUEUED;
			int counter = 1;
			while (!mypq.empty()) {
				distpair p = mypq.top();
				mypq.pop();
				if (localmarked[p.second] == SEEN) continue;

				// check if we are past our distances
				dist_to_vol.push_back({ p.first, counter++ });

				sum += accumulate_function(p.first);

				localmarked[p.second] = SEEN;

				for (unsigned int i = 0; i < m_box_neighbor_offsets.size(); i++) {
					unsigned int nid = p.second + m_box_neighbor_offsets[i];

					if (nid > (m_box_num_verts - 1)) continue;

					// we are looking 

					if (localmarked[nid] != BLANK) continue;
					float tmp_curr_value = m_lowest_definitely_mat_val; // since we are checking a sphere 
					// get the slowdown factor
					float slowdownfactor = slowdown_factor(tmp_curr_value, tmp_curr_value);
					if (slowdownfactor != 0) printf("WHOA sdf = %f not 0\n", slowdownfactor);
					float newdist = p.first + m_distance_to_neighbor[i] * (1.0 + slowdownfactor);

					if (newdist > m_stopping_dist) continue;
					localmarked[nid] = ENQUEUED;
					mypq.push(distpair(newdist, nid));

				}

			} // end while loop
			dist_to_vol.push_back({ m_stopping_dist, counter - 1 });
		}


	public:
		struct DataBox {
			float* box_data;
			unsigned char* localmarked;
			std::priority_queue<distpair, std::vector<distpair>, compare>* mypq;

			DataBox(int num_verts) {
				box_data = new float[num_verts];
				localmarked = new unsigned char[num_verts];
				std::vector<distpair> container;
				container.reserve(num_verts * 6);
				mypq = new std::priority_queue<distpair, std::vector<distpair>, compare>(compare(), container);
			}
			DataBox() {
				box_data = NULL;
				localmarked = NULL;
				mypq = NULL;
			}

			~DataBox() {
				if (box_data != NULL) delete[] box_data;
				if (localmarked != NULL) delete[] localmarked;
				if (mypq != NULL) delete mypq;
			}
		};

	protected:
		// these only work with initialized data boxes

		// the geodesic volume for a ball of specified radius 
		float GeodesicVolume(float travel_time, DataBox* db) {

			int startid = m_box_start_id;
			int current_dist_index = 0;

			double sum = 0;

			db->mypq->push(distpair(0.0f, startid));
			db->localmarked[startid] = ENQUEUED;

			while (!db->mypq->empty()) {
				distpair p = db->mypq->top();
				db->mypq->pop();
				if (db->localmarked[p.second] == SEEN) continue;

				sum += accumulate_function(p.first);

				db->localmarked[p.second] = SEEN;

				for (unsigned int i = 0; i < m_box_neighbor_offsets.size(); i++) {
					unsigned int nid = p.second + m_box_neighbor_offsets[i];

					if (nid > (m_box_num_verts - 1)) continue;

					// check to see if we recurse on neighbor
					// first look at same connected component
					// also skip if it has been processed or enqueued
					float tmp_neg_value = db->box_data[nid];
					if (tmp_neg_value <= m_lowest_passable_value ||
						db->localmarked[nid] != BLANK) continue;	// ---------------- there might be error here - we are not re-enqueueing things that have been enqueued, but due to the slowdown function this might be a bad assumption
					float tmp_curr_value = db->box_data[p.second];
					// get the slowdown factor
					float slowdownfactor = slowdown_factor(tmp_curr_value, tmp_neg_value);

					float newdist = p.first + m_distance_to_neighbor[i] * (1.0 + slowdownfactor);

					if (newdist > travel_time) continue;
					db->localmarked[nid] = ENQUEUED;
					db->mypq->push(distpair(newdist, nid));

				}

			}
			return sum;
		}


		// the geodesic density for a ball reachable in time travel_time, going 1 voxel/unit for
		// definitely material values, with linear ramp to 1/max_time for background (i.e. no voxel is traversed)
		float GeodesicDensity(float travel_time, DataBox* db) {
			float volume = GeodesicVolume(travel_time, db);
			// record last distance seen
			//if (current_dist_index != g_num_outputs - 1) printf("this is wonky: %d, %d\n", current_dist_index, g_num_outputs); // -- actually, near boundary this can be lower
			int last = m_all_dist_to_vol.size() - 1;
			while (m_all_dist_to_vol[last].first > travel_time && last >= 0) {
				last--;
			}
			// the last one to be below the diastance, so the number of verts is good
			return (float)volume / (float)m_all_dist_to_vol[last].second;
		}

		// use the maximum allowable radius to compute the density
		// the geodesic density for a ball of specified radius 
		float GeodesicDensityMaxTime(DataBox* db) {
			float volume = GeodesicVolume(m_stopping_dist, db);
			// the last one to be below the diastance, so the number of verts is good
			return (float)volume / (float)m_all_dist_to_vol[m_all_dist_to_vol.size() - 1].second;
		}

		// returns false if there is no valid box here
		bool initialize_box(GInt::Vec3i position, DataBox* db) {
			INDEX_TYPE i = position[0];
			INDEX_TYPE j = position[1];
			INDEX_TYPE k = position[2];
			INDEX_TYPE X = m_xyz[0];
			INDEX_TYPE Y = m_xyz[1];
			INDEX_TYPE Z = m_xyz[2];
			INDEX_TYPE startid = i + j *X + k *X*Y;
			//newdat[startid] = -1;
			if (m_global_data[startid] <= m_lowest_passable_value) return false;
			memset(db->localmarked, BLANK, m_box_num_verts);

			int maxz = (k + m_max_steps >= Z ? (Z - 1) - k : m_max_steps);
			for (int z = (k - m_max_steps < 0 ? -k : -m_max_steps); z <= maxz; z++) {
				int liz = z + m_max_steps;
				int giz = z + k;
				int lzo = liz * m_box_extent * m_box_extent;
				INDEX_TYPE gzo = giz * X * Y;

				int maxy = (j + m_max_steps >= Y ? (Y - 1) - j : m_max_steps);
				for (int y = (j - m_max_steps < 0 ? -j : -m_max_steps); y <= maxy; y++) {
					int liy = y + m_max_steps;
					int giy = y + j;
					int lyo = liy * m_box_extent;
					INDEX_TYPE gyo = giy * X;

					int maxx = (i + m_max_steps >= X ? (X - 1) - i : m_max_steps);
					for (int x = (i - m_max_steps < 0 ? -i : -m_max_steps); x <= maxx; x++) {
						int lix = x + m_max_steps;
						int gix = x + i;

						int lid = lix + lyo + lzo;
						INDEX_TYPE gid = gix + gyo + gzo;
						db->box_data[lid] = m_global_data[gid];

					}
				}
			}
			while (!db->mypq->empty()) db->mypq->pop();
			return true;
		}

		// returns false if there is no valid box here
		bool initialize_box(GInt::Vec3i position, int tdist, DataBox* db) {
			INDEX_TYPE i = position[0];
			INDEX_TYPE j = position[1];
			INDEX_TYPE k = position[2];
			INDEX_TYPE X = m_xyz[0];
			INDEX_TYPE Y = m_xyz[1];
			INDEX_TYPE Z = m_xyz[2];
			INDEX_TYPE startid = i + j *X + k *X*Y;
			//newdat[startid] = -1;
			if (m_global_data[startid] <= m_lowest_passable_value) return false;
			memset(db->localmarked, BLANK, m_box_num_verts);

			// i j k are coordinates in global X Y Z
			// maxx maxy maxz are maximum offsets from i j k
			int maxz = (k + tdist >= Z ? (Z - 1) - k : tdist);
			for (int z = (k - tdist < 0 ? -k : -tdist); z <= maxz; z++) { // z goes from -dist to dist (unless outside global box)
				int liz = z + m_max_steps;	// m_max_steps is midpoint coordinate
				int giz = z + k;			// giz is global coordinate
				int lzo = liz * m_box_extent * m_box_extent; // the part of the index from the z
				INDEX_TYPE gzo = giz * X * Y;

				int maxy = (j + tdist >= Y ? (Y - 1) - j : tdist);
				for (int y = (j - tdist < 0 ? -j : -tdist); y <= maxy; y++) {
					int liy = y + m_max_steps;
					int giy = y + j;
					int lyo = liy * m_box_extent;
					INDEX_TYPE gyo = giy * X;

					int maxx = (i + tdist >= X ? (X - 1) - i : tdist);
					for (int x = (i - tdist < 0 ? -i : -tdist); x <= maxx; x++) {
						int lix = x + m_max_steps;
						int gix = x + i;

						int lid = lix + lyo + lzo;
						INDEX_TYPE gid = gix + gyo + gzo;
						db->box_data[lid] = m_global_data[gid];

					}
				}
			}
			while (!db->mypq->empty()) db->mypq->pop();
			return true;
		}
		// the size of the local computation box
		int box_size() { return m_box_num_verts; }

	public:

		DataBox* make_box() {
			return new DataBox(m_box_num_verts);
		}

		float GeodesicDensityMaxTime(GInt::Vec3i position, DataBox* db) {
			this->initialize_box(position, db);
			return GeodesicDensityMaxTime(db);
		}

		float GeodesicDensity(float travel_time, GInt::Vec3i position, DataBox* db) {
			//printf("-- Geodesic Density: tt=%f, pos=", travel_time);
			//position.PrintInt();
			int itime = ((int)travel_time) + 1;
			if (itime > m_max_steps) itime = m_max_steps;
			//printf("-- initializing with %d\n", itime);
			this->initialize_box(position, itime, db);
			//printf("-- computing geodesic density %f\n", travel_time);
			return GeodesicDensity(travel_time, db);
		}
		float GeodesicVolume(float travel_time, GInt::Vec3i position, DataBox* db) {
			int itime = ((int)travel_time) + 1;
			if (itime > m_max_steps) itime = m_max_steps;
			this->initialize_box(position, itime, db);
			return GeodesicVolume(travel_time, db);
		}






		// filename, x,y,z , sx, sy, sz, ex, ey, ez, sx, sy, sz
	//	int main(int argc, char** argv) {
	//
	//		if (argc < 9) {
	//			printf("usage: geodesic_weighted_integral <X> <Y> <Z> <i> <j> <k> <input_raw_field>  <free_pass_value> <min_material_value> <num_outputs> [thresholds] <max_dist> \n");
	//
	//			// printf("\tComputes the integral (sum of vertices reachable) within max_dist of each vertex\n");
	//			// printf("\twhere distance between v1 and v2 = euclidean distance * (1.0f - ( (f(v1) + f(v2) / 2)  - min_material_value)/ (free_pass_value - min_material_value)) \n");
	//
	//			return 1;
	//		}
	//
	//		unsigned long long X, Y, Z;
	//		int i, j, k;
	//
	//		printf("input data file: %s\n", argv[4]);
	//
	//		float aval;
	//
	//		sscanf(argv[1], "%llu", &X);
	//		sscanf(argv[2], "%llu", &Y);
	//		sscanf(argv[3], "%llu", &Z);
	//
	//		sscanf(argv[4], "%d", &i);
	//		sscanf(argv[5], "%d", &j);
	//		sscanf(argv[6], "%d", &k);
	//
	//
	//		sscanf(argv[7], "%f", &m_lowest_definitely_mat_val);
	//		sscanf(argv[8], "%f", &m_lowest_passable_value);
	//		sscanf(argv[9], "%f", &m_stopping_dist);
	//
	//		//m_stopping_dist = g_distvec[g_num_outputs - 1];
	//		m_max_steps = (int)m_stopping_dist + 1;
	//		m_oneover_max_minus_min = 1.0f / (m_lowest_definitely_mat_val - m_lowest_passable_value);
	//
	//
	//
	//		printf("using m_max_steps %d size %d %d %d, mv=[%f--%f]\n", m_max_steps, X, Y, Z, m_lowest_definitely_mat_val, m_lowest_passable_value);
	//
	//		long long global_size = X*Y*Z;
	//		float* valfield = new float[X*Y*Z];
	//		FILE* fin = fopen(argv[4], "rb");
	//		fread(valfield, sizeof(float), X*Y*Z, fin);
	//		fclose(fin);
	//
	//		printf("read data\n");
	//		float* newdat = new float[X*Y*Z];
	//		memset(newdat, 0, sizeof(float)*X*Y*Z);
	//		printf("made data structures\n");
	//		// set all neighbor info
	//		m_box_extent = (m_max_steps * 2) + 1;
	//		m_box_middle = m_max_steps;
	//		m_box_num_verts = m_box_extent*m_box_extent*m_box_extent;
	//		for (int i = -1; i <= 1; i++) {
	//			for (int j = -1; j <= 1; j++) {
	//				for (int k = -1; k <= 1; k++) {
	//					if (i == 0 && j == 0 && k == 0) continue;
	//					GInt::Vec3i v(i, j, k);
	//					int noffset = to_box_index(v);
	//					m_box_neighbors.push_back(v);
	//					m_box_neighbor_offsets.push_back(noffset);
	//					m_distance_to_neighbor.push_back(((GInt::Vec3d) v).Mag());
	//				}
	//			}
	//		}
	//		m_box_start_id = to_box_index(GInt::Vec3i(m_box_middle, m_box_middle, m_box_middle)); // initialize m_box_start_id with address of local middle
	//		printf("set up neighborhoods\n");
	//
	//		// now compute volumes an unrestricted sphere would ahve
	//		float* start_local_view = new float[m_box_num_verts];
	//		unsigned char* start_local_marked = new unsigned char[m_box_num_verts];
	//		memset(start_local_marked, BLANK, m_box_num_verts);
	//		for (int i = 0; i < m_box_num_verts; i++) {
	//			start_local_view[i] = m_lowest_definitely_mat_val;
	//		}
	//		g_dist_sphere_volumes_inv_vec = new float[g_num_outputs];
	//		for (int i = 0; i < g_num_outputs; i++) {
	//			g_dist_sphere_volumes_inv_vec[i] = 1.0; // set this to one -> this is a hack! makes so BFSIntegral returns simply the count
	//		}
	//		std::vector<distpair> avec;
	//		avec.reserve(m_box_num_verts);
	//		std::priority_queue<distpair, std::vector<distpair>, compare> apq(compare(), avec);
	//		float* tmpvec = new float[g_num_outputs];
	//		for (int i = 0; i < g_num_outputs; i++) {
	//			printf("distance %f volume = %f\n", g_distvec[i], g_dist_sphere_volumes_inv_vec[i]);
	//			//g_dist_sphere_volumes_inv_vec[i] = 1.0 / tmpvec[i];
	//		}
	//		BFSIntegral(tmpvec, start_local_view, start_local_marked, apq);
	//		for (int i = 0; i < g_num_outputs; i++) {
	//			printf("distance %f volume = %f\n", g_distvec[i], tmpvec[i]);
	//			g_dist_sphere_volumes_inv_vec[i] = 1.0 / tmpvec[i];
	//		}
	//		// prefill the volumes values for a fully filled blob
	//
	//
	//		//for (int i = 0; i < distances.size(); i++)  {
	//		//	m_box_neighbors[i].PrintInt();
	//		//	printf("%d %f\n", m_box_neighbor_offsets[i], distances[i]);
	//		//}
	//		printf("starting computation\n");
	//
	//
	//
	//
	//		//omp_set_num_threads(1);
	//#pragma omp parallel
	//		{
	//			int num_threads = omp_get_num_threads();
	//#pragma omp single
	//			{
	//				printf("using %d threads\n", num_threads);
	//			}
	//			int thread_num = omp_get_thread_num();
	//			float* local_view = new float[m_box_num_verts];
	//			unsigned char* local_marked = new unsigned char[m_box_num_verts];
	//			float* local_res_vec = new float[g_num_outputs];
	//
	//			std::vector<distpair> lvec;
	//			lvec.reserve(m_box_num_verts);
	//			std::priority_queue<distpair, std::vector<distpair>, compare> lpq(compare(), lvec);
	//			for (int k = thread_num; k < Z; k += num_threads) {
	//				printf("doing slice %d\n", k);
	//				for (int j = 0; j < Y; j++) {
	//					for (int i = 0; i < X; i++) {
	//
	//						//if (i != 202 || j != 115 || k != 188) continue;
	//
	//
	//						INDEX_TYPE startid = i + j *X + k *X*Y;
	//						//newdat[startid] = -1;
	//						if (valfield[startid] <= m_lowest_passable_value) continue;
	//						memset(local_marked, BLANK, m_box_num_verts);
	//
	//						int maxz = (k + m_max_steps >= Z ? (Z - 1) - k : m_max_steps);
	//						for (int z = (k - m_max_steps < 0 ? -k : -m_max_steps); z <= maxz; z++) {
	//							int liz = z + m_max_steps;
	//							int giz = z + k;
	//							int lzo = liz * m_box_extent * m_box_extent;
	//							INDEX_TYPE gzo = giz * X * Y;
	//
	//							int maxy = (j + m_max_steps >= Y ? (Y - 1) - j : m_max_steps);
	//							for (int y = (j - m_max_steps < 0 ? -j : -m_max_steps); y <= maxy; y++) {
	//								int liy = y + m_max_steps;
	//								int giy = y + j;
	//								int lyo = liy * m_box_extent;
	//								INDEX_TYPE gyo = giy * X;
	//
	//								int maxx = (i + m_max_steps >= X ? (X - 1) - i : m_max_steps);
	//								for (int x = (i - m_max_steps < 0 ? -i : -m_max_steps); x <= maxx; x++) {
	//									int lix = x + m_max_steps;
	//									int gix = x + i;
	//
	//									int lid = lix + lyo + lzo;
	//									INDEX_TYPE gid = gix + gyo + gzo;
	//									local_view[lid] = valfield[gid];
	//
	//								}
	//							}
	//						}
	//						//FILE* blockf = fopen("TESTBLOCK0.raw", "wb");
	//						//fwrite(local_view, sizeof(unsigned char), m_box_num_verts, blockf);
	//						//fclose(blockf);
	//						//printf("doing slize %d,%d,%d\n",i,j, k);
	//
	//						BFSIntegral(local_res_vec, local_view, local_marked, lpq);
	//
	//						for (int on = 0; on < g_num_outputs; on++) newdat[startid + on * global_size] = local_res_vec[on];
	//
	//						//blockf = fopen("TESTBLOCK.raw", "wb");
	//						//fwrite(local_view, sizeof(unsigned char), m_box_num_verts, blockf);
	//						//fclose(blockf);
	//
	//					}
	//				}
	//			}
	//		}
	//		//return 1;
	//		printf("saving output\n");
	//		for (int on = 0; on < g_num_outputs; on++) {
	//
	//
	//			FILE* fout;
	//			std::string outfile = std::string(argv[4]) + "." + std::to_string(g_distvec[on]) + "_" + std::to_string(X) + "x" + std::to_string(Y) + "x" + std::to_string(Z) + ".pancakes.raw";
	//			fout = fopen(outfile.c_str(), "wb");
	//			fwrite(&(newdat[on*global_size]), sizeof(float), X*Y*Z, fout);
	//			fclose(fout);
	//		}
	//		return 1;
	//	}
	};
}