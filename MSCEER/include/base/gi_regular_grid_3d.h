#ifndef REGULAR_GRID_3D_H
#define REGULAR_GRID_3D_H

#include "gi_basic_types.h"
#include "gi_vectors.h"


namespace GInt {
	// just a class to do some index stuff on a grid
	// e.g. get the 6 neighbors of a point
	// or get the 8 vertices surrounding a sample location.
	// x moves fastest then y then z

	class RegularGrid3D {
	protected:



		const Vec3l m_xyz;			// the extents of the regular grid
		const Vec3d m_xyz_d;		//  extents in doulbe formatting to avoid unnecessary type conversion during bounds check
		const Vec3b m_periodic;
		static const Vec3l kNeighborOffsets6[6];
		static const Vec3l kNeighborOffsets26[26];

	public:

		RegularGrid3D(Vec3l xyz, Vec3b p) : m_xyz(xyz), m_xyz_d(xyz), m_periodic(p) {
            printf(" -- Created RegularGrid3D [%d %d %d] with periodicity [%d %d %d]\n",
                   (int) m_xyz[0], (int) m_xyz[1], (int) m_xyz[2], m_periodic[0], m_periodic[1], m_periodic[2]);
        }

		inline const Vec3l& XYZ() const { return m_xyz; }
		inline const Vec3b& Periodic() const { return m_periodic; }
		inline INDEX_TYPE Index3d(const Vec3l& v) const {
			return v[0] + v[1] * m_xyz[0] + v[2] * m_xyz[0] * m_xyz[1];
		}
		inline Vec3l XYZ3d(INDEX_TYPE id) const {
			Vec3l res(id % m_xyz[0], (id / m_xyz[0]) % m_xyz[1], id / (m_xyz[0] * m_xyz[1]));
			return res;
		}
		static inline INDEX_TYPE PositiveModulo(INDEX_TYPE a, INDEX_TYPE b) {
			return (a % b + b) % b;
		}
		static inline Vec3l PositiveModulo(const Vec3l& a, const Vec3l& b) {
			return (a % b + b) % b;
		}
		// return Inbounds version of the vertex - in case of negative indices
		inline Vec3l Inbounds(const Vec3l& v) const {
			Vec3l res = v;
			for (int i = 0; i < 3; i++) {
				if (m_periodic[i]) {
					res[i] = PositiveModulo(v[i], m_xyz[i]);
				}
				else {
					if (res[i] < 0) res[i] = 0;
					if (res[i] >= m_xyz[i]) res[i] = m_xyz[i] - 1;
				}
			}
			return res;
		}

	

		// return Inbounds version of the vector - find the base point - towards -infinity, 
		inline Vec3d Inbounds(const Vec3d& v) const {


			Vec3d res = v;
			for (int i = 0; i < 3; i++) {
				if (m_periodic[i]) {
					if (res[i] < 0) res[i] += m_xyz_d[i];
					if (res[i] >= m_xyz_d[i]) res[i] -= m_xyz_d[i];
				}
				else {
					if (res[i] < 0) res[i] = 0;
					if (res[i] >= m_xyz_d[i]) res[i] = m_xyz_d[i] - 1;
				}
			}
			return res;

			//Vec3l fl = v.IntFloor();
			//Vec3d diff = fl;
			//diff = v - fl;
			//Vec3l res = Inbounds(fl);
			//Vec3d resd = res;
			//resd = resd + diff;
			//return resd;
		}

		INDEX_TYPE NumElements() const {
			return m_xyz[0] * m_xyz[1] * m_xyz[2];
		}

		static void PrintVector(const Vec3d& v) {
			printf("(%f, %f, %f)\n", v[0], v[1], v[2]);
		}
		static void PrintVector(const Vec3l& v) {
			printf("(%d, %d, %d)\n", v[0], v[1], v[2]);
		}

		//gather all existing neighbors, not restricted to the same boundary type
		int GatherExistingNeighborsAll6(const Vec3l& s, Vec3l* results) const {
			int numsofar = 0;
			if (m_periodic[0]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[0], m_xyz);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[1], m_xyz);
			}
			else {
				if (s[0] < m_xyz[0] - 1) results[numsofar++] = s + kNeighborOffsets6[0];
				if (s[0] > 0) results[numsofar++] = s + kNeighborOffsets6[1];
			}
			if (m_periodic[1]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[2], m_xyz);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[3], m_xyz);
			}
			else {
				if (s[1] < m_xyz[1] - 1) results[numsofar++] = s + kNeighborOffsets6[2];
				if (s[1] > 0) results[numsofar++] = s + kNeighborOffsets6[3];
			}
			if (m_periodic[2]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[4], m_xyz);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[5], m_xyz);
			}
			else {
				if (s[2] < m_xyz[2] - 1) results[numsofar++] = s + kNeighborOffsets6[4];
				if (s[2] > 0) results[numsofar++] = s + kNeighborOffsets6[5];
			}
			return numsofar;
		}

		//gather all existing neighbors, not restricted to the same boundary type
		int GatherExistingNeighborsAll26(const Vec3l& s, Vec3l* results) const {
			int numsofar = 0;
			if (m_periodic[0]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[0], m_xyz);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[1], m_xyz);
			}
			else {
				for (int i = 0; i < 26; i++) {
					Vec3l t = s + kNeighborOffsets26[i];
					if (t[0] < m_xyz[0] - 1 &&
						t[0] > 0 &&
						t[1] < m_xyz[1] - 1 &&
						t[1] > 0 &&
						t[2] < m_xyz[2] - 1 &&
						t[2] > 0)
					{
						results[numsofar++] = t;
					}

				}
			}

			return numsofar;
		}
		// gather existing neighbors, but omit smaller boundary values. e.g. corners
		// of box will have no existing neighbors, but voxels 1 away from boundary will
		// have all six. voxels on boundary plane will have 4 neighbors... etc.
		int GatherExistingNeighborsSameBdry6(const Vec3l& s, Vec3l* results) const {
			int numsofar = 0;
			if (m_periodic[0]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[0], m_xyz);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[1], m_xyz);
			}
			else {
				if (s[0] < m_xyz[0] - 1 && s[0] > 0){
					results[numsofar++] = s + kNeighborOffsets6[0];
					results[numsofar++] = s + kNeighborOffsets6[1];
				}
			}
			if (m_periodic[1]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[2], m_xyz);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[3], m_xyz);
			}
			else {
				if (s[1] < m_xyz[1] - 1 && s[1] > 0) {
					results[numsofar++] = s + kNeighborOffsets6[2];
					results[numsofar++] = s + kNeighborOffsets6[3];
				}
			}
			if (m_periodic[2]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[4], m_xyz);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets6[5], m_xyz);
			}
			else {
				if (s[2] < m_xyz[2] - 1 && s[2] > 0) {
					results[numsofar++] = s + kNeighborOffsets6[4];
					results[numsofar++] = s + kNeighborOffsets6[5];
				}
			}
			return numsofar;
		}

		// get neighbors along chosen axis, where kernelsize is number in each direction
		// results must be big enough to contain, = 2 * kernelsize + 1 (we include original point)
		// for non-m_periodic grids, the boundary point is replicated, e.g. the 2-kernel around point 0 is 0 0 0 1 2
		int Gather1DNeighborhood(const Vec3l& s, int axis, int kernelsize, Vec3l* results) const {
			int current = 0;
			if (m_periodic[axis]) {
				for (int i = -kernelsize; i <= kernelsize; i++) {
					Vec3l offset = s;
					offset[axis] = PositiveModulo(offset[axis] + i, m_xyz[axis]);
					results[i + kernelsize] = offset;
				}
				return kernelsize;
			}
			else {
				INT_TYPE coord = s[axis];
				int dist = DistToBoundary(coord, axis);
				int newkernelsize = my_min(dist, kernelsize);
				for (int i = newkernelsize; i > 0; i--) {
					Vec3l offset = s;
					offset[axis] =offset[axis] - i ;
					if (offset[axis] < 0) {
						offset[axis] = 0;
					}
					results[newkernelsize - i] = offset;
				}
				results[newkernelsize] = s;
				for (int i = 0; i < newkernelsize; i++) {
					Vec3l offset = s;
					offset[axis] = i + 1 + offset[axis];
					if (offset[axis] >= m_xyz[axis]) {
						offset[axis] = m_xyz[axis] - 1;
					}
					results[newkernelsize + i + 1] = offset;
				}
				return newkernelsize;
			}
		}


        int GatherSurrounding(const Vec3d& s, Vec3l* results) const {
			Vec3l basep = s.IntFloor();
			return GatherSurrounding(basep, results);
		}

		int my_min(const int& lhs, const int& rhs) const {
			return (lhs < rhs ? lhs : rhs);
		}
		int DistToBoundary(int coord, int axis) const {
			int b = (m_xyz[axis] - 1) - coord;
			return my_min(coord, b);
		}
		// this gives minimum manhattan distance to a boundary 
		int DistToBoundary(const Vec3l& a) const {
			Vec3l b = (m_xyz - 1) - a; // b stores distance to extents boundary
			// get closest distance to 0,0,0 planes
			INT_TYPE min01a = (a[0] < a[1] ? a[0] : a[1]);
			INT_TYPE min12a = (a[1] < a[2] ? a[1] : a[2]);
			// get closest distance to extents boundary
			INT_TYPE min01b = (b[0] < b[1] ? b[0] : b[1]);
			INT_TYPE min12b = (b[1] < b[2] ? b[1] : b[2]);
			// are we closer to which planes?
			INT_TYPE mina = (min01a < min12a ? min01a : min12a);
			INT_TYPE minb = (min01b < min12b ? min01b : min12b);
			return (mina < minb ? mina : minb);
		}

		// this function should only be used with extreme care! 
		// all bets are off if the queried point lies withing 1 cell of the boundary
        int GatherSurroundingNoBoundaryCheck(const Vec3d& s, Vec3l* results) const {
			Vec3l basep = s;
			return GatherSurroundingNoBoundaryCheck(basep, results);
		}
        int GatherSurroundingNoBoundaryCheck(const Vec3l& basep, Vec3l* results) const {

			//printf("GatherSurrounding::input = "); PrintVector(s);
			//PrintVector(basep);

			int xvecs[2];
			xvecs[0] = basep[0];
			xvecs[1] = basep[0] + 1;
			int yvecs[2];
			yvecs[0] = basep[1];
			yvecs[1] = basep[1] + 1;

			int zvecs[2];
			zvecs[0] = basep[2];
			zvecs[1] = basep[2] + 1;

			int numsofar = 0;
			for (int k = 0; k < 2; k++)
				for (int j = 0; j < 2; j++)
					for (int i = 0; i < 2; i++) {
				results[numsofar++] = Vec3l(xvecs[i], yvecs[j], zvecs[k]);
					}
			return 8;
		}


		// THIS ASSUMES that we are REASONABLE, i.e. basep is not out of bounds
        INT_TYPE GatherSurrounding(const Vec3l& basep, Vec3l* results) const {

			//printf("GatherSurrounding::input = "); PrintVector(s);
			//PrintVector(basep);
			for (int i = 0; i < 3; i++) {
				if (basep[i] < 0 || basep[i] > m_xyz[i]-1) {
					printf("basep="); basep.PrintInt();;
				}
			}

			INT_TYPE xvecs[2];
			if (m_periodic[0]) {
				xvecs[0] = basep[0];
				xvecs[1] = (basep[0] == m_xyz[0]-1 ? 0 : basep[0]+1 ); // logic here: if we are on periodic boundary, wrap around, else increment
			}
			else {
				xvecs[0] = basep[0]; xvecs[1] = 1 + basep[0];
				if (basep[0] < 0) xvecs[0] = 0;
				if (basep[0] >= m_xyz[0] - 1) { xvecs[0] = xvecs[1] = m_xyz[0] - 1; }
				else if (basep[0] >= m_xyz[0]) { xvecs[1] = m_xyz[0] - 1; }
			}
			INT_TYPE yvecs[2];
			if (m_periodic[1]) {
				yvecs[0] = basep[1];
				yvecs[1] = (basep[1] == m_xyz[1]-1 ? 0 : basep[1] + 1); 
			}
			else {
				yvecs[0] = basep[1]; yvecs[1] = 1 + basep[1];
				if (basep[1] < 0) yvecs[0] = 0;
				if (basep[1] >= m_xyz[1] - 1) { yvecs[0] = yvecs[1] = m_xyz[1] - 1; }
				else if (basep[1] >= m_xyz[1]) { yvecs[1] = m_xyz[1] - 1; }
			}

			INT_TYPE zvecs[2];
			if (m_periodic[2]) {
				zvecs[0] = basep[2];
				zvecs[1] = (basep[2] == m_xyz[2]-1 ? 0 : basep[2] + 1); 
			}
			else {
				zvecs[0] = basep[2]; zvecs[1] = 1 + basep[2];
				if (basep[2] < 0) zvecs[0] = 0;
				if (basep[2] >= m_xyz[2] - 1) { zvecs[0] = zvecs[1] = m_xyz[2] - 1; }
				else if (basep[2] >= m_xyz[2]) { zvecs[1] = m_xyz[2] - 1; }
			}

			INT_TYPE numsofar = 0;
			for (int k = 0; k < 2; k++)
				for (int j = 0; j < 2; j++)
					for (int i = 0; i < 2; i++) {
				results[numsofar++] = Vec3l(xvecs[i], yvecs[j], zvecs[k]);
					}
			return 8;
		}
		////int GatherSurrounding(const Vec3d& const s, INDEX_TYPE* results, double* factors) {
		////	Vec3l surrounding[8];
		////	GatherSurrounding(s, surrounding);

		////	return 8;
		////}
	};


}

#endif
