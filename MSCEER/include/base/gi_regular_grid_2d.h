#ifndef REGULAR_GRID_H
#define REGULAR_GRID_H

#include "gi_basic_types.h"
#include "gi_vectors.h"


namespace GInt {
	// just a class to do some index stuff on a grid
	// e.g. get the 6 neighbors of a point
	// or get the 8 vertices surrounding a sample location.
	// x moves fastest then y then z

	class RegularGrid2D {
	protected:



		const Vec2l m_xy;			// the extents of the regular grid
		const Vec2d m_xy_d;		//  extents in doulbe formatting to avoid unnecessary type conversion during bounds check
		const Vec2b m_periodic;
		static const Vec2l kNeighborOffsets4[4];
		static const Vec2l kNeighborOffsets8[8];

	public:

		RegularGrid2D(Vec2l xy, Vec2b p) : m_xy(xy), m_xy_d(xy), m_periodic(p) {
			printf(" -- Created RegularGrid2D [%d %d] with periodicity [%d %d]\n",
				(int)m_xy[0], (int)m_xy[1], m_periodic[0], m_periodic[1] );
		}

		inline const Vec2l& XY() const { return m_xy; }
		inline const Vec2b& Periodic() const { return m_periodic; }
		inline INDEX_TYPE Index2d(const Vec2l& v) const {
			return v[0] + v[1] * m_xy[0];
		}
		inline Vec2l XY2d(INDEX_TYPE id) const {
			Vec2l res(id % m_xy[0], (id / m_xy[0]) );
			return res;
		}
		static inline INDEX_TYPE PositiveModulo(INDEX_TYPE a, INDEX_TYPE b) {
			return (a % b + b) % b;
		}
		static inline Vec2l PositiveModulo(const Vec2l& a, const Vec2l& b) {
			return (a % b + b) % b;
		}
		// return Inbounds version of the vertex - in case of negative indices
		inline Vec2l Inbounds(const Vec2l& v) const {
			Vec2l res = v;
			for (int i = 0; i < 2; i++) {
				if (m_periodic[i]) {
					res[i] = PositiveModulo(v[i], m_xy[i]);
				}
				else {
					if (res[i] < 0) res[i] = 0;
					if (res[i] >= m_xy[i]) res[i] = m_xy[i] - 1;
				}
			}
			return res;
		}



		// return Inbounds version of the vector - find the base point - towards -infinity, 
		inline Vec2d Inbounds(const Vec2d& v) const {


			Vec2d res = v;
			for (int i = 0; i < 2; i++) {
				if (m_periodic[i]) {
					if (res[i] < 0) res[i] += m_xy_d[i];
					if (res[i] >= m_xy_d[i]) res[i] -= m_xy_d[i];
				}
				else {
					if (res[i] < 0) res[i] = 0;
					if (res[i] >= m_xy_d[i]) res[i] = m_xy_d[i] - 1;
				}
			}
			return res;

			//Vec2l fl = v.IntFloor();
			//Vec3d diff = fl;
			//diff = v - fl;
			//Vec2l res = Inbounds(fl);
			//Vec3d resd = res;
			//resd = resd + diff;
			//return resd;
		}

		INDEX_TYPE NumElements() const {
			return m_xy[0] * m_xy[1] ;
		}

		static void PrintVector(const Vec3d& v) {
			printf("(%f, %f)\n", v[0], v[1]);
		}
		static void PrintVector(const Vec2l& v) {
			printf("(%d, %d)\n", v[0], v[1]);
		}

		//gather all existing neighbors, not restricted to the same boundary type
		int GatherExistingNeighborsAll4(const Vec2l& s, Vec2l* results) const {
			int numsofar = 0;
			if (m_periodic[0]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets4[0], m_xy);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets4[1], m_xy);
			}
			else {
				if (s[0] < m_xy[0] - 1) results[numsofar++] = s + kNeighborOffsets4[0];
				if (s[0] > 0) results[numsofar++] = s + kNeighborOffsets4[1];
			}
			if (m_periodic[1]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets4[2], m_xy);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets4[3], m_xy);
			}
			else {
				if (s[1] < m_xy[1] - 1) results[numsofar++] = s + kNeighborOffsets4[2];
				if (s[1] > 0) results[numsofar++] = s + kNeighborOffsets4[3];
			}
			return numsofar;
		}

		//gather all existing neighbors, not restricted to the same boundary type
		int GatherExistingNeighborsAll8(const Vec2l& s, Vec2l* results) const {
			int numsofar = 0;
			if (m_periodic[0]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets8[0], m_xy);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets8[1], m_xy);
			}
			else {
				for (int i = 0; i < 8; i++) {
					Vec2l t = s + kNeighborOffsets8[i];
					if (t[0] < m_xy[0] - 1 &&
						t[0] > 0 &&
						t[1] < m_xy[1] - 1 &&
						t[1] > 0 )
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
		int GatherExistingNeighborsSameBdry4(const Vec2l& s, Vec2l* results) const {
			int numsofar = 0;
			if (m_periodic[0]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets4[0], m_xy);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets4[1], m_xy);
			}
			else {
				if (s[0] < m_xy[0] - 1 && s[0] > 0){
					results[numsofar++] = s + kNeighborOffsets4[0];
					results[numsofar++] = s + kNeighborOffsets4[1];
				}
			}
			if (m_periodic[1]) {
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets4[2], m_xy);
				results[numsofar++] = PositiveModulo(s + kNeighborOffsets4[3], m_xy);
			}
			else {
				if (s[1] < m_xy[1] - 1 && s[1] > 0) {
					results[numsofar++] = s + kNeighborOffsets4[2];
					results[numsofar++] = s + kNeighborOffsets4[3];
				}
			}
			
			return numsofar;
		}

		// get neighbors along chosen axis, where kernelsize is number in each direction
		// results must be big enough to contain, = 2 * kernelsize + 1 (we include original point)
		// for non-m_periodic grids, the boundary point is replicated, e.g. the 2-kernel around point 0 is 0 0 0 1 2
		int Gather1DNeighborhood(const Vec2l& s, int axis, int kernelsize, Vec2l* results) const {
			int current = 0;
			if (m_periodic[axis]) {
				for (int i = -kernelsize; i <= kernelsize; i++) {
					Vec2l offset = s;
					offset[axis] = PositiveModulo(offset[axis] + i, m_xy[axis]);
					results[i + kernelsize] = offset;
				}
				return kernelsize;
			}
			else {
				INT_TYPE coord = s[axis];
				int dist = DistToBoundary(coord, axis);
				int newkernelsize = my_min(dist, kernelsize);
				for (int i = newkernelsize; i > 0; i--) {
					Vec2l offset = s;
					offset[axis] = offset[axis] - i;
					if (offset[axis] < 0) {
						offset[axis] = 0;
					}
					results[newkernelsize - i] = offset;
				}
				results[newkernelsize] = s;
				for (int i = 0; i < newkernelsize; i++) {
					Vec2l offset = s;
					offset[axis] = i + 1 + offset[axis];
					if (offset[axis] >= m_xy[axis]) {
						offset[axis] = m_xy[axis] - 1;
					}
					results[newkernelsize + i + 1] = offset;
				}
				return newkernelsize;
			}
		}


		int GatherSurrounding(const Vec2d& s, Vec2l* results) const {
			Vec2l basep = s.IntFloor();
			return GatherSurrounding(basep, results);
		}

		int my_min(const int& lhs, const int& rhs) const {
			return (lhs < rhs ? lhs : rhs);
		}
		int DistToBoundary(int coord, int axis) const {
			int b = (m_xy[axis] - 1) - coord;
			return my_min(coord, b);
		}
		// this gives minimum manhattan distance to a boundary 
		int DistToBoundary(const Vec2l& a) const {
			Vec2l b = (m_xy - 1) - a; // b stores distance to extents boundary
			// get closest distance to 0,0,0 planes
			INT_TYPE mina = (a[0] < a[1] ? a[0] : a[1]);

			// get closest distance to extents boundary
			INT_TYPE minb = (b[0] < b[1] ? b[0] : b[1]);

			// are we closer to which planes?
			return (mina < minb ? mina : minb);
		}

		// this function should only be used with extreme care! 
		// all bets are off if the queried point lies withing 1 cell of the boundary
		int GatherSurroundingNoBoundaryCheck(const Vec2d& s, Vec2l* results) const {
			Vec2l basep = s;
			return GatherSurroundingNoBoundaryCheck(basep, results);
		}
		int GatherSurroundingNoBoundaryCheck(const Vec2l& basep, Vec2l* results) const {

			//printf("GatherSurrounding::input = "); PrintVector(s);
			//PrintVector(basep);

			int xvecs[2];
			xvecs[0] = basep[0];
			xvecs[1] = basep[0] + 1;
			int yvecs[2];
			yvecs[0] = basep[1];
			yvecs[1] = basep[1] + 1;

			int numsofar = 0;

				for (int j = 0; j < 2; j++)
					for (int i = 0; i < 2; i++) {
				results[numsofar++] = Vec2l(xvecs[i], yvecs[j]);
					}
			return 4;
		}


		// THIS ASSUMES that we are REASONABLE, i.e. basep is not out of bounds
		INT_TYPE GatherSurrounding(const Vec2l& basep, Vec2l* results) const {

			//printf("GatherSurrounding::input = "); PrintVector(s);
			//PrintVector(basep);
			for (int i = 0; i < 2; i++) {
				if (basep[i] < 0 || basep[i] > m_xy[i] - 1) {
					printf("basep="); basep.PrintInt();;
				}
			}

			INT_TYPE xvecs[2];
			if (m_periodic[0]) {
				xvecs[0] = basep[0];
				xvecs[1] = (basep[0] == m_xy[0] - 1 ? 0 : basep[0] + 1); // logic here: if we are on periodic boundary, wrap around, else increment
			}
			else {
				xvecs[0] = basep[0]; xvecs[1] = 1 + basep[0];
				if (basep[0] < 0) xvecs[0] = 0;
				if (basep[0] >= m_xy[0] - 1) { xvecs[0] = xvecs[1] = m_xy[0] - 1; }
				else if (basep[0] >= m_xy[0]) { xvecs[1] = m_xy[0] - 1; }
			}
			INT_TYPE yvecs[2];
			if (m_periodic[1]) {
				yvecs[0] = basep[1];
				yvecs[1] = (basep[1] == m_xy[1] - 1 ? 0 : basep[1] + 1);
			}
			else {
				yvecs[0] = basep[1]; yvecs[1] = 1 + basep[1];
				if (basep[1] < 0) yvecs[0] = 0;
				if (basep[1] >= m_xy[1] - 1) { yvecs[0] = yvecs[1] = m_xy[1] - 1; }
				else if (basep[1] >= m_xy[1]) { yvecs[1] = m_xy[1] - 1; }
			}

	
			INT_TYPE numsofar = 0;
				for (int j = 0; j < 2; j++)
					for (int i = 0; i < 2; i++) {
				results[numsofar++] = Vec2l(xvecs[i], yvecs[j]);
					}
			return 4;
		}
		////int GatherSurrounding(const Vec3d& const s, INDEX_TYPE* results, double* factors) {
		////	Vec2l surrounding[8];
		////	GatherSurrounding(s, surrounding);

		////	return 8;
		////}
	};
}

#endif
