/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef REGULAR_GRID_H
#define REGULAR_GRID_H

#include <assert.h>

#include "gi_basic_types.h"
#include "gi_vectors.h"


namespace GInt {
	// just a class to do some index stuff on a grid
	// e.g. get the 6 neighbors of a point
	// or get the 8 vertices surrounding a sample location.
	// x moves fastest then y then z

	class RegularGrid3D {
	protected:


		const Vec3l m_offsets;
		const Vec3l m_xyz;			// the extents of the regular grid
		const Vec3d m_xyz_d;		//  extents in doulbe formatting to avoid unnecessary type conversion during bounds check
		const Vec3b m_periodic;
		static const Vec3l kNeighborOffsets6[6];
		static const Vec3l kNeighborOffsets26[26];

	public:

		RegularGrid3D(Vec3l xyz, Vec3b p, bool verbose=false) : m_xyz(xyz), m_xyz_d(xyz), m_periodic(p), m_offsets(1, xyz[0], xyz[0] * xyz[1]) {
            if (verbose) printf(" -- Created RegularGrid3D [%d %d %d] with periodicity [%d %d %d]\n",
                   (int) m_xyz[0], (int) m_xyz[1], (int) m_xyz[2], m_periodic[0], m_periodic[1], m_periodic[2]);
        }

		inline const Vec3l& XYZ() const { return m_xyz; }
		inline const Vec3b& Periodic() const { return m_periodic; }
		inline INDEX_TYPE Index3d(const Vec3l& v) const {
			return v[0] + v[1] * m_xyz[0] + v[2] * m_xyz[0] * m_xyz[1];
		}

		inline INDEX_TYPE Offset_z() const {
			return m_offsets[2];
		}
		inline INDEX_TYPE Offset_y() const {
			return m_offsets[1];
		}
		inline INDEX_TYPE Offset_x() const {
			return 1;
		}

		inline Vec3l Coords(INDEX_TYPE id) const {
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

		// return true if v is inside the data box
		inline bool Inside(const Vec3d& v) const {

			for (int i = 0; i < 3; i++) {
				if (! m_periodic[i]) {
					if (v[i] < 0) return false;
					if (v[i] > m_xyz_d[i]-1) return false;
				}
			}
			return true;
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

		struct bdr_iterator {
			const Vec3l xyz;
			Vec3l pos;
			INDEX_TYPE stick_left;
			INDEX_TYPE internal_id; // this is the position in the yz plane
			const RegularGrid3D* mgrid;
			bdr_iterator(RegularGrid3D* grid) : xyz(grid->XYZ()), mgrid(grid) {}
			void begin() {
				pos = { 0,0,0 };
				stick_left = xyz[0];
				internal_id = 0;
			}
			void advance() {

			}
		};




		void GatherNeighborsNoBdryCheck6(INDEX_TYPE id, INDEX_TYPE* results) const {
			results[0] = id + m_offsets[0];
			results[1] = id - m_offsets[0];
			results[2] = id + m_offsets[1];
			results[3] = id - m_offsets[1];
			results[4] = id + m_offsets[2];
			results[5] = id - m_offsets[2];
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
		// NOTE: periodic domain is not supported
		int GatherExistingNeighborsAll26(const Vec3l& s, Vec3l* results) const {
			assert(!m_periodic[0] && !m_periodic[1] && !m_periodic[2]);
			int numsofar = 0;
			for (int i = 0; i < 26; i++) {
				Vec3l t = s + kNeighborOffsets26[i];

				if (t[0] < m_xyz[0] &&
					t[0] >= 0 &&
					t[1] < m_xyz[1] &&
					t[1] >= 0 &&
					t[2] < m_xyz[2] &&
					t[2] >= 0)
				{
					results[numsofar++] = t;
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
	}; // end RegularGrid2D



	
	//class RegularGrid3DDecomposition {
	//public:

	//	struct localExtents {
	//		localExtents(Vec3l start, Vec3l size, RegularGrid3D* grid) : starts(start), sizes(size), grid(grid) {}
	//		localExtents()  {}
	//		Vec3l starts;
	//		Vec3l sizes;
	//		//bool loaded;
	//		//dtype* local_data;
	//		const RegularGrid3D* grid;

	//		struct navigator {
	//			Vec3l start;
	//			INDEX_TYPE currentid;
	//			localExtents* le;
	//			bool navigate(Vec3l move) {
	//				currentid += move[0] + move[1] * le->grid->Offset_y() + move[2] * le->grid->Offset_z();
	//			}
	//			bool next_z() {
	//				currentid += le->grid->Offset_z();
	//			}

	//			INDEX_TYPE value() const {
	//				return currentid;
	//			}

	//		};

	//		navigator nav_begin(Vec3l start) {
	//			return{ start, this->grid->Index3d(start), this };
	//		}
	//	};

	//	localExtents m_global_grid;
	//	Vec3l m_desired_block_sizes;
	//	Vec3uc m_blocksize_pow2;
	//	localExtents* m_blocks;
	//	INDEX_TYPE m_num_blocks;
	//	Vec3l m_blocks_per_axis;
	//	INDEX_TYPE NumBlocks() const { return m_num_blocks; }

	//	bool isGlobalBoundary(int nx, int ny, int nz, localExtents& l) {
	//		int x = nx + l.starts[0] * 2;
	//		int y = ny + l.starts[1] * 2;
	//		int z = nz + l.starts[2] * 2;
	//		return x == 0 || x == ((XYZ()[0] * 2 - 1) - 1) ||
	//			y == 0 || y == ((XYZ()[1] * 2 - 1) - 1) ||
	//			z == 0 || z == ((XYZ()[2] * 2 - 1) - 1);
	//	}

	//	//virtual dtype* loadSubBlockFromFile(INT_TYPE block_id) {
	//	//	localExtents& l = blocks[block_id];
	//	//	// load values from file
	//	//	CELL_INDEX_TYPE num_data = l.sizes[0] * l.sizes[1] * l.sizes[2];
	//	//	//printf("reading %d data\n", num_data);
	//	//	dtype* local_data = new dtype[num_data];

	//	//	FILE* fin = fopen(mFilename, "rb");
	//	//	// seek to start
	//	//	CELL_INDEX_TYPE t_file_start = l.starts[2] * globalSizes.sizes[0] *
	//	//		globalSizes.sizes[1];

	//	//	fseek(fin, sizeof(dtype) * t_file_start, SEEK_SET);
	//	//	dtype* curr_loc = local_data;
	//	//	printf("get here 1 : %d\n", t_file_start);

	//	//	int readcount = 0;
	//	//	for (int z = 0; z < l.sizes[2]; z++) {
	//	//		fseek(fin, sizeof(dtype) * l.starts[1] * globalSizes.sizes[0], SEEK_CUR);
	//	//		for (int y = 0; y < l.sizes[1]; y++) {
	//	//			//skip empty y space

	//	//			//for(int x = 0; x < l.sizes[0]; x++) {
	//	//			// skip empty x space
	//	//			fseek(fin, sizeof(dtype) * l.starts[0], SEEK_CUR);
	//	//			fread(curr_loc, sizeof(dtype), l.sizes[0], fin);
	//	//			readcount += l.sizes[0];

	//	//			curr_loc += /*sizeof(dtype) **/ l.sizes[0];
	//	//			fseek(fin, sizeof(dtype) *
	//	//				(globalSizes.sizes[0] - (l.starts[0] + l.sizes[0])), SEEK_CUR);
	//	//			//}

	//	//		}
	//	//		fseek(fin, sizeof(dtype) * globalSizes.sizes[0] *
	//	//			(globalSizes.sizes[1] - (l.starts[1] + l.sizes[1])), SEEK_CUR);
	//	//	}
	//	//	fclose(fin);
	//	//	printf("load returned: read %d items\n", readcount);
	//	//	return local_data;
	//	//}


	//public:

	//	RegularGrid3DDecomposition(
	//		INDEX_TYPE global_X,
	//		INDEX_TYPE global_Y,
	//		INDEX_TYPE global_Z,
	//		INDEX_TYPE pow_2_size_X,
	//		INDEX_TYPE pow_2_size_Y,
	//		INDEX_TYPE pow_2_size_Z,
	//		RegularGrid3D* global_grid) :
	//		m_ghost(0, 0, 0),
	//		m_global_grid(Vec3l(0, 0, 0), Vec3l(global_X, global_Y, global_Z), global_grid)
	// {
	//		//m_global_grid.starts[0] = 0;
	//		//m_global_grid.starts[1] = 0;
	//		//m_global_grid.starts[2] = 0;
	//		//XYZ()[0] = global_X;
	//		//XYZ()[1] = global_Y;
	//		//XYZ()[2] = global_Z;
	//		SetDesiredSubBlockSizesPow2({ pow_2_size_X, pow_2_size_Y, pow_2_size_Z });
	//	}
	//	virtual ~RegularGrid3DDecomposition() {
	//		delete[] m_blocks;
	//		// memory leak since we don't deallocate 8 regular grids
	//	}

	//	void printBlockInfo(INT_TYPE block_id) {
	//		printf("%d:(%d, %d, %d)->(%d, %d, %d)\n",
	//			block_id,
	//			m_blocks[block_id].starts[0],
	//			m_blocks[block_id].starts[1],
	//			m_blocks[block_id].starts[2],
	//			m_blocks[block_id].starts[0] + m_blocks[block_id].sizes[0],
	//			m_blocks[block_id].starts[1] + m_blocks[block_id].sizes[1],
	//			m_blocks[block_id].starts[2] + m_blocks[block_id].sizes[2]
	//		);
	//	}


	//	void SetDesiredSubBlockSizesPow2(Vec3l blocksize) {
	//		m_blocksize_pow2 = blocksize;
	//		m_desired_block_sizes = { 1 << blocksize[0], 1 << blocksize[1], 1 << blocksize[2] };
	//	}

	//	// for block number at i,j,k (in terms of block indices) 
	//	inline INDEX_TYPE block_id_from_block_pos(Vec3l ijk) const {
	//		return ijk[0] + ijk[1] * m_blocks_per_axis[0] + ijk[2]  * m_blocks_per_axis[0] * m_blocks_per_axis[1];
	//	}
	//	// block ijk from global coords xyz
	//	inline Vec3l block_ijk_from_xyz(Vec3l xyz) const {
	//		return{ xyz[0] >> m_blocksize_pow2[0] ,xyz[1] >> m_blocksize_pow2[1], xyz[2] >> m_blocksize_pow2[2] };
	//	}
	//	// get the i of the ijk block position from x coordinate
	//	inline INDEX_TYPE block_i_from_x(INDEX_TYPE x) const {
	//		return x >> m_blocksize_pow2[0];
	//	}
	//	// get the j of the ijk block position from y coordinate
	//	inline INDEX_TYPE block_j_from_y(INDEX_TYPE y) const {
	//		return y >> m_blocksize_pow2[1];
	//	}
	//	// get the k of the ijk block position from z coordinate
	//	inline INDEX_TYPE block_k_from_z(INDEX_TYPE z) const {
	//		return z >> m_blocksize_pow2[2];
	//	}
	//	// the number of blocks per axis
	//	inline Vec3l blocks_per_axis() const {
	//		return m_blocks_per_axis;
	//	}

	//	// get the localextents for a block number 
	//	inline const localExtents& get_block(INDEX_TYPE block_num) const {
	//		return m_blocks[block_num];
	//	}

	//	// get block number conatining global position x,y,z
	//	inline INDEX_TYPE block_num(Vec3l xyz) const {
	//		return
	//			(xyz[0] >> m_blocksize_pow2[0]) +
	//			(xyz[1] >> m_blocksize_pow2[1]) * m_blocks_per_axis[0] +
	//			(xyz[2] >> m_blocksize_pow2[2]) * m_blocks_per_axis[0] * m_blocks_per_axis[1];
	//	}



	//	virtual void decompose() {
	//		// figure out number in each dimension
	//		for (int i = 0; i < 3; i++) {
	//			m_blocks_per_axis[i] = XYZ()[i] / m_desired_block_sizes[i] +
	//				(XYZ()[i] % m_desired_block_sizes[i] > 0 ? 1 : 0);
	//		}

	//		m_num_blocks = m_blocks_per_axis[0] * m_blocks_per_axis[1] * m_blocks_per_axis[2];
	//		m_blocks = new localExtents[m_num_blocks];


	//		// now simply create the list
	//		for (int z = 0; z < m_blocks_per_axis[2]; z++) {
	//			for (int y = 0; y < m_blocks_per_axis[1]; y++) {
	//				for (int x = 0; x < m_blocks_per_axis[0]; x++) {
	//					localExtents& l = m_blocks[block_id_from_block_pos({ x,y,z })];
	//					//l.loaded = false;

	//					l.starts[0] = x * m_desired_block_sizes[0];
	//					l.starts[1] = y * m_desired_block_sizes[1];
	//					l.starts[2] = z * m_desired_block_sizes[2];
	//					l.sizes[0] = m_desired_block_sizes[0] /*+ 1*/;
	//					if (x == m_blocks_per_axis[0] - 1)
	//						l.sizes[0] = XYZ()[0] - l.starts[0];

	//					l.sizes[1] = m_desired_block_sizes[1] /*+ 1*/;
	//					if (y == m_blocks_per_axis[1] - 1)
	//						l.sizes[1] = XYZ()[1] - l.starts[1];

	//					l.sizes[2] = m_desired_block_sizes[2] /*+ 1*/;
	//					if (z == m_blocks_per_axis[2] - 1)
	//						l.sizes[2] = XYZ()[2] - l.starts[2];
	//					//l.grid = new RegularGrid3D(l.sizes, { 0,0,0 });
	//				}
	//			}
	//		}
	//		// now make regular grid instances- for most we will reuse one constant grid size
	//		RegularGrid3D* interior_block = new RegularGrid3D(m_desired_block_sizes, { 0,0,0 });
	//		for (int z = 0; z < m_blocks_per_axis[2]-1; z++) {
	//			for (int y = 0; y < m_blocks_per_axis[1]-1; y++) {
	//				for (int x = 0; x < m_blocks_per_axis[0]-1; x++) {
	//					localExtents& l = m_blocks[block_id_from_block_pos({ x,y,z })];
	//					l.grid = interior_block;
	//				}
	//			}
	//		}
	//		// last x plane
	//		INDEX_TYPE last_x = m_blocks_per_axis[0] - 1;
	//		RegularGrid3D* last_x_plane = 
	//			new RegularGrid3D(m_blocks[block_id_from_block_pos({ last_x,0,0})].sizes, { 0,0,0 });
	//		for (int z = 0; z < m_blocks_per_axis[2] - 1; z++) {
	//			for (int y = 0; y < m_blocks_per_axis[1] - 1; y++) {
	//				localExtents& l = m_blocks[block_id_from_block_pos({ last_x,y,z })];
	//				l.grid = last_x_plane;
	//			}
	//		}
	//		// last y plane
	//		INDEX_TYPE last_y = m_blocks_per_axis[1] - 1;
	//		RegularGrid3D* last_y_plane =
	//			new RegularGrid3D(m_blocks[block_id_from_block_pos({ 0,last_y,0 })].sizes, { 0,0,0 });
	//		for (int z = 0; z < m_blocks_per_axis[2] - 1; z++) {
	//			for (int x = 0; x < m_blocks_per_axis[0] - 1; x++) {
	//				localExtents& l = m_blocks[block_id_from_block_pos({ x,last_y,z })];
	//				l.grid = last_y_plane;
	//			}
	//		}
	//		// last z plane
	//		INDEX_TYPE last_z = m_blocks_per_axis[2] - 1;
	//		RegularGrid3D* last_z_plane =
	//			new RegularGrid3D(m_blocks[block_id_from_block_pos({ 0,0,last_z })].sizes, { 0,0,0 });
	//		for (int y = 0; y < m_blocks_per_axis[1] - 1; y++) {
	//			for (int x = 0; x < m_blocks_per_axis[0] - 1; x++) {
	//				localExtents& l = m_blocks[block_id_from_block_pos({ x,y,last_z })];
	//				l.grid = last_z_plane;
	//			}
	//		}
	//		// last xy stick
	//		RegularGrid3D* last_xy_stick =
	//			new RegularGrid3D(m_blocks[block_id_from_block_pos({ last_x,last_y,0 })].sizes, { 0,0,0 });
	//		for (int z = 0; z < m_blocks_per_axis[2] - 1; z++) {
	//			localExtents& l = m_blocks[block_id_from_block_pos({ last_x,last_y,z })];
	//			l.grid = last_xy_stick;
	//		}
	//		// last yz stick
	//		RegularGrid3D* last_yz_stick =
	//			new RegularGrid3D(m_blocks[block_id_from_block_pos({ 0,last_y,last_z })].sizes, { 0,0,0 });
	//		for (int x = 0; x < m_blocks_per_axis[0] - 1; x++) {
	//			localExtents& l = m_blocks[block_id_from_block_pos({ x,last_y,last_z })];
	//			l.grid = last_yz_stick;
	//		}
	//		// last xz stick
	//		RegularGrid3D* last_xz_stick =
	//			new RegularGrid3D(m_blocks[block_id_from_block_pos({ last_x,0,last_z })].sizes, { 0,0,0 });
	//		for (int y = 0; y < m_blocks_per_axis[1] - 1; y++) {
	//			localExtents& l = m_blocks[block_id_from_block_pos({last_x,y,last_z })];
	//			l.grid = last_xz_stick;
	//		}
	//		// last xyz corner
	//		localExtents& l = m_blocks[block_id_from_block_pos({ last_x,last_y,last_z })];
	//		l.grid = new RegularGrid3D(m_blocks[block_id_from_block_pos({ last_x,last_y,last_z })].sizes, { 0,0,0 });
	//	}
	//	Vec3l m_ghost;
	//	virtual void AddFullGhostLayer() {

	//		m_ghost += {1, 1, 1};
	//		for (INDEX_TYPE i = 0; i < m_num_blocks; i++) {
	//			auto& b = m_blocks[i];
	//			// extend each block unless it is global boundary!
	//			for (int j = 0; j < 3; j++) {
	//				if (b.starts[j] != 0) {
	//					b.starts[j]--;
	//					b.sizes[j]++;
	//				}
	//				if (b.starts[j] + b.sizes[j] != XYZ()[j]) {
	//					b.sizes[j]++;
	//				}
	//			}
	//		}
	//	}

	//	//void setInputFileName(char* filename) {
	//	//	mFilename = filename;
	//	//}
	//	//void setOutputFileNames(char* grad_name, char* dat_name) {
	//	//	mOutGradFileName = grad_name;
	//	//	mOutDatFileName = dat_name;
	//	//}






	//	//void loadBlock(INT_TYPE block_id) {

	//	//	localExtents& l = blocks[block_id];
	//	//	if (l.loaded) return;

	//	//	// load values from file
	//	//	CELL_INDEX_TYPE num_data = l.sizes[0] * l.sizes[1] * l.sizes[2];
	//	//	dtype* local_data = this->loadSubBlockFromFile(block_id);
	//	//	mscSimplePointerDataHandler<dtype>* data_handler =
	//	//		new mscSimplePointerDataHandler<dtype>();
	//	//	data_handler->set_data(local_data);
	//	//	//data_handler->logify();

	//	//	printf("SIZES: %d %d %d\n", l.sizes[0], l.sizes[1], l.sizes[2]);
	//	//	mscRegularGrid3DImplicitMeshHandler* mesh_handler =
	//	//		new mscRegularGrid3DImplicitMeshHandler(l.sizes[0], l.sizes[1], l.sizes[2]);

	//	//	mscArrayFactory* array_factory = new mscArrayFactory(REGULAR_ARRAY);
	//	//	mscRegularGrid3DMeshFunction<dtype>* mesh_function =
	//	//		new mscRegularGrid3DMeshFunction<dtype>(data_handler, mesh_handler, array_factory);
	//	//	mesh_function->initialize();

	//	//	mscRegularGrid3DGradientField* gradient =
	//	//		new mscRegularGrid3DGradientField(mesh_handler, array_factory);


	//	//	l.data_handler = data_handler;
	//	//	l.mesh_function = mesh_function;
	//	//	l.mesh_handler = mesh_handler;
	//	//	l.grad = gradient;

	//	//	l.loaded = true;
	//	//	// now we read all sub-blocks
	//	//	// create the data function

	//	//}
	//	//void unloadBlock(INT_TYPE block_id) {
	//	//	printf("deleteing %d\n", block_id);
	//	//	localExtents& l = blocks[block_id];
	//	//	if (!l.loaded) return;


	//	//	delete l.data_handler;
	//	//	delete l.mesh_function;
	//	//	delete l.mesh_handler;
	//	//	delete l.grad;
	//	//	l.loaded = false;
	//	//	printf("done\n");

	//	//}


	//	//void writeOutputs(INT_TYPE block_id) {
	//	//	struct bitfield {
	//	//		unsigned char assigned : 1;
	//	//		unsigned char flag : 1;
	//	//		//unsigned char critical : 1;
	//	//		//unsigned char insorter : 1;
	//	//		//unsigned char dimA : 3;
	//	//		unsigned char pair : 3;
	//	//		unsigned char ldir : 3;
	//	//		//unsigned char empty_flag : 1;
	//	//	};
	//	//	localExtents& l = blocks[block_id];
	//	//	// load values from file
	//	//	CELL_INDEX_TYPE num_data =
	//	//		(l.sizes[0] * 2 - 1) *
	//	//		(l.sizes[1] * 2 - 1) *
	//	//		(l.sizes[2] * 2 - 1);

	//	//	FILE* fgrad = fopen(mOutGradFileName, "rb+");
	//	//	FILE* fdat = fopen(mOutDatFileName, "rb+");

	//	//	// seek to start
	//	//	CELL_INDEX_TYPE t_file_start =
	//	//		(l.starts[2] * 2) *
	//	//		(globalSizes.sizes[0] * 2 - 1) *
	//	//		(globalSizes.sizes[1] * 2 - 1);
	//	//	printf("z seek+ :%d\n", t_file_start);
	//	//	fseek(fgrad, sizeof(bitfield) * t_file_start, SEEK_SET);
	//	//	fseek(fdat, sizeof(dtype) * t_file_start, SEEK_SET);
	//	//	CELL_INDEX_TYPE current_cell = 0;

	//	//	CELL_INDEX_TYPE counter = 0;
	//	//	CELL_INDEX_TYPE offsetcounter = t_file_start;
	//	//	for (int z = 0; z < (l.sizes[2] * 2 - 1); z++) {
	//	//		fseek(fgrad, sizeof(bitfield) * (l.starts[1] * 2) * (globalSizes.sizes[0] * 2 - 1), SEEK_CUR);
	//	//		fseek(fdat, sizeof(dtype) * (l.starts[1] * 2) * (globalSizes.sizes[0] * 2 - 1), SEEK_CUR);
	//	//		offsetcounter += sizeof(dtype) * (l.starts[1] * 2) * (globalSizes.sizes[0] * 2 - 1);

	//	//		for (int y = 0; y < (l.sizes[1] * 2 - 1); y++) {
	//	//			//printf("y seek+ :%d\n",(l.starts[1]*2) * (globalSizes.sizes[0]*2-1));


	//	//			//for(int x = 0; x < l.sizes[0]; x++) {

	//	//			//printf("x seek+ :%d",(l.starts[0]*2));
	//	//			fseek(fgrad, sizeof(bitfield) * (l.starts[0] * 2), SEEK_CUR);
	//	//			fseek(fdat, sizeof(dtype) * (l.starts[0] * 2), SEEK_CUR);
	//	//			offsetcounter += sizeof(dtype) * (l.starts[0] * 2);


	//	//			//printf("writing %d\n", sizeof(bitfield) * (l.sizes[0]*2-1));
	//	//			for (int k = 0; k < (l.sizes[0] * 2 - 1); k++) {

	//	//				if (l.grad->get_assigned(current_cell) == 0) {
	//	//					printf("ERROR: not assigned\n");
	//	//				}

	//	//				bitfield temp;
	//	//				temp.assigned = l.grad->get_assigned(current_cell);
	//	//				temp.flag = isGlobalBoundary(k, y, z, l);

	//	//				if (l.grad->get_critical(current_cell)) {
	//	//					temp.pair = 7;
	//	//				}
	//	//				else {
	//	//					temp.pair =
	//	//						(unsigned char)((mscRegularGrid3DImplicitMeshHandler*)l.mesh_handler)->pair_2_offset(
	//	//							l.grad->get_pair(current_cell) - current_cell);
	//	//				}
	//	//				temp.ldir = (unsigned char)l.grad->get_dim_asc_man(current_cell);
	//	//				dtype val = l.mesh_function->cell_value(current_cell);

	//	//				//   	if (current_cell % 1000 == 0)
	//	//				//printf("cell[%d]:(%d, %d, %d, %d)\n", (int) current_cell,
	//	//				//temp.pair, 
	//	//				//temp.ldir,
	//	//				//temp.flag,
	//	//				//temp.assigned);


	//	//				fwrite(&temp, sizeof(bitfield), 1, fgrad);
	//	//				fwrite(&val, sizeof(dtype), 1, fdat);
	//	//				counter++;
	//	//				current_cell++;
	//	//				if (ferror(fdat)) printf("ERROR\n");
	//	//			}
	//	//			//printf("x seek- :%d",((globalSizes.sizes[0]*2-1) - (l.starts[0]*2 + (l.sizes[0]*2-1))));

	//	//			fseek(fgrad, sizeof(bitfield) *
	//	//				((globalSizes.sizes[0] * 2 - 1) - (l.starts[0] * 2 + (l.sizes[0] * 2 - 1))), SEEK_CUR);
	//	//			fseek(fdat, sizeof(dtype) *
	//	//				((globalSizes.sizes[0] * 2 - 1) - (l.starts[0] * 2 + (l.sizes[0] * 2 - 1))), SEEK_CUR);
	//	//			offsetcounter += sizeof(dtype) *
	//	//				((globalSizes.sizes[0] * 2 - 1) - (l.starts[0] * 2 + (l.sizes[0] * 2 - 1)));
	//	//			//}
	//	//			//printf("y seek- :%d\n",(globalSizes.sizes[0]*2-1) *
	//	//			//	((globalSizes.sizes[1]*2-1)- (l.starts[1]*2 + (l.sizes[1]*2-1))));


	//	//		}
	//	//		fseek(fgrad, sizeof(bitfield) * (globalSizes.sizes[0] * 2 - 1) *
	//	//			((globalSizes.sizes[1] * 2 - 1) - (l.starts[1] * 2 + (l.sizes[1] * 2 - 1))), SEEK_CUR);

	//	//		fseek(fdat, sizeof(dtype) * (globalSizes.sizes[0] * 2 - 1) *
	//	//			((globalSizes.sizes[1] * 2 - 1) - (l.starts[1] * 2 + (l.sizes[1] * 2 - 1))), SEEK_CUR);


	//	//		offsetcounter += sizeof(dtype) * (globalSizes.sizes[0] * 2 - 1) *
	//	//			((globalSizes.sizes[1] * 2 - 1) - (l.starts[1] * 2 + (l.sizes[1] * 2 - 1)));
	//	//	}
	//	//	printf("counter=%d    offsetcounter=%d\n", counter, offsetcounter);
	//	//	fclose(fdat);
	//	//	fclose(fgrad);
	//	//}

	//	//mscBasicMeshFunction<dtype>* getMeshFunction(INT_TYPE block_id) {
	//	//	return blocks[block_id].mesh_function;
	//	//}
	//	//mscBasicMeshHandler* getMeshHandler(INT_TYPE block_id) {
	//	//	return blocks[block_id].mesh_handler;
	//	//}
	//	//mscBasicGradientField* getGradientField(INT_TYPE block_id) {
	//	//	return blocks[block_id].grad;
	//	//}
	//	//mscBasicDataHandler<dtype>* getDataHandler(INT_TYPE block_id) {
	//	//	return blocks[block_id].data_handler;
	//	//}
	//	//void writeToDisk(char* filename) {}
	//};

	class RegularGrid3DDecompositionImplicit : public RegularGrid3D {
	public:



		typedef INDEX_TYPE BN_ID_PAIR;


		void PrintBNIDInfo(BN_ID_PAIR bnid) {
			INDEX_TYPE bn, id;
			GetBNAndIDFromPair(bnid, bn, id);
			printf("\n%lld: bn=%lld id=%lld\n",bnid, bn, id);
			printf("\tbn_pos:"); get_grid_of_blocks()->Coords(bn).PrintInt();
			printf("\tid_pos:"); get_block_grid(bn)->Coords(id).PrintInt();
		}
		inline INDEX_TYPE MakeBNPair(INDEX_TYPE block_num, INDEX_TYPE in_block_id) const {
			return (block_num << 32) | in_block_id;
		}

		inline void GetBNAndIDFromPair(BN_ID_PAIR bnid, INDEX_TYPE& bn, INDEX_TYPE& id) const {
			id = bnid & 0x00000000ffffffff;
			bn = bnid >> 32;
		}

		inline INDEX_TYPE GetGlobalIdFromBNID(BN_ID_PAIR bnid) {
			INDEX_TYPE bn, id; 
			GetBNAndIDFromPair(bnid, bn, id);
			Vec3l ijk = block_ijk_from_id(bn);
			Vec3l s;
			block_start(ijk, s);
			Vec3l c = get_block_grid(ijk)->Coords(id);
			auto xyz = s + c;
			return Index3d(xyz);
		}

		//struct InBlockNavigator {
		//	INDEX_TYPE block_num;
		//	INDEX_TYPE id;
		//	Vec3l b_ijk;
		//	Vec3l ijk;

		//	BN_ID_PAIR 
		//};

		struct localExtents {
			localExtents(Vec3l start, Vec3l size, RegularGrid3D* grid) : starts(start), sizes(size), grid(grid) {}
			localExtents() {}
			Vec3l starts;
			Vec3l sizes;
			//bool loaded;
			//dtype* local_data;
			const RegularGrid3D* grid;

			struct navigator {
				Vec3l start;
				INDEX_TYPE currentid;
				localExtents* le;
				void navigate(Vec3l move) {
					currentid += move[0] + move[1] * le->grid->Offset_y() + move[2] * le->grid->Offset_z();
				}
				void next_z() {
					currentid += le->grid->Offset_z();
				}

				INDEX_TYPE value() const {
					return currentid;
				}

			};

			navigator nav_begin(Vec3l start) {
				return{ start, this->grid->Index3d(start), this };
			}
		};
		struct in_block_iterator {
			const Vec3l xyz;	// dimensions of the local block
			const Vec3l b_pos;	// block coordinate position
			Vec3l pos;			// coordinate position in local block
			bool is_valid;
			const INDEX_TYPE block_num; // just block number
			RegularGrid3DDecompositionImplicit::BN_ID_PAIR in_block_id; // actually includes the block mask
									//INDEX_TYPE stick_left;
									//INDEX_TYPE internal_id; // this is the position in the yz plane
			const RegularGrid3D* mgrid;
			const RegularGrid3D* mgridofblocks;

			in_block_iterator(INDEX_TYPE bn, const RegularGrid3DDecompositionImplicit* decomposition) : 
				in_block_iterator(bn, decomposition->get_grid_of_blocks(), decomposition->get_block_grid(bn)) {	
			}

			in_block_iterator(INDEX_TYPE bn, const RegularGrid3D* grid_of_blocks, const RegularGrid3D* local_grid) :
				xyz(local_grid->XYZ()), b_pos(grid_of_blocks->Coords(bn)), block_num(bn), 
				mgrid(local_grid), mgridofblocks(grid_of_blocks) {}
			//in_block_iterator(INDEX_TYPE bn, const RegularGrid3DDecompositionImplicit* decomp) :
			void begin() {
				in_block_id = block_num << 32;
				pos = { 0,0,0 };
				is_valid = true;
			}
			void begin(BN_ID_PAIR start_bnidp) {
				in_block_id = start_bnidp;
				pos = mgrid->Coords(start_bnidp & 0x00000000ffffffff);
				is_valid = true;
			}
			void advance() {
				pos[0]++;
				in_block_id++;
				if (pos[0] == xyz[0]) {
					pos[0] = 0;
					pos[1]++;
					if (pos[1] == xyz[1]) {
						pos[1] = 0;
						pos[2]++;
						is_valid = pos[2] < xyz[2];
					}
				}
			}
			RegularGrid3DDecompositionImplicit::BN_ID_PAIR value() const {
				return in_block_id;
			}
			bool valid() const { return is_valid; }

		};

		struct neighborhood6_iterator {
			int num_neg;
			int m_id;
			RegularGrid3DDecompositionImplicit::BN_ID_PAIR results[6];
			const RegularGrid3DDecompositionImplicit* decomposition;
			neighborhood6_iterator(const RegularGrid3DDecompositionImplicit* decomp) : decomposition(decomp) {}

			void begin(BN_ID_PAIR start_bnidp) {
				INDEX_TYPE bn, id;
				decomposition->GetBNAndIDFromPair(start_bnidp, bn, id);
				in_block_iterator it(bn, decomposition);
				it.begin(start_bnidp);
				this->begin(it);
			}

			void begin(const in_block_iterator& base) {
				num_neg = 0;

				// check local block boundaries
				//// ------- X ---------
				auto gob = decomposition->get_grid_of_blocks();
				// upper boundaries X
				if (base.pos[0] < base.xyz[0] - 1) {
					// we are "internal" with respect to the upper boundary
					results[num_neg++] = base.in_block_id + 1; // base.mgrid->Offset_x();
				}
				else if (base.b_pos[0] < base.mgridofblocks->XYZ()[0] - 1) {					
					Vec3l neg_block_pos( base.b_pos[0] + 1, base.b_pos[1], base.b_pos[2] );
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(0, base.pos[1], base.pos[2]);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}
				// now minus direction
				if (base.pos[0] > 0) {
					// we are "internal" with respect to the lower boundary
					results[num_neg++] = base.in_block_id - 1;// base.mgrid->Offset_x();
				}
				else if (base.b_pos[0] > 0) {
					Vec3l neg_block_pos(base.b_pos[0] - 1, base.b_pos[1], base.b_pos[2]);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(bg->XYZ()[0]-1, base.pos[1], base.pos[2]);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos),	bg->Index3d(neg_id_pos));
				}
								

				//// ------- Y ---------
				// upper boundaries Y
				if (base.pos[1] < base.xyz[1] - 1) {
					// we are "internal" with respect to the upper boundary
					results[num_neg++] = base.in_block_id + base.mgrid->Offset_y();
				}
				else if (base.b_pos[1] < base.mgridofblocks->XYZ()[1] - 1) {
					Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1]+1, base.b_pos[2]);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(base.pos[0], 0, base.pos[2]);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}
				// now minus direction
				if (base.pos[1] > 0) {
					// we are "internal" with respect to the lower boundary
					results[num_neg++] = base.in_block_id - base.mgrid->Offset_y();
				}
				else if (base.b_pos[1] > 0) {
					Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1] - 1, base.b_pos[2]);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(base.pos[0], bg->XYZ()[1] - 1, base.pos[2]);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}


				//// ------- Z ---------
				// upper boundaries Z
				if (base.pos[2] < base.xyz[2] - 1) {
					// we are "internal" with respect to the upper boundary
					results[num_neg++] = base.in_block_id + base.mgrid->Offset_z();
				}
				else if (base.b_pos[2] < base.mgridofblocks->XYZ()[2] - 1) {
					Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1], base.b_pos[2] + 1);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(base.pos[0], base.pos[1], 0);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}
				// now minus direction
				if (base.pos[2] > 0) {
					// we are "internal" with respect to the lower boundary
					results[num_neg++] = base.in_block_id - base.mgrid->Offset_z();
				}
				else if (base.b_pos[2] > 0) {
					Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1], base.b_pos[2] - 1);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(base.pos[0], base.pos[1], bg->XYZ()[2] - 1);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}
				m_id = 0;
			}

			bool valid() {
				return m_id < num_neg;
			}

			void advance() {
				m_id++;
			}
			RegularGrid3DDecompositionImplicit::BN_ID_PAIR value() {
				return results[m_id];
			}

		};
		struct neighborhood27_iterator {
			int num_neg;
			int m_id;
			const RegularGrid3DDecompositionImplicit* decomposition;
			RegularGrid3DDecompositionImplicit::BN_ID_PAIR results[27];
			neighborhood27_iterator(const RegularGrid3DDecompositionImplicit* decomp) : decomposition(decomp) {}

			void begin(BN_ID_PAIR start_bnidp) {
				INDEX_TYPE bn, id;
				decomposition->GetBNAndIDFromPair(start_bnidp, bn, id);
				in_block_iterator it(bn, decomposition);
				it.begin(start_bnidp);
				this->begin(it);
			}

			void begin(const in_block_iterator& base) {

				if (base.pos[0] > 0 && base.pos[0] < base.xyz[0] - 1 &&
					base.pos[1] > 0 && base.pos[1] < base.xyz[1] - 1 &&
					base.pos[2] > 0 && base.pos[2] < base.xyz[2] - 1) {
					// internal
					base.mgrid->kNeighborOffsets26[0];
				}




				num_neg = 0;

				int offsets_count_x = 0;
				BN_ID_PAIR offsets_x[3];
				int offsets_count_y = 0;
				BN_ID_PAIR offsets_y[3];
				int offsets_count_z = 0;
				BN_ID_PAIR offsets_z[3];


				// check local block boundaries
				//// ------- X ---------
				auto gob = decomposition->get_grid_of_blocks();
				
				/// X OFFSETS
				{
					// upper boundaries X
					if (base.pos[0] < base.xyz[0] - 1) {
						// we are "internal" with respect to the upper boundary
						offsets_x[offsets_count_x++] = 1; // base.mgrid->Offset_x();
					}
					else if (base.b_pos[0] < base.mgridofblocks->XYZ()[0] - 1) {
						Vec3l neg_block_pos(base.b_pos[0] + 1, base.b_pos[1], base.b_pos[2]);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(0, base.pos[1], base.pos[2]);
						offsets_x[offsets_count_x++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos))- base.in_block_id;
					}
					// now zero direction
					offsets_x[offsets_count_x++] = 0;
					// now minus direction
					if (base.pos[0] > 0) {
						// we are "internal" with respect to the lower boundary
						offsets_x[offsets_count_x++] = - 1;// base.mgrid->Offset_x();
					}
					else if (base.b_pos[0] > 0) {
						Vec3l neg_block_pos(base.b_pos[0] - 1, base.b_pos[1], base.b_pos[2]);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(bg->XYZ()[0] - 1, base.pos[1], base.pos[2]);
						offsets_x[offsets_count_x++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
				}


				//// ------- Y ---------
				{
					// upper boundaries Y
					if (base.pos[1] < base.xyz[1] - 1) {
						// we are "internal" with respect to the upper boundary
						offsets_y[offsets_count_y++] =  base.mgrid->Offset_y();
					}
					else if (base.b_pos[1] < base.mgridofblocks->XYZ()[1] - 1) {
						Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1] + 1, base.b_pos[2]);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(base.pos[0], 0, base.pos[2]);
						offsets_y[offsets_count_y++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
					// no offset direction
					offsets_y[offsets_count_y++] = 0;

					// now minus direction
					if (base.pos[1] > 0) {
						// we are "internal" with respect to the lower boundary
						offsets_y[offsets_count_y++] = -base.mgrid->Offset_y();
					}
					else if (base.b_pos[1] > 0) {
						Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1] - 1, base.b_pos[2]);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(base.pos[0], bg->XYZ()[1] - 1, base.pos[2]);
						offsets_y[offsets_count_y++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
				}

				//// ------- Z ---------
				{
					// upper boundaries Z
					if (base.pos[2] < base.xyz[2] - 1) {
						// we are "internal" with respect to the upper boundary
						offsets_z[offsets_count_z++] = base.mgrid->Offset_z();
					}
					else if (base.b_pos[2] < base.mgridofblocks->XYZ()[2] - 1) {
						Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1], base.b_pos[2] + 1);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(base.pos[0], base.pos[1], 0);
						offsets_z[offsets_count_z++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
					// no z offset
					offsets_z[offsets_count_z++] = 0;
					// now minus direction
					if (base.pos[2] > 0) {
						// we are "internal" with respect to the lower boundary
						offsets_z[offsets_count_z++] =  -base.mgrid->Offset_z();
					}
					else if (base.b_pos[2] > 0) {
						Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1], base.b_pos[2] - 1);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(base.pos[0], base.pos[1], bg->XYZ()[2] - 1);
						offsets_z[offsets_count_z++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
				}

				num_neg = 0;
				for (int z = 0; z < offsets_count_z; z++)
					for (int y = 0; y < offsets_count_y; y++)
						for (int x = 0; x < offsets_count_x; x++) {
							results[num_neg++] = base.in_block_id + offsets_x[x] + offsets_y[y] + offsets_z[z];
				}
					
				m_id = 0;
			}

			bool valid() {
				return m_id < num_neg;
			}

			void advance() {
				m_id++;
			}
			RegularGrid3DDecompositionImplicit::BN_ID_PAIR value() {
				return results[m_id];
			}

		};
		// this is a comment on the code stuff
		// this is another line of code
		//localExtents m_global_grid;
		Vec3l m_basic_block_dims;
		Vec3uc m_blocksize_pow2;
		RegularGrid3D* m_block_grids[8];
		RegularGrid3D* m_grid_of_blocks;
		INDEX_TYPE m_num_blocks;
		Vec3l m_blocks_per_axis;
		INDEX_TYPE NumBlocks() const { return m_num_blocks; }

		bool isGlobalBoundary(int nx, int ny, int nz, localExtents& l) {
			int x = nx + l.starts[0] * 2;
			int y = ny + l.starts[1] * 2;
			int z = nz + l.starts[2] * 2;
			return x == 0 || x == ((XYZ()[0] * 2 - 1) - 1) ||
				y == 0 || y == ((XYZ()[1] * 2 - 1) - 1) ||
				z == 0 || z == ((XYZ()[2] * 2 - 1) - 1);
		}

	public:

		RegularGrid3DDecompositionImplicit(
			INDEX_TYPE global_X,
			INDEX_TYPE global_Y,
			INDEX_TYPE global_Z,
			INDEX_TYPE pow_2_size_X,
			INDEX_TYPE pow_2_size_Y,
			INDEX_TYPE pow_2_size_Z) :
			RegularGrid3D(Vec3l( global_X, global_Y, global_Z ), Vec3b(0, 0, 0))
		{
			SetDesiredSubBlockSizesPow2({ pow_2_size_X, pow_2_size_Y, pow_2_size_Z });
		}

		virtual ~RegularGrid3DDecompositionImplicit() {
			for (int i = 0; i < 8; i++) delete m_block_grids[i];
		}

		void printBlockInfo(INT_TYPE block_id) {
		}


		void SetDesiredSubBlockSizesPow2(Vec3l blocksize) {
			m_blocksize_pow2 = blocksize;
			m_basic_block_dims = { 1 << blocksize[0], 1 << blocksize[1], 1 << blocksize[2] };
		}

		// for block number at i,j,k (in terms of block indices) 
		inline INDEX_TYPE block_id_from_block_pos(Vec3l ijk) const {
			return ijk[0] + ijk[1] * m_blocks_per_axis[0] + ijk[2] * m_blocks_per_axis[0] * m_blocks_per_axis[1];
		}
		// block ijk from global coords xyz
		inline Vec3l block_ijk_from_xyz(Vec3l xyz) const {
			return{ xyz[0] >> m_blocksize_pow2[0] ,xyz[1] >> m_blocksize_pow2[1], xyz[2] >> m_blocksize_pow2[2] };
		}
		// get the i of the ijk block position from x coordinate
		inline INDEX_TYPE block_i_from_x(INDEX_TYPE x) const {
			return x >> m_blocksize_pow2[0];
		}
		// get the j of the ijk block position from y coordinate
		inline INDEX_TYPE block_j_from_y(INDEX_TYPE y) const {
			return y >> m_blocksize_pow2[1];
		}
		// get the k of the ijk block position from z coordinate
		inline INDEX_TYPE block_k_from_z(INDEX_TYPE z) const {
			return z >> m_blocksize_pow2[2];
		}
		// the number of blocks per axis
		inline Vec3l blocks_per_axis() const {
			return m_blocks_per_axis;
		}

protected:

		inline unsigned char get_ijk_block_type(const Vec3l& ijk) const {
			unsigned char val_i = ijk[0] == (m_blocks_per_axis[0] - 1);
			unsigned char val_j = (ijk[1] == (m_blocks_per_axis[1] - 1)) << 1;
			unsigned char val_k = (ijk[2] == (m_blocks_per_axis[2] - 1)) << 2;
			return val_i | val_j | val_k;
		}
public:
		// get the localextents for a block number 
		inline  void get_block(const Vec3l& ijk, localExtents& res) const {
			block_start(ijk, res.starts);
			block_size(ijk, res.starts, res.sizes);
			res.grid = m_block_grids[get_ijk_block_type(ijk)];
		}
		// the grid associated with block at block pos ijk
		inline  const RegularGrid3D* get_block_grid(const Vec3l& ijk) const {
			return m_block_grids[get_ijk_block_type(ijk)];
		}

		inline const RegularGrid3D* get_grid_of_blocks() const {
			return m_grid_of_blocks;
		}
		// the grid associated with block number
		inline  const RegularGrid3D* get_block_grid(INDEX_TYPE block_num) const {
			return get_block_grid(this->block_ijk_from_id(block_num));
		}
		// block ijk from BLOCK id - not global id
		inline Vec3l block_ijk_from_id(INDEX_TYPE block_num) const {
			INDEX_TYPE slab = (m_blocks_per_axis[0] * m_blocks_per_axis[1]);
			return Vec3l(block_num % m_blocks_per_axis[0],
				(block_num % slab) / m_blocks_per_axis[0],
				block_num / slab);
		}

		// global id to bnid: this is slow
		inline BN_ID_PAIR get_bnid_from_id(INDEX_TYPE id) {
			auto xyz = this->Coords(id);
			auto ijk = block_ijk_from_xyz(xyz);
			Vec3l bs;
			block_start(ijk, bs);
			auto in_block_xyz = xyz - bs;
			INDEX_TYPE inblockid = get_block_grid(ijk)->Index3d(in_block_xyz);			
			return MakeBNPair(block_id_from_block_pos(ijk), inblockid);
		}

		inline void block_start(const Vec3l& ijk, Vec3l& starts) const {
			starts[0] = ijk[0] * m_basic_block_dims[0];
			starts[1] = ijk[1] * m_basic_block_dims[1];
			starts[2] = ijk[2] * m_basic_block_dims[2];
		}
		inline void block_size(const Vec3l& ijk, const Vec3l& starts, Vec3l& sizes) const {
			sizes = m_basic_block_dims;
			if (ijk[0] == m_blocks_per_axis[0] - 1)
				sizes[0] = XYZ()[0] - starts[0];
			sizes[1] = m_basic_block_dims[1] /*+ 1*/;
			if (ijk[1] == m_blocks_per_axis[1] - 1)
				sizes[1] = XYZ()[1] - starts[1];
			sizes[2] = m_basic_block_dims[2] /*+ 1*/;
			if (ijk[2] == m_blocks_per_axis[2] - 1)
				sizes[2] = XYZ()[2] - starts[2];
		}
		inline void block_size(const Vec3l& ijk, Vec3l& sizes) const {
			Vec3l starts;
			block_start(ijk, starts);
			block_size(ijk, starts, sizes);
		}

		// get block number conatining global position x,y,z
		inline INDEX_TYPE block_num(Vec3l xyz) const {
			return
				(xyz[0] >> m_blocksize_pow2[0]) +
				(xyz[1] >> m_blocksize_pow2[1]) * m_blocks_per_axis[0] +
				(xyz[2] >> m_blocksize_pow2[2]) * m_blocks_per_axis[0] * m_blocks_per_axis[1];
		}



		virtual void decompose() {
			// figure out number in each dimension
			for (int i = 0; i < 3; i++) {
				m_blocks_per_axis[i] = XYZ()[i] / m_basic_block_dims[i] +
					(XYZ()[i] % m_basic_block_dims[i] > 0 ? 1 : 0);
			}
			m_grid_of_blocks = new RegularGrid3D(m_blocks_per_axis, { 0,0,0 });
			m_num_blocks = m_grid_of_blocks->NumElements(); 

			// now make regular grid instances- for most we will reuse one constant grid size
			Vec3l last_sizes;		
			// last x plane
			for (INDEX_TYPE i = 0; i < m_blocks_per_axis[0]; i += m_blocks_per_axis[0] - 1) {
				for (INDEX_TYPE j = 0; j < m_blocks_per_axis[1]; j += m_blocks_per_axis[1] - 1) {
					for (INDEX_TYPE k = 0; k < m_blocks_per_axis[2]; k += m_blocks_per_axis[2] - 1) {
						Vec3l pos = { i,j,k };
						block_size(pos, last_sizes);
						m_block_grids[get_ijk_block_type(pos)] = new RegularGrid3D(last_sizes, { 0,0,0 });
					}
				}
			}

		}


	};


	// block size power of 2
//#define BSP2_X 3
//#define BSP2_Y 3
//#define BSP2_Z 3
	template<unsigned char BSP2_X, unsigned char BSP2_Y, unsigned char BSP2_Z>
	class RegularGrid3DDecompositionFixed : public RegularGrid3D {
	public:



		typedef INDEX_TYPE BN_ID_PAIR;


		void PrintBNIDInfo(BN_ID_PAIR bnid) {
			INDEX_TYPE bn, id;
			GetBNAndIDFromPair(bnid, bn, id);
			printf("\n%lld: bn=%lld id=%lld\n", bnid, bn, id);
			printf("\tbn_pos:"); get_grid_of_blocks()->Coords(bn).PrintInt();
			printf("\tid_pos:"); get_block_grid(bn)->Coords(id).PrintInt();
		}
		inline INDEX_TYPE MakeBNPair(INDEX_TYPE block_num, INDEX_TYPE in_block_id) const {
			return (block_num << (BSP2_X+BSP2_Y+BSP2_Z)) | in_block_id;
		}

		inline INDEX_TYPE GetBNFromPair(BN_ID_PAIR bnid) const {
			return bnid >> (BSP2_X + BSP2_Y + BSP2_Z);
		}

		inline INDEX_TYPE GetIDFromPair(BN_ID_PAIR bnid) const {
			return bnid & ((1 << (BSP2_X + BSP2_Y + BSP2_Z)) - 1);
		}

		inline void GetBNAndIDFromPair(BN_ID_PAIR bnid, INDEX_TYPE& bn, INDEX_TYPE& id) const {
			id = GetIDFromPair(bnid);
			bn = GetBNFromPair(bnid);
		}

		inline INDEX_TYPE GetGlobalIdFromBNID(BN_ID_PAIR bnid) {
			INDEX_TYPE bn, id;
			GetBNAndIDFromPair(bnid, bn, id);
			Vec3l ijk = block_ijk_from_id(bn);
			Vec3l s;
			block_start(ijk, s);
			Vec3l c = get_block_grid(ijk)->Coords(id);
			auto b_extents = s + c;
			return Index3d(b_extents);
		}


		struct localExtents {
			localExtents(Vec3l start, Vec3l size, RegularGrid3D* grid) : starts(start), sizes(size), grid(grid) {}
			localExtents() {}
			Vec3l starts;
			Vec3l sizes;
			//bool loaded;
			//dtype* local_data;
			const RegularGrid3D* grid;

			struct navigator {
				Vec3l start;
				INDEX_TYPE currentid;
				localExtents* le;
				void navigate(Vec3l move) {
					currentid += move[0] + move[1] * le->grid->Offset_y() + move[2] * le->grid->Offset_z();
				}
				void next_z() {
					currentid += le->grid->Offset_z();
				}

				INDEX_TYPE value() const {
					return currentid;
				}

			};

			navigator nav_begin(Vec3l start) {
				return{ start, this->grid->Index3d(start), this };
			}
		};

	//protected:
	//	template<bool last_x, bool last_y, bool last_z>
	//	struct _in_block_iterator {

	//	};


	public:


		struct in_block_iterator {
			// ORDER MATTERS HERE FOR INITIALIZERS
			// CONTEXT ABOUT THE BLOCK
			const INDEX_TYPE block_num; // just block number
			const Vec3l b_pos;	// block coordinate position
			const Vec3l b_extents;	// extents of valid indices within local block - NOT size of block
			const Vec3l b_jumps; // the difference between the fixed block size and actula extent of the block
			const BN_ID_PAIR last_id; // the last in BNID to be valid
			// VARIABLES FOR ITERATION
			Vec3l pos;			// coordinate position in local block
			bool is_valid;
			BN_ID_PAIR in_block_id; // actually includes the block mask
																		//INDEX_TYPE stick_left;
			const RegularGrid3DDecompositionFixed<BSP2_X, BSP2_Y, BSP2_Z>* mdecomp;														//INDEX_TYPE internal_id; // this is the position in the yz plane
			const RegularGrid3D* mgrid;
			const RegularGrid3D* mgridofblocks;

			in_block_iterator(INDEX_TYPE bn, const RegularGrid3DDecompositionFixed<BSP2_X, BSP2_Y, BSP2_Z>* decomposition) :
				block_num(bn),
				b_pos(decomposition->get_grid_of_blocks()->Coords(bn)),
				b_extents(decomposition->get_block_size(b_pos)),
				b_jumps(decomposition->get_block_grid()->XYZ() - b_extents),
				last_id(decomposition->MakeBNPair(bn, decomposition->get_block_grid()->Index3d(b_extents - 1))),
				mgridofblocks(decomposition->get_grid_of_blocks()),
				mgrid(decomposition->get_block_grid()),
				mdecomp(decomposition) {}
			//in_block_iterator(INDEX_TYPE bn, const RegularGrid3DDecompositionFixed* decomp) :
			void begin() {
				in_block_id = mdecomp->MakeBNPair(block_num, 0);
				pos = { 0,0,0 };
				is_valid = true;
			}
			void begin(BN_ID_PAIR start_bnidp) {
				in_block_id = start_bnidp;
				pos = mdecomp->get_block_grid()->Coords(mdecomp->GetIDFromPair(start_bnidp));
				is_valid = true;
			}
			void advance() {
				pos[0]++;
				in_block_id++;
				if (pos[0] == b_extents[0]) {
					pos[0] = 0;
					in_block_id += b_jumps[0]; // the gap
					pos[1]++;
					if (pos[1] == b_extents[1]) {
						pos[1] = 0;
						pos[2]++;
						in_block_id += b_extents[0] * b_jumps[1];
						is_valid = in_block_id < last_id;
					}
				}
			}
			BN_ID_PAIR value() const {
				return in_block_id;
			}
			bool valid() const { return is_valid; }

		};

		struct neighborhood6_iterator {
			int num_neg;
			int m_id;
			BN_ID_PAIR results[6];
			const RegularGrid3DDecompositionFixed<BSP2_X, BSP2_Y, BSP2_Z>* decomposition;
			neighborhood6_iterator(const RegularGrid3DDecompositionFixed<BSP2_X, BSP2_Y, BSP2_Z>* decomp) : decomposition(decomp) {}

			void begin(BN_ID_PAIR start_bnidp) {
				INDEX_TYPE bn, id;
				decomposition->GetBNAndIDFromPair(start_bnidp, bn, id);
				in_block_iterator it(bn, decomposition);
				it.begin(start_bnidp);
				this->begin(it);
			}

			void begin(const in_block_iterator& base) {
				num_neg = 0;

				// check local block boundaries
				//// ------- X ---------
				auto gob = decomposition->get_grid_of_blocks();
				// upper boundaries X
				if (base.pos[0] < base.b_extents[0] - 1) {
					// we are "internal" with respect to the upper boundary
					results[num_neg++] = base.in_block_id + 1; // base.mgrid->Offset_x();
				}
				else if (base.b_pos[0] < base.mgridofblocks->XYZ()[0] - 1) {
					Vec3l neg_block_pos(base.b_pos[0] + 1, base.b_pos[1], base.b_pos[2]);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(0, base.pos[1], base.pos[2]);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}
				// now minus direction
				if (base.pos[0] > 0) {
					// we are "internal" with respect to the lower boundary
					results[num_neg++] = base.in_block_id - 1;// base.mgrid->Offset_x();
				}
				else if (base.b_pos[0] > 0) {
					Vec3l neg_block_pos(base.b_pos[0] - 1, base.b_pos[1], base.b_pos[2]);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(bg->XYZ()[0] - 1, base.pos[1], base.pos[2]);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}


				//// ------- Y ---------
				// upper boundaries Y
				if (base.pos[1] < base.b_extents[1] - 1) {
					// we are "internal" with respect to the upper boundary
					results[num_neg++] = base.in_block_id + base.mgrid->Offset_y();
				}
				else if (base.b_pos[1] < base.mgridofblocks->XYZ()[1] - 1) {
					Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1] + 1, base.b_pos[2]);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(base.pos[0], 0, base.pos[2]);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}
				// now minus direction
				if (base.pos[1] > 0) {
					// we are "internal" with respect to the lower boundary
					results[num_neg++] = base.in_block_id - base.mgrid->Offset_y();
				}
				else if (base.b_pos[1] > 0) {
					Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1] - 1, base.b_pos[2]);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(base.pos[0], bg->XYZ()[1] - 1, base.pos[2]);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}


				//// ------- Z ---------
				// upper boundaries Z
				if (base.pos[2] < base.b_extents[2] - 1) {
					// we are "internal" with respect to the upper boundary
					results[num_neg++] = base.in_block_id + base.mgrid->Offset_z();
				}
				else if (base.b_pos[2] < base.mgridofblocks->XYZ()[2] - 1) {
					Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1], base.b_pos[2] + 1);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(base.pos[0], base.pos[1], 0);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}
				// now minus direction
				if (base.pos[2] > 0) {
					// we are "internal" with respect to the lower boundary
					results[num_neg++] = base.in_block_id - base.mgrid->Offset_z();
				}
				else if (base.b_pos[2] > 0) {
					Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1], base.b_pos[2] - 1);
					auto bg = decomposition->get_block_grid(neg_block_pos);
					Vec3l neg_id_pos(base.pos[0], base.pos[1], bg->XYZ()[2] - 1);
					results[num_neg++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos));
				}
				m_id = 0;
			}

			bool valid() {
				return m_id < num_neg;
			}

			void advance() {
				m_id++;
			}
			BN_ID_PAIR value() {
				return results[m_id];
			}

		};
		struct neighborhood27_iterator {
			int num_neg;
			int m_id;
			const RegularGrid3DDecompositionFixed<BSP2_X, BSP2_Y, BSP2_Z>* decomposition;
			BN_ID_PAIR results[27];
			neighborhood27_iterator(const RegularGrid3DDecompositionFixed<BSP2_X, BSP2_Y, BSP2_Z>* decomp) : decomposition(decomp) {}

			void begin(BN_ID_PAIR start_bnidp) {
				INDEX_TYPE bn, id;
				decomposition->GetBNAndIDFromPair(start_bnidp, bn, id);
				in_block_iterator it(bn, decomposition);
				it.begin(start_bnidp);
				this->begin(it);
			}

			void begin(const in_block_iterator& base) {

				//if (base.pos[0] > 0 && base.pos[0] < base.b_extents[0] - 1 &&
				//	base.pos[1] > 0 && base.pos[1] < base.b_extents[1] - 1 &&
				//	base.pos[2] > 0 && base.pos[2] < base.b_extents[2] - 1) {
				//	// internal
				//	base.mgrid->kNeighborOffsets26[0];
				//}




				num_neg = 0;

				int offsets_count_x = 0;
				BN_ID_PAIR offsets_x[3];
				int offsets_count_y = 0;
				BN_ID_PAIR offsets_y[3];
				int offsets_count_z = 0;
				BN_ID_PAIR offsets_z[3];


				// check local block boundaries
				//// ------- X ---------
				auto gob = decomposition->get_grid_of_blocks();

				/// X OFFSETS
				{
					// upper boundaries X
					if (base.pos[0] < base.b_extents[0] - 1) {
						// we are "internal" with respect to the upper boundary
						offsets_x[offsets_count_x++] = 1; // base.mgrid->Offset_x();
					}
					else if (base.b_pos[0] < base.mgridofblocks->XYZ()[0] - 1) {
						Vec3l neg_block_pos(base.b_pos[0] + 1, base.b_pos[1], base.b_pos[2]);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(0, base.pos[1], base.pos[2]);
						offsets_x[offsets_count_x++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
					// now zero direction
					offsets_x[offsets_count_x++] = 0;
					// now minus direction
					if (base.pos[0] > 0) {
						// we are "internal" with respect to the lower boundary
						offsets_x[offsets_count_x++] = -1;// base.mgrid->Offset_x();
					}
					else if (base.b_pos[0] > 0) {
						Vec3l neg_block_pos(base.b_pos[0] - 1, base.b_pos[1], base.b_pos[2]);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(bg->XYZ()[0] - 1, base.pos[1], base.pos[2]);
						offsets_x[offsets_count_x++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
				}


				//// ------- Y ---------
				{
					// upper boundaries Y
					if (base.pos[1] < base.b_extents[1] - 1) {
						// we are "internal" with respect to the upper boundary
						offsets_y[offsets_count_y++] = base.mgrid->Offset_y();
					}
					else if (base.b_pos[1] < base.mgridofblocks->XYZ()[1] - 1) {
						Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1] + 1, base.b_pos[2]);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(base.pos[0], 0, base.pos[2]);
						offsets_y[offsets_count_y++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
					// no offset direction
					offsets_y[offsets_count_y++] = 0;

					// now minus direction
					if (base.pos[1] > 0) {
						// we are "internal" with respect to the lower boundary
						offsets_y[offsets_count_y++] = -base.mgrid->Offset_y();
					}
					else if (base.b_pos[1] > 0) {
						Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1] - 1, base.b_pos[2]);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(base.pos[0], bg->XYZ()[1] - 1, base.pos[2]);
						offsets_y[offsets_count_y++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
				}

				//// ------- Z ---------
				{
					// upper boundaries Z
					if (base.pos[2] < base.b_extents[2] - 1) {
						// we are "internal" with respect to the upper boundary
						offsets_z[offsets_count_z++] = base.mgrid->Offset_z();
					}
					else if (base.b_pos[2] < base.mgridofblocks->XYZ()[2] - 1) {
						Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1], base.b_pos[2] + 1);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(base.pos[0], base.pos[1], 0);
						offsets_z[offsets_count_z++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
					// no z offset
					offsets_z[offsets_count_z++] = 0;
					// now minus direction
					if (base.pos[2] > 0) {
						// we are "internal" with respect to the lower boundary
						offsets_z[offsets_count_z++] = -base.mgrid->Offset_z();
					}
					else if (base.b_pos[2] > 0) {
						Vec3l neg_block_pos(base.b_pos[0], base.b_pos[1], base.b_pos[2] - 1);
						auto bg = decomposition->get_block_grid(neg_block_pos);
						Vec3l neg_id_pos(base.pos[0], base.pos[1], bg->XYZ()[2] - 1);
						offsets_z[offsets_count_z++] = decomposition->MakeBNPair(gob->Index3d(neg_block_pos), bg->Index3d(neg_id_pos)) - base.in_block_id;
					}
				}

				num_neg = 0;
				for (int z = 0; z < offsets_count_z; z++)
					for (int y = 0; y < offsets_count_y; y++)
						for (int x = 0; x < offsets_count_x; x++) {
							results[num_neg++] = base.in_block_id + offsets_x[x] + offsets_y[y] + offsets_z[z];
						}

				m_id = 0;
			}

			bool valid() {
				return m_id < num_neg;
			}

			void advance() {
				m_id++;
			}
			BN_ID_PAIR value() {
				return results[m_id];
			}

		};
		// this is a comment on the code stuff
		// this is another line of code
		//localExtents m_global_grid;

		RegularGrid3D* m_block_grid;
		RegularGrid3D* m_grid_of_blocks;
		INDEX_TYPE m_num_blocks;
		Vec3l m_blocks_per_axis;
		INDEX_TYPE NumBlocks() const { return m_num_blocks; }

		//bool isGlobalBoundary(int nx, int ny, int nz, localExtents& l) {
		//	int x = nx + l.starts[0] * 2;
		//	int y = ny + l.starts[1] * 2;
		//	int z = nz + l.starts[2] * 2;
		//	return x == 0 || x == ((XYZ()[0] * 2 - 1) - 1) ||
		//		y == 0 || y == ((XYZ()[1] * 2 - 1) - 1) ||
		//		z == 0 || z == ((XYZ()[2] * 2 - 1) - 1);
		//}

	public:

		RegularGrid3DDecompositionFixed(
			INDEX_TYPE global_X,
			INDEX_TYPE global_Y,
			INDEX_TYPE global_Z) :
			RegularGrid3D(Vec3l(global_X, global_Y, global_Z), Vec3b(0, 0, 0))
		{}

		virtual ~RegularGrid3DDecompositionFixed() {
			delete m_block_grid;
			delete m_grid_of_blocks;
		}

		//void printBlockInfo(INT_TYPE block_id) {
		//}

		// for block number at i,j,k (in terms of block indices) 
		inline INDEX_TYPE block_id_from_block_pos(Vec3l ijk) const {
			return ijk[0] + ijk[1] * m_blocks_per_axis[0] + ijk[2] * m_blocks_per_axis[0] * m_blocks_per_axis[1];
		}
		// block ijk from global coords b_extents
		inline Vec3l block_ijk_from_xyz(Vec3l b_extents) const {
			return{ b_extents[0] >> BSP2_X ,b_extents[1] >> BSP2_Y, b_extents[2] >> BSP2_Z };
		}
		// get the i of the ijk block position from x coordinate
		inline INDEX_TYPE block_i_from_x(INDEX_TYPE x) const {
			return x >> BSP2_X;
		}
		// get the j of the ijk block position from y coordinate
		inline INDEX_TYPE block_j_from_y(INDEX_TYPE y) const {
			return y >> BSP2_Y;
		}
		// get the k of the ijk block position from z coordinate
		inline INDEX_TYPE block_k_from_z(INDEX_TYPE z) const {
			return z >> BSP2_Z;
		}
		// the number of blocks per axis
		inline Vec3l blocks_per_axis() const {
			return m_blocks_per_axis;
		}

	protected:

		inline unsigned char get_ijk_block_type(const Vec3l& ijk) const {
			unsigned char val_i = ijk[0] == (m_blocks_per_axis[0] - 1);
			unsigned char val_j = (ijk[1] == (m_blocks_per_axis[1] - 1)) << 1;
			unsigned char val_k = (ijk[2] == (m_blocks_per_axis[2] - 1)) << 2;
			return val_i | val_j | val_k;
		}
	public:
		// get the localextents for a block number 
		inline  void get_block(const Vec3l& ijk, localExtents& res) const {
			block_start(ijk, res.starts);
			block_size(ijk, res.starts, res.sizes);
			res.grid = m_block_grid;
		}
		// the grid associated with block at block pos ijk
		inline  const RegularGrid3D* get_block_grid(const Vec3l& ijk) const {
			return m_block_grid;
		}

		inline const RegularGrid3D* get_grid_of_blocks() const {
			return m_grid_of_blocks;
		}
		// the grid associated with block number
		inline  const RegularGrid3D* get_block_grid(INDEX_TYPE block_num) const {
			return m_block_grid;
		}
		inline  const RegularGrid3D* get_block_grid() const {
			return m_block_grid;
		}
		// block ijk from BLOCK id - not global id
		inline Vec3l block_ijk_from_id(INDEX_TYPE block_num) const {
			return m_grid_of_blocks->Coords(block_num);
			//INDEX_TYPE slab = (m_blocks_per_axis[0] * m_blocks_per_axis[1]);
			//return Vec3l(block_num % m_blocks_per_axis[0],
			//	(block_num % slab) / m_blocks_per_axis[0],
			//	block_num / slab);
		}

		// global id to bnid: this is slow
		inline BN_ID_PAIR get_bnid_from_id(INDEX_TYPE id) {
			auto b_extents = this->Coords(id);
			auto ijk = block_ijk_from_xyz(b_extents);
			Vec3l bs;
			block_start(ijk, bs);
			auto in_block_xyz = b_extents - bs;
			INDEX_TYPE inblockid = get_block_grid(ijk)->Index3d(in_block_xyz);
			return MakeBNPair(block_id_from_block_pos(ijk), inblockid);
		}

		inline void block_start(const Vec3l& ijk, Vec3l& starts) const {
			starts[0] = ijk[0] * (1 << BSP2_X);
			starts[1] = ijk[1] * (1 << BSP2_Y);
			starts[2] = ijk[2] * (1 << BSP2_Z);
		}
		inline void block_size(const Vec3l& ijk, const Vec3l& starts, Vec3l& sizes) const {
			sizes = { (1 << BSP2_X), (1 << BSP2_Y), (1 << BSP2_Z) };
			if (ijk[0] == m_blocks_per_axis[0] - 1)
				sizes[0] = XYZ()[0] - starts[0];
			if (ijk[1] == m_blocks_per_axis[1] - 1)
				sizes[1] = XYZ()[1] - starts[1];
			if (ijk[2] == m_blocks_per_axis[2] - 1)
				sizes[2] = XYZ()[2] - starts[2];
		}
		inline void block_size(const Vec3l& ijk, Vec3l& sizes) const {
			Vec3l starts;
			block_start(ijk, starts);
			block_size(ijk, starts, sizes);
		}
		inline Vec3l get_block_start_xyz(const Vec3l& ijk) const {
			return{ ijk[0] * (1 << BSP2_X),
				ijk[1] * (1 << BSP2_Y),
				ijk[2] * (1 << BSP2_Z) };
		}
		inline Vec3l get_block_size(const Vec3l& ijk) const {
			Vec3l sizes;
			block_size(ijk, sizes);
			return sizes;
		}
		//// get block number conatining global position x,y,z
		//inline INDEX_TYPE block_num(Vec3l b_extents) const {
		//	return
		//		(b_extents[0] >> m_blocksize_pow2[0]) +
		//		(b_extents[1] >> m_blocksize_pow2[1]) * m_blocks_per_axis[0] +
		//		(b_extents[2] >> m_blocksize_pow2[2]) * m_blocks_per_axis[0] * m_blocks_per_axis[1];
		//}



		virtual void decompose() {
			// figure out number in each dimension
			m_blocks_per_axis[0] = XYZ()[0] / (1 << BSP2_X) + (XYZ()[0] % (1 << BSP2_X) > 0 ? 1 : 0);
			m_blocks_per_axis[1] = XYZ()[1] / (1 << BSP2_Y) + (XYZ()[1] % (1 << BSP2_Y) > 0 ? 1 : 0);
			m_blocks_per_axis[2] = XYZ()[2] / (1 << BSP2_Z) + (XYZ()[2] % (1 << BSP2_Z) > 0 ? 1 : 0);
			m_grid_of_blocks = new RegularGrid3D(m_blocks_per_axis, { 0,0,0 });
			m_num_blocks = m_grid_of_blocks->NumElements();
			m_block_grid = new RegularGrid3D({ (1 << BSP2_X) , (1 << BSP2_Y) , (1 << BSP2_Z) }, { 0,0,0 });

		}


	};


} // end GInt namespace

#endif
