/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/


#include <vector>
#include <set>
#include <queue>
#include <time.h>
#include "gi_timing.h"
#include "gi_topological_explicit_mesh_function.h"

#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"

#include "gi_topological_utility_functions.h"
#include "gi_numeric_integrator_expanding_region_stop.h" // not where comparer should be
#include "gi_topological_regular_masked_grid.h"
#include "gi_topological_regular_masked_restricted_grid.h"

#include "gi_max_vertex_labeling.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_bifiltration_pairing.h"
#include "gi_fast_robins_noalloc.h"



#include <Visus/IdxDataset.h>
#include <Visus/ApplicationInfo.h>

#include <fstream>

#ifdef WIN32
#pragma warning(disable:4996)
#endif

using namespace Visus;



#define USEMAXV
//using namespace GInt;
//typedef IndexCompareLessThan Comparer;
//typedef RegularGrid3D GridType;
//typedef TopologicalRegularGrid3D MeshType;
//typedef RegularGridTrilinearFunction GridFuncType;
//#ifndef USEMAXV
//typedef TopologicalExplicitDenseMeshFunction<MeshType, float> TopoFuncType;
//#else
//typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
//typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType , GridFuncType, float> TopoFuncType;
//#endif
//typedef DiscreteGradientLabeling<MeshType> GradType;
//typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType> RobinsType;
//typedef MyRobins<MeshType, MaxVLType, GradType> RobinsTypeO;
using namespace GInt;
////typedef int my_item;
////int my_shift = 4;
////class SparseBlockedGridLabeling {
////protected:
////	RegularGrid3D* mGrid;
////	my_item** mElements;
////	INDEX_TYPE mNumBlocks;
////
////	INDEX_TYPE BlockNum(Vec3l c) {
////		
////	}
////public:
////	SparseBlockedGridLabeling(RegularGrid3D* grid) : mGrid(grid) {
////		Vec3l extents = grid->XYZ() >> my_shift;
////		mNumBlocks = extents[0] * extents[1] * extents[2];
////		mElements = new my_item*[mNumBlocks];
////		memset(mElements, NULL, sizeof(my_item*) * mNumBlocks);
////	}
////
////	~SparseBlockedGridLabeling() {
////		for (INDEX_TYPE id = 0; id < mNumBlocks; id++)
////			if (mElements[id] != NULL) delete[] mElements[id];
////		delete[] mElements;
////	}
////
////	void SetLabel(INDEX_TYPE id, my_item label) {
////		m_labels[id] = label;
////	}
////
////	LABEL_TYPE GetLabel(INDEX_TYPE id) const {
////		return m_labels[id];
////	}
////
////	INDEX_TYPE GetNumLabels() const {
////		return m_num_labels;
////	}
////
////	LABEL_TYPE& operator[](const INDEX_TYPE id) { return m_labels[id]; }
////	const LABEL_TYPE& operator[](const INDEX_TYPE id) const { return m_labels[id]; }
////
////	void SetAll(LABEL_TYPE label){
////#pragma omp parallel for schedule(static)
////		for (int i = 0; i < m_num_labels; i++) {
////			m_labels[i] = label;
////		}
////	}
////
////	void CopyValues(const DenseLabeling<LABEL_TYPE>* other) {
////#pragma omp parallel for schedule(static)
////		for (int i = 0; i < m_num_labels; i++) {
////			m_labels[i] = other->m_labels[i];
////		}
////	}
////	//template<typename T>
////	//void ReMapIds(T* output){
////	//  std::unordered_map<INDEX_TYPE, T> unique_ids;
////	//  
////	//  T new_id=0;
////	//  
////	//  for(INDEX_TYPE i=0; i < m_num_labels; i++){
////	//    if(unique_ids.find(m_labels[i]) == unique_ids.end()){
////	//      T set_id = new_id;
////	//      
////	//      if(m_labels[i] < 0)
////	//        set_id = -1;
////	//      else
////	//        set_id = new_id++;
////	//      
////	//      unique_ids[m_labels[i]] = set_id;
////	//      printf("mapping id %lld to %d\n", m_labels[i], set_id);
////	//    }
////	//    
////	//    output[i] = unique_ids.at(m_labels[i]);
////	//  }
////	//  
////	//  printf("remapped %d ids\n", new_id);
////	//}
////
////
////	void ReadFromFile(const char* filename) {
////		FILE* fout = fopen(filename, "rb");
////		fread(m_labels, sizeof(LABEL_TYPE), m_num_labels, fout);
////		fclose(fout);
////	}
////
////	void OutputToFile(const char* filename) const {
////		printf("writing file %s \n", filename);
////		FILE* fout = fopen(filename, "wb");
////		//printf("Sizeof label type: %d\n", sizeof(LABEL_TYPE));
////		fwrite(m_labels, sizeof(LABEL_TYPE), m_num_labels, fout);
////		fclose(fout);
////	}
////	void OutputToIntFile(const char* filename) const {
////		printf("writing file %s \n", filename);
////		FILE* fout = fopen(filename, "wb");
////		//printf("Sizeof int type: %d\n", sizeof(int));
////		for (INDEX_TYPE i = 0; i < m_num_labels; i++) {
////			int tval = (int)m_labels[i];
////			fwrite(&tval, sizeof(int), 1, fout);
////		}
////		fclose(fout);
////	}
////	void OutputToFloatFile(const char* filename) const {
////		printf("writing file %s \n", filename);
////		FILE* fout = fopen(filename, "wb");
////		//printf("Sizeof float type: %d\n", sizeof(float));
////		for (INDEX_TYPE i = 0; i < m_num_labels; i++) {
////			float tval = (float)m_labels[i];
////			fwrite(&tval, sizeof(float), 1, fout);
////		}
////		fclose(fout);
////	}
////};




typedef RegularGrid3D GridType;
typedef RegularGridTrilinearFunction GridFuncType;
//typedef UncachedRegularGridTrilinearFunction GridFuncType;
typedef TopologicalRegularGrid3D MeshType;
typedef DiscreteGradientLabeling<MeshType> GradType;
//typedef UncachedMaximumVertexLabeling<GridType, GridFuncType> MaxVLType;
//typedef MaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
//typedef RegularGridMaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
typedef RegularGridMaxMinVertexLabeling3D<MeshType, GridFuncType> MaxVLType;
typedef MyRobinsNoalloc<MeshType, MaxVLType, GradType, 5, 4> RobinsType;
typedef TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> TopoFuncType;
typedef SlidingWindowRobinsNoalloc < GridType, GridFuncType, MeshType, MaxVLType, GradType> NewRobinsType;

float* GetSubBlockRaw(Vec3l global_size, Vec3l start, Vec3l size, const char* fname) {



	// load values from file
	INDEX_TYPE num_data = size[0] * size[1] * size[2];
	float* local_data = new float[num_data];

	FILE* fin = fopen(fname, "rb");
	// seek to start
	INDEX_TYPE t_file_start = start[2] * global_size[0] *
		global_size[1];

	fseek(fin, sizeof(float) * t_file_start, SEEK_SET);
	float* curr_loc = local_data;
	printf("get here 1 : %d\n", t_file_start);

	int readcount = 0;
	for (int z = 0; z < size[2]; z++) {
		fseek(fin, sizeof(float) * start[1] * global_size[0], SEEK_CUR);
		for (int y = 0; y < size[1]; y++) {
			//skip empty y space

			//for(int x = 0; x < l.sizes[0]; x++) {
			// skip empty x space
			fseek(fin, sizeof(float) * start[0], SEEK_CUR);
			fread(curr_loc, sizeof(float), size[0], fin);
			readcount += size[0];

			curr_loc += /*sizeof(dtype) **/ size[0];
			fseek(fin, sizeof(float) *
				(global_size[0] - (start[0] + size[0])), SEEK_CUR);
			//}

		}
		fseek(fin, sizeof(float) * global_size[0] *
			(global_size[1] - (start[1] + size[1])), SEEK_CUR);
	}
	fclose(fin);
	printf("load returned: read %d items\n", readcount);
	return local_data;

}

float* GetSubBlockIdx(Visus::SharedPtr<Visus::Dataset> dataset, Vec3l start, Vec3l size) {

	//any time you need to read/write data from/to a Dataset create an Access
	auto access = dataset->createAccess();

	BoxNi logic_box(Point3i(start[0], start[1], start[2]), Point3i(start[0] + size[0], start[1] + size[1], start[2] + size[2]));

	PrintInfo("Box query", logic_box.p1, "p2", logic_box.p2);

	auto field = dataset->getDefaultField();
	double timestate = 0;
	auto query = std::make_shared<BoxQuery>(dataset.get(), field, timestate, 'r');
	query->logic_box = logic_box;

	// Set resolution levels
	// You can add multiple resolutions values to end_resolutions
	query->setResolutionRange(dataset->getMaxResolution(), dataset->getMaxResolution());
	dataset->beginQuery(query);
	dataset->executeQuery(access, query);
	Array data = query->buffer;
	return (float*)data.c_ptr();
}

//template<typename dtype>
#define dtype float
//class SparseBlockedGridFunction : public RegularGrid3DDecomposition {
//protected:
//	dtype** m_data_blocks;
//public:
//
//	SparseBlockedGridFunction(INDEX_TYPE global_X,
//		INDEX_TYPE global_Y,
//		INDEX_TYPE global_Z,
//		INDEX_TYPE pow_2_size_X,
//		INDEX_TYPE pow_2_size_Y,
//		INDEX_TYPE pow_2_size_Z,
//		RegularGrid3D* global_grid)  : 
//		RegularGrid3DDecomposition(global_X, global_Y, global_Z, 
//			pow_2_size_X, pow_2_size_Y, pow_2_size_Z, global_grid) {
//	}
//	virtual ~SparseBlockedGridFunction() {
//		for (INDEX_TYPE i = 0; i < m_num_blocks; i++) {
//			if (this->has_block(i)) delete[] m_data_blocks[i];
//		}
//		delete[] m_data_blocks;
//	}
//
//	virtual void decompose() {
//		RegularGrid3DDecomposition::decompose();
//		m_data_blocks = new dtype*[this->NumBlocks()](); // should value-initialize to zero
//		
//	}
//
//	inline dtype SampleData(Vec3l point) const {
//		auto bn = this->block_num(point);
//		const auto& block =  m_blocks[bn];
//		Vec3l ind = point - block.starts;
//		return m_data_blocks[bn][block.grid->Index3d(ind)];
//	}
//
//	bool ReadFromRAWDense(const char* fname) {
//
//		// load values from file
//		INDEX_TYPE slab_size = this->XYZ()[0] * this->XYZ()[1];
//		INDEX_TYPE slab_width = 5;
//#pragma omp parallel
//		{
//			slab_width = omp_get_num_threads() * 1;
//		}
//		printf("slab width = %d\n", slab_width);
//		printf("slab size = (%llu x %llu)\n", this->XYZ()[0], this->XYZ()[1]);
//		dtype* slab_data = new dtype[slab_size * slab_width];
//
//
//		FILE* fin = fopen(fname, "rb");
//		if (fin == NULL) {
//			printf("Error reading %s\n", fname);
//			return false;
//		}
//		printf("allocating data blocks\n");
//		// allocate data blocks
//		for (INDEX_TYPE i = 0; i < NumBlocks(); i++) {
//			//printf("  -- %llu elements\n", m_blocks[i].grid->NumElements());
//			m_data_blocks[i] = new dtype[m_blocks[i].grid->NumElements()];
//		}
//		printf("done allocating\n");
//
//		// outer file reading loop
//		INDEX_TYPE slab_current_z = 0;
//		INDEX_TYPE block_z_size = this->m_desired_block_sizes[2];
//
//		while (slab_current_z < this->XYZ()[2]) {
////#define DEBUG_FIRST_SLICE
//#ifdef DEBUG_FIRST_SLICE
//			if (slab_current_z > 0) break;
//#endif
//			// read next slab
//			if (slab_current_z + slab_width > this->XYZ()[2]) {
//				slab_width = this->XYZ()[2] - slab_current_z;
//			}
//			INDEX_TYPE num_to_read = slab_size * slab_width;
//			fread(slab_data, sizeof(dtype), num_to_read, fin);
//			//printf(" -- read %llu elements\n", num_to_read);
//			// now have read slab_width *X*Y values
//
//			// copy into blocks
//#pragma omp parallel for schedule(dynamic)
//			for (int z = 0; z < slab_width; z++) {
//				INDEX_TYPE t_current_z = slab_current_z + z;
//#ifdef DEBUG_FIRST_SLICE
//				if (z > 0) {
//					slab_current_z++;
//					break;
//				}
//#endif
//				//printf("   -- doing z= %d, slab_current_z = %d\n", z, slab_current_z);
//				// foreach [x, x + blocksize], y stick, copy into block
//				for (int y = 0; y < XYZ()[1]; y++) {
//					//printf("   -- doing y= %d\n", y);
//					// get block number and just add 1 after inner loop
//					INDEX_TYPE bn = this->block_num({ 0, y, t_current_z });
//					//printf("   -- bn = %d\n", bn);
//					// inner loop with memory copy
//					for (int x = 0; x < XYZ()[0]; x += m_desired_block_sizes[0]) {
//						INDEX_TYPE in_block_y = y - m_blocks[bn].starts[1];
//						INDEX_TYPE in_block_z = t_current_z -  m_blocks[bn].starts[2] ;
//						INDEX_TYPE copy_size = m_desired_block_sizes[0];
//						if (x + m_desired_block_sizes[0] > XYZ()[0]) copy_size = XYZ()[0] - x;
//						//printf("   -- x = [%d, %d], in_block=(0, %d, %d)\n", x, x + copy_size, in_block_y, in_block_z);
//						//m_blocks[bn].starts.PrintInt();
//						// address in slab is z sliceinslab + y* slab_x + x
//						dtype* source = &(slab_data[x + this->XYZ()[0] * y + slab_size*z]);
//						
//						// address to copy to gets block at 
//						dtype* dest = &(m_data_blocks[bn][m_blocks[bn].grid->Index3d({ 0, in_block_y, in_block_z })]);
//						
//						//COPY
//						memcpy(dest, source, sizeof(dtype) * copy_size);
//						
//						//for (int t = 0; t < copy_size; t++) dest[t] = z;
//						bn++; // increment block number
//					} // end x loop
//				} // end y loop
//			} // end z loop
//			slab_current_z += slab_width;
//		} //end loop over slabs
//		fclose(fin);
//		return true;
//	}
//
//	unsigned char* blocks_to_include;
//	dtype m_background_value;
//
//	void SetBackgroundValue(dtype val) { m_background_value = val; }
//	const RegularGrid3D* get_global_grid() {
//		return this->m_global_grid.grid;
//	}
//
//	bool has_block(INDEX_TYPE block_num) const {
//		if (blocks_to_include == NULL) return true; // becasue we read dense
//		return blocks_to_include[block_num];
//	}
//	bool has_xyz(Vec3l xyz) const {
//		INDEX_TYPE block_num = this->block_num(xyz);
//		if (blocks_to_include == NULL) return true; // becasue we read dense
//		return blocks_to_include[block_num];
//	}
//
//	bool ReadFromRAWSparse(const char* fname, const char* maskname) {
//
//		// load values from file
//		INDEX_TYPE slab_size = this->XYZ()[0] * this->XYZ()[1];
//		INDEX_TYPE slab_width = 5;
//#pragma omp parallel
//		{
//			slab_width = omp_get_num_threads() * 1;
//		}
//		printf("slab width = %d\n", slab_width);
//		printf("slab size = (%llu x %llu)\n", this->XYZ()[0], this->XYZ()[1]);
//		dtype* slab_data = new dtype[slab_size * slab_width];
//		unsigned char* slab_mask = new unsigned char[slab_size * slab_width];
//
//
//		FILE* fin = fopen(fname, "rb");
//		if (fin == NULL) {
//			printf("Error reading %s\n", fname);
//			return false;
//		}
//
//		FILE* fin_mask = fopen(maskname, "rb");
//		if (fin == NULL) {
//			printf("Error reading %s\n", maskname);
//			fclose(fin);
//			return false;
//		}
//
//		//printf("allocating data blocks\n");
//		//// allocate data blocks
//		//for (INDEX_TYPE i = 0; i < NumBlocks(); i++) {
//		//	//printf("  -- %llu elements\n", m_blocks[i].grid->NumElements());
//		//	m_data_blocks[i] = new dtype[m_blocks[i].grid->NumElements()];
//		//}
//		//printf("done allocating\n");
//
//		// outer file reading loop
//		INDEX_TYPE slab_current_z = 0;
//		INDEX_TYPE block_z_size = this->m_desired_block_sizes[2];
//
//		blocks_to_include = new unsigned char[this->NumBlocks()](); // initialized to 0
//
//		while (slab_current_z < this->XYZ()[2]) {
//			//#define DEBUG_FIRST_SLICE
//#ifdef DEBUG_FIRST_SLICE
//			if (slab_current_z > 0) break;
//#endif
//			// read next slab
//			if (slab_current_z + slab_width > this->XYZ()[2]) {
//				slab_width = this->XYZ()[2] - slab_current_z;
//			}
//			INDEX_TYPE num_to_read = slab_size * slab_width;
//			fread(slab_data, sizeof(dtype), num_to_read, fin);
//			fread(slab_mask, sizeof(unsigned char), num_to_read, fin_mask);
//			// now have read slab_width *X*Y values
//			//printf("read slabs %d\n", slab_current_z);
//			// now make sure all NEEDED blocks are allocated
//			// - for each z in slab, check if data blocks intersecting z have non-zero in mask slab
//			// - will allocate needed data blocks
//			// - can check either non-zero data block pointer or 1 on blocks_to_include -- redundant? 
//			for (int z = 0; z < slab_width; z++) {
//				INDEX_TYPE t_current_z = slab_current_z + z;
//
//				INDEX_TYPE k = this->block_k_from_z(t_current_z);
//#pragma omp for schedule(dynamic)
//				for (auto j = 0; j < this->blocks_per_axis()[1]; j++) {
//					for (auto i = 0; i < this->blocks_per_axis()[0]; i++) {
//						INDEX_TYPE t_block_num = this->block_id_from_block_pos({ i,j,k });
//						// check if ijk block exists - if so, skip
//						if (blocks_to_include[t_block_num] != 0) continue;
//
//						//check z plane in block for non-zero
//						const Vec3l& starts = get_block(t_block_num).starts;
//						const Vec3l& sizes = get_block(t_block_num).sizes;
//
//						for (INDEX_TYPE y = 0; y < sizes[1]; y++) {
//							auto gy = y + starts[1];
//							for (INDEX_TYPE x = 0; x < sizes[0]; x++) {
//								auto gx = x + starts[0];
//								auto val = slab_mask[gx + this->XYZ()[0] * gy + slab_size*z];
//								if (val != 0) {
//									blocks_to_include[t_block_num] = 1; // mark the block
//									m_data_blocks[t_block_num] = new dtype[m_blocks[t_block_num].grid->NumElements()]; // allocate the space
//									y = sizes[1]; // skip out of y loop
//									break; // skip out of x loop
//								}
//							}
//						}
//
//
//					}
//				} // end parallel over j blocks
//			} // end iteration over z slices
//
//			//printf("read masks and allocatd\n");
//			// copy into blocks
//#pragma omp parallel for schedule(dynamic)
//			for (int z = 0; z < slab_width; z++) {
//				INDEX_TYPE t_current_z = slab_current_z + z;
//
//#ifdef DEBUG_FIRST_SLICE
//				if (z > 0) {
//					slab_current_z++;
//					break;
//				}
//#endif
//				//printf("   -- doing z= %d, slab_current_z = %d\n", z, slab_current_z);
//				// foreach [x, x + blocksize], y stick, copy into block
//				for (int y = 0; y < XYZ()[1]; y++) {
//					//printf("   -- doing y= %d\n", y);
//					// get block number and just add 1 after inner loop
//					INDEX_TYPE bn = this->block_num({ 0, y, t_current_z });
//					//printf("   -- bn = %d\n", bn);
//					// inner loop with memory copy
//					for (int x = 0; x < XYZ()[0]; x += m_desired_block_sizes[0]) {
//						if (blocks_to_include[bn] == 0) {
//							bn++;
//							continue;
//						}
//
//						INDEX_TYPE in_block_y = y - m_blocks[bn].starts[1];
//						INDEX_TYPE in_block_z = t_current_z - m_blocks[bn].starts[2];
//						INDEX_TYPE copy_size = m_desired_block_sizes[0];
//						if (x + m_desired_block_sizes[0] > XYZ()[0]) copy_size = XYZ()[0] - x;
//						//printf("   -- x = [%d, %d], in_block=(0, %d, %d)\n", x, x + copy_size, in_block_y, in_block_z);
//						//m_blocks[bn].starts.PrintInt();
//						// address in slab is z sliceinslab + y* slab_x + x
//						dtype* source = &(slab_data[x + this->XYZ()[0] * y + slab_size*z]);
//						unsigned char* mask = &(slab_mask[x + this->XYZ()[0] * y + slab_size*z]);
//
//						// address to copy to gets block at 
//						dtype* dest = &(m_data_blocks[bn][m_blocks[bn].grid->Index3d({ 0, in_block_y, in_block_z })]);
//						//COPY
//						//memcpy(dest, source, sizeof(dtype) * copy_size);
//
//						for (int t = 0; t < copy_size; t++) dest[t] = (mask[t] ? source[t] : m_background_value); // masked values
//						bn++; // increment block number
//					} // end x loop
//				} // end y loop
//			} // end z loop
//			slab_current_z += slab_width;
//		} //end loop over slabs
//		fclose(fin);
//		return true;
//	}
//
//
//
////	bool ReadFromRAWSparseImplicit(const char* fname) {
////
////		// load values from file
////		INDEX_TYPE slab_size = this->XYZ()[0] * this->XYZ()[1];
////		INDEX_TYPE slab_width = 5;
////#pragma omp parallel
////		{
////			slab_width = omp_get_num_threads() * 1;
////		}
////		printf("slab width = %d\n", slab_width);
////		printf("slab size = (%llu x %llu)\n", this->XYZ()[0], this->XYZ()[1]);
////		dtype* slab_data = new dtype[slab_size * slab_width];
////
////		FILE* fin = fopen(fname, "rb");
////		if (fin == NULL) {
////			printf("Error reading %s\n", fname);
////			return false;
////		}
////
////		// outer file reading loop
////		INDEX_TYPE slab_current_z = 0;
////		INDEX_TYPE block_z_size = this->m_desired_block_sizes[2];
////
////		blocks_to_include = new unsigned char[this->NumBlocks()](); // initialized to 0
////
////		while (slab_current_z < this->XYZ()[2]) {
////			//#define DEBUG_FIRST_SLICE
////#ifdef DEBUG_FIRST_SLICE
////			if (slab_current_z > 0) break;
////#endif
////			// read next slab
////			if (slab_current_z + slab_width > this->XYZ()[2]) {
////				slab_width = this->XYZ()[2] - slab_current_z;
////			}
////			INDEX_TYPE num_to_read = slab_size * slab_width;
////			fread(slab_data, sizeof(dtype), num_to_read, fin);
////			//fread(slab_mask, sizeof(unsigned char), num_to_read, fin_mask);
////			// now have read slab_width *X*Y values
////			//printf("read slabs %d\n", slab_current_z);
////			// now make sure all NEEDED blocks are allocated
////			// - for each z in slab, check if data blocks intersecting z have non-zero in mask slab
////			// - will allocate needed data blocks
////			// - can check either non-zero data block pointer or 1 on blocks_to_include -- redundant? 
////			for (int z = 0; z < slab_width; z++) {
////				INDEX_TYPE t_current_z = slab_current_z + z;
////
////				INDEX_TYPE k = this->block_k_from_z(t_current_z);
////#pragma omp for schedule(dynamic)
////				for (auto j = 0; j < this->blocks_per_axis()[2]; j++) {
////					for (auto i = 0; i < this->blocks_per_axis()[0]; i++) {
////						INDEX_TYPE t_block_num = this->block_id_from_block_pos({ i,j,k });
////						// check if ijk block exists - if so, skip
////						if (blocks_to_include[t_block_num] != 0) continue;
////
////						//check z plane in block for non-zero
////						const Vec3l& starts = get_block(t_block_num).starts;
////						const Vec3l& sizes = get_block(t_block_num).sizes;
////
////						for (INDEX_TYPE y = 0; y < sizes[1]; y++) {
////							auto gy = y + starts[1];
////							for (INDEX_TYPE x = 0; x < sizes[0]; x++) {
////								auto gx = x + starts[0];
////								auto val = slab_mask[gx + this->XYZ()[0] * gy + slab_size*z];
////								if (val != 0) {
////									blocks_to_include[t_block_num] = 1; // mark the block
////									m_data_blocks[t_block_num] = new dtype[m_blocks[t_block_num].grid->NumElements()]; // allocate the space
////									y = sizes[1]; // skip out of y loop
////									break; // skip out of x loop
////								}
////							}
////						}
////
////
////					}
////				} // end parallel over j blocks
////			} // end iteration over z slices
////
////			  //printf("read masks and allocatd\n");
////			  // copy into blocks
////#pragma omp parallel for schedule(dynamic)
////			for (int z = 0; z < slab_width; z++) {
////				INDEX_TYPE t_current_z = slab_current_z + z;
////
////#ifdef DEBUG_FIRST_SLICE
////				if (z > 0) {
////					slab_current_z++;
////					break;
////				}
////#endif
////				//printf("   -- doing z= %d, slab_current_z = %d\n", z, slab_current_z);
////				// foreach [x, x + blocksize], y stick, copy into block
////				for (int y = 0; y < XYZ()[1]; y++) {
////					//printf("   -- doing y= %d\n", y);
////					// get block number and just add 1 after inner loop
////					INDEX_TYPE bn = this->block_num({ 0, y, t_current_z });
////					//printf("   -- bn = %d\n", bn);
////					// inner loop with memory copy
////					for (int x = 0; x < XYZ()[0]; x += m_desired_block_sizes[0]) {
////						if (blocks_to_include[bn] == 0) {
////							bn++;
////							continue;
////						}
////
////						INDEX_TYPE in_block_y = y - m_blocks[bn].starts[1];
////						INDEX_TYPE in_block_z = t_current_z - m_blocks[bn].starts[2];
////						INDEX_TYPE copy_size = m_desired_block_sizes[0];
////						if (x + m_desired_block_sizes[0] > XYZ()[0]) copy_size = XYZ()[0] - x;
////						//printf("   -- x = [%d, %d], in_block=(0, %d, %d)\n", x, x + copy_size, in_block_y, in_block_z);
////						//m_blocks[bn].starts.PrintInt();
////						// address in slab is z sliceinslab + y* slab_x + x
////						dtype* source = &(slab_data[x + this->XYZ()[0] * y + slab_size*z]);
////
////						// address to copy to gets block at 
////						dtype* dest = &(m_data_blocks[bn][m_blocks[bn].grid->Index3d({ 0, in_block_y, in_block_z })]);
////
////						//COPY
////						memcpy(dest, source, sizeof(dtype) * copy_size);
////
////						//for (int t = 0; t < copy_size; t++) dest[t] = z;
////						bn++; // increment block number
////					} // end x loop
////				} // end y loop
////			} // end z loop
////			slab_current_z += slab_width;
////		} //end loop over slabs
////		fclose(fin);
////		return true;
////	}
//
//
//};

//class SparseBlockedGridFunctionImplicit : public GridDecompType {
//protected:
//	dtype** m_data_blocks;
//public:
//
//	SparseBlockedGridFunctionImplicit(INDEX_TYPE global_X,
//		INDEX_TYPE global_Y,
//		INDEX_TYPE global_Z,
//		INDEX_TYPE pow_2_size_X,
//		INDEX_TYPE pow_2_size_Y,
//		INDEX_TYPE pow_2_size_Z) :
//		GridDecompType(global_X, global_Y, global_Z,
//			pow_2_size_X, pow_2_size_Y, pow_2_size_Z) {
//	}
//	virtual ~SparseBlockedGridFunctionImplicit() {
//		for (INDEX_TYPE i = 0; i < m_num_blocks; i++) {
//			if (this->has_block(i)) delete[] m_data_blocks[i];
//		}
//		delete[] m_data_blocks;
//	}
//	virtual void decompose() {
//		GridDecompType::decompose();
//		m_data_blocks = new dtype*[this->NumBlocks()](); // should value-initialize to zero
//
//	}
//
//	inline dtype SampleData(Vec3l point) const {
//		auto ijk = this->block_ijk_from_xyz(point);
//		Vec3l starts;
//		block_start(ijk, starts);
//		Vec3l ind = point - starts;
//		return m_data_blocks[block_id_from_block_pos(ijk)][get_block_grid(ijk)->Index3d(ind)];
//	}
//
//	bool ReadFromRAWDense(const char* fname) {
//
//		// load values from file
//		INDEX_TYPE slab_size = this->XYZ()[0] * this->XYZ()[1];
//		INDEX_TYPE slab_width = 5;
//#pragma omp parallel
//		{
//			slab_width = omp_get_num_threads() * 1;
//		}
//		printf("slab width = %d\n", slab_width);
//		printf("slab size = (%llu x %llu)\n", this->XYZ()[0], this->XYZ()[1]);
//		dtype* slab_data = new dtype[slab_size * slab_width];
//
//
//		FILE* fin = fopen(fname, "rb");
//		if (fin == NULL) {
//			printf("Error reading %s\n", fname);
//			return false;
//		}
//		printf("allocating data blocks\n");
//		// allocate data blocks
//		for (INDEX_TYPE i = 0; i < NumBlocks(); i++) {
//			//printf("  -- %llu elements\n", m_blocks[i].grid->NumElements());
//			Vec3l ijk = block_ijk_from_id(i);
//			m_data_blocks[i] = new dtype[get_block_grid(ijk)->NumElements()];
//		}
//		printf("done allocating\n");
//
//		// outer file reading loop
//		INDEX_TYPE slab_current_z = 0;
//		INDEX_TYPE block_z_size = this->m_desired_block_sizes[2];
//
//		while (slab_current_z < this->XYZ()[2]) {
//			//#define DEBUG_FIRST_SLICE
//#ifdef DEBUG_FIRST_SLICE
//			if (slab_current_z > 0) break;
//#endif
//			// read next slab
//			if (slab_current_z + slab_width > this->XYZ()[2]) {
//				slab_width = this->XYZ()[2] - slab_current_z;
//			}
//			INDEX_TYPE num_to_read = slab_size * slab_width;
//			fread(slab_data, sizeof(dtype), num_to_read, fin);
//			//printf(" -- read %llu elements\n", num_to_read);
//			// now have read slab_width *X*Y values
//
//			// copy into blocks
//#pragma omp parallel for schedule(dynamic)
//			for (int z = 0; z < slab_width; z++) {
//				INDEX_TYPE t_current_z = slab_current_z + z;
//#ifdef DEBUG_FIRST_SLICE
//				if (z > 0) {
//					slab_current_z++;
//					break;
//				}
//#endif
//				//printf("   -- doing z= %d, slab_current_z = %d\n", z, slab_current_z);
//				// foreach [x, x + blocksize], y stick, copy into block
//				for (int y = 0; y < XYZ()[1]; y++) {
//					//printf("   -- doing y= %d\n", y);
//					// get block number and just add 1 after inner loop
//					Vec3l ijk= this->block_ijk_from_xyz({ 0, y, t_current_z });
//					//printf("   -- bn = %d\n", bn);
//					// inner loop with memory copy
//					for (int x = 0; x < XYZ()[0]; x += m_desired_block_sizes[0]) {
//						localExtents le;
//						this->get_block(ijk, le);
//						INDEX_TYPE bn = this->block_id_from_block_pos(ijk);
//						INDEX_TYPE in_block_y = y - le.starts[1];
//						INDEX_TYPE in_block_z = t_current_z - le.starts[2];
//						INDEX_TYPE copy_size = m_desired_block_sizes[0];
//						if (x + m_desired_block_sizes[0] > XYZ()[0]) copy_size = XYZ()[0] - x;
//						//printf("   -- x = [%d, %d], in_block=(0, %d, %d)\n", x, x + copy_size, in_block_y, in_block_z);
//						//m_blocks[bn].starts.PrintInt();
//						// address in slab is z sliceinslab + y* slab_x + x
//						dtype* source = &(slab_data[x + this->XYZ()[0] * y + slab_size*z]);
//
//						// address to copy to gets block at 
//						dtype* dest = &(m_data_blocks[bn][le.grid->Index3d({ 0, in_block_y, in_block_z })]);
//
//						//COPY
//						memcpy(dest, source, sizeof(dtype) * copy_size);
//
//						//for (int t = 0; t < copy_size; t++) dest[t] = z;
//						ijk[0]++; // increment block number
//					} // end x loop
//				} // end y loop
//			} // end z loop
//			slab_current_z += slab_width;
//		} //end loop over slabs
//		fclose(fin);
//		return true;
//	}
//
//	unsigned char* blocks_to_include;
//	dtype m_background_value;
//
//	void SetBackgroundValue(dtype val) { m_background_value = val; }
//	const RegularGrid3D* get_global_grid() {
//		return this->m_global_grid.grid;
//	}
//
//	bool has_block(INDEX_TYPE block_num) const {
//		if (blocks_to_include == NULL) return true; // becasue we read dense
//		return blocks_to_include[block_num];
//	}
//	bool has_xyz(Vec3l xyz) const {
//		INDEX_TYPE block_num = this->block_num(xyz);
//		if (blocks_to_include == NULL) return true; // becasue we read dense
//		return blocks_to_include[block_num];
//	}
//
//	bool ReadFromRAWSparse(const char* fname, const char* maskname) {
//
//		// load values from file
//		INDEX_TYPE slab_size = this->XYZ()[0] * this->XYZ()[1];
//		INDEX_TYPE slab_width = 5;
//#pragma omp parallel
//		{
//			slab_width = omp_get_num_threads() * 1;
//		}
//		printf("slab width = %d\n", slab_width);
//		printf("slab size = (%llu x %llu)\n", this->XYZ()[0], this->XYZ()[1]);
//		dtype* slab_data = new dtype[slab_size * slab_width];
//		unsigned char* slab_mask = new unsigned char[slab_size * slab_width];
//
//
//		FILE* fin = fopen(fname, "rb");
//		if (fin == NULL) {
//			printf("Error reading %s\n", fname);
//			return false;
//		}
//
//		FILE* fin_mask = fopen(maskname, "rb");
//		if (fin == NULL) {
//			printf("Error reading %s\n", maskname);
//			fclose(fin);
//			return false;
//		}
//
//		//printf("allocating data blocks\n");
//		//// allocate data blocks
//		//for (INDEX_TYPE i = 0; i < NumBlocks(); i++) {
//		//	//printf("  -- %llu elements\n", m_blocks[i].grid->NumElements());
//		//	m_data_blocks[i] = new dtype[m_blocks[i].grid->NumElements()];
//		//}
//		//printf("done allocating\n");
//
//		// outer file reading loop
//		INDEX_TYPE slab_current_z = 0;
//		INDEX_TYPE block_z_size = this->m_desired_block_sizes[2];
//
//		blocks_to_include = new unsigned char[this->NumBlocks()](); // initialized to 0
//
//		while (slab_current_z < this->XYZ()[2]) {
//			//#define DEBUG_FIRST_SLICE
//#ifdef DEBUG_FIRST_SLICE
//			if (slab_current_z > 0) break;
//#endif
//			// read next slab
//			if (slab_current_z + slab_width > this->XYZ()[2]) {
//				slab_width = this->XYZ()[2] - slab_current_z;
//			}
//			INDEX_TYPE num_to_read = slab_size * slab_width;
//			fread(slab_data, sizeof(dtype), num_to_read, fin);
//			fread(slab_mask, sizeof(unsigned char), num_to_read, fin_mask);
//			// now have read slab_width *X*Y values
//			//printf("read slabs %d\n", slab_current_z);
//			// now make sure all NEEDED blocks are allocated
//			// - for each z in slab, check if data blocks intersecting z have non-zero in mask slab
//			// - will allocate needed data blocks
//			// - can check either non-zero data block pointer or 1 on blocks_to_include -- redundant? 
//			for (int z = 0; z < slab_width; z++) {
//				INDEX_TYPE t_current_z = slab_current_z + z;
//
//				INDEX_TYPE k = this->block_k_from_z(t_current_z);
//#pragma omp for schedule(dynamic)
//				for (auto j = 0; j < this->blocks_per_axis()[1]; j++) {
//					for (auto i = 0; i < this->blocks_per_axis()[0]; i++) {
//						Vec3l ijk(i, j, k);
//						INDEX_TYPE t_block_num = this->block_id_from_block_pos(ijk);
//						// check if ijk block exists - if so, skip
//						if (blocks_to_include[t_block_num] != 0) continue;
//						
//						localExtents le;
//						get_block(ijk, le);
//
//						for (INDEX_TYPE y = 0; y < le.sizes[1]; y++) {
//							auto gy = y + le.starts[1];
//							for (INDEX_TYPE x = 0; x < le.sizes[0]; x++) {
//								auto gx = x + le.starts[0];
//								auto val = slab_mask[gx + this->XYZ()[0] * gy + slab_size*z];
//								if (val != 0) {
//									blocks_to_include[t_block_num] = 1; // mark the block
//									m_data_blocks[t_block_num] = new dtype[le.grid->NumElements()]; // allocate the space
//									y = le.sizes[1]; // skip out of y loop
//									break; // skip out of x loop
//								}
//							}
//						}
//
//
//					}
//				} // end parallel over j blocks
//			} // end iteration over z slices
//
//			  //printf("read masks and allocatd\n");
//			  // copy into blocks
//#pragma omp parallel for schedule(dynamic)
//			for (int z = 0; z < slab_width; z++) {
//				INDEX_TYPE t_current_z = slab_current_z + z;
//
//#ifdef DEBUG_FIRST_SLICE
//				if (z > 0) {
//					slab_current_z++;
//					break;
//				}
//#endif
//				//printf("   -- doing z= %d, slab_current_z = %d\n", z, slab_current_z);
//				// foreach [x, x + blocksize], y stick, copy into block
//				for (int y = 0; y < XYZ()[1]; y++) {
//					//printf("   -- doing y= %d\n", y);
//					// get block number and just add 1 after inner loop
//					Vec3l ijk = this->block_ijk_from_xyz({ 0, y, t_current_z });
//					INDEX_TYPE bn = this->block_id_from_block_pos(ijk);
//					//printf("   -- bn = %d\n", bn);
//					// inner loop with memory copy
//					for (int x = 0; x < XYZ()[0]; x += m_desired_block_sizes[0]) {
//						if (blocks_to_include[bn] == 0) {
//							bn++;
//							ijk[0]++;
//							continue;
//						}
//						localExtents le;
//						get_block(ijk, le);
//
//						INDEX_TYPE in_block_y = y - le.starts[1];
//						INDEX_TYPE in_block_z = t_current_z - le.starts[2];
//						INDEX_TYPE copy_size = m_desired_block_sizes[0];
//						if (x + m_desired_block_sizes[0] > XYZ()[0]) copy_size = XYZ()[0] - x;
//						//printf("   -- x = [%d, %d], in_block=(0, %d, %d)\n", x, x + copy_size, in_block_y, in_block_z);
//						//m_blocks[bn].starts.PrintInt();
//						// address in slab is z sliceinslab + y* slab_x + x
//						dtype* source = &(slab_data[x + this->XYZ()[0] * y + slab_size*z]);
//						unsigned char* mask = &(slab_mask[x + this->XYZ()[0] * y + slab_size*z]);
//
//						// address to copy to gets block at 
//						dtype* dest = &(m_data_blocks[bn][le.grid->Index3d({ 0, in_block_y, in_block_z })]);
//						//COPY
//						//memcpy(dest, source, sizeof(dtype) * copy_size);
//
//						for (int t = 0; t < copy_size; t++) dest[t] = (mask[t] ? source[t] : m_background_value); // masked values
//						bn++; // increment block number
//						ijk[0]++;
//					} // end x loop
//				} // end y loop
//			} // end z loop
//			slab_current_z += slab_width;
//		} //end loop over slabs
//		fclose(fin);
//		return true;
//	}
//
//};


//class SparseBlockedGridFunctionImplicitNotSubclass   {
//protected:
//	dtype** m_data_blocks;
//	const GridDecompType* m_sparse_grid;
//public:
//
//	SparseBlockedGridFunctionImplicitNotSubclass(RegularGrid3D* global_grid, GridDecompType* grid) :
//		m_sparse_grid(grid) {
//	}
//	virtual ~SparseBlockedGridFunctionImplicitNotSubclass() {
//		for (INDEX_TYPE i = 0; i < m_sparse_grid->NumBlocks(); i++) {
//			if (this->has_block(i)) delete[] m_data_blocks[i];
//		}
//		delete[] m_data_blocks;
//	}
//	virtual void initialize() {
//		m_data_blocks = new dtype*[m_sparse_grid->NumBlocks()](); // should value-initialize to zero
//	}
//
//	inline dtype SampleData(Vec3l point) const {
//		auto ijk = m_sparse_grid->block_ijk_from_xyz(point);
//		Vec3l starts;
//		m_sparse_grid->block_start(ijk, starts);
//		Vec3l ind = point - starts;
//		return m_data_blocks[m_sparse_grid->block_id_from_block_pos(ijk)][m_sparse_grid->get_block_grid(ijk)->Index3d(ind)];
//	}
//
//	bool ReadFromRAWDense(const char* fname) {
//
//		// load values from file
//		INDEX_TYPE slab_size = m_sparse_grid->XYZ()[0] * m_sparse_grid->XYZ()[1];
//		INDEX_TYPE slab_width = 5;
//#pragma omp parallel
//		{
//			slab_width = omp_get_num_threads() * 1;
//		}
//		printf("slab width = %d\n", slab_width);
//		printf("slab size = (%llu x %llu)\n", m_sparse_grid->XYZ()[0], m_sparse_grid->XYZ()[1]);
//		dtype* slab_data = new dtype[slab_size * slab_width];
//
//
//		FILE* fin = fopen(fname, "rb");
//		if (fin == NULL) {
//			printf("Error reading %s\n", fname);
//			return false;
//		}
//		printf("allocating data blocks\n");
//		// allocate data blocks
//		for (INDEX_TYPE i = 0; i < m_sparse_grid->NumBlocks(); i++) {
//			//printf("  -- %llu elements\n", m_blocks[i].grid->NumElements());
//			Vec3l ijk = m_sparse_grid->block_ijk_from_id(i);
//			m_data_blocks[i] = new dtype[m_sparse_grid->get_block_grid(ijk)->NumElements()];
//		}
//		printf("done allocating\n");
//
//		// outer file reading loop
//		INDEX_TYPE slab_current_z = 0;
//		INDEX_TYPE block_z_size = m_sparse_grid->m_desired_block_sizes[2];
//
//		while (slab_current_z < m_sparse_grid->XYZ()[2]) {
//			//#define DEBUG_FIRST_SLICE
//#ifdef DEBUG_FIRST_SLICE
//			if (slab_current_z > 0) break;
//#endif
//			// read next slab
//			if (slab_current_z + slab_width > m_sparse_grid->XYZ()[2]) {
//				slab_width = m_sparse_grid->XYZ()[2] - slab_current_z;
//			}
//			INDEX_TYPE num_to_read = slab_size * slab_width;
//			fread(slab_data, sizeof(dtype), num_to_read, fin);
//			//printf(" -- read %llu elements\n", num_to_read);
//			// now have read slab_width *X*Y values
//
//			// copy into blocks
//#pragma omp parallel for schedule(dynamic)
//			for (int z = 0; z < slab_width; z++) {
//				INDEX_TYPE t_current_z = slab_current_z + z;
//#ifdef DEBUG_FIRST_SLICE
//				if (z > 0) {
//					slab_current_z++;
//					break;
//				}
//#endif
//				//printf("   -- doing z= %d, slab_current_z = %d\n", z, slab_current_z);
//				// foreach [x, x + blocksize], y stick, copy into block
//				for (int y = 0; y < m_sparse_grid->XYZ()[1]; y++) {
//					//printf("   -- doing y= %d\n", y);
//					// get block number and just add 1 after inner loop
//					Vec3l ijk = m_sparse_grid->block_ijk_from_xyz({ 0, y, t_current_z });
//					//printf("   -- bn = %d\n", bn);
//					// inner loop with memory copy
//					for (int x = 0; x < m_sparse_grid->XYZ()[0]; x += m_sparse_grid->m_desired_block_sizes[0]) {
//						GridDecompType::localExtents le;
//						m_sparse_grid->get_block(ijk, le);
//						INDEX_TYPE bn = m_sparse_grid->block_id_from_block_pos(ijk);
//						INDEX_TYPE in_block_y = y - le.starts[1];
//						INDEX_TYPE in_block_z = t_current_z - le.starts[2];
//						INDEX_TYPE copy_size = m_sparse_grid->m_desired_block_sizes[0];
//						if (x + m_sparse_grid->m_desired_block_sizes[0] > m_sparse_grid->XYZ()[0]) 
//							copy_size = m_sparse_grid->XYZ()[0] - x;
//						//printf("   -- x = [%d, %d], in_block=(0, %d, %d)\n", x, x + copy_size, in_block_y, in_block_z);
//						//m_blocks[bn].starts.PrintInt();
//						// address in slab is z sliceinslab + y* slab_x + x
//						dtype* source = &(slab_data[x + m_sparse_grid->XYZ()[0] * y + slab_size*z]);
//
//						// address to copy to gets block at 
//						dtype* dest = &(m_data_blocks[bn][le.grid->Index3d({ 0, in_block_y, in_block_z })]);
//
//						//COPY
//						memcpy(dest, source, sizeof(dtype) * copy_size);
//
//						//for (int t = 0; t < copy_size; t++) dest[t] = z;
//						ijk[0]++; // increment block number
//					} // end x loop
//				} // end y loop
//			} // end z loop
//			slab_current_z += slab_width;
//		} //end loop over slabs
//		fclose(fin);
//		return true;
//	}
//
//	unsigned char* blocks_to_include;
//	dtype m_background_value;
//
//	void SetBackgroundValue(dtype val) { m_background_value = val; }
//	const RegularGrid3D* get_global_grid() {
//		return m_sparse_grid->m_global_grid.grid;
//	}
//
//	bool has_block(INDEX_TYPE block_num) const {
//		if (blocks_to_include == NULL) return true; // becasue we read dense
//		return blocks_to_include[block_num];
//	}
//	bool has_xyz(Vec3l xyz) const {
//		INDEX_TYPE block_num = m_sparse_grid->block_num(xyz);
//		if (blocks_to_include == NULL) return true; // becasue we read dense
//		return blocks_to_include[block_num];
//	}
//
//	bool ReadFromRAWSparse(const char* fname, const char* maskname) {
//
//		// load values from file
//		INDEX_TYPE slab_size = m_sparse_grid->XYZ()[0] * m_sparse_grid->XYZ()[1];
//		INDEX_TYPE slab_width = 5;
//#pragma omp parallel
//		{
//			slab_width = omp_get_num_threads() * 1;
//		}
//		printf("slab width = %d\n", slab_width);
//		printf("slab size = (%llu x %llu)\n", m_sparse_grid->XYZ()[0], m_sparse_grid->XYZ()[1]);
//		dtype* slab_data = new dtype[slab_size * slab_width];
//		unsigned char* slab_mask = new unsigned char[slab_size * slab_width];
//
//
//		FILE* fin = fopen(fname, "rb");
//		if (fin == NULL) {
//			printf("Error reading %s\n", fname);
//			return false;
//		}
//
//		FILE* fin_mask = fopen(maskname, "rb");
//		if (fin == NULL) {
//			printf("Error reading %s\n", maskname);
//			fclose(fin);
//			return false;
//		}
//
//		//printf("allocating data blocks\n");
//		//// allocate data blocks
//		//for (INDEX_TYPE i = 0; i < NumBlocks(); i++) {
//		//	//printf("  -- %llu elements\n", m_blocks[i].grid->NumElements());
//		//	m_data_blocks[i] = new dtype[m_blocks[i].grid->NumElements()];
//		//}
//		//printf("done allocating\n");
//
//		// outer file reading loop
//		INDEX_TYPE slab_current_z = 0;
//		INDEX_TYPE block_z_size = m_sparse_grid->m_desired_block_sizes[2];
//
//		blocks_to_include = new unsigned char[m_sparse_grid->NumBlocks()](); // initialized to 0
//
//		while (slab_current_z < m_sparse_grid->XYZ()[2]) {
//			//#define DEBUG_FIRST_SLICE
//#ifdef DEBUG_FIRST_SLICE
//			if (slab_current_z > 0) break;
//#endif
//			// read next slab
//			if (slab_current_z + slab_width > m_sparse_grid->XYZ()[2]) {
//				slab_width = m_sparse_grid->XYZ()[2] - slab_current_z;
//			}
//			INDEX_TYPE num_to_read = slab_size * slab_width;
//			fread(slab_data, sizeof(dtype), num_to_read, fin);
//			fread(slab_mask, sizeof(unsigned char), num_to_read, fin_mask);
//			// now have read slab_width *X*Y values
//			//printf("read slabs %d\n", slab_current_z);
//			// now make sure all NEEDED blocks are allocated
//			// - for each z in slab, check if data blocks intersecting z have non-zero in mask slab
//			// - will allocate needed data blocks
//			// - can check either non-zero data block pointer or 1 on blocks_to_include -- redundant? 
//			for (int z = 0; z < slab_width; z++) {
//				INDEX_TYPE t_current_z = slab_current_z + z;
//
//				INDEX_TYPE k = m_sparse_grid->block_k_from_z(t_current_z);
//#pragma omp for schedule(dynamic)
//				for (auto j = 0; j < m_sparse_grid->blocks_per_axis()[1]; j++) {
//					for (auto i = 0; i < m_sparse_grid->blocks_per_axis()[0]; i++) {
//						Vec3l ijk(i, j, k);
//						INDEX_TYPE t_block_num = m_sparse_grid->block_id_from_block_pos(ijk);
//						// check if ijk block exists - if so, skip
//						if (blocks_to_include[t_block_num] != 0) continue;
//
//						GridDecompType::localExtents le;
//						m_sparse_grid->get_block(ijk, le);
//
//						for (INDEX_TYPE y = 0; y < le.sizes[1]; y++) {
//							auto gy = y + le.starts[1];
//							for (INDEX_TYPE x = 0; x < le.sizes[0]; x++) {
//								auto gx = x + le.starts[0];
//								auto val = slab_mask[gx + m_sparse_grid->XYZ()[0] * gy + slab_size*z];
//								if (val != 0) {
//									blocks_to_include[t_block_num] = 1; // mark the block
//									m_data_blocks[t_block_num] = new dtype[le.grid->NumElements()]; // allocate the space
//									y = le.sizes[1]; // skip out of y loop
//									break; // skip out of x loop
//								}
//							}
//						}
//
//
//					}
//				} // end parallel over j blocks
//			} // end iteration over z slices
//
//			  //printf("read masks and allocatd\n");
//			  // copy into blocks
//#pragma omp parallel for schedule(dynamic)
//			for (int z = 0; z < slab_width; z++) {
//				INDEX_TYPE t_current_z = slab_current_z + z;
//
//#ifdef DEBUG_FIRST_SLICE
//				if (z > 0) {
//					slab_current_z++;
//					break;
//				}
//#endif
//				//printf("   -- doing z= %d, slab_current_z = %d\n", z, slab_current_z);
//				// foreach [x, x + blocksize], y stick, copy into block
//				for (int y = 0; y < m_sparse_grid->XYZ()[1]; y++) {
//					//printf("   -- doing y= %d\n", y);
//					// get block number and just add 1 after inner loop
//					Vec3l ijk = m_sparse_grid->block_ijk_from_xyz({ 0, y, t_current_z });
//					INDEX_TYPE bn = m_sparse_grid->block_id_from_block_pos(ijk);
//					//printf("   -- bn = %d\n", bn);
//					// inner loop with memory copy
//					for (int x = 0; x < m_sparse_grid->XYZ()[0]; x += m_sparse_grid->m_desired_block_sizes[0]) {
//						if (blocks_to_include[bn] == 0) {
//							bn++;
//							ijk[0]++;
//							continue;
//						}
//						GridDecompType::localExtents le;
//						m_sparse_grid->get_block(ijk, le);
//
//						INDEX_TYPE in_block_y = y - le.starts[1];
//						INDEX_TYPE in_block_z = t_current_z - le.starts[2];
//						INDEX_TYPE copy_size = m_sparse_grid->m_desired_block_sizes[0];
//						if (x + m_sparse_grid->m_desired_block_sizes[0] > m_sparse_grid->XYZ()[0]) 
//							copy_size = m_sparse_grid->XYZ()[0] - x;
//						//printf("   -- x = [%d, %d], in_block=(0, %d, %d)\n", x, x + copy_size, in_block_y, in_block_z);
//						//m_blocks[bn].starts.PrintInt();
//						// address in slab is z sliceinslab + y* slab_x + x
//						dtype* source = &(slab_data[x + m_sparse_grid->XYZ()[0] * y + slab_size*z]);
//						unsigned char* mask = &(slab_mask[x + m_sparse_grid->XYZ()[0] * y + slab_size*z]);
//
//						// address to copy to gets block at 
//						dtype* dest = &(m_data_blocks[bn][le.grid->Index3d({ 0, in_block_y, in_block_z })]);
//						//COPY
//						//memcpy(dest, source, sizeof(dtype) * copy_size);
//
//						for (int t = 0; t < copy_size; t++) dest[t] = (mask[t] ? source[t] : m_background_value); // masked values
//						bn++; // increment block number
//						ijk[0]++;
//					} // end x loop
//				} // end y loop
//			} // end z loop
//			slab_current_z += slab_width;
//		} //end loop over slabs
//		fclose(fin);
//		return true;
//	}
//
//};



std::chrono::steady_clock::time_point g_start_time;


//#define USE_IDX

#ifdef USE_IDX

int main(int argc, char** argv) {

	int block_x, block_y, block_z;

	if (argc < 5) { printf("Usage: filename.idx per_x per_y per_z  [outputdebug=0]\n"); return 0; }

	auto dataset_url = std::string(argv[1]);
	sscanf(argv[2], "%d", &block_x);
	sscanf(argv[3], "%d", &block_y);
	sscanf(argv[4], "%d", &block_z);
	int outputdebug = 0;
	if (argc >= 6)
		sscanf(argv[5], "%d", &outputdebug);
	omp_set_num_threads(outputdebug);

	DbModule::attach();

	DoAtExit do_at_exit([] {
		DbModule::detach();
	});

	auto dataset = LoadDataset(dataset_url);
	VisusReleaseAssert(dataset);

	auto world_box = dataset->getLogicBox();
	auto max_resolution = dataset->getMaxResolution();
	PrintInfo("Data size:", world_box.p1.toString(), world_box.p2.toString(), "max res:", max_resolution);


	printf("making decomposed grid\n");
	RegularGrid3D* single_grid = new RegularGrid3D(world_box.size, { 0,0,0 });
	RegularGrid3DDecomposition* grid =
		new RegularGrid3DDecomposition(world_box.p2[0], world_box.p2[1], world_box.p2[2], block_x, block_y, block_z, single_grid);
	printf(" -- decomposing...\n");
	grid->decompose();
	printf(" -- adding ghost....\n");
	//grid->AddFullGhostLayer();
	int numblocks = grid->numBlocks();
	printf(" -- made %d blocks\n", numblocks);
	float** blocks = new float*[numblocks];





	g_start_time = std::chrono::steady_clock::now();

#pragma omp for schedule(dynamic)
	for (int i = 0; i < numblocks; i++) {
		auto& lstarts = grid->m_blocks[i].starts;
		auto& lsizes = grid->m_blocks[i].sizes;
		blocks[i] = GetSubBlockIdx(dataset, lstarts, lsizes);
	}

	printf("TIMING: %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - g_start_time).count());


	return 1;


}

#else 

//void do_block_maxvl(INDEX_TYPE block_num,
//	TopologicalSparseBlockedMesh* blocked_mesh,
//	SparseBlockedLabeling<unsigned char>* maxvl,
//	SparseBlockedLabeling<float>* data) {
//
//	auto block_ijk = blocked_mesh->block_ijk_from_id(block_num);
//	auto sub_mesh = blocked_mesh->get_submesh(block_num);
//	auto sub_grid = data->GetGrid()->get_block_grid(block_num);
//	
//	auto const xyz = sub_grid->XYZ(); // number of data points in grid
//	
//	// DO INTERNAL
//
//	// now check each direction
//	
//	// has in each direction....
//
//	// if has +x do slab
//
//	// if has +y do slab
//
//	// if has +z do slab
//
//	// if has +x +y +xy stick of quads n hexs
//
//	// for each yz position do fastest stick
//	for (INDEX_TYPE y = 0; y < xyz[1]; y++) {
//		for (INDEX_TYPE z = 0; z < xyz[2]; z++) {
//
//			//now do fastest stick
//	
//
//
//		}
//	}
//
//}

//void do_max_vert_label(TopologicalSparseBlockedMesh* blocked_mesh,
//	SparseBlockedLabeling<unsigned char>* maxvl,
//	SparseBlockedLabeling<float>* data) {
//
//	// for each block do fastest stick
//	
//	// gather block nums
//	std::vector<INDEX_TYPE> present_blocks;
//	for (INDEX_TYPE i = 0; i < blocked_mesh->NumBlocks(); i++) {
//		if (data->has_block(i)) present_blocks.push_back(i);
//	}
//	INDEX_TYPE num_blocks = present_blocks.size();
//	printf("about to do maxvl over %lld blocks...\n", num_blocks);
//
//	// first do fastest stick on each block
//#pragma omp for schedule(dynamic)
//	for (INDEX_TYPE bn = 0; bn < num_blocks; bn++) {
//		auto block_num = present_blocks[bn];
//		do_block_maxvl(block_num, blocked_mesh, maxvl, data);
//	}
//
//
//
//
//}


template<class GridDecompType>
void watershed_phase_1(const GridDecompType* decomp,
	const SparseBlockedLabeling<float, GridDecompType>* scalar_function,
	SparseBlockedLabeling<typename GridDecompType::BN_ID_PAIR, GridDecompType>* watershed) {
	typedef GridDecompType::neighborhood6_iterator GridDecompNeighborIteratorType;
	printf("pass 1, finding highest neighbors\n");
	INDEX_TYPE data_counter = 0;
	INDEX_TYPE crit_counter = 0;
	int NUM_BLOCKS = decomp->NumBlocks();
#pragma omp parallel for schedule(dynamic)
	for (INDEX_TYPE bn = 0; bn < NUM_BLOCKS; bn++) {
		// skip nonexistest blocks
		if (!scalar_function->has_block(bn)) continue;
		auto in_b_iter = GridDecompType::in_block_iterator(bn, decomp);
		auto n_iter = GridDecompNeighborIteratorType(decomp);
		for (in_b_iter.begin(); in_b_iter.valid(); in_b_iter.advance()) {
			auto id_p = in_b_iter.value();
			// skip nonexistent voxels
			if (!scalar_function->has_bnidpair(id_p)) continue;

			auto highest_idp = id_p; // set highest to itself
			float highest_val = scalar_function->GetLabelBNID(id_p);
			for (n_iter.begin(in_b_iter); n_iter.valid(); n_iter.advance()) {
				auto n_id_p = n_iter.value();
				// skip non-existent neighbors
				if (!scalar_function->has_bnidpair(n_id_p)) continue;
				// we HAVE this neighbor
				auto tmp_val2 = scalar_function->GetLabelBNID(n_id_p);
				if (tmp_val2 > highest_val ||
					(tmp_val2 == highest_val && n_id_p > highest_idp)) {
					highest_idp = n_id_p;
					highest_val = tmp_val2;
				}

				// pointer to highest neighbor
			}
			data_counter++;
			crit_counter += id_p == highest_idp;
			watershed->SetLabelBNID(id_p, highest_idp);
		}
	}
	printf("found %lld crits from %lld data points\n", crit_counter, data_counter);
}

template<class GridDecompType>
void watershed_phase_2(const GridDecompType* decomp,
	SparseBlockedLabeling<typename GridDecompType::BN_ID_PAIR, GridDecompType>* watershed) {
	typedef GridDecompType::neighborhood6_iterator GridDecompNeighborIteratorType;
	printf("pass 2, path compress ids\n");
	int NUM_BLOCKS = decomp->NumBlocks();
#pragma omp parallel for schedule(dynamic)
	for (INDEX_TYPE bn = 0; bn < NUM_BLOCKS; bn++) {		// skip nonexistest blocks
		if (!watershed->has_block(bn)) continue;
		std::vector<typename GridDecompType::BN_ID_PAIR> path;
		path.reserve(1024);
		auto in_b_iter = GridDecompType::in_block_iterator(bn, decomp);
		for (in_b_iter.begin(); in_b_iter.valid(); in_b_iter.advance()) {
			auto id_p = in_b_iter.value();
			// skip nonexistent voxels
			if (!watershed->has_bnidpair(id_p)) continue;

			path.clear();
			auto start = id_p;
			auto next = watershed->GetLabelBNID(start);

			while (true) {
				if (start == next) {
					// this is critical
					for (auto p : path) {
						watershed->SetLabelBNID(p, start);
					}
					break;
				}
				path.push_back(start);
				//INDEX_TYPE debug_bn, debug_id;
				//decomp->GetBNAndIDFromPair(next, debug_bn, debug_id);
				start = next;
				next = watershed->GetLabelBNID(start);
			}
		}
	}
}

template<unsigned char BSX, unsigned char BSY, unsigned char BSZ> int do_work(int X, int  Y, int Z, std::string filename, std::string maskname);



int(*func_array[125])(int, int, int, std::string, std::string) = {
	do_work<1,1,1>,
	do_work<2,1,1>,
	do_work<3,1,1>,
	do_work<4,1,1>,
	do_work<5,1,1>,

	do_work<1,2,1>,
	do_work<2,2,1>,
	do_work<3,2,1>,
	do_work<4,2,1>,
	do_work<5,2,1>,

	do_work<1,3,1>,
	do_work<2,3,1>,
	do_work<3,3,1>,
	do_work<4,3,1>,
	do_work<5,3,1>,

	do_work<1,4,1>,
	do_work<2,4,1>,
	do_work<3,4,1>,
	do_work<4,4,1>,
	do_work<5,4,1>,

	do_work<1,5,1>,
	do_work<2,5,1>,
	do_work<3,5,1>,
	do_work<4,5,1>,
	do_work<5,5,1>,

	do_work<1,1,2>,
	do_work<2,1,2>,
	do_work<3,1,2>,
	do_work<4,1,2>,
	do_work<5,1,2>,

	do_work<1,2,2>,
	do_work<2,2,2>,
	do_work<3,2,2>,
	do_work<4,2,2>,
	do_work<5,2,2>,

	do_work<1,3,2>,
	do_work<2,3,2>,
	do_work<3,3,2>,
	do_work<4,3,2>,
	do_work<5,3,2>,

	do_work<1,4,2>,
	do_work<2,4,2>,
	do_work<3,4,2>,
	do_work<4,4,2>,
	do_work<5,4,2>,

	do_work<1,5,2>,
	do_work<2,5,2>,
	do_work<3,5,2>,
	do_work<4,5,2>,
	do_work<5,5,2>,

	do_work<1,1,3>,
	do_work<2,1,3>,
	do_work<3,1,3>,
	do_work<4,1,3>,
	do_work<5,1,3>,

	do_work<1,2,3>,
	do_work<2,2,3>,
	do_work<3,2,3>,
	do_work<4,2,3>,
	do_work<5,2,3>,

	do_work<1,3,3>,
	do_work<2,3,3>,
	do_work<3,3,3>,
	do_work<4,3,3>,
	do_work<5,3,3>,

	do_work<1,4,3>,
	do_work<2,4,3>,
	do_work<3,4,3>,
	do_work<4,4,3>,
	do_work<5,4,3>,

	do_work<1,5,3>,
	do_work<2,5,3>,
	do_work<3,5,3>,
	do_work<4,5,3>,
	do_work<5,5,3>,

	do_work<1,1,4>,
	do_work<2,1,4>,
	do_work<3,1,4>,
	do_work<4,1,4>,
	do_work<5,1,4>,

	do_work<1,2,4>,
	do_work<2,2,4>,
	do_work<3,2,4>,
	do_work<4,2,4>,
	do_work<5,2,4>,

	do_work<1,3,4>,
	do_work<2,3,4>,
	do_work<3,3,4>,
	do_work<4,3,4>,
	do_work<5,3,4>,

	do_work<1,4,4>,
	do_work<2,4,4>,
	do_work<3,4,4>,
	do_work<4,4,4>,
	do_work<5,4,4>,

	do_work<1,5,4>,
	do_work<2,5,4>,
	do_work<3,5,4>,
	do_work<4,5,4>,
	do_work<5,5,4>,

		do_work<1, 1, 5>,
		do_work<2, 1, 5>,
		do_work<3, 1, 5>,
		do_work<4, 1, 5>,
		do_work<5, 1, 5>,

		do_work<1, 2, 5>,
		do_work<2, 2, 5>,
		do_work<3, 2, 5>,
		do_work<4, 2, 5>,
		do_work<5, 2, 5>,

		do_work<1, 3, 5>,
		do_work<2, 3, 5>,
		do_work<3, 3, 5>,
		do_work<4, 3, 5>,
		do_work<5, 3, 5>,

		do_work<1, 4, 5>,
		do_work<2, 4, 5>,
		do_work<3, 4, 5>,
		do_work<4, 4, 5>,
		do_work<5, 4, 5>,

		do_work<1, 5, 5>,
		do_work<2, 5, 5>,
		do_work<3, 5, 5>,
		do_work<4, 5, 5>,
		do_work<5, 5, 5>
};

int main(int argc, char** argv) {


#if 0
	printf("a\n");

	FILE* fout = fopen("index_5x5x5.raw", "wb");
	if (fout != NULL) {
		printf("b\n");
		float* vals = new float[5 * 5 * 5];
		printf("c\n");
		for (int i = 0; i < 5 * 5 * 5; i++) {
			vals[i] = i;
		}
		printf("d\n");
		fwrite(vals, sizeof(float), 5 * 5 * 5, fout);
		fclose(fout);
	}
	printf("done\n");
	return 1;
#endif

	ThreadedTimer timer(1);
	timer.StartGlobal();

	int X, Y, Z;
	int block_x, block_y, block_z;
	std::string filename;
	std::string maskname;

#define USE_MASK
#ifdef USE_MASK
	if (argc < 8) { printf("Usage: X Y Z filename per_x per_y per_z  [outputdebug=0]\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
	maskname = std::string(argv[5]);
	sscanf(argv[6], "%d", &block_x);
	sscanf(argv[7], "%d", &block_y);
	sscanf(argv[8], "%d", &block_z);
	int outputdebug = 0;
	if (argc >= 10)
		sscanf(argv[9], "%d", &outputdebug);
	omp_set_num_threads(outputdebug);

	printf("%d, %d, %d, %d, %d, %d, %d\n", X, Y, Z, block_x, block_y, block_z, outputdebug);

	int func_id = (block_x - 1) + 5 * (block_y - 1) + 25 * (block_z - 1);
	printf("doing function #%d\n", func_id);
	auto work_func = func_array[func_id](X, Y, Z, filename, maskname);

};

template<unsigned char BSX, unsigned char BSY, unsigned char BSZ> int do_work(int X, int  Y, int Z, std::string filename, std::string maskname) {
	//printf("making decomposed grid\n");
	//RegularGrid3D* single_grid = new RegularGrid3D({ X, Y, Z }, { 0,0,0 });
	//SparseBlockedGridFunction<float>* data_grid = new SparseBlockedGridFunction<float>(X, Y, Z, block_x, block_y, block_z, single_grid);

	typedef GInt::RegularGrid3DDecompositionFixed<BSX, BSY, BSZ> GridDecompType;
	typedef GridDecompType::neighborhood6_iterator GridDecompNeighborIteratorType;
#else
	if (argc < 8) { printf("Usage: X Y Z filename per_x per_y per_z  [outputdebug=0]\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
	sscanf(argv[5], "%d", &block_x);
	sscanf(argv[6], "%d", &block_y);
	sscanf(argv[7], "%d", &block_z);
	int outputdebug = 0;
	if (argc >= 9)
		sscanf(argv[8], "%d", &outputdebug);
	omp_set_num_threads(outputdebug);

	printf("%d, %d, %d, %d, %d, %d, %d\n", X, Y, Z, block_x, block_y, block_z, outputdebug);
	g_start_time = std::chrono::steady_clock::now();

	printf("making decomposed grid\n");
	RegularGrid3D* single_grid = new RegularGrid3D({ X, Y, Z }, { 0,0,0 });
	//SparseBlockedGridFunction<float>* data_grid = new SparseBlockedGridFunction<float>(X, Y, Z, block_x, block_y, block_z, single_grid);
	SparseBlockedGridFunction* data_grid = new SparseBlockedGridFunction(X, Y, Z, block_x, block_y, block_z, single_grid);

	printf("got here after new data_grid\n");

	data_grid->decompose();
	printf("got here after decompose\n");
	data_grid->ReadFromRAWDense(filename.c_str());
#endif

#if 0
	// IMPLICIT GRID LOADING
	for (int i = 0; i < 2; i++) {
		g_start_time = std::chrono::steady_clock::now();

		SparseBlockedGridFunctionImplicit* scalar_function = new SparseBlockedGridFunctionImplicit(X, Y, Z, block_x, block_y, block_z, single_grid);
		scalar_function->SetBackgroundValue(NumericLimits<float>::highest());
		printf("got here after new data_grid\n");
		scalar_function->decompose();
		printf("got here after decompose\n");
		scalar_function->ReadFromRAWSparse(filename.c_str(), maskname.c_str());
		printf("========> explicit %d %d %d\n", i, block_x,
			std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - g_start_time).count());
		delete scalar_function;
	}
#endif
	// IMPLICIT NONSUBGRID GRID LOADING
	for (int i = 0; i < 1; i++) {
		g_start_time = std::chrono::steady_clock::now();
		GridDecompType* decomp = new GridDecompType(X, Y, Z/*, block_x, block_y, block_z*/);
		decomp->decompose();
		printf("created decomposition into %dx%dx%d cubes:\n",
			decomp->get_block_grid()->XYZ()[0],
			decomp->get_block_grid()->XYZ()[1],
			decomp->get_block_grid()->XYZ()[2]);
		printf("num blocks per dim: ");
		decomp->blocks_per_axis().PrintInt();

		auto bs = decomp->get_block_grid()->XYZ();

		//printf("initializing data grid...\n");
		SparseBlockedLabeling<float, GridDecompType>* scalar_function =
			new SparseBlockedLabeling<float, GridDecompType>(decomp);
		scalar_function->SetMaskValue(NumericLimits<float>::highest());
		scalar_function->initialize();
		//printf("reading data and mask...\n");
		scalar_function->ReadFromRAWSparse(filename.c_str(), maskname.c_str());
		//scalar_function->ReadFromRAWDense(filename.c_str());
		printf("TIMING bs %d %d %d data_read %d\n", bs[0], bs[1], bs[2],
			std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - g_start_time).count());


		// debug block iterator
		if (false) {
			auto iter_time = std::chrono::steady_clock::now();
			int count = 0;
			int block_count = 0;
			printf("iterating existing neighbors of existing items\n");
			for (INDEX_TYPE bn = 0; bn < decomp->NumBlocks(); bn++) {
				if (!scalar_function->has_block(bn)) continue;
				block_count++;
				auto in_b_iter = GridDecompType::in_block_iterator(bn, decomp);
				auto n_iter = GridDecompNeighborIteratorType(decomp);
				for (in_b_iter.begin(); in_b_iter.valid(); in_b_iter.advance()) {
					auto id_p = in_b_iter.value();
					if (scalar_function->has_bnidpair(id_p)) {
						count++;
						for (n_iter.begin(in_b_iter); n_iter.valid(); n_iter.advance()) {
							auto n_id_p = n_iter.value();
							if (!scalar_function->has_bnidpair(n_id_p)) continue; // skip nonexistent neighbors
						}
					}
				}
			}
			printf("has %d values int %d blocks\n", count, block_count);
			printf("========> iteration time: %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - iter_time).count());
		}


		// sanity checking index computation
		if (false) {
			printf("beginning sanity check:\n");
			INDEX_TYPE hard_id = 0;
			for (INDEX_TYPE z = 0; z < Z; z++) {
				for (INDEX_TYPE y = 0; y < Y; y++) {
					for (INDEX_TYPE x = 0; x < X; x++) {
						INDEX_TYPE comp_id = decomp->Index3d({ x,y,z });
						if (comp_id != hard_id) {
							printf("ERROR: hard %lld != %lld comp\n", hard_id, comp_id);
						}
						Vec3l xyz = decomp->Coords(hard_id);
						if (xyz != Vec3l(x, y, z)) {
							printf("ERROR: hard (%lld, %lld, %lld) != (%lld, %lld, %lld) comp\n",
								x, y, z, xyz[0], xyz[1], xyz[2]);
						}

						auto bnid = decomp->get_bnid_from_id(hard_id);
						INDEX_TYPE ind_id = decomp->GetGlobalIdFromBNID(bnid);

						if (hard_id != ind_id) {
							printf("ERROR: hard %lld != %lld ind\n", hard_id, ind_id);
						}

						INDEX_TYPE bn, id;
						decomp->GetBNAndIDFromPair(bnid, bn, id);
						auto ijk = decomp->block_ijk_from_xyz(xyz);
						Vec3l starts; decomp->block_start(ijk, starts);
						auto in_b_xyz = xyz - starts;
						INDEX_TYPE BN = decomp->block_id_from_block_pos(ijk);
						INDEX_TYPE INBID = decomp->get_block_grid(ijk)->Index3d(in_b_xyz);
						INDEX_TYPE tbn, tid;
						decomp->GetBNAndIDFromPair(decomp->MakeBNPair(BN, INBID), tbn, tid);
						//printf("(%lld, %lld) => (%lld, %lld)\n", BN, INBID, tbn, tid);
						if (bn != decomp->block_id_from_block_pos(ijk)) {
							printf("WHOA: %lld != %lld block num\n", bn, decomp->block_id_from_block_pos(ijk));
							printf("\thardid=%lld, inbid=%lld, bn=%lld, (%d, %d)\n", hard_id, id, bn, ((int*)&bnid)[0], ((int*)&bnid)[1]);
							auto bnid2 = decomp->MakeBNPair(decomp->block_id_from_block_pos(ijk), decomp->get_block_grid(ijk)->Index3d(in_b_xyz));
							printf("\t(%d, %d)\n", ((int*)&bnid2)[0], ((int*)&bnid2)[1]);
							xyz.PrintInt();
							ijk.PrintInt();
							starts.PrintInt();
							in_b_xyz.PrintInt();
						}

						if (id != decomp->get_block_grid(ijk)->Index3d(in_b_xyz)) {
							printf("WHOA: inid %lld != %lld\n", id, decomp->get_block_grid(ijk)->Index3d(in_b_xyz));
							xyz.PrintInt();
							ijk.PrintInt();
							starts.PrintInt();
							in_b_xyz.PrintInt();
						}
						//printf("dong id %lld:\n", hard_id);
						//xyz.PrintInt();
						//ijk.PrintInt();
						//starts.PrintInt();
						//in_b_xyz.PrintInt();

						hard_id++;
					}
				}
			}
			printf("end sanity check\n");
		}


		if (false) {
			char* debugbuffer = new char[X*Y*Z];

			auto fill_time = std::chrono::steady_clock::now();
			for (long long i = 0; i < X*Y*Z; i++) {
				debugbuffer[i] = 0;
				auto bnid = decomp->get_bnid_from_id(i);
				INDEX_TYPE bn, in_b_id;
				decomp->GetBNAndIDFromPair(bnid, bn, in_b_id);
				if (scalar_function->has_block(bn)) {
					debugbuffer[i] = 1;
				}
				if (!scalar_function->has_block(bn) && scalar_function->has_bnid(bn, in_b_id)) {
					printf("WHWHWHWHEWHWHWOA: %lld, %lld\n", bn, in_b_id);
				}
				if (scalar_function->has_bnid(bn, in_b_id)) {
					debugbuffer[i] = 2;
				}
			}
			printf("========> fill time: %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - fill_time).count());

			auto biter_time = std::chrono::steady_clock::now();
			char* debugbuffer2 = new char[X*Y*Z];
			memset(debugbuffer2, 0, X*Y*Z * sizeof(char));
			for (INDEX_TYPE bn = 0; bn < decomp->NumBlocks(); bn++) {
				if (!scalar_function->has_block(bn)) continue;


				auto in_b_iter = GridDecompType::in_block_iterator(bn, decomp);
				for (in_b_iter.begin(); in_b_iter.valid(); in_b_iter.advance()) {
					auto id_p = in_b_iter.value();
					auto gid = decomp->GetGlobalIdFromBNID(id_p);
					if (scalar_function->has_bnidpair(id_p)) {
						debugbuffer2[gid] = 2;
					}
					else {
						debugbuffer2[gid] = 1;
					}
				}
			}
			printf("========> biter time: %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - biter_time).count());

			printf("comparing buffers:\n");
			for (long long i = 0; i < X*Y*Z; i++) {
				if (false && debugbuffer[i] != debugbuffer2[i]) {
					printf("debugbuffer %llu:%d != %llu:%d debugbuffer2\n", i, debugbuffer[i], i, debugbuffer2[i]);
					auto bnid = decomp->get_bnid_from_id(i);
					decomp->PrintBNIDInfo(bnid);
				}

			}
			printf("done comparing buffers\n");
			char foutname[2048];
			sprintf(foutname, "debug1_%d_%dx%dx%d_%dx%dx%d.raw",
				i, decomp->get_block_grid()->XYZ()[0],
				decomp->get_block_grid()->XYZ()[1],
				decomp->get_block_grid()->XYZ()[2],
				X, Y, Z);
			FILE* fout = fopen(foutname, "wb");
			fwrite(debugbuffer, sizeof(char), X*Y*Z, fout);
			fclose(fout);
			sprintf(foutname, "debug2_%d_%dx%dx%d_%dx%dx%d.raw",
				i, decomp->get_block_grid()->XYZ()[0],
				decomp->get_block_grid()->XYZ()[1],
				decomp->get_block_grid()->XYZ()[2],
				X, Y, Z);
			fout = fopen(foutname, "wb");
			fwrite(debugbuffer2, sizeof(char), X*Y*Z, fout);
			fclose(fout);
			delete[] debugbuffer;
			delete[] debugbuffer2;
		}

		// test reciprocality in neighbor iteration (does my neighbor have me as neighbor)
		if (false) {
			printf("testing neighbor iteration reciprical\n");
			auto ncount_time = std::chrono::steady_clock::now();
			for (INDEX_TYPE bn = 0; bn < decomp->NumBlocks(); bn++) {
				if (!scalar_function->has_block(bn)) continue;
				auto in_b_iter = GridDecompType::in_block_iterator(bn, decomp);
				auto n_iter = GridDecompNeighborIteratorType(decomp);
				auto n_n_iter = GridDecompNeighborIteratorType(decomp);
				for (in_b_iter.begin(); in_b_iter.valid(); in_b_iter.advance()) {
					auto id_p = in_b_iter.value();
					if (!scalar_function->has_bnidpair(id_p)) continue;
					// iterate over neighbors
					for (n_iter.begin(in_b_iter); n_iter.valid(); n_iter.advance()) {
						auto n_id_p = n_iter.value();
						if (!scalar_function->has_bnidpair(n_id_p)) continue;
						int has_reciprical = 0;
						for (n_n_iter.begin(n_id_p); n_n_iter.valid(); n_n_iter.advance()) {
							auto n_n_id_p = n_n_iter.value();
							if (n_n_id_p == id_p) {
								// the neighbor of my neighbor is me!!
								has_reciprical++;
							}
						}
						if (has_reciprical != 1) {
							// print some debug info
							printf("Error: %lld neighbor %lld does not have reciprical\n", id_p, n_id_p);
							printf("base_bnid: ");
							decomp->PrintBNIDInfo(id_p);
							printf("neig_bnid: ");
							decomp->PrintBNIDInfo(n_id_p);

							printf("\nbase_neighbors:\n");
							for (n_iter.begin(in_b_iter); n_iter.valid(); n_iter.advance()) {
								auto tn_id_p = n_iter.value();
								decomp->PrintBNIDInfo(tn_id_p);
							}

							printf("\nneig_neighbors:\n");
							for (n_n_iter.begin(n_id_p); n_n_iter.valid(); n_n_iter.advance()) {
								auto tn_n_id_p = n_n_iter.value();
								decomp->PrintBNIDInfo(tn_n_id_p);
							}
							printf("endcompare\n");
						}
					}
				}
			}
			printf("========> reciprical time: %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - ncount_time).count());

		}
		// test neighbor iteration
		if (false) {
			printf("allocating structure for neg count:\n");
			auto alloc_time = std::chrono::steady_clock::now();
			typedef GridDecompType::BN_ID_PAIR BNIDPair;
			SparseBlockedLabeling<char, GridDecompType>* ncount =
				new SparseBlockedLabeling<char, GridDecompType>(decomp);
			ncount->SetMaskValue(0);
			ncount->initialize<SparseBlockedLabeling<float, GridDecompType>>(scalar_function);
			printf("========> alloclowest time: %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - alloc_time).count());
			auto ncount_time = std::chrono::steady_clock::now();
			for (INDEX_TYPE bn = 0; bn < decomp->NumBlocks(); bn++) {
				if (!scalar_function->has_block(bn)) continue;
				auto in_b_iter = GridDecompType::in_block_iterator(bn, decomp);
				auto n_iter = GridDecompNeighborIteratorType(decomp);
				for (in_b_iter.begin(); in_b_iter.valid(); in_b_iter.advance()) {
					auto id_p = in_b_iter.value();
					if (scalar_function->has_bnidpair(id_p)) {
						//n_iter.begin(in_b_iter);
						//ncount->SetLabelBNID(id_p, n_iter.num_neg);
						char ncounter = 0;
						for (n_iter.begin(in_b_iter); n_iter.valid(); n_iter.advance()) {
							auto n_id_p = n_iter.value();
							if (scalar_function->has_bnidpair(n_id_p)) {
								// we HAVE this neighbor
								//ncount->LabelRefBNID(n_id_p)++;
								ncounter++;
							}
						}
						ncount->SetLabelBNID(id_p, ncounter);
					}
				}
			}
			printf("========> ncount time: %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - ncount_time).count());

			auto biter_time = std::chrono::steady_clock::now();
			char* debugbuffer2 = new char[X*Y*Z];
			memset(debugbuffer2, 0, X*Y*Z * sizeof(char));
			printf("Flattening Result\n");
			for (INDEX_TYPE bn = 0; bn < decomp->NumBlocks(); bn++) {
				if (!ncount->has_block(bn)) continue;
				//printf("doing block %llu, f:%d, n:%d\n", bn, scalar_function->has_block(bn), ncount->has_block(bn));
				auto in_b_iter = GridDecompType::in_block_iterator(bn, decomp);
				for (in_b_iter.begin(); in_b_iter.valid(); in_b_iter.advance()) {
					auto id_p = in_b_iter.value();
					auto gid = decomp->GetGlobalIdFromBNID(id_p);
					debugbuffer2[gid] = 1 + ncount->GetLabelBNID(id_p);
				}
			}
			printf("========> biter time: %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - biter_time).count());

			char foutname[2048];
			sprintf(foutname, "debugnegc_%d_%dx%dx%d_%dx%dx%d.raw",
				i, decomp->get_block_grid()->XYZ()[0],
				decomp->get_block_grid()->XYZ()[1],
				decomp->get_block_grid()->XYZ()[2],
				X, Y, Z);
			FILE* fout = fopen(foutname, "wb");
			fwrite(debugbuffer2, sizeof(char), X*Y*Z, fout);
			fclose(fout);
			delete[] debugbuffer2;
		}

		if (true) {
			printf("allocating structure watershed:\n");
			auto alloc_time = std::chrono::steady_clock::now();
			typedef GridDecompType::BN_ID_PAIR BNIDPair;
			SparseBlockedLabeling<BNIDPair, GridDecompType>* watershed =
				new SparseBlockedLabeling<BNIDPair, GridDecompType>(decomp);
			watershed->SetMaskValue(99999);
			watershed->initialize<SparseBlockedLabeling<float, GridDecompType>>(scalar_function);
			printf("TIMING bs %d %d %d wat_alloc %d\n", bs[0], bs[1], bs[2], std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - alloc_time).count());
			auto watershed_time = std::chrono::steady_clock::now();
			watershed_phase_1<GridDecompType>(decomp, scalar_function, watershed);
			auto watershed2_time = std::chrono::steady_clock::now();
			printf("TIMING bs %d %d %d wat_p_1 %d\n", bs[0], bs[1], bs[2], std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - watershed_time).count());
			watershed_phase_2<GridDecompType>(decomp, watershed);
			printf("TIMING bs %d %d %d wat_p_2 %d\n", bs[0], bs[1], bs[2], std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - watershed2_time).count());
			printf("TIMING bs %d %d %d wat_tot %d\n", bs[0], bs[1], bs[2], std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - watershed_time).count());
			printf("TIMING bs %d %d %d overall %d\n", bs[0], bs[1], bs[2], std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - g_start_time).count());
			printf("\n");
			if (false) {
				auto biter_time = std::chrono::steady_clock::now();
				unsigned char* debugbuffer2 = new unsigned char[X*Y*Z];
				memset(debugbuffer2, 0, X*Y*Z * sizeof(unsigned char));
				printf("Flattening Result\n");
				std::unordered_map<BNIDPair, unsigned char> id_2_col;
				for (INDEX_TYPE bn = 0; bn < decomp->NumBlocks(); bn++) {
					if (!watershed->has_block(bn)) continue;
					//printf("doing block %llu, f:%d, n:%d\n", bn, scalar_function->has_block(bn), ncount->has_block(bn));
					auto in_b_iter = GridDecompType::in_block_iterator(bn, decomp);
					for (in_b_iter.begin(); in_b_iter.valid(); in_b_iter.advance()) {
						auto id_p = in_b_iter.value();
						if (!watershed->has_bnidpair(id_p)) continue;
						auto rep_bnid = watershed->GetLabelBNID(id_p);
						unsigned char val;
						if (id_2_col.count(rep_bnid) == 0) {
							id_2_col[rep_bnid] = (decomp->GetGlobalIdFromBNID(rep_bnid) % 128) + 127;
						}
						val = id_2_col[rep_bnid];
						auto gid = decomp->GetGlobalIdFromBNID(id_p);
						debugbuffer2[gid] = val;
					}
				}
				printf("========> biter time: %dms\n", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - biter_time).count());

				char foutname[2048];
				sprintf(foutname, "debugwatershed_%d_%dx%dx%d_%dx%dx%d.raw",
					i, decomp->get_block_grid()->XYZ()[0],
					decomp->get_block_grid()->XYZ()[1],
					decomp->get_block_grid()->XYZ()[2],
					X, Y, Z);
				FILE* fout = fopen(foutname, "wb");
				fwrite(debugbuffer2, sizeof(char), X*Y*Z, fout);
				fclose(fout);
				delete[] debugbuffer2;
			}
			delete watershed;
		}
		//#pragma omp for schedule(dynamic)
		//		for (INDEX_TYPE bn = 0; bn < decomp->NumBlocks(); bn++) {
		//			if (!scalar_function->has_block(bn)) continue;
		//			auto lgrid = decomp->get_block_grid(bn);
		//			auto xyz = lgrid->XYZ();
		//			// do interior
		//			for (int z = 1; z < xyz[2] - 1; z++) {
		//				for (int y = 1; y < xyz[1] - 1; y++) {
		//					for (int x = 1; x < xyz[0] - 1; x++) {
		//						// 
		//						Vec3l res[6];
		//						lgrid->GatherSurroundingNoBoundaryCheck(Vec3l( x,y,z ), res);
		//						//
		//					}
		//				}
		//			}
		//		}


				//printf(" -- about to decompose mesh\n");
				//TopologicalSparseBlockedMesh* blocked_mesh =
				//	new TopologicalSparseBlockedMesh(decomp);
				//printf("   -- intermediate\n");
				////blocked_mesh->decompose();
				//printf(" -- done\n about to make maxvl\n");
				//SparseBlockedLabeling<unsigned char>* maxvl =
				//	new SparseBlockedLabeling<unsigned char>(blocked_mesh);
				//maxvl->SetMaskValue(255);
				//printf("    --> initializing\n");
				//maxvl->initialize<float>(scalar_function);
				//printf(" -- done\n");

				//maxvl->PrintBlockComparison<float>(scalar_function);
				//TopologicalRegularGrid3D* mesh = new TopologicalRegularGrid3D(single_grid);
				//auto mesh_xyz = mesh->XYZ();
				//RegularGrid3D* mesh_grid = new RegularGrid3D(mesh->XYZ(), { 0,0,0 });
				//GridDecompType* decomp_grid = new GridDecompType(mesh_xyz[0], mesh_xyz[1], mesh_xyz[2], block_x*2, block_y*2, block_z*2, single_grid);
				//decomp_grid->decompose();
				//SparseBlockedLabeling<float>* mesh_grid2 =
				//	new SparseBlockedLabeling<float>(single_grid, decomp);

		delete scalar_function;
		delete decomp;
		//delete maxvl;
		//delete blocked_mesh;
	} // end for loop for testing speed

	//// Explicit grid loading
	//for (int i = 0; i < 2; i++) {
	//	g_start_time = std::chrono::steady_clock::now();
	//	SparseBlockedGridFunction* data_grid = new SparseBlockedGridFunction(X, Y, Z, block_x, block_y, block_z, single_grid);
	//	data_grid->SetBackgroundValue(NumericLimits<float>::highest());
	//	printf("got here after new data_grid\n");
	//	data_grid->decompose();
	//	printf("got here after decompose\n");
	//	data_grid->ReadFromRAWSparse(filename.c_str(), maskname.c_str());
	//	printf("========> implicit %d %d %d\n", i, block_x,
	//		std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - g_start_time).count());
	//	delete data_grid;
	//}
	return 1;
#if 0
	printf("BLOCK SIZE: "); scalar_function->m_desired_block_sizes.PrintInt();
	int loaded = 0;
	INDEX_TYPE numloadedvals = 0;
	for (int i = 0; i < data_grid->NumBlocks(); i++) {
		loaded += data_grid->has_block(i);
		if (data_grid->has_block(i)) numloadedvals += data_grid->get_block(i).grid->NumElements();
	}
	printf("Num blocks loaded: %d of %d\n", loaded, data_grid->NumBlocks());
	printf(" -- block ratio: %f\n", ((float)loaded) / data_grid->NumBlocks());
	printf(" -- num floats: %llu of %llu\n", numloadedvals, data_grid->NumElements());
	printf(" -- float ratio: %f\n", ((float)numloadedvals) / data_grid->NumElements());
	INDEX_TYPE overhead = (data_grid->NumBlocks() * (1 + sizeof(float*) + sizeof(RegularGrid3DDecomposition::localExtents))); //localextents for each block + data block pointer + block_bool
	INDEX_TYPE data_bytes = numloadedvals * sizeof(float);	// size of floating point values
	INDEX_TYPE total_bytes = overhead + data_bytes;

	printf(" -- TOTAL SIZE: %llu bytes\n", total_bytes);
	printf("    -- overhead: %llu bytes\n", overhead);
	printf("    -- datasize: %llu bytes\n", data_bytes);
	printf("\n");


	printf("TIMING ACCESS SPARSE GRID......\n");

	float* data;
	data = new float[data_grid->get_global_grid()->NumElements()];
	DenseLabeling<unsigned char>* mask = new DenseLabeling<unsigned char>(single_grid->NumElements());
	mask->ReadFromFile(maskname.c_str());

	g_start_time = std::chrono::steady_clock::now();

	printf("copy values...\n");
	for (INDEX_TYPE z = 0; z < Z; z++) {
		for (INDEX_TYPE y = 0; y < Y; y++) {
			for (INDEX_TYPE x = 0; x < X; x++) {
				INDEX_TYPE bn = data_grid->block_num({ x,y,z });
				if (!data_grid->has_block(bn)) continue;
				INDEX_TYPE id = data_grid->Index3d({ x,y,z });
				if (!mask->GetLabel(id)) continue;
				data[id] = data_grid->SampleData({ x,y,z });
			}
		}
	}
	printf("========> SPARSE ALL DATA ACCESS TIMING: %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - g_start_time).count());
	g_start_time = std::chrono::steady_clock::now();
	printf("copy values...\n");
	for (INDEX_TYPE z = 0; z < Z; z++) {
		for (INDEX_TYPE y = 0; y < Y; y++) {
			for (INDEX_TYPE x = 0; x < X; x++) {
				INDEX_TYPE bn = scalar_function->block_num({ x,y,z });
				if (!scalar_function->has_block(bn)) continue;
				INDEX_TYPE id = scalar_function->Index3d({ x,y,z });
				if (!mask->GetLabel(id)) continue;
				data[id] = scalar_function->SampleData({ x,y,z });
			}
		}
	}
	printf("========> SPARSE IMPLICIT ALL DATA ACCESS TIMING: %d\n",
		std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - g_start_time).count());

	return 1;
#endif

	//	g_start_time = std::chrono::steady_clock::now();
	//	printf("single block basic reading\n");
	//	RegularGridTrilinearFunction* reg_Func = new RegularGridTrilinearFunction(single_grid);
	//	reg_Func->LoadImageFromFloatFile(filename.c_str());
	//	printf("========> FULL LOAD TIMING: %d\n",
	//		std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - g_start_time).count());
	//#ifdef USE_MASK
	//	DenseLabeling<unsigned char>* mask = new DenseLabeling<unsigned char>(single_grid->NumElements());
	//	mask->ReadFromFile(maskname.c_str());
		//printf("validating...\n");
		//for (INDEX_TYPE z = 0; z < Z; z++) {
		//	for (INDEX_TYPE y = 0; y < Y; y++) {
		//		for (INDEX_TYPE x = 0; x < X; x++) {
		//			INDEX_TYPE bn = data_grid->block_num({ x,y,z });
		//			if (! data_grid->has_block(bn)) continue;
		//			if (!mask->GetLabel(data_grid->Index3d({ x,y,z }))) continue;
		//			if (reg_Func->SampleImage({ x, y, z }) != data_grid->SampleData({ x,y,z })) {
		//				printf("images dont match %f, %f at (%d, %d, %d)\n", reg_Func->SampleImage({ x, y, z }),
		//					data_grid->SampleData({ x,y,z }), x, y, z);
		//			}
		//		}
		//	}
		//}
	//	printf("done\n");
	//
	//
	//#endif
#if 0
	printf("==============BLOCKED VALUES===================\n");
	for (INDEX_TYPE z = 0; z < Z; z++) {
		printf("z=%d\n", z);
		for (INDEX_TYPE y = 0; y < Y; y++) {
			for (INDEX_TYPE x = 0; x < X; x++) {
				printf("%.4f ", data_grid->SampleData({ x,y,z }));
			}
			printf("\n");
		}
	}

	printf("==================RAW VALUES===================\n");
	for (INDEX_TYPE z = 0; z < Z; z++) {
		printf("z=%d\n", z);
		for (INDEX_TYPE y = 0; y < Y; y++) {
			for (INDEX_TYPE x = 0; x < X; x++) {
				printf("%.4f ", reg_Func->SampleImage({ x,y,z }));
			}
			printf("\n");
		}
	}
	printf("==============BLOCKED NUMS=====================\n");
	for (INDEX_TYPE z = 0; z < Z; z++) {
		printf("z=%d\n", z);
		for (INDEX_TYPE y = 0; y < Y; y++) {
			for (INDEX_TYPE x = 0; x < X; x++) {
				printf("%d ", data_grid->block_num({ x,y,z }));
			}
			printf("\n");
		}
	}
	//printf("==============BLOCK DIMS=======================\n");
	//for (int i = 0; i < data_grid->NumBlocks(); i++) {
	//	printf("block %d:\n", i);
	//	data_grid->m_blocks[i].starts.PrintInt();
	//	data_grid->m_blocks[i].sizes.PrintInt();
	//}
#endif

#ifndef USE_MASK
	printf("validating...\n");
	for (INDEX_TYPE z = 0; z < Z; z++) {
		for (INDEX_TYPE y = 0; y < Y; y++) {
			for (INDEX_TYPE x = 0; x < X; x++) {
				if (reg_Func->SampleImage({ x, y, z }) != data_grid->SampleData({ x,y,z })) {
					printf("images dont match %f, %f at (%d, %d, %d)\n", reg_Func->SampleImage({ x, y, z }),
						data_grid->SampleData({ x,y,z }), x, y, z);
				}
			}
		}
	}
	printf("done\n");
#endif
	return 1;


}
#endif

