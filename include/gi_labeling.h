#ifndef VERTEX_LABELING_H
#define VERTEX_LABELING_H

#include <unordered_map>
#include "gi_basic_types.h"
#include "gi_regular_grid.h"

namespace GInt {
    template<typename LABEL_TYPE>
    class DenseLabeling {
    protected:
        LABEL_TYPE* m_labels;
        INDEX_TYPE m_num_labels;
    public:

        DenseLabeling(INDEX_TYPE num_labels) : m_num_labels(num_labels) {
            m_labels = new LABEL_TYPE[num_labels];
            //printf(" allocated space for %d values of size %d\n", num_labels, sizeof(LABEL_TYPE));
        }
        ~DenseLabeling() {
            delete[] m_labels;
        }

        void SetLabel(INDEX_TYPE id, LABEL_TYPE label) {
            m_labels[id] = label;
        }

        LABEL_TYPE GetLabel(INDEX_TYPE id) const {
            return m_labels[id];
        }

        INDEX_TYPE GetNumLabels() const {
            return m_num_labels;
        }

        LABEL_TYPE& operator[](const INDEX_TYPE id) { return m_labels[id]; }
        const LABEL_TYPE& operator[](const INDEX_TYPE id) const { return m_labels[id]; }

        void SetAll(LABEL_TYPE label){
#pragma omp parallel for schedule(static)
            for (int i = 0; i < m_num_labels; i++) {
                m_labels[i] = label;
            }
        }

		void CopyValues(const DenseLabeling<LABEL_TYPE>* other) {
#pragma omp parallel for schedule(static)
			for (int i = 0; i < m_num_labels; i++) {
				m_labels[i] = other->m_labels[i];
			}
		}
      template<typename T>
      void ReMapIds(T* output){
        std::unordered_map<INDEX_TYPE, T> unique_ids;

        T new_id=0;

        for(INDEX_TYPE i=0; i < m_num_labels; i++){
          if(unique_ids.find(m_labels[i]) == unique_ids.end()){
            T set_id = new_id;

            if(m_labels[i] < 0)
              set_id = -1;
            else
              set_id = new_id++;

            unique_ids[m_labels[i]] = set_id;
            printf("mapping id %lld to %d\n", m_labels[i], set_id);
          }

          output[i] = unique_ids.at(m_labels[i]);
        }

        printf("remapped %d ids\n", new_id);
      }


        void ReadFromFile(const char* filename) {
            FILE* fout = fopen(filename, "rb");
            fread(m_labels, sizeof(LABEL_TYPE), m_num_labels, fout);
            fclose(fout);
        }

        void OutputToFile(const char* filename) const {
			printf("writing file %s \n", filename);
            FILE* fout = fopen(filename, "wb");
            //printf("Sizeof label type: %d\n", sizeof(LABEL_TYPE));
            fwrite(m_labels, sizeof(LABEL_TYPE), m_num_labels, fout);
            fclose(fout);
        }
        void OutputToIntFile(const char* filename) const {
			printf("writing file %s \n", filename);
			FILE* fout = fopen(filename, "wb");
            //printf("Sizeof int type: %d\n", sizeof(int));
            for (INDEX_TYPE i = 0; i < m_num_labels; i++) {
                int tval = (int)m_labels[i];
                fwrite(&tval, sizeof(int), 1, fout);
            }
            fclose(fout);
        }
        void OutputToFloatFile(const char* filename) const {
			printf("writing file %s \n", filename);
			FILE* fout = fopen(filename, "wb");
            //printf("Sizeof float type: %d\n", sizeof(float));
            for (INDEX_TYPE i = 0; i < m_num_labels; i++) {
                float tval = (float)m_labels[i];
                fwrite(&tval, sizeof(float), 1, fout);
            }
            fclose(fout);
        }

        LABEL_TYPE* LabelArray() { return m_labels; }
    };



	template<typename dtype, class DecompType>
	class SparseBlockedLabeling {
	protected:
		dtype** m_data_blocks;
		const DecompType* m_sparse_grid;
	public:

		SparseBlockedLabeling(DecompType* grid) :
			m_sparse_grid(grid) {
		}
		virtual ~SparseBlockedLabeling() {
			for (INDEX_TYPE i = 0; i < m_sparse_grid->NumBlocks(); i++) {
				if (this->has_block(i)) delete[] m_data_blocks[i];
			}
			delete[] m_data_blocks;
		}
		virtual void initialize() {
			m_data_blocks = new dtype*[m_sparse_grid->NumBlocks()](); // should value-initialize to zero
		}
		const DecompType* GetGrid() const {
			return m_sparse_grid;
		}
		// this version allocates the same blocks as other sparse label - but with the
		// own underlying sparse grid. e.g. used to create same mesh label structure as grid but with larger blocks
		// values are initialized to zero?
		template<class OtherSparseBlockedLabelingType>
		void initialize(OtherSparseBlockedLabelingType* other) {
			if (m_sparse_grid->NumBlocks() != other->GetGrid()->NumBlocks()) {
				printf("ERROR: SparseBlockedLabeling other has different number of blocks\n");
				return;
			}
			INDEX_TYPE num_blocks = m_sparse_grid->NumBlocks();
			blocks_to_include = new unsigned char[num_blocks]();
			m_data_blocks = new dtype*[num_blocks](); // should value-initialize to zero
#pragma omp for schedule(dynamic)
			for (INDEX_TYPE bn = 0; bn < num_blocks; bn++) {
				if (other->has_block(bn)) {
					blocks_to_include[bn] = 1;
					INDEX_TYPE block_size = m_sparse_grid->get_block_grid(bn)->NumElements();
					m_data_blocks[bn] = new dtype[block_size]();
					for (INDEX_TYPE i = 0; i < block_size; i++) m_data_blocks[bn][i] = m_mask_value;
				}
			}
		}

		template<class OtherSparseBlockedLabelingType>
		void PrintBlockComparison(OtherSparseBlockedLabelingType* other) const {
			printf("Comparing blocks between labelings:\n");
			bool hasdif = false;
			for (int i = 0; i < this->m_sparse_grid->NumBlocks(); i++) {
				if (this->has_block(i) != other->has_block(i)) {
					printf("%d: this=%d != other=%d\n", i, this->has_block(i), other->has_block(i));
					hasdif = true;
				}
			}
			if (!hasdif) {
				printf(" -- same blocks present in both\n");
			}
		}

		inline dtype GetLabelBNID(typename DecompType::BN_ID_PAIR bnid) const {
			INDEX_TYPE id, bn;
			m_sparse_grid->GetBNAndIDFromPair(bnid, bn, id);
			return m_data_blocks[bn][id];
		}
		inline dtype GetLabel(Vec3l point) const {
			auto ijk = m_sparse_grid->block_ijk_from_xyz(point);
			Vec3l starts;
			m_sparse_grid->block_start(ijk, starts);
			Vec3l ind = point - starts;
			return m_data_blocks[m_sparse_grid->block_id_from_block_pos(ijk)][m_sparse_grid->get_block_grid(ijk)->Index3d(ind)];
		}
		inline dtype GetLabel(INDEX_TYPE id) const {
			Vec3l point = m_sparse_grid->Coords(id);
			return GetLabel(point);
		}

		inline void SetLabel(Vec3l point, dtype label) {
			auto ijk = m_sparse_grid->block_ijk_from_xyz(point);
			Vec3l starts;
			m_sparse_grid->block_start(ijk, starts);
			Vec3l ind = point - starts;
			m_data_blocks[m_sparse_grid->block_id_from_block_pos(ijk)][m_sparse_grid->get_block_grid(ijk)->Index3d(ind)] = label;
		}
		inline void SetLabel(INDEX_TYPE id, dtype label) {
			Vec3l point = m_sparse_grid->Coords(id);
			SetLabel(point, label);
		}
		inline void SetLabelBNID(typename DecompType::BN_ID_PAIR bnid, dtype label) {
			INDEX_TYPE id, bn;
			m_sparse_grid->GetBNAndIDFromPair(bnid, bn, id);
			m_data_blocks[bn][id] = label;
		}
		inline dtype& LabelRefBNID(typename DecompType::BN_ID_PAIR bnid) {
			INDEX_TYPE id, bn;
			m_sparse_grid->GetBNAndIDFromPair(bnid, bn, id);
			return m_data_blocks[bn][id];
		}
		bool ReadFromRAWDense(const char* fname) {

			// load values from file
			INDEX_TYPE slab_size = m_sparse_grid->XYZ()[0] * m_sparse_grid->XYZ()[1];
			INDEX_TYPE slab_width = 5;
#pragma omp parallel
			{
				slab_width = omp_get_num_threads() * 1;
			}
			printf("slab width = %d\n", slab_width);
			printf("slab size = (%llu x %llu)\n", m_sparse_grid->XYZ()[0], m_sparse_grid->XYZ()[1]);
			dtype* slab_data = new dtype[slab_size * slab_width];


			FILE* fin = fopen(fname, "rb");
			if (fin == NULL) {
				printf("Error reading %s\n", fname);
				return false;
			}
			blocks_to_include = new unsigned char[m_sparse_grid->NumBlocks()]; // initialized to 0
			printf("allocating data blocks\n");
			// allocate data blocks
			for (INDEX_TYPE i = 0; i < m_sparse_grid->NumBlocks(); i++) {
				blocks_to_include[i] = 1; // include every block
				//printf("  -- %llu elements\n", m_blocks[i].grid->NumElements());
				Vec3l ijk = m_sparse_grid->block_ijk_from_id(i);
				m_data_blocks[i] = new dtype[m_sparse_grid->get_block_grid(ijk)->NumElements()];
			}
			printf("done allocating\n");

			// outer file reading loop
			INDEX_TYPE slab_current_z = 0;
			//INDEX_TYPE block_z_size = m_sparse_grid->m_basic_block_dims[2];

			while (slab_current_z < m_sparse_grid->XYZ()[2]) {
				//#define DEBUG_FIRST_SLICE
#ifdef DEBUG_FIRST_SLICE
				if (slab_current_z > 0) break;
#endif
				// read next slab
				if (slab_current_z + slab_width > m_sparse_grid->XYZ()[2]) {
					slab_width = m_sparse_grid->XYZ()[2] - slab_current_z;
				}
				INDEX_TYPE num_to_read = slab_size * slab_width;
				fread(slab_data, sizeof(dtype), num_to_read, fin);
				//printf(" -- read %llu elements\n", num_to_read);
				// now have read slab_width *X*Y values

				// copy into blocks
#pragma omp parallel for schedule(dynamic)
				for (int z = 0; z < slab_width; z++) {
					INDEX_TYPE t_current_z = slab_current_z + z;
#ifdef DEBUG_FIRST_SLICE
					if (z > 0) {
						slab_current_z++;
						break;
					}
#endif
					//printf("   -- doing z= %d, slab_current_z = %d\n", z, slab_current_z);
					// foreach [x, x + blocksize], y stick, copy into block
					for (int y = 0; y < m_sparse_grid->XYZ()[1]; y++) {
						//printf("   -- doing y= %d\n", y);
						// get block number and just add 1 after inner loop
						Vec3l ijk = m_sparse_grid->block_ijk_from_xyz({ 0, y, t_current_z });
						//printf("   -- bn = %d\n", bn);
						// inner loop with memory copy
						for (int x = 0; x < m_sparse_grid->XYZ()[0]; x += m_sparse_grid->m_basic_block_dims[0]) {
							typename DecompType::localExtents le;
							m_sparse_grid->get_block(ijk, le);
							INDEX_TYPE bn = m_sparse_grid->block_id_from_block_pos(ijk);
							INDEX_TYPE in_block_y = y - le.starts[1];
							INDEX_TYPE in_block_z = t_current_z - le.starts[2];
							INDEX_TYPE copy_size = m_sparse_grid->m_basic_block_dims[0];
							if (x + m_sparse_grid->m_basic_block_dims[0] > m_sparse_grid->XYZ()[0])
								copy_size = m_sparse_grid->XYZ()[0] - x;
							//printf("   -- x = [%d, %d], in_block=(0, %d, %d)\n", x, x + copy_size, in_block_y, in_block_z);
							//m_blocks[bn].starts.PrintInt();
							// address in slab is z sliceinslab + y* slab_x + x
							dtype* source = &(slab_data[x + m_sparse_grid->XYZ()[0] * y + slab_size*z]);

							// address to copy to gets block at 
							dtype* dest = &(m_data_blocks[bn][le.grid->Index3d({ 0, in_block_y, in_block_z })]);

							//COPY
							memcpy(dest, source, sizeof(dtype) * copy_size);

							//for (int t = 0; t < copy_size; t++) dest[t] = z;
							ijk[0]++; // increment block number
						} // end x loop
					} // end y loop
				} // end z loop
				slab_current_z += slab_width;
			} //end loop over slabs
			fclose(fin);
			return true;
		}

		unsigned char* blocks_to_include;
		dtype m_mask_value;

		// Extreme caution needs to be used with this. THe intent is to 
		// avoid holding a separate "mask" field - rather bake into 
		// the label value a special "background" value. e.g. with floats
		// you could use NumericLimits<float>::highest(). 
		void SetMaskValue(dtype val) { m_mask_value = val; }

		bool has_block(INDEX_TYPE block_num) const {
			return blocks_to_include[block_num];
		}
		bool has_xyz(const Vec3l& xyz) const {
			INDEX_TYPE block_num = m_sparse_grid->block_num(xyz);
			return blocks_to_include[block_num];
		}
		bool has_bnidpair(typename DecompType::BN_ID_PAIR bnid) const {
			INDEX_TYPE block_num, in_block_id;
			m_sparse_grid->GetBNAndIDFromPair(bnid, block_num, in_block_id);
			return has_bnid(block_num, in_block_id);
		}
		// block_num, local_id are direct lookup
		bool has_bnid(INDEX_TYPE block_num, INDEX_TYPE local_id) const {
			return has_block(block_num) && this->m_data_blocks[block_num][local_id] != m_mask_value;
		}
		// ijk block position, local_id is id inside block
		bool has_bpos_id(const Vec3l& ijk, INDEX_TYPE local_id) const {
			return has_bnid(m_sparse_grid->block_id_from_block_pos(ijk), local_id);
		}
	

		bool ReadFromRAWSparse(const char* fname, const char* maskname) {

			// load values from file
			INDEX_TYPE slab_size = m_sparse_grid->XYZ()[0] * m_sparse_grid->XYZ()[1];
			INDEX_TYPE slab_width = 5;
#pragma omp parallel
			{
				slab_width = omp_get_num_threads() * 1;
			}
			printf("slab width = %d\n", slab_width);
			printf("slab size = (%llu x %llu)\n", m_sparse_grid->XYZ()[0], m_sparse_grid->XYZ()[1]);
			dtype* slab_data = new dtype[slab_size * slab_width];
			unsigned char* slab_mask = new unsigned char[slab_size * slab_width];


			FILE* fin = fopen(fname, "rb");
			if (fin == NULL) {
				printf("Error reading %s\n", fname);
				return false;
			}

			FILE* fin_mask = fopen(maskname, "rb");
			if (fin == NULL) {
				printf("Error reading %s\n", maskname);
				fclose(fin);
				return false;
			}

			//printf("allocating data blocks\n");
			//// allocate data blocks
			//for (INDEX_TYPE i = 0; i < NumBlocks(); i++) {
			//	//printf("  -- %llu elements\n", m_blocks[i].grid->NumElements());
			//	m_data_blocks[i] = new dtype[m_blocks[i].grid->NumElements()];
			//}
			//printf("done allocating\n");

			// outer file reading loop
			INDEX_TYPE slab_current_z = 0;
			//INDEX_TYPE block_z_size = m_sparse_grid->m_basic_block_dims[2];

			blocks_to_include = new unsigned char[m_sparse_grid->NumBlocks()](); // initialized to 0

			while (slab_current_z < m_sparse_grid->XYZ()[2]) {
				//#define DEBUG_FIRST_SLICE
#ifdef DEBUG_FIRST_SLICE
				if (slab_current_z > 0) break;
#endif
				// read next slab
				if (slab_current_z + slab_width > m_sparse_grid->XYZ()[2]) {
					slab_width = m_sparse_grid->XYZ()[2] - slab_current_z;
				}
				INDEX_TYPE num_to_read = slab_size * slab_width;
				fread(slab_data, sizeof(dtype), num_to_read, fin);
				fread(slab_mask, sizeof(unsigned char), num_to_read, fin_mask);
				// now have read slab_width *X*Y values
				//printf("read slabs %d\n", slab_current_z);
				// now make sure all NEEDED blocks are allocated
				// - for each z in slab, check if data blocks intersecting z have non-zero in mask slab
				// - will allocate needed data blocks
				// - can check either non-zero data block pointer or 1 on blocks_to_include -- redundant? 
				for (int z = 0; z < slab_width; z++) {
					INDEX_TYPE t_current_z = slab_current_z + z;

					INDEX_TYPE k = m_sparse_grid->block_k_from_z(t_current_z);
#pragma omp for schedule(dynamic)
					for (auto j = 0; j < m_sparse_grid->blocks_per_axis()[1]; j++) {
						for (auto i = 0; i < m_sparse_grid->blocks_per_axis()[0]; i++) {
							Vec3l ijk(i, j, k);
							INDEX_TYPE t_block_num = m_sparse_grid->block_id_from_block_pos(ijk);
							// check if ijk block exists - if so, skip
							if (blocks_to_include[t_block_num] != 0) continue;

							typename DecompType::localExtents le;
							m_sparse_grid->get_block(ijk, le);

							for (INDEX_TYPE y = 0; y < le.sizes[1]; y++) {
								auto gy = y + le.starts[1];
								for (INDEX_TYPE x = 0; x < le.sizes[0]; x++) {
									auto gx = x + le.starts[0];
									auto val = slab_mask[gx + m_sparse_grid->XYZ()[0] * gy + slab_size*z];
									if (val != 0) {
										blocks_to_include[t_block_num] = 1; // mark the block
										m_data_blocks[t_block_num] = new dtype[le.grid->NumElements()]; // allocate the space
										auto lim = le.grid->NumElements();
										auto* ptr = m_data_blocks[t_block_num];
										for (int i = 0; i < lim; i++) ptr[i] = m_mask_value;
										y = le.sizes[1]; // skip out of y loop
										break; // skip out of x loop
									}
								}
							}


						}
					} // end parallel over j blocks
				} // end iteration over z slices

				  //printf("read masks and allocatd\n");
				  // copy into blocks
#pragma omp parallel for schedule(dynamic)
				for (int z = 0; z < slab_width; z++) {
					INDEX_TYPE t_current_z = slab_current_z + z;

#ifdef DEBUG_FIRST_SLICE
					if (z > 0) {
						slab_current_z++;
						break;
					}
#endif
					//printf("   -- doing z= %d, slab_current_z = %d\n", z, slab_current_z);
					// foreach [x, x + blocksize], y stick, copy into block
					for (int y = 0; y < m_sparse_grid->XYZ()[1]; y++) {
						//printf("   -- doing y= %d\n", y);
						// get block number and just add 1 after inner loop
						Vec3l ijk = m_sparse_grid->block_ijk_from_xyz({ 0, y, t_current_z });
						INDEX_TYPE bn = m_sparse_grid->block_id_from_block_pos(ijk);
						//printf("   -- bn = %d\n", bn);
						// inner loop with memory copy
						for (int x = 0; x < m_sparse_grid->XYZ()[0]; x += m_sparse_grid->get_block_grid()->XYZ()[0]) {
							if (blocks_to_include[bn] == 0) {
								bn++;
								ijk[0]++;
								continue;
							}
							typename DecompType::localExtents le;
							m_sparse_grid->get_block(ijk, le);

							INDEX_TYPE in_block_y = y - le.starts[1];
							INDEX_TYPE in_block_z = t_current_z - le.starts[2];
							INDEX_TYPE copy_size = m_sparse_grid->get_block_grid()->XYZ()[0];
							if (x + m_sparse_grid->get_block_grid()->XYZ()[0] > m_sparse_grid->XYZ()[0])
								copy_size = m_sparse_grid->XYZ()[0] - x;
							//printf("   -- x = [%d, %d], in_block=(0, %d, %d)\n", x, x + copy_size, in_block_y, in_block_z);
							//m_blocks[bn].starts.PrintInt();
							// address in slab is z sliceinslab + y* slab_x + x
							dtype* source = &(slab_data[x + m_sparse_grid->XYZ()[0] * y + slab_size*z]);
							unsigned char* mask = &(slab_mask[x + m_sparse_grid->XYZ()[0] * y + slab_size*z]);

							// address to copy to gets block at 
							dtype* dest = &(m_data_blocks[bn][le.grid->Index3d({ 0, in_block_y, in_block_z })]);
							//COPY
							//memcpy(dest, source, sizeof(dtype) * copy_size);

							for (int t = 0; t < copy_size; t++) dest[t] = (mask[t] ? source[t] : m_mask_value); // masked values
							bn++; // increment block number
							ijk[0]++;
						} // end x loop
					} // end y loop
				} // end z loop
				slab_current_z += slab_width;
			} //end loop over slabs
			fclose(fin);
			return true;
		}

	};




	template<typename LABEL_TYPE>
	class SparseActuallyDenseLabeling : public DenseLabeling<LABEL_TYPE> {
	public: 
		bool Has(INDEX_TYPE id) const {
			return true;
		}
	};

	template<typename LABEL_TYPE>
	class SparseLabeling {
	protected:
		std::unordered_map<INDEX_TYPE, LABEL_TYPE> m_labels;
		
	public:
		typedef typename std::unordered_map<INDEX_TYPE, LABEL_TYPE>::iterator SparseIterator;
		SparseIterator begin() { return m_labels.begin(); }
		SparseIterator end() { return m_labels.end(); }

		class SparseKeyIterator : public SparseIterator {
		public:
			SparseKeyIterator() : SparseIterator() {};
			SparseKeyIterator(SparseIterator s) : SparseIterator(s) {};
			INDEX_TYPE* operator->() { return (INDEX_TYPE* const)&(SparseIterator::operator->()->first); }
			INDEX_TYPE operator*() { return SparseIterator::operator*().first; }
		};

		SparseLabeling()  {
		}
		~SparseLabeling() {
		}
		bool Has(INDEX_TYPE id) const {
			return m_labels.count(id) == 1;
		}

		void SetLabel(INDEX_TYPE id, LABEL_TYPE label) {
			m_labels[id] = label;
		}

		LABEL_TYPE GetLabel(INDEX_TYPE id) const {
			return m_labels.find(id)->second;
		}

		INDEX_TYPE GetNumLabels() const {
			return m_labels.size();
		}

		LABEL_TYPE& operator[](const INDEX_TYPE id) { return m_labels[id]; }
		const LABEL_TYPE& operator[](const INDEX_TYPE id) const { 
			return m_labels.find(id)->second;
		}

		void SetAll(LABEL_TYPE label) {
			for (auto& val_pair : m_labels) {
				val_pair.second = label;
			}
		}

	};



	template<typename LABEL_TYPE>
	class IndirectLabeling {
	protected:
		DenseLabeling<LABEL_TYPE>* m_direct_labeling;
		SparseLabeling<int>* m_remap;
	public:
		IndirectLabeling(INDEX_TYPE count) { printf("ERROR: DONT do this\n"); }
		IndirectLabeling(SparseLabeling<int>* sparse_labeling) : m_remap(sparse_labeling) {
			m_direct_labeling = new DenseLabeling<LABEL_TYPE>(m_remap->GetNumLabels());
		}
		~IndirectLabeling() {
			delete m_direct_labeling;
		}

		void SetLabel(INDEX_TYPE id, LABEL_TYPE label) {
			m_direct_labeling->SetLabel(m_remap->GetLabel(id), label);
		}

		LABEL_TYPE GetLabel(INDEX_TYPE id) const {
			return m_direct_labeling->GetLabel(m_remap->GetLabel(id));

		}

		INDEX_TYPE GetNumLabels() const {
			return m_remap->GetNumLabels();
		}

		LABEL_TYPE& operator[](const INDEX_TYPE id) { return m_direct_labeling->operator[](m_remap->GetLabel(id)); }
		const LABEL_TYPE& operator[](const INDEX_TYPE id) const { return m_direct_labeling->operator[](m_remap->GetLabel(id)); }

		void SetAll(LABEL_TYPE label) {
			m_direct_labeling->SetAll(label);
		}

		void CopyValues(const DenseLabeling<LABEL_TYPE>* other) {
			m_direct_labeling->CopyValues(other);
		}



		//void ReadFromFile(const char* filename) {
		//	FILE* fout = fopen(filename, "rb");
		//	fread(m_labels, sizeof(LABEL_TYPE), m_num_labels, fout);
		//	fclose(fout);
		//}

		//void OutputToFile(const char* filename) const {
		//	printf("writing file %s \n", filename);
		//	FILE* fout = fopen(filename, "wb");
		//	//printf("Sizeof label type: %d\n", sizeof(LABEL_TYPE));
		//	fwrite(m_labels, sizeof(LABEL_TYPE), m_num_labels, fout);
		//	fclose(fout);
		//}
		//void OutputToIntFile(const char* filename) const {
		//	printf("writing file %s \n", filename);
		//	FILE* fout = fopen(filename, "wb");
		//	//printf("Sizeof int type: %d\n", sizeof(int));
		//	for (INDEX_TYPE i = 0; i < m_num_labels; i++) {
		//		int tval = (int)m_labels[i];
		//		fwrite(&tval, sizeof(int), 1, fout);
		//	}
		//	fclose(fout);
		//}
		//void OutputToFloatFile(const char* filename) const {
		//	printf("writing file %s \n", filename);
		//	FILE* fout = fopen(filename, "wb");
		//	//printf("Sizeof float type: %d\n", sizeof(float));
		//	for (INDEX_TYPE i = 0; i < m_num_labels; i++) {
		//		float tval = (float)m_labels[i];
		//		fwrite(&tval, sizeof(float), 1, fout);
		//	}
		//	fclose(fout);
		//}

		LABEL_TYPE* LabelArray() { return m_direct_labeling->LabelArray(); }
	};

}

#endif
