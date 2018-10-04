#ifndef VERTEX_LABELING_H
#define VERTEX_LABELING_H

#include "gi_basic_types.h"

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



}

#endif
