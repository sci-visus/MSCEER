/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef DISCRETE_GRADIENT_LABELING_H
#define DISCRETE_GRADIENT_LABELING_H

#include "gi_basic_types.h"
//#include "gi_regular_grid.h"
#include "gi_labeling.h"
//#include "gi_topological_regular_grid.h"

namespace GInt {
	template <class TopoMeshType>
	class DiscreteGradientLabeling {
	public:
		struct GradBitfield {
			unsigned char assigned : 1;
			unsigned char flag : 1;
			//unsigned char critical : 1;
			//unsigned char insorter : 1;
			//unsigned char dimA : 3;
			unsigned char pair : 3;
			unsigned char ldir : 3;
		};
	protected:
		TopoMeshType* m_mesh;
		

		DenseLabeling<GradBitfield>* m_dgrad;

	public:
		DiscreteGradientLabeling(TopoMeshType* mesh) : m_mesh(mesh) {
			m_dgrad = new DenseLabeling<GradBitfield>(m_mesh->numCells());
		}
		////GradBitfield& GetBitField(INDEX_TYPE id) {
		////	return m_dgrad->GetLabel[id];
		////}
		void copyValues(const DiscreteGradientLabeling<TopoMeshType>* other) {
			if (m_dgrad->GetNumLabels() != other->m_dgrad->GetNumLabels()) {
				printf("Error: copy needs same sizes\n");
			}
			m_dgrad->CopyValues(other->m_dgrad);
		}
		
		GradBitfield& getBitfield(INDEX_TYPE cellid) { return (*m_dgrad)[cellid]; }
		char getAsChar(INDEX_TYPE cellid) { 
			char* val;
			GradBitfield vv = (*m_dgrad)[cellid];
			val = (char*) &vv;
			return *val;
		}
		void clearGrad(INDEX_TYPE cellid) {
			(*m_dgrad)[cellid] = GradBitfield{ 0 };
		}
		// this is a super dangerous call - use at your own risk!
		// basically sometimes it is useful to change the mesh while 
		// keeping the same underlying data structure - e.g. when you 
		// want to switch to the logic of a complement mesh, but keep
		// the computation you did on the original mesh
		void resetMeshHandler(TopoMeshType* mesh) {
			m_mesh = mesh;
		}
		void ClearAllGradient() {
			memset(m_dgrad->LabelArray(), 0, sizeof(GradBitfield)*m_mesh->numCells());
		}
		// ALWAYS VALID 
        ASSIGNED_TYPE getAssigned(INDEX_TYPE cellid) const {
			return (*m_dgrad)[cellid].assigned;
		}

		void setAssigned(INDEX_TYPE cellid, ASSIGNED_TYPE value) {
			(*m_dgrad)[cellid].assigned = value;
		}

		// VALID ONLY AFTER FIRST ASSIGNMENT
        INDEX_TYPE getPair(INDEX_TYPE cellid) const{
			return m_mesh->UncompressByteTo6NeighborOffset(cellid, (*m_dgrad)[cellid].pair);
		}

		void setPair(INDEX_TYPE cellid, INDEX_TYPE value) {
			(*m_dgrad)[cellid].pair = m_mesh->Compress6NeighborOffsetToByte(cellid, value);
		}

        bool getCritical(INDEX_TYPE cellid) const {
			return (*m_dgrad)[cellid].pair == 7;
		}
		void setCritical(INDEX_TYPE cellid, bool value) {
			if (value) (*m_dgrad)[cellid].pair = 7;
		}

		DIM_TYPE getDimAscMan(INDEX_TYPE cellid){
			return (*m_dgrad)[cellid].ldir;
		}
		void setDimAscMan(INDEX_TYPE cellid, DIM_TYPE value){
			(*m_dgrad)[cellid].ldir = value;
		}

		// VALID ONLY BEFORE FIRST ASSIGNEMT
        INDEX_TYPE getNumUnpairedFacets(INDEX_TYPE cellid) {
            return (*m_dgrad)[cellid].ldir;
        }
        void setNumUnpairedFacets(INDEX_TYPE cellid, INDEX_TYPE value) {
            (*m_dgrad)[cellid].ldir = value;
        }

        unsigned char getNumUnpairedFacets_new(INDEX_TYPE cellid) {
            return (*m_dgrad)[cellid].ldir;
        }
        unsigned char decrementNumUnpairedFacets_new(INDEX_TYPE cellid) {
            return --((*m_dgrad)[cellid].ldir);
        }
        void setNumUnpairedFacets_new(INDEX_TYPE cellid, unsigned char value) {
            (*m_dgrad)[cellid].ldir = value;
        }

		DIM_TYPE getMark(INDEX_TYPE cellid){
			return (*m_dgrad)[cellid].pair;
		}
		void setMark(INDEX_TYPE cellid, DIM_TYPE value){
			(*m_dgrad)[cellid].pair = value;
		}



		void outputToRenderer(const char* filenamebase) {

			char grad_name[1024];
			sprintf(grad_name, "%s.grad", filenamebase);

			FILE* fgrad = fopen(grad_name, "wb");

			for (int i = 0; i < m_mesh->numCells(); i++) {

				(*m_dgrad)[i].flag = (bool)m_mesh->boundaryValue(i);
				//if (i % 1000 == 0)
				//	printf("cell[%d]:(%d, %d, %d, %d)\n", i,
				//	m_dgrad[i].pair,
				//	m_dgrad[i].ldir,
				//	m_dgrad[i].flag,
				//	m_dgrad[i].assigned);


				fwrite(&(*m_dgrad)[i], sizeof(GradBitfield), 1, fgrad);
			}
			fclose(fgrad);
		}

		void outputToFile(const char* filenamebase) {


			FILE* fgrad = fopen(filenamebase, "wb");
			printf("Going to output %lu cells in %s...\n", m_mesh->numCells(), filenamebase);
			for (size_t i = 0; i < m_mesh->numCells(); i++) {

                (*m_dgrad)[i].flag = (bool)m_mesh->boundaryValue(i);

                //(*m_dgrad)[i].flag = 0;
                //(*m_dgrad)[i].assigned = 0;
                //(*m_dgrad)[i].pair = 0;
                //(*m_dgrad)[i].ldir = 0;
		                fwrite(&((*m_dgrad)[i]), sizeof(GradBitfield), 1, fgrad);
			}

            printf("   -- Output to %s.. %lu values of size %d\n", filenamebase, m_mesh->numCells (), sizeof(GradBitfield));
			fclose(fgrad);
		}


		bool load_from_file(const char* filename) {

			FILE* fdat = fopen(filename, "rb");
			if (fdat == NULL) {
				return false;

			}

			fread(m_dgrad->LabelArray(), sizeof(GradBitfield), m_mesh->numCells(), fdat);

			//for (int i = 0; i < m_mesh->numCells(); i++) {
			//	
			//	bitfield val;
			//	fread(&val, sizeof(bitfield), 1, fdat);
			//	m_dgrad[i]=val;

			//}

			fclose(fdat);
			return true;
		}

	};

}
#endif
