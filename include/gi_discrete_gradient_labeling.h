#ifndef DISCRETE_GRADIENT_LABELING_H
#define DISCRETE_GRADIENT_LABELING_H

#include "gi_basic_types.h"
#include "gi_regular_grid.h"
#include "gi_labeling.h"

namespace GInt {

	struct GradBitfield {
		unsigned char assigned : 1;
		unsigned char flag : 1;
		//unsigned char critical : 1;
		//unsigned char insorter : 1;
		//unsigned char dimA : 3;
		unsigned char pair : 3;
		unsigned char ldir : 3;
	};

	template <class TopoMeshType, class LabelingType = DenseLabeling<GradBitfield> >
	class DiscreteGradientLabeling {
	public:
		TopoMeshType* m_mesh;
		LabelingType* m_dgrad;
	
		DiscreteGradientLabeling(TopoMeshType* mesh) : m_mesh(mesh) {
			m_dgrad = new LabelingType(m_mesh->numCells());
		}
		~DiscreteGradientLabeling() {
			delete m_dgrad;
		}
		DiscreteGradientLabeling(TopoMeshType* mesh, LabelingType* dgrad) : m_mesh(mesh) {
			m_dgrad = dgrad;
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

		//GradBitfield& getBitfield(INDEX_TYPE cellid) { return (*m_dgrad)[cellid]; }
		char getAsChar(INDEX_TYPE cellid) {
			char* val;
			GradBitfield vv = m_dgrad->GetLabel(cellid);//(*m_dgrad)[cellid];
			val = (char*)&vv;
			return *val;
		}
		void clearGrad(INDEX_TYPE cellid) {
			m_dgrad->SetLabel(cellid, GradBitfield{ 0 });              // store it back atomically

			//(*m_dgrad)[cellid] = GradBitfield{ 0 };
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
			//m_dgrad->SetAll(0);
			memset(m_dgrad->LabelArray(), 0, sizeof(GradBitfield)*m_mesh->numCells());
		}
		// ALWAYS VALID 
		ASSIGNED_TYPE getAssigned(INDEX_TYPE cellid) const {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.assigned;
			//return (*m_dgrad)[cellid].assigned;
		}

		void setAssigned(INDEX_TYPE cellid, ASSIGNED_TYPE value) {

			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.assigned = value; //m_mesh->Compress6NeighborOffsetToByte(cellid, value);                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically

			//(*m_dgrad)[cellid].assigned = value;
		}

		// VALID ONLY AFTER FIRST ASSIGNMENT
		INDEX_TYPE getPair(INDEX_TYPE cellid) const {
			return m_mesh->UncompressByteTo6NeighborOffset(cellid, (*m_dgrad)[cellid].pair);
		}

		void setPair(INDEX_TYPE cellid, INDEX_TYPE value) {
			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.pair = m_mesh->Compress6NeighborOffsetToByte(cellid, value);                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically

			//(*m_dgrad)[cellid].pair = m_mesh->Compress6NeighborOffsetToByte(cellid, value);
		}

		bool getCritical(INDEX_TYPE cellid) const {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.pair == 7;
			//return (*m_dgrad)[cellid].pair == 7;
		}
		void setCritical(INDEX_TYPE cellid, bool value) {

			if (value) {
				GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
				temp.pair = 7;                            // modify the bitfield
				m_dgrad->SetLabel(cellid, temp);              // store it back atomically
			}
			//(*m_dgrad)[cellid].pair = 7;
		}

		DIM_TYPE getDimAscMan(INDEX_TYPE cellid) {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.ldir;
			//return (*m_dgrad)[cellid].ldir;
		}
		//void setDimAscMan(INDEX_TYPE cellid, DIM_TYPE value){
		//	(*m_dgrad)[cellid].ldir = value;
		//}
		void setDimAscMan(INDEX_TYPE cellid, DIM_TYPE value) {
			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.ldir = value;                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically
		}
		// VALID ONLY BEFORE FIRST ASSIGNEMT
		INDEX_TYPE getNumUnpairedFacets(INDEX_TYPE cellid) {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.ldir;
			//return (*m_dgrad)[cellid].ldir;
		}
		void setNumUnpairedFacets(INDEX_TYPE cellid, INDEX_TYPE value) {
			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.ldir = value;                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically

			//(*m_dgrad)[cellid].ldir = value;
		}

		unsigned char getNumUnpairedFacets_new(INDEX_TYPE cellid) {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.ldir;
			//return (*m_dgrad)[cellid].ldir;
		}
		unsigned char decrementNumUnpairedFacets_new(INDEX_TYPE cellid) {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			--label.ldir;                                  // modify the field
			m_dgrad->SetLabel(cellid, label);              // atomic store
			return label.ldir;
			//return --((*m_dgrad)[cellid].ldir);
		}
		void setNumUnpairedFacets_new(INDEX_TYPE cellid, unsigned char value) {
			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.ldir = value;                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically
			//(*m_dgrad)[cellid].ldir = value;
		}

		DIM_TYPE getMark(INDEX_TYPE cellid) {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.pair;
			//return (*m_dgrad)[cellid].pair;
		}
		void setMark(INDEX_TYPE cellid, DIM_TYPE value) {
			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.pair = value;                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically

			//(*m_dgrad)[cellid].pair = value;
		}



		void outputToRenderer(const char* filenamebase) {

			char grad_name[2048];
			sprintf(grad_name, "%s.grad", filenamebase);

			FILE* fgrad = fopen(grad_name, "wb");

			for (int i = 0; i < m_mesh->numCells(); i++) {
				GradBitfield label = m_dgrad->GetLabel(i);
				label.flag = (bool)m_mesh->boundaryValue(i);
				//(*m_dgrad)[i].flag = (bool)m_mesh->boundaryValue(i);
				//if (i % 1000 == 0)
				//	printf("cell[%d]:(%d, %d, %d, %d)\n", i,
				//	m_dgrad[i].pair,
				//	m_dgrad[i].ldir,
				//	m_dgrad[i].flag,
				//	m_dgrad[i].assigned);


				fwrite(&label, sizeof(GradBitfield), 1, fgrad);
				//fwrite(&(*m_dgrad)[i], sizeof(GradBitfield), 1, fgrad);
			}
			fclose(fgrad);
		}

		void outputToFile(const char* filenamebase) {


			FILE* fgrad = fopen(filenamebase, "wb");
			printf("Going to output %lu cells in %s...\n", m_mesh->numCells(), filenamebase);
			INDEX_TYPE num_cells = m_mesh->numCells();
#pragma omp parallel for
			for (INDEX_TYPE i = 0; i < num_cells; i++) {
				GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
				temp.flag = (bool)m_mesh->boundaryValue(i);   // modify the bitfield
				m_dgrad->SetLabel(cellid, temp);
				//(*m_dgrad)[i].flag = (bool)m_mesh->boundaryValue(i);
			}

			//(*m_dgrad)[i].flag = 0;
			//(*m_dgrad)[i].assigned = 0;
			//(*m_dgrad)[i].pair = 0;
			//(*m_dgrad)[i].ldir = 0;
		fwrite(m_dgrad->LabelArray(), sizeof(GradBitfield), m_mesh->numCells(), fgrad);
			//for (INDEX_TYPE i = 0; i < m_mesh->numCells(); ++i) {
			//	GradBitfield val = m_dgrad->GetLabel(i);  // atomic load
			//	fwrite(&val, sizeof(GradBitfield), 1, fgrad);
			//}

			printf("   -- Output to %s.. %lu values of size %d\n", filenamebase, m_mesh->numCells(), sizeof(GradBitfield));
			fclose(fgrad);
		}


		bool load_from_file(const char* filename) {

			FILE* fdat = fopen(filename, "rb");
			if (fdat == NULL) {
				return false;

			}
			//for (INDEX_TYPE i = 0; i < m_mesh->numCells(); ++i) {
			//	GradBitfield val;
			//	fread(&val, sizeof(GradBitfield), 1, fdat);
			//	m_dgrad->SetLabel(i, val);  // atomic store
			//}
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
//#endif
//
//
//#ifndef DISCRETE_GRADIENT_LABELING_H
//#define DISCRETE_GRADIENT_LABELING_H

//#include "gi_basic_types.h"
//#include "gi_regular_grid.h"
//#include "gi_labeling.h"

namespace GInt {

	//struct GradBitfield {
	//	unsigned char assigned : 1;
	//	unsigned char flag : 1;
	//	//unsigned char critical : 1;
	//	//unsigned char insorter : 1;
	//	//unsigned char dimA : 3;
	//	unsigned char pair : 3;
	//	unsigned char ldir : 3;
	//};

	template <class TopoMeshType, class LabelingType = DenseLabelingAtomic<GradBitfield> >
	class DiscreteGradientLabelingAtomic {
	public:
		TopoMeshType* m_mesh;
		LabelingType* m_dgrad;

		DiscreteGradientLabelingAtomic(TopoMeshType* mesh) : m_mesh(mesh) {
			m_dgrad = new LabelingType(m_mesh->numCells());
		}
		~DiscreteGradientLabelingAtomic() {
			delete m_dgrad;
		}
		DiscreteGradientLabelingAtomic(TopoMeshType* mesh, LabelingType* dgrad) : m_mesh(mesh) {
			m_dgrad = dgrad;
		}
		////GradBitfield& GetBitField(INDEX_TYPE id) {
		////	return m_dgrad->GetLabel[id];
		////}
		void copyValues(const DiscreteGradientLabelingAtomic<TopoMeshType>* other) {
			if (m_dgrad->GetNumLabels() != other->m_dgrad->GetNumLabels()) {
				printf("Error: copy needs same sizes\n");
			}
			m_dgrad->CopyValues(other->m_dgrad);
		}

		GradBitfield& getBitfield(INDEX_TYPE cellid) { return (*m_dgrad)[cellid]; }
		char getAsChar(INDEX_TYPE cellid) {
			char* val;
			GradBitfield vv = m_dgrad->GetLabel(cellid);//(*m_dgrad)[cellid];
			val = (char*)&vv;
			return *val;
		}
		void clearGrad(INDEX_TYPE cellid) {
			m_dgrad->SetLabel(cellid, GradBitfield{ 0 });              // store it back atomically

			//(*m_dgrad)[cellid] = GradBitfield{ 0 };
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
			m_dgrad->SetAll(0);
			//memset(m_dgrad->LabelArray(), 0, sizeof(GradBitfield)*m_mesh->numCells());
		}
		// ALWAYS VALID 
		ASSIGNED_TYPE getAssigned(INDEX_TYPE cellid) const {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.assigned;
			//return (*m_dgrad)[cellid].assigned;
		}

		void setAssigned(INDEX_TYPE cellid, ASSIGNED_TYPE value) {

			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.assigned = m_mesh->Compress6NeighborOffsetToByte(cellid, value);                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically

			//(*m_dgrad)[cellid].assigned = value;
		}

		// VALID ONLY AFTER FIRST ASSIGNMENT
		INDEX_TYPE getPair(INDEX_TYPE cellid) const {
			return m_mesh->UncompressByteTo6NeighborOffset(cellid, (*m_dgrad)[cellid].pair);
		}

		void setPair(INDEX_TYPE cellid, INDEX_TYPE value) {
			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.pair = m_mesh->Compress6NeighborOffsetToByte(cellid, value);                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically

			//(*m_dgrad)[cellid].pair = m_mesh->Compress6NeighborOffsetToByte(cellid, value);
		}

		bool getCritical(INDEX_TYPE cellid) const {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.pair == 7;
			//return (*m_dgrad)[cellid].pair == 7;
		}
		void setCritical(INDEX_TYPE cellid, bool value) {

			if (value) {
				GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
				temp.pair = 7;                            // modify the bitfield
				m_dgrad->SetLabel(cellid, temp);              // store it back atomically
			}
			//(*m_dgrad)[cellid].pair = 7;
		}

		DIM_TYPE getDimAscMan(INDEX_TYPE cellid) {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.ldir;
			//return (*m_dgrad)[cellid].ldir;
		}
		//void setDimAscMan(INDEX_TYPE cellid, DIM_TYPE value){
		//	(*m_dgrad)[cellid].ldir = value;
		//}
		void setDimAscMan(INDEX_TYPE cellid, DIM_TYPE value) {
			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.ldir = value;                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically
		}
		// VALID ONLY BEFORE FIRST ASSIGNEMT
		INDEX_TYPE getNumUnpairedFacets(INDEX_TYPE cellid) {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.ldir;
			//return (*m_dgrad)[cellid].ldir;
		}
		void setNumUnpairedFacets(INDEX_TYPE cellid, INDEX_TYPE value) {
			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.ldir = value;                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically

			//(*m_dgrad)[cellid].ldir = value;
		}

		unsigned char getNumUnpairedFacets_new(INDEX_TYPE cellid) {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.ldir;
			//return (*m_dgrad)[cellid].ldir;
		}
		unsigned char decrementNumUnpairedFacets_new(INDEX_TYPE cellid) {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			--label.ldir;                                  // modify the field
			m_dgrad->SetLabel(cellid, label);              // atomic store
			return label.ldir;
			//return --((*m_dgrad)[cellid].ldir);
		}
		void setNumUnpairedFacets_new(INDEX_TYPE cellid, unsigned char value) {
			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.ldir = value;                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically
			//(*m_dgrad)[cellid].ldir = value;
		}

		DIM_TYPE getMark(INDEX_TYPE cellid) {
			GradBitfield label = m_dgrad->GetLabel(cellid);  // atomic load
			return label.pair;
			//return (*m_dgrad)[cellid].pair;
		}
		void setMark(INDEX_TYPE cellid, DIM_TYPE value) {
			GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
			temp.pair = value;                            // modify the bitfield
			m_dgrad->SetLabel(cellid, temp);              // store it back atomically

			//(*m_dgrad)[cellid].pair = value;
		}



		void outputToRenderer(const char* filenamebase) {

			char grad_name[2048];
			sprintf(grad_name, "%s.grad", filenamebase);

			FILE* fgrad = fopen(grad_name, "wb");

			for (int i = 0; i < m_mesh->numCells(); i++) {
				GradBitfield label = m_dgrad->GetLabel(i);
				label.flag = (bool)m_mesh->boundaryValue(i);
				//(*m_dgrad)[i].flag = (bool)m_mesh->boundaryValue(i);
				//if (i % 1000 == 0)
				//	printf("cell[%d]:(%d, %d, %d, %d)\n", i,
				//	m_dgrad[i].pair,
				//	m_dgrad[i].ldir,
				//	m_dgrad[i].flag,
				//	m_dgrad[i].assigned);


				fwrite(&label, sizeof(GradBitfield), 1, fgrad);
				//fwrite(&(*m_dgrad)[i], sizeof(GradBitfield), 1, fgrad);
			}
			fclose(fgrad);
		}

		void outputToFile(const char* filenamebase) {


			FILE* fgrad = fopen(filenamebase, "wb");
			printf("Going to output %lu cells in %s...\n", m_mesh->numCells(), filenamebase);
			INDEX_TYPE num_cells = m_mesh->numCells();
#pragma omp parallel for
			for (INDEX_TYPE i = 0; i < num_cells; i++) {
				GradBitfield temp = m_dgrad->GetLabel(cellid);  // load the current value
				temp.flag = (bool)m_mesh->boundaryValue(i);   // modify the bitfield
				m_dgrad->SetLabel(cellid, temp);
				//(*m_dgrad)[i].flag = (bool)m_mesh->boundaryValue(i);
			}

			//(*m_dgrad)[i].flag = 0;
			//(*m_dgrad)[i].assigned = 0;
			//(*m_dgrad)[i].pair = 0;
			//(*m_dgrad)[i].ldir = 0;
		//fwrite(m_dgrad->LabelArray(), sizeof(GradBitfield), m_mesh->numCells(), fgrad);
			for (INDEX_TYPE i = 0; i < m_mesh->numCells(); ++i) {
				GradBitfield val = m_dgrad->GetLabel(i);  // atomic load
				fwrite(&val, sizeof(GradBitfield), 1, fgrad);
			}

			printf("   -- Output to %s.. %lu values of size %d\n", filenamebase, m_mesh->numCells(), sizeof(GradBitfield));
			fclose(fgrad);
		}


		bool load_from_file(const char* filename) {

			FILE* fdat = fopen(filename, "rb");
			if (fdat == NULL) {
				return false;

			}
			for (INDEX_TYPE i = 0; i < m_mesh->numCells(); ++i) {
				GradBitfield val;
				fread(&val, sizeof(GradBitfield), 1, fdat);
				m_dgrad->SetLabel(i, val);  // atomic store
			}
			//fread(m_dgrad->LabelArray(), sizeof(GradBitfield), m_mesh->numCells(), fdat);

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
