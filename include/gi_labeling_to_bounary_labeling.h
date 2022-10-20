#ifndef VERTEX_LABELING_TO_BOUNDARY_LABELING_H
#define VERTEX_LABELING_TO_BOUNDARY_LABELING_H

#include "gi_basic_types.h"
#include "gi_regular_grid.h"
#include "gi_topological_regular_grid.h"
#include "gi_labeling.h"
#include "gi_array_index_partition.h"

namespace GInt {
	



	template<typename INPUT_LABEL_TYPE, class MaxVLType>
	class VertexLabelingToBoundaryLabeling {
	protected:
		DenseLabeling<char>* m_output_labels;	
		TopologicalRegularGrid3D* m_topological_grid;
	public:

		VertexLabelingToBoundaryLabeling(TopologicalRegularGrid3D* topological_grid, DenseLabeling<char>* output_labels) :
		m_topological_grid(topological_grid) {
			m_output_labels = output_labels;
		}

		VertexLabelingToBoundaryLabeling(TopologicalRegularGrid3D* topological_grid) :
			m_topological_grid(topological_grid) {
			m_output_labels = new DenseLabeling<char>(m_topological_grid->numCells());
		}
        void OutputEdgesToFile(const char* filename) {
            FILE* fout = fopen(filename, "wb");

            TopologicalRegularGrid3D::DCellsIterator edges(m_topological_grid, 1);
            for (edges.begin(); edges.valid(); edges.advance()) {
                INDEX_TYPE edge = edges.value();
                if (m_output_labels->GetLabel(edge) == 1)
                    fwrite(&edge, sizeof(INDEX_TYPE), 1, fout);
            }
            fclose(fout);

        }

		DenseLabeling<char>* GetOutputLabels() { return m_output_labels; }
		void InitializeFirst() {
			m_output_labels->SetAll(0);
		}
		DenseLabeling<char>* ComputeMINBoundary(DenseLabeling<INPUT_LABEL_TYPE>* input_labels) {

			//input_labels->OutputToFile("INPUTE_LABELS.raw");
#pragma omp parallel
			{
				int num_threads = omp_get_num_threads();
				int thread_num = omp_get_thread_num();
				std::vector<INDEX_TYPE> partition;
				ArrayIndexPartitioner::EvenChunkSplit(m_topological_grid->numCells(), num_threads, partition);

				TopologicalRegularGrid3D::DCellsIterator edges(m_topological_grid, 1, partition[thread_num], partition[thread_num + 1]);
				for (edges.begin(); edges.valid(); edges.advance()) {
					TopologicalRegularGrid3D::FacetsIterator vertices(m_topological_grid);
					INDEX_TYPE edge = edges.value();
					vertices.begin(edge);
					INDEX_TYPE vertex1 = vertices.value();
					vertices.advance();
					INDEX_TYPE vertex2 = vertices.value();

					INDEX_TYPE vertex_number1 = m_topological_grid->VertexNumberFromCellID(vertex1);
					INDEX_TYPE vertex_number2 = m_topological_grid->VertexNumberFromCellID(vertex2);
					
					int lab1 = input_labels->GetLabel(vertex_number1);
					int lab2 = input_labels->GetLabel(vertex_number2);
					if (lab1 == -1 || lab2 == -1) {
						//printf("GOT HERE: lab1 %d lab2 %d\n", lab1, lab2);
						continue;
					}
					if (lab1 != lab2) {
						(*m_output_labels)[edge] = 1;
					}

				}
#pragma omp barrier
				TopologicalRegularGrid3D::DCellsIterator quads(m_topological_grid, 2, partition[thread_num], partition[thread_num + 1]);
				for (quads.begin(); quads.valid(); quads.advance()) {
					INDEX_TYPE quad = quads.value();
					TopologicalRegularGrid3D::FacetsIterator quadedges(m_topological_grid);
					
					for (quadedges.begin(quad); quadedges.valid(); quadedges.advance()) {
						if ((*m_output_labels)[quadedges.value()] == 1) {
							(*m_output_labels)[quad] = 1;
							break;
						}
					}
					

				}
#pragma omp barrier
				TopologicalRegularGrid3D::DCellsIterator voxels(m_topological_grid, 3, partition[thread_num], partition[thread_num + 1]);
				for (voxels.begin(); voxels.valid(); voxels.advance()) {
					INDEX_TYPE voxel = voxels.value();
					TopologicalRegularGrid3D::FacetsIterator voxelquads(m_topological_grid);
					for (voxelquads.begin(voxel); voxelquads.valid(); voxelquads.advance()) {
						if ((*m_output_labels)[voxelquads.value()] == 1) {
							(*m_output_labels)[voxel] = 1;
							break;
						}
					}

				}

					
			}
			return m_output_labels;
		}

		void ComputeMAXBoundary(DenseLabeling<INPUT_LABEL_TYPE>* input_labels, MaxVLType* maxv_labeling) {
#pragma omp parallel
			{
				int num_threads = omp_get_num_threads();
				int thread_num = omp_get_thread_num();
				std::vector<INDEX_TYPE> partition;
				ArrayIndexPartitioner::EvenChunkSplit(m_topological_grid->numCells(), num_threads, partition);

				TopologicalRegularGrid3D::DCellsIterator quads(m_topological_grid, 2, partition[thread_num], partition[thread_num + 1]);
				for (quads.begin(); quads.valid(); quads.advance()) {
					TopologicalRegularGrid3D::CofacetsIterator hexs(m_topological_grid);
					INDEX_TYPE quad = quads.value();
					if (m_topological_grid->boundaryValue(quad) != 0) continue;
					hexs.begin(quad);
					INDEX_TYPE hex1 = hexs.value();
					hexs.advance();
					INDEX_TYPE hex2 = hexs.value();

					INDEX_TYPE v1gid = maxv_labeling->Cell2HighestVertex(hex1);
					INDEX_TYPE v2gid = maxv_labeling->Cell2HighestVertex(hex2);
					if (v1gid == v2gid) continue; // they are part of same lower star so no worries here

					INDEX_TYPE vertex_number1 = m_topological_grid->VertexNumberFromCellID(v1gid);
					INDEX_TYPE vertex_number2 = m_topological_grid->VertexNumberFromCellID(v2gid);

					auto lab1 = input_labels->GetLabel(vertex_number1);
					auto lab2 = input_labels->GetLabel(vertex_number2);
					if (lab1 == -1 || lab2 == -1) continue;
					
					if (lab1 != lab2) {
						(*m_output_labels)[quad] += 2;
					}

				}
#pragma omp barrier
				TopologicalRegularGrid3D::DCellsIterator edges(m_topological_grid, 1, partition[thread_num], partition[thread_num + 1]);
				for (edges.begin(); edges.valid(); edges.advance()) {
					INDEX_TYPE edge = edges.value();
					if (m_topological_grid->boundaryValue(edge) != 0) continue;
					TopologicalRegularGrid3D::CofacetsIterator quadedges(m_topological_grid);

					for (quadedges.begin(edge); quadedges.valid(); quadedges.advance()) {
						if ((*m_output_labels)[quadedges.value()] > 1) {
							(*m_output_labels)[edge] += 2;
							break;
						}
					}


				}
#pragma omp barrier
				TopologicalRegularGrid3D::DCellsIterator vertices(m_topological_grid, 0, partition[thread_num], partition[thread_num + 1]);
				for (vertices.begin(); vertices.valid(); vertices.advance()) {
					INDEX_TYPE vertex = vertices.value();
					if (m_topological_grid->boundaryValue(vertex) != 0) continue;
					TopologicalRegularGrid3D::CofacetsIterator voxelquads(m_topological_grid);
					for (voxelquads.begin(vertex); voxelquads.valid(); voxelquads.advance()) {
						if ((*m_output_labels)[voxelquads.value()] > 1) {
							(*m_output_labels)[vertex] += 2;
							break;
						}
					}

				}


			}
			return;// m_output_labels;
		}

		DenseLabeling<char>* ComputeBoundaryHACK() {
#pragma omp parallel
			{
				int num_threads = omp_get_num_threads();
				int thread_num = omp_get_thread_num();
				std::vector<INDEX_TYPE> partition;
				ArrayIndexPartitioner::EvenChunkSplit(m_topological_grid->numCells(), num_threads, partition);

				TopologicalRegularGrid3D::DCellsIterator edges(m_topological_grid, 1, partition[thread_num], partition[thread_num + 1]);
				for (edges.begin(); edges.valid(); edges.advance()) {
					TopologicalRegularGrid3D::FacetsIterator vertices(m_topological_grid);
					INDEX_TYPE edge = edges.value();
					vertices.begin(edge);
					INDEX_TYPE vertex1 = vertices.value();
					vertices.advance();
					INDEX_TYPE vertex2 = vertices.value();

					INDEX_TYPE vertex_number1 = m_topological_grid->VertexNumberFromCellID(vertex1);
					INDEX_TYPE vertex_number2 = m_topological_grid->VertexNumberFromCellID(vertex2);

					if ((*(this->m_input_labels))[vertex_number1] != (*(this->m_input_labels))[vertex_number2]) {
						INDEX_TYPE lv = ((*(this->m_input_labels))[vertex_number1] > (*(this->m_input_labels))[vertex_number2] ? (*(this->m_input_labels))[vertex_number1] : (*(this->m_input_labels))[vertex_number2]);
						(*(this->m_output_labels))[edge] = (lv % 126 + 1);
					}

				}
#pragma omp barrier
				TopologicalRegularGrid3D::DCellsIterator quads(m_topological_grid, 2, partition[thread_num], partition[thread_num + 1]);
				for (quads.begin(); quads.valid(); quads.advance()) {
					INDEX_TYPE quad = quads.value();
					TopologicalRegularGrid3D::FacetsIterator quadedges(m_topological_grid);

					for (quadedges.begin(quad); quadedges.valid(); quadedges.advance()) {
						if ((*m_output_labels)[quadedges.value()] == 1) {
							(*m_output_labels)[quad] = 1;
							break;
						}
					}


				}
#pragma omp barrier
				TopologicalRegularGrid3D::DCellsIterator voxels(m_topological_grid, 3, partition[thread_num], partition[thread_num + 1]);
				for (voxels.begin(); voxels.valid(); voxels.advance()) {
					INDEX_TYPE voxel = voxels.value();
					TopologicalRegularGrid3D::FacetsIterator voxelquads(m_topological_grid);
					for (voxelquads.begin(voxel); voxelquads.valid(); voxelquads.advance()) {
						if ((*m_output_labels)[voxelquads.value()] == 1) {
							(*m_output_labels)[voxel] = 1;
							break;
						}
					}

				}


			}
			return m_output_labels;
		}



		// takes a labeling from each vertex of a GRID, and then makes a MESH labeling where 
		// cells with 2 differently labeled vertices are type 1, 
		// 3 differently labeled are type 2, and > 3 is min (cell dimension, count) (so up to 3)
		DenseLabeling<char>* ComputeRegionBoundaryKind(DenseLabeling<INPUT_LABEL_TYPE>* input_labels) {
#pragma omp parallel
			{
				int num_threads = omp_get_num_threads();
				int thread_num = omp_get_thread_num();
				std::vector<INDEX_TYPE> partition;
				ArrayIndexPartitioner::EvenChunkSplit(m_topological_grid->numCells(), num_threads, partition);

				TopologicalRegularGrid3D::DCellsIterator edges(m_topological_grid, 1, partition[thread_num], partition[thread_num + 1]);
				for (edges.begin(); edges.valid(); edges.advance()) {
					TopologicalRegularGrid3D::FacetsIterator vertices(m_topological_grid);
					INDEX_TYPE edge = edges.value();
					vertices.begin(edge);
					INDEX_TYPE vertex1 = vertices.value();
					vertices.advance();
					INDEX_TYPE vertex2 = vertices.value();

					INDEX_TYPE vertex_number1 = m_topological_grid->VertexNumberFromCellID(vertex1);
					INDEX_TYPE vertex_number2 = m_topological_grid->VertexNumberFromCellID(vertex2);

					auto lab1 = input_labels->GetLabel(vertex_number1);
					auto lab2 = input_labels->GetLabel(vertex_number2);
					if (lab1 == -1 || lab2 == -1) continue;

					if (lab1 != lab2) {
						(*m_output_labels)[edge] = 1;
					}
					else {
						(*m_output_labels)[edge] = 0;
					}

				}
#pragma omp barrier
				TopologicalRegularGrid3D::DCellsIterator quads(m_topological_grid, 2, partition[thread_num], partition[thread_num + 1]);
				for (quads.begin(); quads.valid(); quads.advance()) {
					INDEX_TYPE quad = quads.value();
					TopologicalRegularGrid3D::FacetsIterator quadedges(m_topological_grid);
					int sum = 0;
					for (quadedges.begin(quad); quadedges.valid(); quadedges.advance()) {
						if ((*m_output_labels)[quadedges.value()] == 1) {
							sum++;
						}
					}
					if (sum == 0) {
						(*m_output_labels)[quad] = 0;
					}
					else if (sum == 1 || sum == 2) {
						(*m_output_labels)[quad] = 1;
					}
					else {
						// now we have to check if it *REALLY* is a triple bounday
						INPUT_LABEL_TYPE vert_types[4];
						int num_unique = 0;
						TopologicalRegularGrid3D::CellVerticesIterator vertices(m_topological_grid);
						for (vertices.begin(quad); vertices.valid(); vertices.advance()) {
							auto v_cellid = vertices.value();
							INDEX_TYPE vertex_number = m_topological_grid->VertexNumberFromCellID(v_cellid);
							INPUT_LABEL_TYPE lab1 = input_labels->GetLabel(vertex_number);

							bool is_unique = true;
							for (int i = 0; i < num_unique; i++) {
								if (vert_types[i] == lab1) {
									is_unique = false;
									break;
								}
							}
							if (is_unique) {
								vert_types[num_unique] = lab1;
								num_unique++;
								if (num_unique > 2) break;
							}
						}

						if (num_unique > 2) {
							(*m_output_labels)[quad] = 2;
						}
						else {
							(*m_output_labels)[quad] = 1;
						}
					}
				}
#pragma omp barrier
				TopologicalRegularGrid3D::DCellsIterator voxels(m_topological_grid, 3, partition[thread_num], partition[thread_num + 1]);
				for (voxels.begin(); voxels.valid(); voxels.advance()) {
					INDEX_TYPE voxel = voxels.value();
					TopologicalRegularGrid3D::FacetsIterator voxelquads(m_topological_grid);
					int sum = 0;
					int maxv = 0;
					for (voxelquads.begin(voxel); voxelquads.valid(); voxelquads.advance()) {
						auto tmp = (*m_output_labels)[voxelquads.value()];
						if (tmp > maxv) {
							maxv = tmp;
						}
						if (tmp == 2)
							sum++;
					}
					if (sum > 2) {
						// now we have to check if it *REALLY* is a triple bounday
						INPUT_LABEL_TYPE vert_types[4];
						int num_unique = 0;
						TopologicalRegularGrid3D::CellVerticesIterator vertices(m_topological_grid);
						for (vertices.begin(voxel); vertices.valid(); vertices.advance()) {
							auto v_cellid = vertices.value();
							INDEX_TYPE vertex_number = m_topological_grid->VertexNumberFromCellID(v_cellid);
							INPUT_LABEL_TYPE lab1 = input_labels->GetLabel(vertex_number);

							bool is_unique = true;
							for (int i = 0; i < num_unique; i++) {
								if (vert_types[i] == lab1) {
									is_unique = false;
									break;
								}
							}
							if (is_unique) {
								vert_types[num_unique] = lab1;
								num_unique++;
								if (num_unique > 3) break;
							}
						}

						if (num_unique > 3) {
							(*m_output_labels)[voxel] = 3;
						}
						else {
							(*m_output_labels)[voxel] = 2;
						}
					}
					else {
						(*m_output_labels)[voxel] = maxv;
					}
				}


			}
			return m_output_labels;
		}
		};
}

#endif
