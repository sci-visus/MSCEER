#ifndef GI_MORSE_SMALE_COMPLEX_PARTITIONED_H
#define GI_MORSE_SMALE_COMPLEX_PARTITIONED_H

#include <algorithm>
#include <chrono>
#include <memory>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "gi_morse_smale_complex_basic.h"
#include "gi_partitioned_topological_regular_grid.h"

namespace GInt {

	template<typename SCALAR_TYPE, class MESH_TYPE, class FUNC_TYPE, class GRAD_TYPE>
	class MorseSmaleComplexPartitioned {
	public:
		typedef MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE> BaseMscType;

		struct DelayedArcRecord {
			INDEX_TYPE lower_cell_id;
			INDEX_TYPE upper_cell_id;
			INT_TYPE lower_partition_id;
			INT_TYPE upper_partition_id;
			INT_TYPE source_partition_id;
			INT_TYPE source_upper_local_node_id;
			SCALAR_TYPE persistence_hint;
		};

		struct LineageTransferRecord {
			INDEX_TYPE representative_cell_id;
			std::vector<INDEX_TYPE> constituent_cell_ids;
		};

		class PartitionLocalMsc : public BaseMscType {
		protected:
			typedef MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE> Base;
			typedef typename Base::ArcCandidate ArcCandidate;

			const PartitionedTopologicalRegularGrid2D* m_partition_grid;
			INT_TYPE m_partition_id;
			bool m_enforce_interior_only;
			std::unordered_set<INT_TYPE> m_frozen_local_nodes;

			void trace_down_candidates_partitioned(
				const INDEX_TYPE& cellid,
				DIM_TYPE& temp_dim,
				INT_TYPE startNodeID,
				INDEX_TYPE startCellID,
				std::vector<ArcCandidate>& local_candidates,
				std::vector<DelayedArcRecord>& delayed_records) {
				if (m_partition_grid == NULL) {
					throw std::runtime_error("PartitionLocalMsc partition grid is null during candidate tracing.");
				}
				typename MESH_TYPE::FacetsIterator facets(this->mMesh);
				for (facets.begin(cellid); facets.valid(); facets.advance()) {
					INDEX_TYPE temp_id = facets.value();
					if (this->mGrad->getCritical(temp_id)) {
						const INT_TYPE lower_partition = m_partition_grid->cell_id_to_partition_num(temp_id);
						if (lower_partition == m_partition_id) {
							INT_TYPE lowerNodeID = this->nodeIdForCell(temp_id);
							if (lowerNodeID != NULLID) {
								ArcCandidate c;
								c.lowerNodeID = lowerNodeID;
								c.upperNodeID = startNodeID;
								local_candidates.push_back(c);
							}
						}
						else {
							DelayedArcRecord dr;
							dr.lower_cell_id = temp_id;
							dr.upper_cell_id = startCellID;
							dr.lower_partition_id = lower_partition;
							dr.upper_partition_id = m_partition_id;
							dr.source_partition_id = m_partition_id;
							dr.source_upper_local_node_id = startNodeID;
							const SCALAR_TYPE lower_val = this->mFunc->cellValue(temp_id);
							const SCALAR_TYPE upper_val = this->nodes[startNodeID].value;
							dr.persistence_hint = upper_val - lower_val;
							delayed_records.push_back(dr);
							// Critical policy: if a node participates in a cross-partition arc, freeze it.
							m_frozen_local_nodes.insert(startNodeID);
						}
					}
					else if (this->mGrad->getDimAscMan(temp_id) == temp_dim) {
						INDEX_TYPE pair = this->mGrad->getPair(temp_id);
						if (pair != cellid && this->mMesh->dimension(pair) == this->mMesh->dimension(cellid)) {
							trace_down_candidates_partitioned(pair, temp_dim, startNodeID, startCellID, local_candidates, delayed_records);
						}
					}
				}
			}

			void add_arcs_for_partition_node(
				INT_TYPE node_id,
				std::vector<ArcCandidate>& local_candidates,
				std::vector<DelayedArcRecord>& delayed_records) {
				node<SCALAR_TYPE>& n = this->getNode(node_id);
				if (n.dim == 0) return;
				if (!this->mBuildArcAtAll[n.dim - 1]) return;
				if (this->mBuildArcGeometry[n.dim - 1]) {
					throw std::runtime_error("Partitioned local arc construction currently requires no-geometry mode.");
				}

				const INDEX_TYPE startCellID = n.cellindex;
				DIM_TYPE temp_dim = this->mGrad->getDimAscMan(startCellID) + 1;
				trace_down_candidates_partitioned(startCellID, temp_dim, node_id, startCellID, local_candidates, delayed_records);
			}

			virtual bool isValid(INT_TYPE a, arc<SCALAR_TYPE>& ap) const {
				if (!Base::isValid(a, ap)) return false;
				if (!m_enforce_interior_only || m_partition_grid == NULL) return true;
				if (m_frozen_local_nodes.count(ap.lower) != 0) return false;
				if (m_frozen_local_nodes.count(ap.upper) != 0) return false;

				const INDEX_TYPE lower_cell = this->nodes[ap.lower].cellindex;
				const INDEX_TYPE upper_cell = this->nodes[ap.upper].cellindex;

				if (m_partition_grid->cell_id_to_partition_num(lower_cell) != m_partition_id) return false;
				if (m_partition_grid->cell_id_to_partition_num(upper_cell) != m_partition_id) return false;
				if (!m_partition_grid->is_partition_interior_cell(lower_cell, m_partition_id)) return false;
				if (!m_partition_grid->is_partition_interior_cell(upper_cell, m_partition_id)) return false;
				return true;
			}

		public:
			PartitionLocalMsc(GRAD_TYPE* grad, MESH_TYPE* mesh, FUNC_TYPE* func) :
				Base(grad, mesh, func),
				m_partition_grid(NULL),
				m_partition_id(-1),
				m_enforce_interior_only(false) {}

			void ConfigurePartitionRestrictions(
				const PartitionedTopologicalRegularGrid2D* partition_grid,
				INT_TYPE partition_id,
				bool enforce_interior_only) {
				m_partition_grid = partition_grid;
				m_partition_id = partition_id;
				m_enforce_interior_only = enforce_interior_only;
			}

			bool TryGetNodeIdForCell(INDEX_TYPE cellID, INT_TYPE& nodeID) const {
				nodeID = this->nodeIdForCell(cellID);
				return nodeID != NULLID;
			}

			bool IsNodeAliveAtLocalEnd(INT_TYPE nodeID) const {
				return this->nodes[nodeID].destroyed > this->num_cancelled;
			}

			bool IsArcAliveAtLocalEnd(INT_TYPE arcID) const {
				return this->isAlive(this->arcs[arcID], this->num_cancelled);
			}

			void CollectLivingNodeCells(std::vector<INDEX_TYPE>& out_cells) const {
				out_cells.clear();
				out_cells.reserve((size_t)this->numNodes());
				for (INT_TYPE nid = 0; nid < this->numNodes(); nid++) {
					if (!IsNodeAliveAtLocalEnd(nid)) continue;
					out_cells.push_back(this->nodes[nid].cellindex);
				}
			}

			void CollectLivingArcCellPairs(std::vector<std::pair<INDEX_TYPE, INDEX_TYPE> >& out_pairs) const {
				out_pairs.clear();
				out_pairs.reserve((size_t)this->numArcs());
				for (INT_TYPE aid = 0; aid < this->numArcs(); aid++) {
					if (!IsArcAliveAtLocalEnd(aid)) continue;
					const arc<SCALAR_TYPE>& a = this->arcs[aid];
					out_pairs.push_back(std::make_pair(this->nodes[a.lower].cellindex, this->nodes[a.upper].cellindex));
				}
			}

			void GatherRepresentativeCellsForBaseNode(INT_TYPE base_node_id, std::vector<INDEX_TYPE>& reps) const {
				reps.clear();
				if (base_node_id < 0 || base_node_id >= this->numNodes()) return;

				const node<SCALAR_TYPE>& base = this->nodes[base_node_id];
				std::unordered_set<INDEX_TYPE> rep_cells;

				if (IsNodeAliveAtLocalEnd(base_node_id)) {
					rep_cells.insert(base.cellindex);
				}

				for (INT_TYPE nid = 0; nid < this->numNodes(); nid++) {
					if (!IsNodeAliveAtLocalEnd(nid)) continue;
					if (this->nodes[nid].dim != base.dim) continue;

					bool matched = false;
					std::set<INT_TYPE> constituents;
					this->GatherNodes(nid, constituents, true);
					if (constituents.count(base_node_id) != 0) matched = true;

					if (!matched) {
						constituents.clear();
						this->GatherNodes(nid, constituents, false);
						if (constituents.count(base_node_id) != 0) matched = true;
					}

					if (matched) {
						rep_cells.insert(this->nodes[nid].cellindex);
					}
				}

				if (rep_cells.empty()) {
					rep_cells.insert(base.cellindex);
				}

				reps.assign(rep_cells.begin(), rep_cells.end());
				std::sort(reps.begin(), reps.end());
			}

			void CollectLivingExtremaLineage(bool ascending, std::vector<LineageTransferRecord>& out_records) const {
				out_records.clear();
				for (INT_TYPE nid = 0; nid < this->numNodes(); nid++) {
					if (!IsNodeAliveAtLocalEnd(nid)) continue;
					const node<SCALAR_TYPE>& n = this->nodes[nid];
					if (ascending) {
						if (n.dim != 0) continue;
					}
					else {
						if (n.dim != 2) continue;
					}

					std::set<INT_TYPE> constituents;
					this->GatherNodes(nid, constituents, ascending);
					std::unordered_set<INDEX_TYPE> cell_ids;
					for (std::set<INT_TYPE>::const_iterator it = constituents.begin(); it != constituents.end(); ++it) {
						cell_ids.insert(this->nodes[*it].cellindex);
					}
					if (cell_ids.empty()) {
						cell_ids.insert(n.cellindex);
					}

					LineageTransferRecord rec;
					rec.representative_cell_id = n.cellindex;
					rec.constituent_cell_ids.assign(cell_ids.begin(), cell_ids.end());
					std::sort(rec.constituent_cell_ids.begin(), rec.constituent_cell_ids.end());
					out_records.push_back(rec);
				}
			}

			void ComputeFromGradInPartition(
				const PartitionedTopologicalRegularGrid2D& partition_grid,
				INT_TYPE partition_id,
				std::vector<DelayedArcRecord>& delayed_records) {
				printf("[partitioned] partition=%d ComputeFromGradInPartition enter\n", (int)partition_id);
				fflush(stdout);
				delayed_records.clear();
				m_frozen_local_nodes.clear();
				printf("[partitioned] partition=%d delayed_records cleared\n", (int)partition_id);
				fflush(stdout);
				typename PartitionedTopologicalRegularGrid2D::PartitionCellsIterator pit(&partition_grid, partition_id);
				printf("[partitioned] partition=%d PartitionCellsIterator created start_id=%lld start_x=%lld start_y=%lld x_range=%lld y_range=%lld\n",
					(int)partition_id,
					(long long)pit.start_id(),
					(long long)pit.start_x(),
					(long long)pit.start_y(),
					(long long)pit.x_range(),
					(long long)pit.y_range());
				fflush(stdout);

				int created_nodes = 0;
				for (pit.begin(); pit.valid(); pit.advance()) {
					const INDEX_TYPE cid = pit.value();
					if (this->mGrad->getCritical(cid)) {
						this->createNode(cid);
						created_nodes++;
					}
				}
				printf("[partitioned] partition=%d node scan created_nodes=%d total_nodes=%d\n",
					(int)partition_id, created_nodes, (int)this->numNodes());
				fflush(stdout);

				std::vector<ArcCandidate> local_candidates;
				printf("[partitioned] partition=%d arc candidate generation total_nodes=%d\n",
					(int)partition_id, (int)this->numNodes());
				fflush(stdout);
				for (INT_TYPE nid = 0; nid < this->numNodes(); nid++) {
					add_arcs_for_partition_node(nid, local_candidates, delayed_records);
				}
				printf("[partitioned] partition=%d arc candidate generation candidates=%llu delayed=%llu\n",
					(int)partition_id,
					(unsigned long long)local_candidates.size(),
					(unsigned long long)delayed_records.size());
				fflush(stdout);

				printf("[partitioned] partition=%d arc materialization\n", (int)partition_id);
				fflush(stdout);
				for (size_t i = 0; i < local_candidates.size(); i++) {
					this->createArcFromNodeIDs(local_candidates[i].lowerNodeID, local_candidates[i].upperNodeID);
				}
				printf("[partitioned] partition=%d arc materialization total_arcs=%d\n",
					(int)partition_id, (int)this->numArcs());
				fflush(stdout);
				printf("[partitioned] partition=%d frozen_cross_boundary_nodes=%llu\n",
					(int)partition_id, (unsigned long long)m_frozen_local_nodes.size());
				fflush(stdout);
				printf("[partitioned] partition=%d ComputeFromGradInPartition exit\n", (int)partition_id);
				fflush(stdout);
			}
		};

		struct PartitionRunResult {
			INT_TYPE partition_id;
			std::unique_ptr<PartitionLocalMsc> msc;
			std::vector<DelayedArcRecord> delayed_arcs;
			std::vector<LineageTransferRecord> ascending_lineage;
			std::vector<LineageTransferRecord> descending_lineage;
		};
		struct TimingBreakdown {
			long long local_stage_total_ms;
			long long local_build_ms;
			long long local_simplify_ms;
			long long reconcile_ms;
			long long global_simplify_ms;
			TimingBreakdown() :
				local_stage_total_ms(0),
				local_build_ms(0),
				local_simplify_ms(0),
				reconcile_ms(0),
				global_simplify_ms(0) {}
		};

	protected:
		GRAD_TYPE* mGrad;
		MESH_TYPE* mMesh;
		FUNC_TYPE* mFunc;
		
		inline long long elapsedMilliseconds(
			const std::chrono::steady_clock::time_point& start,
			const std::chrono::steady_clock::time_point& end) const {
			return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		}
		inline void printAggregateTiming(const char* label,
			const std::chrono::steady_clock::time_point& start,
			const std::chrono::steady_clock::time_point& end) const {
			printf("TIMING: %s ms=%lld\n", label, elapsedMilliseconds(start, end));
		}

		template <typename OutputMscType>
		static void append_arcs_no_dedupe(
			OutputMscType* out_msc,
			const std::vector<std::pair<INDEX_TYPE, INDEX_TYPE> >& arcs) {
			for (size_t i = 0; i < arcs.size(); i++) {
				const INDEX_TYPE lower_cell = arcs[i].first;
				const INDEX_TYPE upper_cell = arcs[i].second;
				if (lower_cell == upper_cell) continue;
				(void)out_msc->CreateArcByCellIDs(lower_cell, upper_cell);
			}
		}

	public:
		class ReconciledGlobalMsc : public BaseMscType {
		protected:
			std::unordered_map<INDEX_TYPE, std::vector<INDEX_TYPE> > m_ascending_lineage_by_rep_cell;
			std::unordered_map<INDEX_TYPE, std::vector<INDEX_TYPE> > m_descending_lineage_by_rep_cell;

			static void merge_lineage_cells(
				std::unordered_map<INDEX_TYPE, std::vector<INDEX_TYPE> >& lineage_map,
				INDEX_TYPE representative_cell_id,
				const std::vector<INDEX_TYPE>& cells) {
				std::vector<INDEX_TYPE>& dst = lineage_map[representative_cell_id];
				dst.insert(dst.end(), cells.begin(), cells.end());
				std::sort(dst.begin(), dst.end());
				dst.erase(std::unique(dst.begin(), dst.end()), dst.end());
			}
		public:
			ReconciledGlobalMsc(GRAD_TYPE* grad, MESH_TYPE* mesh, FUNC_TYPE* func) :
				BaseMscType(grad, mesh, func) {}

			INT_TYPE EnsureNodeByCellID(INDEX_TYPE cellID) {
				INT_TYPE nid = this->nodeIdForCell(cellID);
				if (nid != NULLID) return nid;
				return this->createNode(cellID);
			}

			INT_TYPE CreateArcByCellIDs(INDEX_TYPE lowerCellID, INDEX_TYPE upperCellID) {
				(void)EnsureNodeByCellID(lowerCellID);
				(void)EnsureNodeByCellID(upperCellID);
				return this->createArc(lowerCellID, upperCellID);
			}

			void AddAscendingLineageCells(INDEX_TYPE representative_cell_id, const std::vector<INDEX_TYPE>& cells) {
				merge_lineage_cells(m_ascending_lineage_by_rep_cell, representative_cell_id, cells);
			}

			void AddDescendingLineageCells(INDEX_TYPE representative_cell_id, const std::vector<INDEX_TYPE>& cells) {
				merge_lineage_cells(m_descending_lineage_by_rep_cell, representative_cell_id, cells);
			}

			bool GetAscendingLineageCells(INDEX_TYPE representative_cell_id, std::vector<INDEX_TYPE>& out_cells) const {
				typename std::unordered_map<INDEX_TYPE, std::vector<INDEX_TYPE> >::const_iterator it =
					m_ascending_lineage_by_rep_cell.find(representative_cell_id);
				if (it == m_ascending_lineage_by_rep_cell.end()) return false;
				out_cells = it->second;
				return true;
			}

			bool GetDescendingLineageCells(INDEX_TYPE representative_cell_id, std::vector<INDEX_TYPE>& out_cells) const {
				typename std::unordered_map<INDEX_TYPE, std::vector<INDEX_TYPE> >::const_iterator it =
					m_descending_lineage_by_rep_cell.find(representative_cell_id);
				if (it == m_descending_lineage_by_rep_cell.end()) return false;
				out_cells = it->second;
				return true;
			}
		};

		MorseSmaleComplexPartitioned(GRAD_TYPE* grad, MESH_TYPE* mesh, FUNC_TYPE* func) :
			mGrad(grad), mMesh(mesh), mFunc(func) {}

		std::vector<PartitionRunResult> BuildPartitionLocalMSCs(
			INT_TYPE num_partitions,
			SCALAR_TYPE local_persistence_abs,
			TimingBreakdown* timings = NULL) {
			if (timings != NULL) {
				timings->local_stage_total_ms = 0;
				timings->local_build_ms = 0;
				timings->local_simplify_ms = 0;
			}
			printf("[partitioned] BuildPartitionLocalMSCs enter num_partitions=%d local_persistence_abs=%f\n",
				(int)num_partitions, (double)local_persistence_abs);
			fflush(stdout);
			const std::chrono::steady_clock::time_point localStart = std::chrono::steady_clock::now();
			printf("[partitioned] BuildPartitionLocalMSCs dynamic_cast mesh\n");
			fflush(stdout);
			const TopologicalRegularGrid2D* mesh2d = dynamic_cast<const TopologicalRegularGrid2D*>(mMesh);
			if (mesh2d == NULL) {
				printf("[partitioned] BuildPartitionLocalMSCs dynamic_cast mesh FAILED\n");
				fflush(stdout);
				throw std::runtime_error("Partitioned local MSC build currently supports TopologicalRegularGrid2D meshes.");
			}
			printf("[partitioned] BuildPartitionLocalMSCs dynamic_cast mesh ok x=%lld y=%lld\n",
				(long long)mesh2d->numCellsAxis(0), (long long)mesh2d->numCellsAxis(1));
			fflush(stdout);
			printf("[partitioned] BuildPartitionLocalMSCs constructing partition grid\n");
			fflush(stdout);
			PartitionedTopologicalRegularGrid2D partition_grid(mesh2d, num_partitions);
			printf("[partitioned] BuildPartitionLocalMSCs partition grid ready rows=%lld cols=%lld parts=%d\n",
				(long long)partition_grid.tile_rows(), (long long)partition_grid.tile_cols(), (int)partition_grid.num_partitions());
			fflush(stdout);

			const INT_TYPE part_count = partition_grid.num_partitions();
			std::vector<PartitionRunResult> results((size_t)part_count);
			std::vector<long long> per_partition_build_ms((size_t)part_count, 0);
			std::vector<long long> per_partition_simplify_ms((size_t)part_count, 0);
			printf("[partitioned] BuildPartitionLocalMSCs prepared results=%d\n", (int)part_count);
			fflush(stdout);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (INT_TYPE pid = 0; pid < part_count; pid++) {
				printf("[partitioned] BuildPartitionLocalMSCs partition loop pid=%d\n", (int)pid);
				fflush(stdout);
				const std::chrono::steady_clock::time_point partitionStart = std::chrono::steady_clock::now();
				PartitionRunResult run;
				run.partition_id = pid;
				printf("[partitioned] partition=%d allocate local MSC\n", (int)pid);
				fflush(stdout);
				run.msc.reset(new PartitionLocalMsc(mGrad, mMesh, mFunc));
				// Configure partition context before local arc tracing so cell->partition lookups are valid.
				run.msc->ConfigurePartitionRestrictions(&partition_grid, pid, false);
				printf("[partitioned] partition=%d SetBuildArcGeometry\n", (int)pid);
				fflush(stdout);
				run.msc->SetBuildArcGeometry(Vec3b(false, false, false));
				printf("[partitioned] partition=%d ComputeFromGradInPartition\n", (int)pid);
				fflush(stdout);
				const std::chrono::steady_clock::time_point localBuildStart = std::chrono::steady_clock::now();
				run.msc->ComputeFromGradInPartition(partition_grid, pid, run.delayed_arcs);
				const std::chrono::steady_clock::time_point localBuildEnd = std::chrono::steady_clock::now();
				printf("[partitioned] partition=%d ComputeFromGradInPartition nodes=%d arcs=%d delayed=%d\n",
					(int)pid, (int)run.msc->numNodes(), (int)run.msc->numArcs(), (int)run.delayed_arcs.size());
				fflush(stdout);
				printf("   -- Partition %d local graph: nodes=%d arcs=%d delayed=%d\n",
					(int)pid, (int)run.msc->numNodes(), (int)run.msc->numArcs(), (int)run.delayed_arcs.size());
				fflush(stdout);
				printf("[partitioned] partition=%d ConfigurePartitionRestrictions\n", (int)pid);
				fflush(stdout);
				run.msc->ConfigurePartitionRestrictions(&partition_grid, pid, true);
				printf("[partitioned] partition=%d ComputeHierarchy\n", (int)pid);
				fflush(stdout);
				const std::chrono::steady_clock::time_point localSimplifyStart = std::chrono::steady_clock::now();
				run.msc->ComputeHierarchy(local_persistence_abs);
				printf("[partitioned] partition=%d SetSelectPersAbs\n", (int)pid);
				fflush(stdout);
				run.msc->SetSelectPersAbs(local_persistence_abs);
				run.msc->CollectLivingExtremaLineage(true, run.ascending_lineage);
				run.msc->CollectLivingExtremaLineage(false, run.descending_lineage);
				const std::chrono::steady_clock::time_point localSimplifyEnd = std::chrono::steady_clock::now();
				per_partition_build_ms[(size_t)pid] = elapsedMilliseconds(localBuildStart, localBuildEnd);
				per_partition_simplify_ms[(size_t)pid] = elapsedMilliseconds(localSimplifyStart, localSimplifyEnd);
				const std::chrono::steady_clock::time_point partitionEnd = std::chrono::steady_clock::now();
				printAggregateTiming("Partition local build+cancel", partitionStart, partitionEnd);
				fflush(stdout);
				results[(size_t)pid] = std::move(run);
			}
			const std::chrono::steady_clock::time_point localEnd = std::chrono::steady_clock::now();
			printAggregateTiming("Partition stage total", localStart, localEnd);
			if (timings != NULL) {
				long long total_build = 0;
				long long total_simplify = 0;
				for (INT_TYPE pid = 0; pid < part_count; pid++) {
					total_build += per_partition_build_ms[(size_t)pid];
					total_simplify += per_partition_simplify_ms[(size_t)pid];
				}
				timings->local_build_ms = total_build;
				timings->local_simplify_ms = total_simplify;
				timings->local_stage_total_ms = elapsedMilliseconds(localStart, localEnd);
			}
			fflush(stdout);
			printf("[partitioned] BuildPartitionLocalMSCs exit results=%llu\n", (unsigned long long)results.size());
			fflush(stdout);
			return results;
		}

		std::unique_ptr<ReconciledGlobalMsc> BuildReconciledGlobalBase(
			const PartitionedTopologicalRegularGrid2D& partition_grid,
			const std::vector<PartitionRunResult>& partition_results,
			TimingBreakdown* timings = NULL) const {
			const std::chrono::steady_clock::time_point reconcileStart = std::chrono::steady_clock::now();
			(void)partition_grid;
			std::unique_ptr<ReconciledGlobalMsc> out(new ReconciledGlobalMsc(mGrad, mMesh, mFunc));

			for (size_t i = 0; i < partition_results.size(); i++) {
				const PartitionLocalMsc* local = partition_results[i].msc.get();
				std::vector<INDEX_TYPE> living_node_cells;
				local->CollectLivingNodeCells(living_node_cells);
				for (size_t j = 0; j < living_node_cells.size(); j++) {
					(void)out->EnsureNodeByCellID(living_node_cells[j]);
				}

				std::vector<std::pair<INDEX_TYPE, INDEX_TYPE> > living_arcs;
				local->CollectLivingArcCellPairs(living_arcs);
				append_arcs_no_dedupe(out.get(), living_arcs);

				const std::vector<LineageTransferRecord>& asc = partition_results[i].ascending_lineage;
				for (size_t j = 0; j < asc.size(); j++) {
					out->AddAscendingLineageCells(asc[j].representative_cell_id, asc[j].constituent_cell_ids);
				}
				const std::vector<LineageTransferRecord>& dsc = partition_results[i].descending_lineage;
				for (size_t j = 0; j < dsc.size(); j++) {
					out->AddDescendingLineageCells(dsc[j].representative_cell_id, dsc[j].constituent_cell_ids);
				}
			}

			int delayed_records_processed = 0;
			for (size_t i = 0; i < partition_results.size(); i++) {
				const std::vector<DelayedArcRecord>& delayed = partition_results[i].delayed_arcs;
				for (size_t j = 0; j < delayed.size(); j++) {
					delayed_records_processed++;
					const DelayedArcRecord& dr = delayed[j];
					if (dr.lower_cell_id == dr.upper_cell_id) continue;
					(void)out->CreateArcByCellIDs(dr.lower_cell_id, dr.upper_cell_id);
				}
			}
			printf("   -- Reconcile summary: partitions=%d delayed_records=%d merged_nodes=%d merged_arcs=%d\n",
				(int)partition_results.size(), delayed_records_processed, (int)out->numNodes(), (int)out->numArcs());
			const std::chrono::steady_clock::time_point reconcileEnd = std::chrono::steady_clock::now();
			printAggregateTiming("Partition reconcile", reconcileStart, reconcileEnd);
			if (timings != NULL) {
				timings->reconcile_ms = elapsedMilliseconds(reconcileStart, reconcileEnd);
			}

			return out;
		}

		std::unique_ptr<ReconciledGlobalMsc> BuildPartitionedThenContinueSerial(
			INT_TYPE num_partitions,
			SCALAR_TYPE local_persistence_abs,
			SCALAR_TYPE final_persistence_abs) {
			const std::chrono::steady_clock::time_point pipelineStart = std::chrono::steady_clock::now();
			const TopologicalRegularGrid2D* mesh2d = dynamic_cast<const TopologicalRegularGrid2D*>(mMesh);
			if (mesh2d == NULL) {
				throw std::runtime_error("Partitioned pipeline currently supports TopologicalRegularGrid2D meshes.");
			}

			PartitionedTopologicalRegularGrid2D partition_grid(mesh2d, num_partitions);
			std::vector<PartitionRunResult> local_results = BuildPartitionLocalMSCs(num_partitions, local_persistence_abs);
			std::unique_ptr<ReconciledGlobalMsc> reconciled = BuildReconciledGlobalBase(partition_grid, local_results);
			const std::chrono::steady_clock::time_point serialStart = std::chrono::steady_clock::now();
			reconciled->ComputeHierarchy(final_persistence_abs);
			reconciled->SetSelectPersAbs(final_persistence_abs);
			const std::chrono::steady_clock::time_point serialEnd = std::chrono::steady_clock::now();
			printAggregateTiming("Partition serial continuation", serialStart, serialEnd);
			const std::chrono::steady_clock::time_point pipelineEnd = std::chrono::steady_clock::now();
			printAggregateTiming("Partition pipeline total", pipelineStart, pipelineEnd);
			return reconciled;
		}
	};
}

#endif
