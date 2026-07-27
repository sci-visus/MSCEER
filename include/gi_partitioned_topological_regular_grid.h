#ifndef GI_PARTITIONED_TOPOLOGICAL_REGULAR_GRID_H
#define GI_PARTITIONED_TOPOLOGICAL_REGULAR_GRID_H

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "gi_topological_regular_grid.h"

namespace GInt {

	class PartitionedTopologicalRegularGrid2D {
	public:
		struct PartitionBounds {
			INT_TYPE partition_id;
			INDEX_TYPE start_x;
			INDEX_TYPE start_y;
			INDEX_TYPE x_range;
			INDEX_TYPE y_range;
			INDEX_TYPE start_id;
		};

		class PartitionCellsIterator {
		protected:
			const PartitionedTopologicalRegularGrid2D* m_partitioned_mesh;
			PartitionBounds m_bounds;
			INDEX_TYPE m_curr_x;
			INDEX_TYPE m_curr_y;
			INDEX_TYPE m_curr_id;
			bool m_valid;

		public:
			PartitionCellsIterator(const PartitionedTopologicalRegularGrid2D* mesh, INT_TYPE partition_number) :
				m_partitioned_mesh(mesh),
				m_bounds(mesh->partition_bounds(partition_number)),
				m_curr_x(0),
				m_curr_y(0),
				m_curr_id(0),
				m_valid(false) {}

			void begin() {
				m_curr_x = m_bounds.start_x;
				m_curr_y = m_bounds.start_y;
				if (m_bounds.x_range == 0 || m_bounds.y_range == 0) {
					m_valid = false;
					m_curr_id = m_bounds.start_id;
					return;
				}
				m_curr_id = m_bounds.start_id;
				m_valid = true;
			}

			void advance() {
				if (!m_valid) return;
				const INDEX_TYPE end_x = m_bounds.start_x + m_bounds.x_range;
				const INDEX_TYPE end_y = m_bounds.start_y + m_bounds.y_range;
				m_curr_x++;
				if (m_curr_x < end_x) {
					m_curr_id++;
					return;
				}

				m_curr_x = m_bounds.start_x;
				m_curr_y++;
				if (m_curr_y >= end_y) {
					m_valid = false;
					return;
				}
				m_curr_id = m_curr_x + m_curr_y * m_partitioned_mesh->numCellsAxis(0);
			}

			bool valid() const {
				return m_valid;
			}

			INDEX_TYPE value() const {
				return m_curr_id;
			}

			INDEX_TYPE start_id() const { return m_bounds.start_id; }
			INDEX_TYPE start_x() const { return m_bounds.start_x; }
			INDEX_TYPE start_y() const { return m_bounds.start_y; }
			INDEX_TYPE x_range() const { return m_bounds.x_range; }
			INDEX_TYPE y_range() const { return m_bounds.y_range; }
		};

	protected:
		const TopologicalRegularGrid2D* m_mesh;
		INT_TYPE m_num_partitions;
		INDEX_TYPE m_tile_rows;
		INDEX_TYPE m_tile_cols;
		std::vector<PartitionBounds> m_partition_bounds;
		std::vector<INT_TYPE> m_x_to_tile_col;
		std::vector<INT_TYPE> m_y_to_tile_row;

		static void tile_shape_for_partitions(INT_TYPE num_partitions, INDEX_TYPE& tile_rows, INDEX_TYPE& tile_cols) {
			switch (num_partitions) {
			case 1: tile_rows = 1; tile_cols = 1; return;
			case 2: tile_rows = 1; tile_cols = 2; return;
			case 3: tile_rows = 1; tile_cols = 3; return;
			case 4: tile_rows = 2; tile_cols = 2; return;
			case 6: tile_rows = 2; tile_cols = 3; return;
			case 8: tile_rows = 2; tile_cols = 4; return;
			case 9: tile_rows = 3; tile_cols = 3; return;
			case 12: tile_rows = 3; tile_cols = 4; return;
			case 16: tile_rows = 4; tile_cols = 4; return;
			default:
				throw std::invalid_argument("Unsupported partition count. Supported values are {1,2,3,4,6,8,9,12,16}.");
			}
		}

		static void make_axis_partitions(INDEX_TYPE axis_size, INDEX_TYPE num_tiles,
			std::vector<INDEX_TYPE>& starts, std::vector<INDEX_TYPE>& ranges) {
			starts.assign(num_tiles, 0);
			ranges.assign(num_tiles, 0);
			const INDEX_TYPE base = axis_size / num_tiles;
			const INDEX_TYPE rem = axis_size % num_tiles;
			INDEX_TYPE cursor = 0;
			for (INDEX_TYPE i = 0; i < num_tiles; i++) {
				const INDEX_TYPE width = base + (i < rem ? 1 : 0);
				starts[i] = cursor;
				ranges[i] = width;
				cursor += width;
			}
		}

		void initialize_mapping() {
			INDEX_TYPE tile_rows = 0, tile_cols = 0;
			tile_shape_for_partitions(m_num_partitions, tile_rows, tile_cols);
			m_tile_rows = tile_rows;
			m_tile_cols = tile_cols;

			const INDEX_TYPE xdim = numCellsAxis(0);
			const INDEX_TYPE ydim = numCellsAxis(1);

			std::vector<INDEX_TYPE> x_starts, x_ranges, y_starts, y_ranges;
			make_axis_partitions(xdim, m_tile_cols, x_starts, x_ranges);
			make_axis_partitions(ydim, m_tile_rows, y_starts, y_ranges);

			m_partition_bounds.clear();
			m_partition_bounds.reserve((size_t)m_num_partitions);
			for (INDEX_TYPE ty = 0; ty < m_tile_rows; ty++) {
				for (INDEX_TYPE tx = 0; tx < m_tile_cols; tx++) {
					PartitionBounds b;
					b.partition_id = (INT_TYPE)(ty * m_tile_cols + tx);
					b.start_x = x_starts[tx];
					b.start_y = y_starts[ty];
					b.x_range = x_ranges[tx];
					b.y_range = y_ranges[ty];
					b.start_id = b.start_x + b.start_y * xdim;
					m_partition_bounds.push_back(b);
				}
			}

			m_x_to_tile_col.assign((size_t)xdim, 0);
			for (INDEX_TYPE tx = 0; tx < m_tile_cols; tx++) {
				const INDEX_TYPE start_x = x_starts[tx];
				const INDEX_TYPE end_x = start_x + x_ranges[tx];
				for (INDEX_TYPE x = start_x; x < end_x; x++) {
					m_x_to_tile_col[(size_t)x] = (INT_TYPE)tx;
				}
			}
			m_y_to_tile_row.assign((size_t)ydim, 0);
			for (INDEX_TYPE ty = 0; ty < m_tile_rows; ty++) {
				const INDEX_TYPE start_y = y_starts[ty];
				const INDEX_TYPE end_y = start_y + y_ranges[ty];
				for (INDEX_TYPE y = start_y; y < end_y; y++) {
					m_y_to_tile_row[(size_t)y] = (INT_TYPE)ty;
				}
			}
		}

	public:
		PartitionedTopologicalRegularGrid2D(const TopologicalRegularGrid2D* mesh, INT_TYPE num_partitions) :
			m_mesh(mesh), m_num_partitions(num_partitions), m_tile_rows(0), m_tile_cols(0) {
			if (m_mesh == NULL) throw std::invalid_argument("PartitionedTopologicalRegularGrid2D requires a non-null mesh.");
			initialize_mapping();
		}

		static bool IsSupportedPartitionCount(INT_TYPE num_partitions) {
			switch (num_partitions) {
			case 1:
			case 2:
			case 3:
			case 4:
			case 6:
			case 8:
			case 9:
			case 12:
			case 16:
				return true;
			default:
				return false;
			}
		}

		INDEX_TYPE numCellsAxis(int axis) const { return m_mesh->numCellsAxis(axis); }
		INDEX_TYPE numCells() const { return m_mesh->numCells(); }
		INDEX_TYPE numCells(DIM_TYPE dim) const { return m_mesh->numCells(dim); }
		const Vec2l& XY() const { return m_mesh->XY(); }

		void cellid2Coords(INDEX_TYPE cellid, Vec2l& coords) const { m_mesh->cellid2Coords(cellid, coords); }
		INDEX_TYPE coords2Cellid(const Vec2l& coords) const { return m_mesh->coords2Cellid(coords); }
		DIM_TYPE dimension(const INDEX_TYPE& cellid) const { return m_mesh->dimension(cellid); }
		DIM_TYPE dimension(const Vec2l& coords) const { return m_mesh->dimension(coords); }
		BOUNDARY_TYPE boundaryValue(INDEX_TYPE cellid) const { return m_mesh->boundaryValue(cellid); }
		BOUNDARY_TYPE boundaryValue(const Vec2l& coords) const { return m_mesh->boundaryValue(coords); }

		INT_TYPE num_partitions() const { return m_num_partitions; }
		INDEX_TYPE tile_rows() const { return m_tile_rows; }
		INDEX_TYPE tile_cols() const { return m_tile_cols; }

		const PartitionBounds& partition_bounds(INT_TYPE partition_number) const {
			if (partition_number < 0 || partition_number >= (INT_TYPE)m_partition_bounds.size()) {
				throw std::out_of_range("partition_number out of range");
			}
			return m_partition_bounds[(size_t)partition_number];
		}

		INT_TYPE coords_to_partition_num(const Vec2l& coords) const {
			const INT_TYPE tx = m_x_to_tile_col[(size_t)coords[0]];
			const INT_TYPE ty = m_y_to_tile_row[(size_t)coords[1]];
			return (INT_TYPE)(ty * m_tile_cols + tx);
		}

		INT_TYPE cell_id_to_partition_num(INDEX_TYPE cellid) const {
			Vec2l coords;
			cellid2Coords(cellid, coords);
			return coords_to_partition_num(coords);
		}

		bool is_partition_boundary_coords(const Vec2l& coords, INT_TYPE partition_number) const {
			const PartitionBounds& b = partition_bounds(partition_number);
			const INDEX_TYPE end_x = b.start_x + b.x_range;
			const INDEX_TYPE end_y = b.start_y + b.y_range;
			return coords[0] == b.start_x ||
				coords[0] == end_x - 1 ||
				coords[1] == b.start_y ||
				coords[1] == end_y - 1;
		}

		bool is_partition_boundary_cell(INDEX_TYPE cellid, INT_TYPE partition_number) const {
			Vec2l coords;
			cellid2Coords(cellid, coords);
			return is_partition_boundary_coords(coords, partition_number);
		}

		bool is_partition_interior_cell(INDEX_TYPE cellid, INT_TYPE partition_number) const {
			return !is_partition_boundary_cell(cellid, partition_number);
		}
	};

	class PartitionedTopologicalRegularGrid3D {
	public:
		struct PartitionBounds {
			INT_TYPE partition_id;
			INDEX_TYPE start_x;
			INDEX_TYPE start_y;
			INDEX_TYPE start_z;
			INDEX_TYPE x_range;
			INDEX_TYPE y_range;
			INDEX_TYPE z_range;
			INDEX_TYPE start_id;
		};

		class PartitionCellsIterator {
		protected:
			const PartitionedTopologicalRegularGrid3D* m_partitioned_mesh;
			PartitionBounds m_bounds;
			INDEX_TYPE m_curr_x;
			INDEX_TYPE m_curr_y;
			INDEX_TYPE m_curr_z;
			INDEX_TYPE m_curr_id;
			bool m_valid;

		public:
			PartitionCellsIterator(const PartitionedTopologicalRegularGrid3D* mesh, INT_TYPE partition_number) :
				m_partitioned_mesh(mesh),
				m_bounds(mesh->partition_bounds(partition_number)),
				m_curr_x(0),
				m_curr_y(0),
				m_curr_z(0),
				m_curr_id(0),
				m_valid(false) {}

			void begin() {
				m_curr_x = m_bounds.start_x;
				m_curr_y = m_bounds.start_y;
				m_curr_z = m_bounds.start_z;
				if (m_bounds.x_range == 0 || m_bounds.y_range == 0 || m_bounds.z_range == 0) {
					m_valid = false;
					m_curr_id = m_bounds.start_id;
					return;
				}
				m_curr_id = m_bounds.start_id;
				m_valid = true;
			}

			void advance() {
				if (!m_valid) return;
				const INDEX_TYPE end_x = m_bounds.start_x + m_bounds.x_range;
				const INDEX_TYPE end_y = m_bounds.start_y + m_bounds.y_range;
				const INDEX_TYPE end_z = m_bounds.start_z + m_bounds.z_range;
				m_curr_x++;
				if (m_curr_x < end_x) {
					m_curr_id++;
					return;
				}

				m_curr_x = m_bounds.start_x;
				m_curr_y++;
				if (m_curr_y < end_y) {
					m_curr_id = m_curr_x +
						m_curr_y * m_partitioned_mesh->numCellsAxis(0) +
						m_curr_z * m_partitioned_mesh->numCellsAxis(0) * m_partitioned_mesh->numCellsAxis(1);
					return;
				}

				m_curr_y = m_bounds.start_y;
				m_curr_z++;
				if (m_curr_z >= end_z) {
					m_valid = false;
					return;
				}
				m_curr_id = m_curr_x +
					m_curr_y * m_partitioned_mesh->numCellsAxis(0) +
					m_curr_z * m_partitioned_mesh->numCellsAxis(0) * m_partitioned_mesh->numCellsAxis(1);
			}

			bool valid() const {
				return m_valid;
			}

			INDEX_TYPE value() const {
				return m_curr_id;
			}

			INDEX_TYPE start_id() const { return m_bounds.start_id; }
			INDEX_TYPE start_x() const { return m_bounds.start_x; }
			INDEX_TYPE start_y() const { return m_bounds.start_y; }
			INDEX_TYPE start_z() const { return m_bounds.start_z; }
			INDEX_TYPE x_range() const { return m_bounds.x_range; }
			INDEX_TYPE y_range() const { return m_bounds.y_range; }
			INDEX_TYPE z_range() const { return m_bounds.z_range; }
		};

	protected:
		const TopologicalRegularGrid3D* m_mesh;
		INT_TYPE m_num_partitions;
		INDEX_TYPE m_tile_rows;
		INDEX_TYPE m_tile_cols;
		INDEX_TYPE m_tile_layers;
		Vec3l m_partition_splits;
		std::vector<PartitionBounds> m_partition_bounds;
		std::vector<INT_TYPE> m_x_to_tile_col;
		std::vector<INT_TYPE> m_y_to_tile_row;
		std::vector<INT_TYPE> m_z_to_tile_layer;

		static void make_axis_partitions(INDEX_TYPE axis_size, INDEX_TYPE num_tiles,
			std::vector<INDEX_TYPE>& starts, std::vector<INDEX_TYPE>& ranges) {
			starts.assign(num_tiles, 0);
			ranges.assign(num_tiles, 0);
			const INDEX_TYPE base = axis_size / num_tiles;
			const INDEX_TYPE rem = axis_size % num_tiles;
			INDEX_TYPE cursor = 0;
			for (INDEX_TYPE i = 0; i < num_tiles; i++) {
				const INDEX_TYPE width = base + (i < rem ? 1 : 0);
				starts[i] = cursor;
				ranges[i] = width;
				cursor += width;
			}
		}

		static double split_ratio_score(
			INDEX_TYPE sx, INDEX_TYPE sy, INDEX_TYPE sz,
			INDEX_TYPE xdim, INDEX_TYPE ydim, INDEX_TYPE zdim) {
			const double rx = static_cast<double>(sx) / static_cast<double>(xdim);
			const double ry = static_cast<double>(sy) / static_cast<double>(ydim);
			const double rz = static_cast<double>(sz) / static_cast<double>(zdim);
			const double dxy = rx - ry;
			const double dxz = rx - rz;
			const double dyz = ry - rz;
			return dxy * dxy + dxz * dxz + dyz * dyz;
		}

		static Vec3l infer_splits_for_requested_partitions(
			INDEX_TYPE xdim, INDEX_TYPE ydim, INDEX_TYPE zdim, INT_TYPE requested_partitions) {
			if (requested_partitions <= 0) {
				throw std::invalid_argument("requested_partitions must be positive.");
			}
			const INDEX_TYPE p = static_cast<INDEX_TYPE>(requested_partitions);

			bool found = false;
			Vec3l best(1, 1, 1);
			double best_score = std::numeric_limits<double>::infinity();
			double best_max_ratio = std::numeric_limits<double>::infinity();
			INDEX_TYPE best_max_split = std::numeric_limits<INDEX_TYPE>::max();

			for (INDEX_TYPE sx = 1; sx <= p; sx++) {
				if (p % sx != 0) continue;
				const INDEX_TYPE rem_xy = p / sx;
				for (INDEX_TYPE sy = 1; sy <= rem_xy; sy++) {
					if (rem_xy % sy != 0) continue;
					const INDEX_TYPE sz = rem_xy / sy;
					if (sx > xdim || sy > ydim || sz > zdim) continue;

					const double score = split_ratio_score(sx, sy, sz, xdim, ydim, zdim);
					const double rx = static_cast<double>(sx) / static_cast<double>(xdim);
					const double ry = static_cast<double>(sy) / static_cast<double>(ydim);
					const double rz = static_cast<double>(sz) / static_cast<double>(zdim);
					const double max_ratio = std::max(rx, std::max(ry, rz));
					const INDEX_TYPE max_split = std::max(sx, std::max(sy, sz));
					const bool better =
						(!found) ||
						(score < best_score - 1e-12) ||
						(std::abs(score - best_score) <= 1e-12 && max_ratio < best_max_ratio - 1e-12) ||
						(std::abs(score - best_score) <= 1e-12 && std::abs(max_ratio - best_max_ratio) <= 1e-12 && max_split < best_max_split) ||
						(std::abs(score - best_score) <= 1e-12 && std::abs(max_ratio - best_max_ratio) <= 1e-12 && max_split == best_max_split &&
						 (sx < best[0] || (sx == best[0] && (sy < best[1] || (sy == best[1] && sz < best[2])))));
					if (better) {
						found = true;
						best = Vec3l(sx, sy, sz);
						best_score = score;
						best_max_ratio = max_ratio;
						best_max_split = max_split;
					}
				}
			}

			if (!found) {
				throw std::invalid_argument("Unable to infer valid 3D partition splits for requested_partitions with current axis sizes.");
			}
			return best;
		}

		void initialize_mapping() {
			m_tile_cols = m_partition_splits[0];
			m_tile_rows = m_partition_splits[1];
			m_tile_layers = m_partition_splits[2];
			m_num_partitions = static_cast<INT_TYPE>(m_tile_cols * m_tile_rows * m_tile_layers);

			const INDEX_TYPE xdim = numCellsAxis(0);
			const INDEX_TYPE ydim = numCellsAxis(1);
			const INDEX_TYPE zdim = numCellsAxis(2);

			std::vector<INDEX_TYPE> x_starts, x_ranges, y_starts, y_ranges, z_starts, z_ranges;
			make_axis_partitions(xdim, m_tile_cols, x_starts, x_ranges);
			make_axis_partitions(ydim, m_tile_rows, y_starts, y_ranges);
			make_axis_partitions(zdim, m_tile_layers, z_starts, z_ranges);

			m_partition_bounds.clear();
			m_partition_bounds.reserve(static_cast<size_t>(m_num_partitions));
			for (INDEX_TYPE tz = 0; tz < m_tile_layers; tz++) {
				for (INDEX_TYPE ty = 0; ty < m_tile_rows; ty++) {
					for (INDEX_TYPE tx = 0; tx < m_tile_cols; tx++) {
						PartitionBounds b;
						b.partition_id = static_cast<INT_TYPE>((tz * m_tile_rows + ty) * m_tile_cols + tx);
						b.start_x = x_starts[tx];
						b.start_y = y_starts[ty];
						b.start_z = z_starts[tz];
						b.x_range = x_ranges[tx];
						b.y_range = y_ranges[ty];
						b.z_range = z_ranges[tz];
						b.start_id = b.start_x +
							b.start_y * xdim +
							b.start_z * xdim * ydim;
						m_partition_bounds.push_back(b);
					}
				}
			}

			m_x_to_tile_col.assign(static_cast<size_t>(xdim), 0);
			for (INDEX_TYPE tx = 0; tx < m_tile_cols; tx++) {
				const INDEX_TYPE start_x = x_starts[tx];
				const INDEX_TYPE end_x = start_x + x_ranges[tx];
				for (INDEX_TYPE x = start_x; x < end_x; x++) {
					m_x_to_tile_col[static_cast<size_t>(x)] = static_cast<INT_TYPE>(tx);
				}
			}
			m_y_to_tile_row.assign(static_cast<size_t>(ydim), 0);
			for (INDEX_TYPE ty = 0; ty < m_tile_rows; ty++) {
				const INDEX_TYPE start_y = y_starts[ty];
				const INDEX_TYPE end_y = start_y + y_ranges[ty];
				for (INDEX_TYPE y = start_y; y < end_y; y++) {
					m_y_to_tile_row[static_cast<size_t>(y)] = static_cast<INT_TYPE>(ty);
				}
			}
			m_z_to_tile_layer.assign(static_cast<size_t>(zdim), 0);
			for (INDEX_TYPE tz = 0; tz < m_tile_layers; tz++) {
				const INDEX_TYPE start_z = z_starts[tz];
				const INDEX_TYPE end_z = start_z + z_ranges[tz];
				for (INDEX_TYPE z = start_z; z < end_z; z++) {
					m_z_to_tile_layer[static_cast<size_t>(z)] = static_cast<INT_TYPE>(tz);
				}
			}
		}

	public:
		PartitionedTopologicalRegularGrid3D(const TopologicalRegularGrid3D* mesh, const Vec3l& partition_splits) :
			m_mesh(mesh),
			m_num_partitions(0),
			m_tile_rows(0),
			m_tile_cols(0),
			m_tile_layers(0),
			m_partition_splits(partition_splits) {
			if (m_mesh == NULL) throw std::invalid_argument("PartitionedTopologicalRegularGrid3D requires a non-null mesh.");
			if (partition_splits[0] <= 0 || partition_splits[1] <= 0 || partition_splits[2] <= 0) {
				throw std::invalid_argument("3D partition_splits must be positive along each axis.");
			}
			if (partition_splits[0] > numCellsAxis(0) || partition_splits[1] > numCellsAxis(1) || partition_splits[2] > numCellsAxis(2)) {
				throw std::invalid_argument("3D partition_splits cannot exceed axis cell counts.");
			}
			initialize_mapping();
		}

		PartitionedTopologicalRegularGrid3D(const TopologicalRegularGrid3D* mesh, INT_TYPE requested_partitions) :
			m_mesh(mesh),
			m_num_partitions(0),
			m_tile_rows(0),
			m_tile_cols(0),
			m_tile_layers(0),
			m_partition_splits(1, 1, 1) {
			if (m_mesh == NULL) throw std::invalid_argument("PartitionedTopologicalRegularGrid3D requires a non-null mesh.");
			m_partition_splits = infer_splits_for_requested_partitions(
				numCellsAxis(0), numCellsAxis(1), numCellsAxis(2), requested_partitions);
			initialize_mapping();
		}

		static Vec3l InferPartitionSplits(const TopologicalRegularGrid3D* mesh, INT_TYPE requested_partitions) {
			if (mesh == NULL) throw std::invalid_argument("InferPartitionSplits requires a non-null mesh.");
			return infer_splits_for_requested_partitions(
				mesh->numCellsAxis(0), mesh->numCellsAxis(1), mesh->numCellsAxis(2), requested_partitions);
		}

		INDEX_TYPE numCellsAxis(int axis) const { return m_mesh->numCellsAxis(axis); }
		INDEX_TYPE numCells() const { return m_mesh->numCells(); }
		INDEX_TYPE numCells(DIM_TYPE dim) const { return m_mesh->numCells(dim); }
		const Vec3l& XYZ() const { return m_mesh->XYZ(); }

		void cellid2Coords(INDEX_TYPE cellid, Vec3l& coords) const { m_mesh->cellid2Coords(cellid, coords); }
		INDEX_TYPE coords2Cellid(const Vec3l& coords) const { return m_mesh->coords2Cellid(coords); }
		DIM_TYPE dimension(const INDEX_TYPE& cellid) const { return m_mesh->dimension(cellid); }
		DIM_TYPE dimension(const Vec3l& coords) const { return m_mesh->dimension(coords); }
		BOUNDARY_TYPE boundaryValue(INDEX_TYPE cellid) const { return m_mesh->boundaryValue(cellid); }
		BOUNDARY_TYPE boundaryValue(const Vec3l& coords) const { return m_mesh->boundaryValue(coords); }

		INT_TYPE num_partitions() const { return m_num_partitions; }
		INDEX_TYPE tile_rows() const { return m_tile_rows; }
		INDEX_TYPE tile_cols() const { return m_tile_cols; }
		INDEX_TYPE tile_layers() const { return m_tile_layers; }
		const Vec3l& partition_splits() const { return m_partition_splits; }

		const PartitionBounds& partition_bounds(INT_TYPE partition_number) const {
			if (partition_number < 0 || partition_number >= static_cast<INT_TYPE>(m_partition_bounds.size())) {
				throw std::out_of_range("partition_number out of range");
			}
			return m_partition_bounds[static_cast<size_t>(partition_number)];
		}

		INT_TYPE coords_to_partition_num(const Vec3l& coords) const {
			const INT_TYPE tx = m_x_to_tile_col[static_cast<size_t>(coords[0])];
			const INT_TYPE ty = m_y_to_tile_row[static_cast<size_t>(coords[1])];
			const INT_TYPE tz = m_z_to_tile_layer[static_cast<size_t>(coords[2])];
			return static_cast<INT_TYPE>((tz * m_tile_rows + ty) * m_tile_cols + tx);
		}

		INT_TYPE cell_id_to_partition_num(INDEX_TYPE cellid) const {
			Vec3l coords;
			cellid2Coords(cellid, coords);
			return coords_to_partition_num(coords);
		}

		bool is_partition_boundary_coords(const Vec3l& coords, INT_TYPE partition_number) const {
			const PartitionBounds& b = partition_bounds(partition_number);
			const INDEX_TYPE end_x = b.start_x + b.x_range;
			const INDEX_TYPE end_y = b.start_y + b.y_range;
			const INDEX_TYPE end_z = b.start_z + b.z_range;
			return coords[0] == b.start_x ||
				coords[0] == end_x - 1 ||
				coords[1] == b.start_y ||
				coords[1] == end_y - 1 ||
				coords[2] == b.start_z ||
				coords[2] == end_z - 1;
		}

		bool is_partition_boundary_cell(INDEX_TYPE cellid, INT_TYPE partition_number) const {
			Vec3l coords;
			cellid2Coords(cellid, coords);
			return is_partition_boundary_coords(coords, partition_number);
		}

		bool is_partition_interior_cell(INDEX_TYPE cellid, INT_TYPE partition_number) const {
			return !is_partition_boundary_cell(cellid, partition_number);
		}
	};
}

#endif
