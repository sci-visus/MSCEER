#ifndef PY_GRAD_BUILDER_H
#define PY_GRAD_BUILDER_H



#include <vector>
#include <set>
#include <queue>
#include <time.h>

#include "base/gi_timing.h"
#include "base/gi_topological_regular_grid_3d.h"
#include "base/gi_isolated_region_remover.h"
#include "base/gi_isolated_region_remover_masked.h"
#include "base/gi_numeric_integrator_path_compressing.h"
#include "base/gi_numeric_streamline_integrator_digitizing.h"
#include "base/gi_timing.h"
#include "base/gi_adaptive_euler_advector_2d.h"
#include "base/gi_adaptive_euler_advector_3d.h"
#include "base/gi_advection_checkers.h"
#include "base/gi_advection_events.h"
#include "base/gi_index_comparer.h"
#include "base/gi_maxmin_vertex_labeling.h"
#include "base/gi_conforming_discrete_gradient.h"
#include "base/gi_robins_sliding_regular_grid.h"
#include "base/gi_labeling_to_bounary_labeling.h"
#include "base/gi_topological_gradient_using_algorithms.h"
#include "base/gi_topological_gradient_using_algorithms.h"
#include "base/gi_isolated_region_remover.h"
#include "base/gi_bifiltration_pairing.h"
#include "base/gi_topological_max_vertex_mesh_function.h"
#include "base/gi_extrema_region_builder.h"
#include "base/gi_numeric_integrator_path_compressing.h"


class DiscreteGradientBuilder {
public:

	enum ComputeMode {
		STEEPEST_ROBINS,
		STEEPEST_GYU,
		CONFORMING,
		CONVERGENT,
		ONDEMAND_ACCURATE
	};

protected:
	int m_native_dimension;
	int m_data_x;
	int m_data_y;
	int m_data_z;
	int m_periodic_x;
	int m_periodic_y;
	int m_periodic_z;
	ComputeMode m_mode;

	int m_source_mode; // -1 = not set, 1 = file, 2 = array pointer
	bool m_has_filename;
	std::string m_raw_filename;
	float* m_raw_array;

public:
	DiscreteGradientBuilder() {
		m_mode = STEEPEST_ROBINS;
		m_native_dimension = -1;
		m_source_mode = -1;
		m_has_filename = false;
		m_raw_array = NULL;
	}

	void SetDataDims(int x, int y, int z) {
		m_data_x = x; m_data_y = y; m_data_z = z;
		m_native_dimension = 3;
	}

	void SetDataDims(int x, int y) {
		m_data_x = x; m_data_y = y;
		m_native_dimension = 2;
	}

	void SetPeriodicity(int x, int y, int z) {
		m_periodic_x = x; m_periodic_y = y; m_periodic_z = z;
	}

	void SetPeriodicity(int x, int y) {
		m_periodic_x = x; m_periodic_y = y;
	}

	void SetInputRawFileName(std::string name) {
		m_raw_filename = name;
		m_source_mode = 1;
		m_has_filename = true;
	}

	void SetInputRawArray(float* vals) {
		m_raw_array = vals;
		m_source_mode = 2;
	}

	void SetComputeMode(ComputeMode mode) {
		m_mode = mode;
	}

	bool CheckValidSettings() {
		bool valid = true;
		if (m_native_dimension == -1) {
			printf("GradientBuilder: warning => data dimensions not set - use SetDataDims()\n");
			valid = false;
		}
		if (m_source_mode == -1) {
			printf("GradientBuilder: warning => no source set for data - use SetInputRawFileName() or SetInputRawArray()\n");
			valid = false;
		}
		return valid;
	}

	void ComputeGradient() {




	}
}; // class DiscreteGradientBuilder



#endif