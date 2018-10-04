#ifndef REGULAR_GRID_BILINEAR_FUNCTION
#define REGULAR_GRID_BILINEAR_FUNCTION

#include <algorithm>
#include <cmath>

#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "base/gi_regular_grid_2d.h"


namespace GInt {

	class RegularGridBilinearFunction {
	protected:
		RegularGrid2D * m_grid;
		Vec2d* m_grad;

		FLOATTYPE* m_image;

		FLOATTYPE m_min_value;
		FLOATTYPE m_max_value;

		bool m_i_made_gradient;
		bool m_i_made_image;

		void fill_extents() {
			FLOATTYPE t_max_val = m_max_value = m_image[0];
			FLOATTYPE t_min_val = m_min_value = m_image[0];
			INDEX_TYPE num_elements = m_grid->NumElements();
			INDEX_TYPE ii;

#pragma omp parallel shared(num_elements) private(ii) firstprivate(t_max_val,t_min_val)
			{
#pragma omp for nowait
				for (ii = 0; ii<num_elements; ++ii)
				{
					if (m_image[ii] > t_max_val)
					{

						t_max_val = m_image[ii];

					}
					if (m_image[ii] < t_min_val)
					{

						t_min_val = m_image[ii];
					}
				}
#pragma omp critical 
				{
					if (t_max_val > m_max_value) m_max_value = t_max_val;
					if (t_min_val < m_min_value) m_min_value = t_min_val;
				}
			}

		}

	public:
		FLOATTYPE GetMinValue() const { return m_min_value; }
		FLOATTYPE GetMaxValue() const { return m_max_value; }

		RegularGridBilinearFunction(RegularGrid2D* grid, FLOATTYPE *image = 0) : m_grid(grid) {
			m_i_made_image = false;
			m_i_made_gradient = false;
			m_image = NULL;
			m_grad = NULL;
			// use the function if it is passed, otherwise simply allocate memory
			if (image != 0) { m_image = image; }

			//m_grad = new Vec2d[m_grid->NumElements()];
		}
		~RegularGridBilinearFunction() {
			if (m_i_made_gradient) delete[] m_grad;
			if (m_i_made_image) delete[] m_image;
		}

		// return pointer to underlying mesh and function
		const RegularGrid2D* GetGrid() const { return m_grid; }
		FLOATTYPE* GetImage() const { return m_image; }

		// sample the image at integral location
		FLOATTYPE SampleImage(const Vec2l& p) const {
			return m_image[m_grid->Index2d(p)];
		}

		// sample the image at integral location
		FLOATTYPE SampleImage(const INDEX_TYPE id) const {
			return m_image[id];
		}

		// sample the gradient at integral location
		const Vec2d& SampleGrad(const Vec2l& p) const {
			return m_grad[m_grid->Index2d(p)];
		}


		FLOATTYPE BiLinInterpValue(const Vec2d& s) const {
			Vec2l n[4]; // for 8 vertices around s - some may be repeated based on boundary cond.
			m_grid->GatherSurrounding(s, n);
			Vec2d b = n[0];
			//s.print_vf();
			//b.print_vf();
			Vec2d factors = s - b;

			FLOATTYPE x0 = (1 - factors[0]) * SampleImage(n[0]) + SampleImage(n[1]) * factors[0];
			FLOATTYPE x1 = (1 - factors[0]) * SampleImage(n[2]) + SampleImage(n[3]) * factors[0];

			return (1 - factors[1]) *x0 + x1 * factors[1];
		}

		// return trilinearly interpolated value
		Vec2d BiLinInterpGrad(const Vec2d& s) const {
			Vec2l n[4]; // for 8 vertices around s - some may be repeated based on boundary cond.
			m_grid->GatherSurrounding(s, n);
			Vec2d b = n[0];
			//s.print_vf();
			//b.print_vf();
			Vec2d factors = s - b;

			Vec2d x0 = Vec2d::Lerp(SampleGrad(n[0]), SampleGrad(n[1]), factors[0]);
			Vec2d x1 = Vec2d::Lerp(SampleGrad(n[2]), SampleGrad(n[3]), factors[0]);
			
			return  Vec2d::Lerp(x0, x1, factors[1]);
		}

		void SetGradExplicit(INDEX_TYPE id, Vec2d vec) {
			this->m_grad[id] = vec;
		}

		// fill in vals with the 8 values of hte gradient around sample poitn
		void GetGradSurrounding(const Vec2d& s, Vec2d* vals) const {
			Vec2l n[4]; // for 8 vertices around s - some may be repeated based on boundary cond.
			m_grid->GatherSurrounding(s, n);
			for (int i = 0; i < 4; i++) vals[i] = SampleGrad(n[i]);
		}
		void GetGradSurrounding(const Vec2l& s, Vec2d* vals) const {
			Vec2l n[4]; // for 8 vertices around s - some may be repeated based on boundary cond.
			m_grid->GatherSurrounding(s, n);
			for (int i = 0; i < 4; i++) vals[i] = SampleGrad(n[i]);
		}

		// use with extreme care - no boundary checks, only do on really interior poitns
		void GetGradSurroundingNoBoundaryCheck(const Vec2d& s, Vec2d* vals) const {
			Vec2l n[4]; // for 8 vertices around s - some may be repeated based on boundary cond.
			m_grid->GatherSurroundingNoBoundaryCheck(s, n);
			for (int i = 0; i < 4; i++) vals[i] = SampleGrad(n[i]);
		}

		FLOATTYPE InterpolatedValue(const Vec2d& s) const {
			return BiLinInterpValue(s);
		}

		Vec2d InterpolatedGrad(const Vec2d& s) const {
			return BiLinInterpGrad(s);
		}

		// allow reuse of sampled gradient - the assumption that vals has the gradient arrows around s
		Vec2d BiLinInterpGrad(const Vec2d& s, const Vec2l& int_base, Vec2d* vals) const {

			//if (!(s.IntFloor() == int_base)) {
			//	printf("s=");  s.PrintFloat(); printf("d="); int_base.PrintFloat();
			//}
			//
			//Vec2d d = int_base.IntFloor();
			Vec2d factors = s - int_base;

			Vec2d x0 = Vec2d::Lerp(vals[0], vals[1], factors[0]);
			Vec2d x1 = Vec2d::Lerp(vals[2], vals[3], factors[0]);
			
			return Vec2d::Lerp(x0, x1, factors[1]);

		}

		void LoadImageFromFile(const char* fname) {

			size_t image_size = m_grid->NumElements();

			// fill in image
			m_image = new FLOATTYPE[image_size]; m_i_made_image = true;

			FILE* fin = fopen(fname, "rb");
			fread(m_image, sizeof(FLOATTYPE), image_size, fin);
			fclose(fin);

			fill_extents();
			printf("min = %e, max = %e\n", this->m_min_value, this->m_max_value);
		}

		void ShallowCopyImage(FLOATTYPE *image) {

			m_image = image;

			INDEX_TYPE image_size = m_grid->NumElements();
			fill_extents();
			printf("min = %e, max = %e\n", this->m_min_value, this->m_max_value);
		}

		void DeepCopyImage(const FLOATTYPE *image) {

			m_image = new FLOATTYPE[m_grid->NumElements()];  m_i_made_image = true;

			INDEX_TYPE image_size = m_grid->NumElements();
			memcpy(m_image, image, image_size*sizeof(FLOATTYPE));

			fill_extents();
			printf("min = %e, max = %e\n", this->m_min_value, this->m_max_value);

		}

		static const FLOATTYPE kRKCoefficients[5][9];

		Vec2d GradientFromImage(const Vec2l& p, int rklevel) {
			Vec2l negs[9]; // don't support more than 4th order - cmon. would be ridiculous

			double res_x = 0.0;
			int rklevel_x = m_grid->Gather1DNeighborhood(p, 0, rklevel, negs);
			int nume_x = rklevel_x * 2 + 1; // number of entries to average
			for (int i = 0; i < nume_x; i++) {
				res_x += kRKCoefficients[rklevel_x][i] * SampleImage(negs[i]);
			}
			double res_y = 0.0;
			int rklevel_y = m_grid->Gather1DNeighborhood(p, 1, rklevel, negs);
			int nume_y = rklevel_y * 2 + 1; // number of entries to average
			for (int i = 0; i < nume_y; i++) {
				res_y += kRKCoefficients[rklevel_y][i] * SampleImage(negs[i]);
			}

			return Vec2d(res_x, res_y);
		}

		inline bool IsGreater(INDEX_TYPE a, INDEX_TYPE b) const {
			if (m_image[a] > m_image[b]) return true;
			if (m_image[b] > m_image[a]) return false;
			//if (a == b) printf("WHOA THERE NELLY\n");
			return a > b;
		}

		void ComputeGradFromImage(int rklevel) {
			m_grad = new Vec2d[m_grid->NumElements()];
			m_i_made_gradient = true;
#pragma omp parallel for
			for (int i = 0; i < m_grid->XY()[0]; i++) {
				for (int j = 0; j < m_grid->XY()[1]; j++) {

					Vec2l p(i, j);
					m_grad[m_grid->Index2d(p)] = GradientFromImage(p, rklevel);

				}
			}

		}



		void Negate() {
			if (m_grad != NULL) {
#pragma omp parallel for schedule(static)
				for (INDEX_TYPE i = 0; i < m_grid->NumElements(); i++) {
					this->m_image[i] *= -1;
					this->m_grad[i] *= -1.0;
				}

			}
			else {
#pragma omp parallel for schedule(static)
				for (INDEX_TYPE i = 0; i < m_grid->NumElements(); i++) {
					this->m_image[i] *= -1;
				}
			}
		}

	};



};

#endif
