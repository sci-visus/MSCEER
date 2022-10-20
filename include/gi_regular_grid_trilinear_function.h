#ifndef REGULAR_GRID_TRILINEAR_FUNCTION
#define REGULAR_GRID_TRILINEAR_FUNCTION

#include "gi_vectors.h"
#include "gi_regular_grid.h"
#include "gi_labeling.h"
#include <algorithm>
#include <cmath>
#include <limits>


namespace GInt {
    // store image and gradient tied to a 3d grid.
    class RegularGridTrilinearFunction {
    protected:
        RegularGrid3D * m_grid;
        Vec3d* m_grad;

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

		void fill_extentsSerial() {
			FLOATTYPE m_max_value = m_image[0];
			FLOATTYPE m_min_value = m_image[0];
			INDEX_TYPE num_elements = m_grid->NumElements();
			for (INDEX_TYPE ii = 0; ii<num_elements; ++ii)
			{
				if (m_image[ii] > m_max_value)
				{

					m_max_value = m_image[ii];

				}
				if (m_image[ii] < m_min_value)
				{

					m_min_value = m_image[ii];
				}
			}
		}

    public:
      FLOATTYPE GetMinValue() const { return m_min_value; }
      FLOATTYPE GetMaxValue() const { return m_max_value; }
	  typedef float DType;

        RegularGridTrilinearFunction(RegularGrid3D* grid, FLOATTYPE *image = 0) : m_grid(grid) {
            m_i_made_image = false;
            m_i_made_gradient = false;
            m_image = NULL;
            m_grad = NULL;
            // use the function if it is passed, otherwise simply allocate memory
            if(image != 0) {   
				m_image = image;
				fill_extents();
				//printf("min = %e, max = %e\n", this->m_min_value, this->m_max_value);
			}

            //m_grad = new Vec3d[m_grid->NumElements()];
        }
	RegularGridTrilinearFunction() {};
	void RegularGridTrilinearFunctionSerial(RegularGrid3D* grid, FLOATTYPE *image = 0) {
		m_grid = grid;
		m_i_made_image = false;
		m_i_made_gradient = false;
		m_image = NULL;
		m_grad = NULL;
		// use the function if it is passed, otherwise simply allocate memory
		if(image != 0) {   
			m_image = image;
			fill_extentsSerial();
			//printf("min = %e, max = %e\n", this->m_min_value, this->m_max_value);
		}

		//m_grad = new Vec3d[m_grid->NumElements()];
	}
        ~RegularGridTrilinearFunction() {
            if (m_i_made_gradient) delete[] m_grad;
            if (m_i_made_image) delete[] m_image;
        }

        // return pointer to underlying mesh and function
        const RegularGrid3D* GetGrid() const {    return m_grid;  }
        FLOATTYPE* GetImage() const {         return m_image; }

        // sample the image at integral location
        FLOATTYPE SampleImage(const Vec3l& p) const {
            return m_image[m_grid->Index3d(p)];
        }

        // sample the image at integral location
        FLOATTYPE SampleImage(const INDEX_TYPE id) const {
            return m_image[id];
        }

        // sample the gradient at integral location
        const Vec3d& SampleGrad(const Vec3l& p) const {
            return m_grad[m_grid->Index3d(p)];
        }


        FLOATTYPE TriLinInterpValue(const Vec3d& s) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            Vec3d b = n[0];
            //s.print_vf();
            //b.print_vf();
            Vec3d factors = s - b;

            FLOATTYPE x0 = (1 - factors[0]) * SampleImage(n[0])  +  SampleImage(n[1]) * factors[0];
            FLOATTYPE x1 = (1 - factors[0]) * SampleImage(n[2]) + SampleImage(n[3]) * factors[0];
            FLOATTYPE x2 = (1 - factors[0]) * SampleImage(n[4]) + SampleImage(n[5]) * factors[0];
            FLOATTYPE x3 = (1 - factors[0]) * SampleImage(n[6]) + SampleImage(n[7]) * factors[0];

            FLOATTYPE y0 = (1 - factors[1]) *x0 + x1 * factors[1];
            FLOATTYPE y1 = (1 - factors[1]) *x2 + x3 * factors[1];

            return (1 - factors[2]) *y0 + y1 * factors[2];
        }

        // return trilinearly interpolated value
        Vec3d TriLinInterpGrad(const Vec3d& s) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            Vec3d b = n[0];
            //s.print_vf();
            //b.print_vf();
            Vec3d factors = s - b;

            Vec3d x0 = Vec3d::Lerp(SampleGrad(n[0]), SampleGrad(n[1]), factors[0]);
            Vec3d x1 = Vec3d::Lerp(SampleGrad(n[2]), SampleGrad(n[3]), factors[0]);
            Vec3d x2 = Vec3d::Lerp(SampleGrad(n[4]), SampleGrad(n[5]), factors[0]);
            Vec3d x3 = Vec3d::Lerp(SampleGrad(n[6]), SampleGrad(n[7]), factors[0]);

            Vec3d y0 = Vec3d::Lerp(x0, x1, factors[1]);
            Vec3d y1 = Vec3d::Lerp(x2, x3, factors[1]);

            return Vec3d::Lerp(y0, y1, factors[2]);
        }

        void SetGradExplicit(INDEX_TYPE id, Vec3d vec) {
            this->m_grad[id] = vec;
        }

        // fill in vals with the 8 values of hte gradient around sample poitn
        void GetGradSurrounding(const Vec3d& s, Vec3d* vals) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
        }
        void GetGradSurrounding(const Vec3l& s, Vec3d* vals) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
        }

        // use with extreme care - no boundary checks, only do on really interior poitns
        void GetGradSurroundingNoBoundaryCheck(const Vec3d& s, Vec3d* vals) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurroundingNoBoundaryCheck(s, n);
            for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
        }

        FLOATTYPE InterpolatedValue(const Vec3d& s) const {
            return TriLinInterpValue(s);
        }

        Vec3d InterpolatedGrad(const Vec3d& s) const {
            return TriLinInterpGrad(s);
        }

        // allow reuse of sampled gradient - the assumption that vals has the gradient arrows around s
        Vec3d TriLinInterpGrad(const Vec3d& s, const Vec3l& int_base, Vec3d* vals) const {

            //if (!(s.IntFloor() == int_base)) {
            //	printf("s=");  s.PrintFloat(); printf("d="); int_base.PrintFloat();
            //}
            //
            //Vec3d d = int_base.IntFloor();
            Vec3d factors = s - int_base;

            Vec3d x0 = Vec3d::Lerp(vals[0], vals[1], factors[0]);
            Vec3d x1 = Vec3d::Lerp(vals[2], vals[3], factors[0]);
            Vec3d x2 = Vec3d::Lerp(vals[4], vals[5], factors[0]);
            Vec3d x3 = Vec3d::Lerp(vals[6], vals[7], factors[0]);

            Vec3d y0 = Vec3d::Lerp(x0, x1, factors[1]);
            Vec3d y1 = Vec3d::Lerp(x2, x3, factors[1]);

            return Vec3d::Lerp(y0, y1, factors[2]);
        }
        void LoadImageFromFloatFile(const char* fname) {

            size_t image_size = m_grid->NumElements();

            // fill in image
            m_image = new FLOATTYPE[image_size];    m_i_made_image = true;


            FILE* fin = fopen(fname, "rb");
			fread(m_image, sizeof(float), image_size, fin);
		

            fclose(fin);
			fill_extents();
            printf("min = %e, max = %e\n", this->m_min_value, this->m_max_value);
        }
        void LoadImageFromFile(const char* fname) {

            size_t image_size = m_grid->NumElements();

            // fill in image
            m_image = new FLOATTYPE[image_size]; m_i_made_image = true;

            FILE* fin = fopen(fname, "rb");
            fread(m_image, sizeof(FLOATTYPE), image_size, fin);
            fclose(fin);

			fill_extents();
			printf("loaded %s...\n", fname);
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

        Vec3d GradientFromImage(const Vec3l& p, int rklevel) {
            Vec3l negs[9]; // don't support more than 4th order - cmon. would be ridiculous

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
            double res_z = 0.0;
			int rklevel_z = m_grid->Gather1DNeighborhood(p, 2, rklevel, negs);		
			int nume_z = rklevel_z * 2 + 1; // number of entries to average
			for (int i = 0; i < nume_z; i++) {
				res_z += kRKCoefficients[rklevel_z][i] * SampleImage(negs[i]);
            }
            return Vec3d(res_x, res_y, res_z);
        }

        inline bool IsGreater(INDEX_TYPE a, INDEX_TYPE b) const {
            if (m_image[a] > m_image[b]) return true;
            if (m_image[b] > m_image[a]) return false;
            //if (a == b) printf("WHOA THERE NELLY\n");
            return a > b;
        }
        //Vec3d IStep(const Vec3d& p, const Vec3d& grad, const FLOATTYPE h) const {
        //	return m_grid->Inbounds(p + (grad * h));
        //}
        //Vec3d IStepNoBoundaryCheck(const Vec3d& p, const Vec3d& grad, const FLOATTYPE h) const {
        //	return p + (grad * h);
        //}


        // add in block structure




        void ComputeGradFromImage(int rklevel) {
            m_grad = new Vec3d[m_grid->NumElements()];
            m_i_made_gradient = true;
#pragma omp parallel for
            for (int i = 0; i < m_grid->XYZ()[0]; i++) {
                for (int j = 0; j < m_grid->XYZ()[1]; j++) {
                    for (int k = 0; k < m_grid->XYZ()[2]; k++) {
                        Vec3l p(i, j, k);
                        m_grad[m_grid->Index3d(p)] = GradientFromImage(p, rklevel);
                    }
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


	    class UncachedRegularGridTrilinearFunction {
    protected:
        RegularGrid3D * m_grid;


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
      
	  UncachedRegularGridTrilinearFunction(RegularGrid3D* grid, FLOATTYPE *image = 0) : m_grid(grid) {
            m_i_made_image = false;
            m_i_made_gradient = false;
            m_image = NULL;

            // use the function if it is passed, otherwise simply allocate memory
            if(image != 0) {    m_image = image;                                }

            //m_grad = new Vec3d[m_grid->NumElements()];
        }
        ~UncachedRegularGridTrilinearFunction() {
            if (m_i_made_image) delete[] m_image;
        }

        // return pointer to underlying mesh and function
        const RegularGrid3D* GetGrid() const {    return m_grid;  }
        FLOATTYPE* GetImage() const {         return m_image; }

        // sample the image at integral location
        FLOATTYPE SampleImage(const Vec3l& p) const {
            return m_image[m_grid->Index3d(p)];
        }

        // sample the image at integral location
        FLOATTYPE SampleImage(const INDEX_TYPE id) const {
            return m_image[id];
        }
		static const FLOATTYPE kRKCoefficients[5][9];

		Vec3d GradientFromImage(const Vec3l& p, int rklevel) const {
			Vec3l negs[9]; // don't support more than 4th order - cmon. would be ridiculous

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
			double res_z = 0.0;
			int rklevel_z = m_grid->Gather1DNeighborhood(p, 2, rklevel, negs);
			int nume_z = rklevel_z * 2 + 1; // number of entries to average
			for (int i = 0; i < nume_z; i++) {
				res_z += kRKCoefficients[rklevel_z][i] * SampleImage(negs[i]);
			}
			return Vec3d(res_x, res_y, res_z);
		}
		// sample the gradient at integral location
		const Vec3d SampleGrad(const Vec3l& p) const {
			return GradientFromImage(p, 1);
		}


        FLOATTYPE TriLinInterpValue(const Vec3d& s) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            Vec3d b = n[0];
            //s.print_vf();
            //b.print_vf();
            Vec3d factors = s - b;

            FLOATTYPE x0 = (1 - factors[0]) * SampleImage(n[0])  +  SampleImage(n[1]) * factors[0];
            FLOATTYPE x1 = (1 - factors[0]) * SampleImage(n[2]) + SampleImage(n[3]) * factors[0];
            FLOATTYPE x2 = (1 - factors[0]) * SampleImage(n[4]) + SampleImage(n[5]) * factors[0];
            FLOATTYPE x3 = (1 - factors[0]) * SampleImage(n[6]) + SampleImage(n[7]) * factors[0];

            FLOATTYPE y0 = (1 - factors[1]) *x0 + x1 * factors[1];
            FLOATTYPE y1 = (1 - factors[1]) *x2 + x3 * factors[1];

            return (1 - factors[2]) *y0 + y1 * factors[2];
        }

        // return trilinearly interpolated value
        Vec3d TriLinInterpGrad(const Vec3d& s) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            Vec3d b = n[0];
            //s.print_vf();
            //b.print_vf();
            Vec3d factors = s - b;

            Vec3d x0 = Vec3d::Lerp(SampleGrad(n[0]), SampleGrad(n[1]), factors[0]);
            Vec3d x1 = Vec3d::Lerp(SampleGrad(n[2]), SampleGrad(n[3]), factors[0]);
            Vec3d x2 = Vec3d::Lerp(SampleGrad(n[4]), SampleGrad(n[5]), factors[0]);
            Vec3d x3 = Vec3d::Lerp(SampleGrad(n[6]), SampleGrad(n[7]), factors[0]);

            Vec3d y0 = Vec3d::Lerp(x0, x1, factors[1]);
            Vec3d y1 = Vec3d::Lerp(x2, x3, factors[1]);

            return Vec3d::Lerp(y0, y1, factors[2]);
        }

        void SetGradExplicit(INDEX_TYPE id, Vec3d vec) {
            //this->m_grad[id] = vec;
        }

        // fill in vals with the 8 values of hte gradient around sample poitn
        void GetGradSurrounding(const Vec3d& s, Vec3d* vals) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
        }
        void GetGradSurrounding(const Vec3l& s, Vec3d* vals) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurrounding(s, n);
            for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
        }

        // use with extreme care - no boundary checks, only do on really interior poitns
        void GetGradSurroundingNoBoundaryCheck(const Vec3d& s, Vec3d* vals) const {
            Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
            m_grid->GatherSurroundingNoBoundaryCheck(s, n);
            for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
        }

        FLOATTYPE InterpolatedValue(const Vec3d& s) const {
            return TriLinInterpValue(s);
        }

        Vec3d InterpolatedGrad(const Vec3d& s) const {
            return TriLinInterpGrad(s);
        }

        // allow reuse of sampled gradient - the assumption that vals has the gradient arrows around s
        Vec3d TriLinInterpGrad(const Vec3d& s, const Vec3l& int_base, Vec3d* vals) const {

            //if (!(s.IntFloor() == int_base)) {
            //	printf("s=");  s.PrintFloat(); printf("d="); int_base.PrintFloat();
            //}
            //
            //Vec3d d = int_base.IntFloor();
            Vec3d factors = s - int_base;

            Vec3d x0 = Vec3d::Lerp(vals[0], vals[1], factors[0]);
            Vec3d x1 = Vec3d::Lerp(vals[2], vals[3], factors[0]);
            Vec3d x2 = Vec3d::Lerp(vals[4], vals[5], factors[0]);
            Vec3d x3 = Vec3d::Lerp(vals[6], vals[7], factors[0]);

            Vec3d y0 = Vec3d::Lerp(x0, x1, factors[1]);
            Vec3d y1 = Vec3d::Lerp(x2, x3, factors[1]);

            return Vec3d::Lerp(y0, y1, factors[2]);
        }
        void LoadImageFromFloatFile(const char* fname) {

            size_t image_size = m_grid->NumElements();

            // fill in image
            m_image = new FLOATTYPE[image_size];    m_i_made_image = true;


            FILE* fin = fopen(fname, "rb");
            for (size_t i = 0; i < image_size; i++) {

                float tval = 0;
                fread(&tval, sizeof(float), 1, fin);
                m_image[i] = tval;
            }

            fclose(fin);
			fill_extents();
            printf("min = %e, max = %e\n", this->m_min_value, this->m_max_value);
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

 

        inline bool IsGreater(INDEX_TYPE a, INDEX_TYPE b) const {
            if (m_image[a] > m_image[b]) return true;
            if (m_image[b] > m_image[a]) return false;
            //if (a == b) printf("WHOA THERE NELLY\n");
            return a > b;
        }
        //Vec3d IStep(const Vec3d& p, const Vec3d& grad, const FLOATTYPE h) const {
        //	return m_grid->Inbounds(p + (grad * h));
        //}
        //Vec3d IStepNoBoundaryCheck(const Vec3d& p, const Vec3d& grad, const FLOATTYPE h) const {
        //	return p + (grad * h);
        //}


        // add in block structure




        void ComputeGradFromImage(int rklevel) {


        }



        void Negate() {

#pragma omp parallel for schedule(static)
                for (INDEX_TYPE i = 0; i < m_grid->NumElements(); i++) {
                    this->m_image[i] *= -1;
                }
            }
        

    };


#if 1

#endif



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
			if (image != 0) {
				m_image = image; 
				fill_extents();
				printf("min = %e, max = %e\n", this->m_min_value, this->m_max_value);
			}

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


	class SparseRegularGridTrilinearFunction {
	protected:
		RegularGrid3D* m_grid;
		SparseLabeling<Vec3d> m_grad;
		SparseLabeling<FLOATTYPE> m_image;

		FLOATTYPE m_min_value;
		FLOATTYPE m_max_value;

		bool m_i_made_gradient;
		bool m_i_made_image;

		void fill_extents() {
			FLOATTYPE t_max_val = m_max_value = (std::numeric_limits<FLOATTYPE>::min)(); // wrapped with parentheses to avoid macro from Windows.h
			FLOATTYPE t_min_val = m_min_value = (std::numeric_limits<FLOATTYPE>::max)();
			INDEX_TYPE num_elements = m_grid->NumElements();
			INDEX_TYPE ii;

#pragma omp parallel shared(num_elements) private(ii) firstprivate(t_max_val,t_min_val)
			{
#pragma omp for nowait
				for (ii = 0; ii<num_elements; ++ii)
				{
					if (!m_image.Has(ii)) continue;
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

		SparseLabeling<FLOATTYPE>& GetSparseMap() { return m_image; }
		FLOATTYPE GetMinValue() const { return m_min_value; }
		FLOATTYPE GetMaxValue() const { return m_max_value; }
		typedef float DType;

		SparseRegularGridTrilinearFunction(RegularGrid3D* grid) : m_grid(grid) {
			m_i_made_image = false;
			m_i_made_gradient = false;
			//m_grad = new Vec3d[m_grid->NumElements()];
		}
		~SparseRegularGridTrilinearFunction() {
		}

		// return pointer to underlying mesh and function
		const RegularGrid3D* GetGrid() const { return m_grid; }
		//FLOATTYPE* GetImage() const { return m_image; }

		// sample the image at integral location
		FLOATTYPE SampleImage(const Vec3l& p) const {
			return m_image[m_grid->Index3d(p)];
		}

		// sample the image at integral location
		FLOATTYPE SampleImage(const INDEX_TYPE id) const {
			return m_image[id];
		}

		// sample the gradient at integral location
		const Vec3d& SampleGrad(const Vec3l& p) const {
			return m_grad[m_grid->Index3d(p)];
		}


		FLOATTYPE TriLinInterpValue(const Vec3d& s) const {
			Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
			int res = m_grid->GatherSurrounding(s, n);
			for (int i = 0; i < res; i++) {
				if (!m_image.Has(m_grid->Index3d(n[i]))) return 0;
			}
			Vec3d b = n[0];
			//s.print_vf();
			//b.print_vf();
			Vec3d factors = s - b;

			FLOATTYPE x0 = (1 - factors[0]) * SampleImage(n[0]) + SampleImage(n[1]) * factors[0];
			FLOATTYPE x1 = (1 - factors[0]) * SampleImage(n[2]) + SampleImage(n[3]) * factors[0];
			FLOATTYPE x2 = (1 - factors[0]) * SampleImage(n[4]) + SampleImage(n[5]) * factors[0];
			FLOATTYPE x3 = (1 - factors[0]) * SampleImage(n[6]) + SampleImage(n[7]) * factors[0];

			FLOATTYPE y0 = (1 - factors[1]) *x0 + x1 * factors[1];
			FLOATTYPE y1 = (1 - factors[1]) *x2 + x3 * factors[1];

			return (1 - factors[2]) *y0 + y1 * factors[2];
		}

		// return trilinearly interpolated value
		Vec3d TriLinInterpGrad(const Vec3d& s) const {
			Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
			m_grid->GatherSurrounding(s, n);
			Vec3d b = n[0];
			//s.print_vf();
			//b.print_vf();
			Vec3d factors = s - b;

			Vec3d x0 = Vec3d::Lerp(SampleGrad(n[0]), SampleGrad(n[1]), factors[0]);
			Vec3d x1 = Vec3d::Lerp(SampleGrad(n[2]), SampleGrad(n[3]), factors[0]);
			Vec3d x2 = Vec3d::Lerp(SampleGrad(n[4]), SampleGrad(n[5]), factors[0]);
			Vec3d x3 = Vec3d::Lerp(SampleGrad(n[6]), SampleGrad(n[7]), factors[0]);

			Vec3d y0 = Vec3d::Lerp(x0, x1, factors[1]);
			Vec3d y1 = Vec3d::Lerp(x2, x3, factors[1]);

			return Vec3d::Lerp(y0, y1, factors[2]);
		}

		void SetGradExplicit(INDEX_TYPE id, Vec3d vec) {
			this->m_grad[id] = vec;
		}

		// fill in vals with the 8 values of hte gradient around sample poitn
		void GetGradSurrounding(const Vec3d& s, Vec3d* vals) const {
			Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
			m_grid->GatherSurrounding(s, n);
			for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
		}
		void GetGradSurrounding(const Vec3l& s, Vec3d* vals) const {
			Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
			m_grid->GatherSurrounding(s, n);
			for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
		}

		// use with extreme care - no boundary checks, only do on really interior poitns
		void GetGradSurroundingNoBoundaryCheck(const Vec3d& s, Vec3d* vals) const {
			Vec3l n[8]; // for 8 vertices around s - some may be repeated based on boundary cond.
			m_grid->GatherSurroundingNoBoundaryCheck(s, n);
			for (int i = 0; i < 8; i++) vals[i] = SampleGrad(n[i]);
		}

		FLOATTYPE InterpolatedValue(const Vec3d& s) const {
			return TriLinInterpValue(s);
		}

		Vec3d InterpolatedGrad(const Vec3d& s) const {
			return TriLinInterpGrad(s);
		}

		// allow reuse of sampled gradient - the assumption that vals has the gradient arrows around s
		Vec3d TriLinInterpGrad(const Vec3d& s, const Vec3l& int_base, Vec3d* vals) const {

			//if (!(s.IntFloor() == int_base)) {
			//	printf("s=");  s.PrintFloat(); printf("d="); int_base.PrintFloat();
			//}
			//
			//Vec3d d = int_base.IntFloor();
			Vec3d factors = s - int_base;

			Vec3d x0 = Vec3d::Lerp(vals[0], vals[1], factors[0]);
			Vec3d x1 = Vec3d::Lerp(vals[2], vals[3], factors[0]);
			Vec3d x2 = Vec3d::Lerp(vals[4], vals[5], factors[0]);
			Vec3d x3 = Vec3d::Lerp(vals[6], vals[7], factors[0]);

			Vec3d y0 = Vec3d::Lerp(x0, x1, factors[1]);
			Vec3d y1 = Vec3d::Lerp(x2, x3, factors[1]);

			return Vec3d::Lerp(y0, y1, factors[2]);
		}
		void LoadImageFromFloatAndMaskFile(const char* fname, const char* fmask) {

			size_t image_size = m_grid->NumElements();

			FILE* fin = fopen(fname, "rb");
			FILE* fma = fopen(fmask, "rb");
			auto dims = m_grid->XYZ();
			float* x_line_f = new float[dims[0]];
			BYTE_TYPE* x_line_m = new BYTE_TYPE[dims[0]];
			INDEX_TYPE counter = 0;
			for (INDEX_TYPE yz = 0; yz < dims[1] * dims[2]; yz++) {
				fread(x_line_f, sizeof(float), dims[0], fin);
				fread(x_line_m, sizeof(BYTE_TYPE), dims[0], fma);
				for (INDEX_TYPE x = 0; x < dims[0]; x++) {
					if (x_line_m[x] > 0) {
						m_image.SetLabel(counter, x_line_f[x]);
					}
					counter++;
				}
			}

			fclose(fin);
			fclose(fma);
			fill_extents();
			printf("min = %e, max = %e\n", this->m_min_value, this->m_max_value);
		}

		static const FLOATTYPE kRKCoefficients[5][9];


		inline bool IsGreater(INDEX_TYPE a, INDEX_TYPE b) const {
			if (m_image[a] > m_image[b]) return true;
			if (m_image[b] > m_image[a]) return false;
			//if (a == b) printf("WHOA THERE NELLY\n");
			return a > b;
		}



	}; // end SparseRegularGridTrilinearFunction



}; // end namespace





#endif
