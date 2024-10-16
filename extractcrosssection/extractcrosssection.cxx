
#include <vector>
#include <map>
#include <string>
#include <mutex>
#include <functional>
#include <sstream>
#include <random>
#include <algorithm>

#define VTK_ENABLED
#ifdef VTK_ENABLED
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkNamedColors.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkMatrix4x4.h>
#include <vtkUnsignedCharArray.h>

#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>

#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataReader.h>
#include <vtkRenderer.h>
#include <vtkImageReslice.h>
#include <vtkTransform.h>

// for re-do of slice
#include <vtkActor.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>

#endif

// For identifying connected components
#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif

#include <fstream>

#include "gi_basic_types.h"
#include "gi_timing.h"
#include "gi_topological_explicit_mesh_function.h"

#include "gi_topological_gradient_using_algorithms.h"
#include "gi_topological_regular_grid_restricted.h"

#include "gi_topological_utility_functions.h"
#include "gi_numeric_integrator_expanding_region_stop.h" // not where comparer should be
#include "gi_topological_regular_masked_grid.h"
#include "gi_topological_regular_masked_restricted_grid.h"

#include "gi_max_vertex_labeling.h"
#include "gi_topological_max_vertex_mesh_function.h"
#include "gi_modified_robins.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_basic_geometry.h"


#include <list>
#include <algorithm>
#include <iterator>

/*
 * Arguments Passed on Launch:
 * X Y Z - dimentions of volume
 * connections_remapped.vtp - skeleton file
 * volume.raw - cube volume made from stacked images.
 *
 * Output:
 * prints ligament being iterated over
 * prints points mapped to connected component label
 * observed area.
 * dumps file 'slice.vtp' which is the cross sectional slice of
 * the volume orthogonal to ligament currently being iterated over.
 *
 * Note:
 * currently iterates over WINDOW_SIZE points centered at middle of ligament
 * to identify cross sectional slice of interest points in interior
 * of ligament marked in order to identify correct connected component
 * to sum interior points for area.
 * Boundary cases pose issue if marked points not found. results in large
 * incorrect area calculation.
 */
typedef GInt::RegularGrid3D GridType;


using namespace GInt;

using namespace std;

int X, Y, Z;
int lig_boundary_x, lig_boundary_y, lig_boundary_z;
static float* dist_field;
static float* marked_dist_field;
static float* dist_field_slice;

#define WINDOW_SIZE 12

// crosssection atribute list to fill

vtkSmartPointer<vtkDoubleArray> cs_area_values = vtkSmartPointer<vtkDoubleArray>::New();
vtkSmartPointer<vtkDoubleArray> cs_minradius_values = vtkSmartPointer<vtkDoubleArray>::New();
vtkSmartPointer<vtkDoubleArray> cs_perimeter_values = vtkSmartPointer<vtkDoubleArray>::New();
//vtkSmartPointer<vtkDoubleArray> lig_curvature_values = vtkSmartPointer<vtkDoubleArray>::New();

struct xyzi {
  double v[3];
  int i;
  vtkIdType pid;
};

std::map<int, vtkIdType> pointid_map;

xyzi operator +(const xyzi& u, const xyzi& v) {
  xyzi p;
  p.v[0] = u.v[0] + v.v[0];
  p.v[1] = u.v[1] + v.v[1];
  p.v[2] = u.v[2] + v.v[2];
  return p;
}

xyzi operator -(const xyzi& u, const xyzi& v) {
  xyzi p;
  p.v[0] = u.v[0] - v.v[0];
  p.v[1] = u.v[1] - v.v[1];
  p.v[2] = u.v[2] - v.v[2];
  return p;
}

double operator *(const xyzi& u, const xyzi& v) {
  double dot;
  dot = u.v[0]*v.v[0] + u.v[1]*v.v[1] + u.v[2]*v.v[2];
  return dot;
}

xyzi operator *(const xyzi& u, const double& v) {
  xyzi p;
  
  p.v[0] = u.v[0]*v;
  p.v[1] = u.v[1]*v;
  p.v[2] = u.v[2]*v;
  
  return p;
}

xyzi toXYZI(std::vector<double> p, int i){
  xyzi v;
  v.v[0] = p[0];
  v.v[1] = p[1];
  v.v[2] = p[2];
  v.i = i;
  return v;
}

struct ligradius{
  xyzi s;
  xyzi r;
};

struct crosssection{
  float perimeter;
  double area;
  float curvature;
  double maxr;
  xyzi p;
};


void CoordinateFromLinearIndex(int idx, int dim_x, int dim_y, double& x, double& y, double& z){
  x =  idx % (dim_x);
  idx /= (dim_x);
  y = idx % (dim_y);
  idx /= (dim_y);
  z = idx;
}

inline unsigned long LinearIndexFromCoordinate(int x, int y, int z, unsigned long X, unsigned long Y) {
  return x + y*X + z*(X*Y);//i*(N*M) + j*M + k; //index=y*WIDTH + x (2d)
}

template<typename I>
void dumpArray(I* f, int* size, string filename){
  unsigned long outsize = size[0]*size[1]*size[2];
  
  I* outdata = new I[outsize];
  for (int i = 0; i < size[0]; ++i)
    for (int j = 0; j < size[1]; ++j)
      for (int k = 0; k < size[2]; ++k) {
        
        unsigned long idx = LinearIndexFromCoordinate(i,j,k,size[0],size[1]);
        outdata[idx] = f[idx];
      }
  
  ofstream outFile;
  outFile.open(filename, ios::out|ios::binary|ios::ate);
  outFile.write((char*)outdata, outsize*sizeof(I));
  outFile.close();
  
  delete [] outdata;
  printf("data written\n");
}

void print(std::vector<double> v){
  
  for(int i=0;i< v.size();++i) {
    printf("%.2f ", v[i]);
  }
  printf("\n");
}
void print(xyzi v, string name = " " ){
  cout << name << v.v[0] << " " << v.v[1] << " " << v.v[2] << endl;
}

float max_element(float a[]){
  float max = 0;
  
  for(int i=0; i < X*Y*Z; i++)
    if(a[i] > max)
      max = a[i];
  return max;
}

void Norm(xyzi* p){
  float mag = sqrt((*p).v[0]*(*p).v[0] + (*p).v[1]*(*p).v[1] + (*p).v[2]*(*p).v[2]);
  
  (*p).v[0] *= 1.0/mag;
  (*p).v[1] *= 1.0/mag;
  (*p).v[2] *= 1.0/mag;
}

xyzi Norm(xyzi p){
  float mag = sqrt(p*p);
  
  p.v[0] *= 1.0/mag;
  p.v[1] *= 1.0/mag;
  p.v[2] *= 1.0/mag;
  return p;
}


/*
 Rotate a point p by angle theta around an arbitrary axis r
 Return the rotated point.
 Positive angles are anticlockwise looking down the axis
 towards the origin.
 Assume right hand coordinate system.
 */
xyzi ArbitraryRotate(xyzi p,double theta, xyzi r)
{
  xyzi q;
  q.v[0] = 0.0;
  q.v[1] = 0.0;
  q.v[2] = 0.0;
  q.i = 97;
  double costheta,sintheta;
  
  Norm(&r);
  costheta = cos(theta);
  sintheta = sin(theta);
  
  q.v[0] += (costheta + (1 - costheta) * r.v[0] * r.v[0]) * p.v[0];
  q.v[0] += ((1 - costheta) * r.v[0] * r.v[1] - r.v[2] * sintheta) * p.v[1];
  q.v[0] += ((1 - costheta) * r.v[0] * r.v[2] + r.v[1] * sintheta) * p.v[2];
  
  q.v[1] += ((1 - costheta) * r.v[0] * r.v[1] + r.v[2] * sintheta) * p.v[0];
  q.v[1] += (costheta + (1 - costheta) * r.v[1] * r.v[1]) * p.v[1];
  q.v[1] += ((1 - costheta) * r.v[1] * r.v[2] - r.v[0] * sintheta) * p.v[2];
  
  q.v[2] += ((1 - costheta) * r.v[0] * r.v[2] - r.v[1] * sintheta) * p.v[0];
  q.v[2] += ((1 - costheta) * r.v[1] * r.v[2] + r.v[0] * sintheta) * p.v[1];
  q.v[2] += (costheta + (1 - costheta) * r.v[2] * r.v[2]) * p.v[2];
  
  return(q);
}

/*
 Rotate a point p by angle theta around an arbitrary line segment p1-p2
 Return the rotated point.
 Positive angles are anticlockwise looking down the axis
 towards the origin.
 Assume right hand coordinate system.
 */
xyzi ArbitraryRotate2(xyzi p,double theta,xyzi p1,xyzi p2)
{
  xyzi q;
  q.v[0] = 0.0;
  q.v[1] = 0.0;
  q.v[2] = 0.0;
  q.i = 97;
  double costheta,sintheta;
  xyzi r;
  
  r.v[0] = p2.v[0] - p1.v[0];
  r.v[1] = p2.v[1] - p1.v[1];
  r.v[2] = p2.v[2] - p1.v[2];
  p.v[0] -= p1.v[0];
  p.v[1] -= p1.v[1];
  p.v[2] -= p1.v[2];
  Norm(&r);
  
  costheta = cos(theta);
  sintheta = sin(theta);
  
  q.v[0] += (costheta + (1 - costheta) * r.v[0] * r.v[0]) * p.v[0];
  q.v[0] += ((1 - costheta) * r.v[0] * r.v[1] - r.v[2] * sintheta) * p.v[1];
  q.v[0] += ((1 - costheta) * r.v[0] * r.v[2] + r.v[1] * sintheta) * p.v[2];
  
  q.v[1] += ((1 - costheta) * r.v[0] * r.v[1] + r.v[2] * sintheta) * p.v[0];
  q.v[1] += (costheta + (1 - costheta) * r.v[1] * r.v[1]) * p.v[1];
  q.v[1] += ((1 - costheta) * r.v[1] * r.v[2] - r.v[0] * sintheta) * p.v[2];
  
  q.v[2] += ((1 - costheta) * r.v[0] * r.v[2] - r.v[1] * sintheta) * p.v[0];
  q.v[2] += ((1 - costheta) * r.v[1] * r.v[2] + r.v[0] * sintheta) * p.v[0];
  q.v[2] += (costheta + (1 - costheta) * r.v[2] * r.v[2]) * p.v[2];
  
  q.v[0] += p1.v[0];
  q.v[1] += p1.v[1];
  q.v[2] += p1.v[2];
  return(q);
}



/*
 *
 *
 |---------------- Begin Slice / Ligament Attribute Computation ----------------|
 *
 *
 */


void smooth(int times, vector<xyzi>& line) {
  
  for (int k = 0; k < times; k++) {
    vector<xyzi> copyline = line;
    for (int i = 1; i < line.size() - 1; i++) {
      
      double v1[3], v2[3], v3[3];
      for(int d=0;d<3; d++){
        v1[d]=line[i - 1].v[d]*0.25;
        v2[d]=line[i].v[d]*0.5;
        v3[d]=line[i + 1].v[d]*0.25;
        
        copyline[i].v[d]=v1[d]+v2[d]+v3[d];
      }
      //copyline[i] = line[i - 1]* 0.25 + line[i]* 0.5 + line[i + 1] * 0.25;
    }
    line = copyline;
  }
  
}



//Group points in slice of volume orthogonally intersecting ligament
//by connected components based on gradient(iso) values. Sum points in connected
// component previously marked with largest gradient value in the interior.
crosssection ConnectedComponentsAttr(vtkImageData* cross_section_slice, xyzi p, std::vector<double> bbox, bool output_slice=false, int lig_id=0, std::string folder_path=""){
  
  float threshold = 0.0;
  
  int* dims = cross_section_slice->GetDimensions();
  crosssection cs;
  
  GridType* underlying_grid;
  // set up structures to navigate grid, and load the 3d image
  underlying_grid = new GridType(GInt::Vec3i(dims[0], dims[1], dims[2]), GInt::Vec3b(0, 0, 0));
  
  printf("loaded field function\n");
  
  VolumeConnectedComponents cc(underlying_grid);
  
  // for area
  DenseLabeling<char> *maskvol = new DenseLabeling<char>(underlying_grid->NumElements());
  
  float max_field_val = max_element(marked_dist_field);
  
  std::vector<int> lig_idx;
  int lig_cc_idx = 0;
  for(INDEX_TYPE i=0; i < underlying_grid->NumElements();i++){
    float gradient_value = dist_field_slice[i];
    bool in_cc = false;
    if(gradient_value > max_field_val){
      lig_cc_idx = i;
      lig_idx.insert(lig_idx.end(), i);
      in_cc = true;
      maskvol->SetLabel(i, 44);
    }
    else
      maskvol->SetLabel(i, gradient_value >= threshold);
  }

  cc.PerformUnionFind(maskvol);
  
  int* output = new int[underlying_grid->NumElements()];
  //int my_dims[3]={dims[0], dims[1], dims[2]};
  
  //dumpArray(output, my_dims, "maskvol.raw");
  // new connected components given unique IDs
  cc.mIDVol->ReMapIds(output);
  
  int perim = 0;
  int area = 0;
  int num_elm = underlying_grid->NumElements();
  for(int i=0; i<num_elm; i++){
      //int observed_label = maskvol->GetLabel(lig_cc_idx);
      //observed_label //==output[lig_cc_idx]
    if(output[i] == output[lig_cc_idx]){ //>= 0 && dist_field_slice[i] >= 0.0 ){ //max_field_val){
      //in connected component of interest and inside lig
      area += 1;
      // get coordinates to mask image and find border
      double x, y, z;
      CoordinateFromLinearIndex(i, dims[0], dims[1], x, y, z);
      // mask image over area observed
      double* pixel = static_cast<double*>(cross_section_slice->GetScalarPointer(x,y,z));
      *pixel = max_field_val + 1.0;
      //see if on perimeter
      long nbr_pixel_x_up = LinearIndexFromCoordinate(x+1, y, z, dims[0], dims[1]);
      long nbr_pixel_y_up = LinearIndexFromCoordinate(x, y+1, z, dims[0], dims[1]);
      long nbr_pixel_x_dwn = LinearIndexFromCoordinate(x-1, y, z, dims[0], dims[1]);
      long nbr_pixel_y_dwn = LinearIndexFromCoordinate(x, y-1, z, dims[0], dims[1]);
      if(dist_field_slice[nbr_pixel_x_up] <= 0 || dist_field_slice[nbr_pixel_y_up] <= 0
         || dist_field_slice[nbr_pixel_x_dwn] <= 0 || dist_field_slice[nbr_pixel_y_dwn] <= 0){
        perim+=1;
        // mask perim
        *pixel = max_field_val + 2.0;
      }
      
      //long neighbor_pixel_z = LinearIndexFromCoordinate(x, y, z, dims[0], dims[1]);
    }
  }

  ofstream outFile;
  std::stringstream filename;
  filename<<folder_path<<"/cs_seg_"<<lig_id<<"_"<<dims[0]<<"_"<<dims[1]<<"_"<<dims[2]<<".raw";
    
  if(output_slice){
    fprintf(stderr,"saving %s\n", filename.str().c_str());
    outFile.open(filename.str().c_str(), std::ofstream::binary);
    
    outFile.write((char*)output, underlying_grid->NumElements()*sizeof(int));
    outFile.close();
   
  }

  //cout << "AREA FROM CONNECTED COMP " << area <<endl;
  //cout << "PERIMETER FROM CONNECTED COMP " << perim <<endl;
  cs.area = area;
  cs.p = p;
  cs.perimeter = perim;
  return cs;
}

//return the gradient field value at point p
float Gradient(xyzi p){
  float gradient = dist_field[LinearIndexFromCoordinate(p.v[0],p.v[1] ,p.v[2],X,Y)];
  return gradient;
}

void SetMinMax(std::vector<xyzi>& pointset, double* minv, double* maxv) {
  for (int j = 0; j < 3; j++) {
    minv[j] = maxv[j] = pointset[0].v[j];
  }
  
  for (auto p : pointset) {
    for (int j = 0; j < 3; j++) {
      if (minv[j] > p.v[j]) minv[j] = p.v[j];
      if (maxv[j] < p.v[j]) maxv[j] = p.v[j];
    }
  }
}
// Use predefined transformation matrix to translate / rotate point
// based on basis and column in vmatrix. Now, the original input volume
// is in the "input coordinate system" of vtkImageReslice. We can call the
// input coordinate system "x" and the output coordinate system "x'".
// The relationship between these coordinates is as follows:  x = T*M*x' where
// "T" is ResliceTransform and "M" is ResliceMatrix
double* Transform(vtkMatrix4x4 *vmatrix, vtkAbstractTransform* vtktransform, double p[3]){
  
  vtkSmartPointer<vtkTransform> vtkTranslation = vtkSmartPointer<vtkTransform>::New();
  
  vtkTranslation->SetMatrix(vmatrix);
  
  //vtkTranslation->Inverse();
  
  vtkTranslation->Update();
  
  
  double* tp = vtkTranslation->TransformPoint(p[0], p[1], p[2]);
  return p;
}


// Project point onto orthonormal basis defining slice plane perpendicular to
// point on line at ligament center
void Proj(double* p, std::vector<double> x, std::vector<double> y , std::vector<double> z ){
  
  double dot_x = 0.0;
  double dot_y = 0.0;
  double dot_z = 0.0;
  double norm_x = 0.0;
  double norm_y = 0.0;
  double norm_z = 0.0;
  
  for( int i = 0; i < x.size(); i++){
    norm_x += x[i]*x[i];
    norm_y += y[i]*y[i];
    norm_z += z[i]*z[i];
    dot_x += p[i]*x[i];
    dot_y += p[i]*y[i];
    dot_z += p[i]*z[i];
  }
  
  // scale by projection factor
  
  for( int i = 0; i < x.size(); i++){
    p[i] = (dot_x/sqrt(norm_x))*x[i] + (dot_y/sqrt(norm_y))*y[i] + (dot_z/sqrt(norm_z))*z[i];
    //std::cout << " p.v: " << p[i] << std::endl;
  }
  
}

//inner product.
double Inner( double u[3], std::vector<double> v){
  double s_prod = 0;
  s_prod += u[0] * v[0];
  s_prod += u[1] * v[1];
  s_prod += u[2] * v[2];
  return s_prod;
}
// overloaded inner product
double Inner( std::vector<double> u, std::vector<double> v){
  double s_prod = 0;
  s_prod += u[0] * v[0];
  s_prod += u[1] * v[1];
  s_prod += u[2] * v[2];
  return s_prod;
}

// Return orientation from Z axis of point
double orientation_z(xyzi p){
  double z[3] = {0,0,1};
  double inner_prod = Inner(z, {p.v[0], p.v[1], p.v[2]});
  double norm = sqrt(p.v[0]*p.v[0]+p.v[1]*p.v[1]+p.v[2]*p.v[2]);
  // divisor = norm * norm(z);
  double orient = acos(inner_prod/norm)*57.2957;
  
  if(orient >= 90)
    orient = 180-orient;
  
  return orient;
}

//L2 norm
double L2(std::vector<double> u){
  return sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
}

//cosine similarity between vectors u and v
double cosSim(std::vector<double> u, std::vector<double> v){
  double numerator = Inner(u,v);
  double denominator = L2(u)*L2(v);
  double rad = acos(numerator/denominator);
  double pi = 4*atan(1);
  double degree = rad*(180.0/pi);
  return degree;
}

//cosine similarity between vectors u and v as type xyzi
double cosSim(xyzi u, xyzi v){
  double numerator = u*v;
  double denominator = sqrt(u*u)*sqrt(v*v);
  double rad = acos(numerator/denominator);
  double pi = 4*atan(1);
  double degree = rad*(180.0/pi);
  return degree;
}

// increase magnitude of vector v by scale
xyzi Scale(xyzi v, float scale){
  double mag = sqrt(v.v[0]*v.v[0]+v.v[1]*v.v[1]+v.v[2]*v.v[2]);
  v.v[0] *= scale;///mag;
  v.v[1] *= scale;///mag;
  v.v[2] *= scale;///mag;
  return v;
}

//walk up direction of vector p with gradient(iso) value g by stepsize step
ligradius WalkUpToZero(xyzi p, xyzi p_theta, float g, float step){
  float gradient = g;
  ligradius r;
  
  while( gradient  > 0.0){
    p_theta = Scale( p_theta, 1.0+step);
    p = p + p_theta;
    r.r = p_theta;
    r.s = p;
    gradient = Gradient(p);
    
  }
  return r;
}




/*
 given the vtkImageData crosssectional slice of the volume intersecting orthogonly
 at current ligament collect distance field values of slice and compute
 connected components of slice. Values aroung p and along ligament interior marked
 with larger gradient value for identifying in connected component. Count points in
 marked connected component to obtain area.
 */
crosssection cross_section( vtkImageData* imageData, xyzi p, int lig_id, std::vector<double> bbox, bool save_img, std::string folder_path){
  
  // Create an image data
  
  // Specify the size of the image data
#if VTK_MAJOR_VERSION <= 5
  imageData->SetNumberOfScalarComponents(1);
  imageData->SetScalarTypeToDouble();
#else
  imageData->AllocateScalars(VTK_DOUBLE,1);
#endif
  
  int* dims = imageData->GetDimensions();
  
  //std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;
  //std::cout << "Number of points: " << imageData->GetNumberOfPoints() << std::endl;
  //std::cout << "Number of cells: " << imageData->GetNumberOfCells() << std::endl;
  
  
  dist_field_slice = new float[dims[0]*dims[1]*dims[2]];
  
  for (int x = 0; x < dims[0]; x++){
    for (int y = 0; y < dims[1]; y++){
      for (int z = 0; z < dims[2] ; z++){
        double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
        dist_field_slice[LinearIndexFromCoordinate(x, y, z, dims[0], dims[1]) ] = pixel[0];
      }
    }
  }
  
  //compute connected components of slice. Identify connected component with
  //pre-marked gradient points and count points within connected components
  crosssection cs = ConnectedComponentsAttr(imageData, p, bbox, save_img, lig_id, folder_path);
  cout << "Observed Area: " << cs.area <<endl;
  cout << "Observed Perim: " << cs.perimeter <<endl;
  //cs_area_values -> InsertComponent(p.i, 0, cs.area);
  //add found area of crosssectional slice at point p to cs_area_values vtksmartpointer
  
  if(save_img){
    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    
    std::string f_name = folder_path+"/cs_orig_"+std::to_string(lig_id)+".vti";
    printf("saving %s\n",f_name.c_str());
    writer->SetFileName(f_name.c_str());
    writer->SetInputData(imageData);
    writer->Write();
  }
  return cs;
}


#ifdef VTK_ENABLED

vtkSmartPointer<vtkPolyData> ReadVTKLigaments(char* filename, std::vector<xyzi>& pointset,
                                              std::map<int,std::vector<xyzi> >& lines, std::map<int,int>* point_to_line_map=NULL) {
  
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename);
  reader->Update();
  
  vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
  
  //polyData->Print(std::cout);
  
  vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
  
  // Create a cell array to store the lines in and add the lines to it
  vtkSmartPointer<vtkCellArray> cells = polyData->GetLines();
  
  vtkSmartPointer<vtkPointData> pointdata = polyData->GetPointData();
  vtkSmartPointer<vtkIntArray> ids = vtkIntArray::SafeDownCast(pointdata->GetArray("Id"));
  
  //printf("%s has:\t%d points and %d ids\n", filename, points->GetNumberOfPoints(), ids->GetNumberOfValues());
  //std::cout << "There are " << polyData->GetNumberOfLines() << " lines." << std::endl;
  
  unsigned long point_id=0;
  
  polyData->GetLines()->InitTraversal();
  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
  while(polyData->GetLines()->GetNextCell(idList)){
    int npoints = idList->GetNumberOfIds();
    //std::cout << "Line has " << npoints << " points." << std::endl;
    
    if(npoints==0)
      continue;
    
    int id = ids->GetValue(idList->GetId(0)); //line id
    
    if(lines.count(id) > 0)
      continue;
    
    std::vector<xyzi>& l = lines[id];
    
    for(vtkIdType pointId = 0; pointId < npoints; pointId++){
      //std::cout << idList->GetId(pointId) << " ";
      vtkIdType pid = idList->GetId(pointId);
      double *p = points->GetPoint(pid);
      
      int lid = ids->GetValue(pid);
      xyzi t;
      for (int j = 0; j < 3; j++) t.v[j] = p[j];
      t.i = lid;
      t.pid=pid;
      
      pointid_map[lid]=pid;
      
      pointset.push_back(t);
      l.push_back(t);
      
      point_id++;
      
      if(point_to_line_map!=NULL)
        (*point_to_line_map)[idList->GetId(pointId)] = id;
      
    }
    
    if(l.size()!=npoints){
      printf("error lines size %d vtk %d\n", l.size(), npoints);
      assert(false);
    }
    //std::cout << std::endl;
  }
  
  printf("Read complete correctly\n");
  
  return polyData;
}


void updateVTKLigaments(vtkSmartPointer<vtkPolyData> polyData, std::string filename_out){
  
  vtkSmartPointer<vtkPolyData> polyData_out = vtkSmartPointer<vtkPolyData>::New();
  polyData_out->DeepCopy(polyData);
  
  // Remove exisisting arrays for the features we are computing here
  polyData_out->GetPointData()->RemoveArray("Cross Section Area");
  polyData_out->GetPointData()->RemoveArray("Cross Section Perimeter");
  //polyData_out->GetPointData()->RemoveArray("Curvature");
  
  // Add updated arrays to the file
  polyData_out->GetPointData()->AddArray(cs_area_values);
  polyData_out->GetPointData()->AddArray(cs_minradius_values);


//
// need to add perim
//
//polyData_out->GetPointData()->AddArray(cs_perimeter_values);
  //polyData_out->GetPointData()->AddArray(lig_curvature_values);
  
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(polyData_out);
  writer->SetFileName(filename_out.c_str());
  writer->Write();
  
  printf("Ligaments updated\n");
}


#endif

/*
 create vtkImageData cube filled with gradient values in gradient_field.
 mark points around p with fabricated max gradient value to later identify connected
 component of interest in crosssectional slice.
 */
vtkImageData* create_vtkFunction(xyzi p, xyzi p_og, vtkImageData* cube, std::vector<double> bbox) {
  
  
  //printf(" --- creating a vtk vol around ligament point...");
  fflush(stdout);
  
  cube->SetOrigin(0.0, 0.0, 0.0);
  cube->SetSpacing(1.0,1.0,1.0);
  
  cube->SetDimensions((int)bbox[1]-(int)bbox[0]+1, (int)bbox[3]-(int)bbox[2]+1, (int)bbox[5]-(int)bbox[4]+1); //set dim to max-min of bbox
  
  cube->AllocateScalars(VTK_DOUBLE, 1);
  cube->GetPointData()->GetScalars()->SetName("function");
  
  float max_field_val = max_element(dist_field);
  
  for (int x = bbox[0]; x < bbox[1]; x++) { //iterate of min to max all axes of bbox(within original vol)
    for (int y = bbox[2]; y < bbox[3]; y++) {
      for (int z = bbox[4]; z < bbox[5]; z++) {
        
        double* pixel = static_cast<double*>(cube->GetScalarPointer(x-(int)bbox[0],y -(int)bbox[2],z - (int)bbox[4]));
        
        // fill vtkimage with marker around p
        float field_val = marked_dist_field[LinearIndexFromCoordinate(x,y ,z,X,Y)];
        pixel[0]= field_val;
        
        // mark cube around point p original to identify ligament in connected component computation
        // if cube is < 5 the second ligament can be missed. I'm not sure what size value
        // to use to be able to mark all ligaments including boundary cases while not
        //risking marking too  large a cube around point p
        if( abs(x-int(p_og.v[0])) < 1 && abs( y- int(p_og.v[1])) < 1 && abs( z - int(p_og.v[2])) < 1 && field_val > 0)
          pixel[0]= max_field_val+ 0.0001;// 1.0; // ! I worry this will sometimes mark pixels outside interior of ligament giving false area computations. If only interior of ligament marked (in main) then false zero gradient values can occur in ligament. If too small a numer of pixels marked here boundary ligaments pose a problem of not being found.
      }
    }
  }
  /*
   vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
   writer->SetFileName("cube_window.vti");
   writer->SetInputData(cube);
   writer->Write();
   */
  printf(" vtk function done.\n");
  
  
  return cube;
  
}

// Using transormation matrix vmatrix compute slice vtkSlice of volume vtkdata
vtkImageData* compute_slice(vtkImageData *vtkdata, vtkMatrix4x4 *vmatrix, vtkAbstractTransform *vtktransform, vtkImageData *vtkSlice) {
  //#ifndef USE_VTK
  //  print Visit not available!;
  //#else
  vtkSmartPointer<vtkImageReslice> vtkSlicer = vtkSmartPointer<vtkImageReslice>::New();
  vtkSlicer->SetOutputDimensionality(2);
  vtkSlicer->SetTransformInputSampling(1);
  // having auto crop on ensures nothing cropped however adds points of positive value to edge
  //vtkSlicer->SetAutoCropOutput(true);
  vtkSlicer->SetBorder(false);
  vtkSlicer->SetWrap(1);
  vtkSlicer->SetInterpolationModeToCubic();
  
  vtkSlicer->SetResliceAxes(vmatrix);
  
  vtkSlicer->SetInputData(vtkdata);
  vtkSlicer->Update();
  
  //vtkImageData *slice = vtkImageData::New();
  vtkSlice->DeepCopy(vtkSlicer->GetOutput());
  
  //vtkSmartPointer<vtkTransform> vtkTransform = vtkSmartPointer<vtkTransform>::New();
  vtktransform = vtkSlicer -> GetResliceTransform();
  //vtkTransform = vtkSlicer -> ResliceTransform();
  return vtkSlice;
  
}

std::vector<double> bounding_box(xyzi p, int bounds_x, int bounds_y, int bounds_z){

    float px = p.v[0]; float py = p.v[1]; float pz = p.v[2];

    while(px - bounds_x < 0)
        bounds_x += 1;
    while(py - bounds_y < 0)
        bounds_y += 1;
    while(pz - bounds_z < 0)
        bounds_z += 1;
    while(px + bounds_x > X-1)
        bounds_x -= 1;
    while(py + bounds_y > Y-1)
        bounds_y -= 1;
    while(pz + bounds_z > Z-1)
        bounds_z -= 1;

  float min_x = ((px-bounds_x)<0)?0:(px-bounds_x);
  float max_x = ((px+bounds_x)>X)?X-1:(px+bounds_x);
  float min_y = ((py-bounds_y)<0)?0:(py-bounds_y);
  float max_y = ((py+bounds_y)>Y)?Y-1:(py+bounds_y);
  float min_z = ((pz-bounds_z)<0)?0:(pz-bounds_z);
  float max_z = ((pz+bounds_z)>Z)?Z-1:(pz+bounds_z);
  std::vector<double> box(6);
  box[0] = min_x; box[2] = min_y; box[4] = min_z;
  box[1] = max_x; box[3] = max_y; box[5] = max_z;
  return box;
}


std::vector<int> Sum( double u[3], std::vector<double> v){
  std::vector<int> s(3);
  s[0] = u[0] + v[0];
  s[1] = u[1] + v[1];
  s[2] = u[2] + v[2];
  return s;
}

std::vector<int> Sub( double u[3], std::vector<double> v){
  std::vector<int> s(3);
  s[0] = u[0] - v[0];
  s[1] = u[1] - v[1];
  s[2] = u[2] - v[2];
  return s;
}

std::vector<double> Sub( double u[3], double v[3]){
  std::vector<double> s(3);
  s[0] = u[0] - v[0];
  s[1] = u[1] - v[1];
  s[2] = u[2] - v[2];
  return s;
}

// return unit vector of v
std::vector<double> Unit(std::vector<double> v){
  float mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] *= 1.0/mag;
  v[1] *= 1.0/mag;
  v[2] *= 1.0/mag;
  return v;
}
xyzi Unit(xyzi p){
  float mag = sqrt(p.v[0]*p.v[0] + p.v[1]*p.v[1] + p.v[2]*p.v[2]);
  p.v[0] *= 1.0/mag;
  p.v[1] *= 1.0/mag;
  p.v[2] *= 1.0/mag;
  return p;
}

// unit normal vector of two points
std::vector<double> Normal(xyzi p1, xyzi p2){
  std::vector<double> nz(3);
  nz[0] = p1.v[0] - p2.v[0];
  nz[1] = p1.v[1] - p2.v[1];
  nz[2] = p1.v[2] - p2.v[2];

  
  nz = Unit(nz);
  return nz;
}

// return u X v
std::vector<double> Cross(std::vector<double> u, std::vector<double> v){
  std::vector<double> cp(3);
  
  cp[0] = u[1]*v[2] - u[2]*v[1];
  cp[1] = u[2]*v[0] - u[0]*v[2];
  cp[2] = u[0]*v[1] - u[1]*v[0];
  
  return cp;
}

// increase magnitude of vector by scale
std::vector<double> Scale(std::vector<double> v, float scale){
  double mag = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  //scale = 1.0+scale;
  v[0] *= scale/mag;
  v[1] *= scale/mag;
  v[2] *= scale/mag;
  return v;
}

double to360(double x){
  x = fmod(x,360);
  if (x < 0)
    x += 360;
  return x;
}

int main(int argc, char** argv) {
  
  if(argc < 6){
    printf("usage: extractcrosssection <X> <Y> <Z> <lines.vtp> <dist_field.raw> <only center slice>\n");
    return 1;
  }
//  int X=300;
//  int Y=300;
//  int Z=300;
//  argv[4] = "/home/sam/Documents/PhD/Research/GradIntegrator/data/000_dist.raw.lig40.vtp";
//  argv[5] = "/home/sam/Documents/PhD/Research/GradIntegrator/data/000_dist.raw";
//  argv[6] = 0;

  sscanf(argv[1], "%d", &X);
  sscanf(argv[2], "%d", &Y);
  sscanf(argv[3], "%d", &Z);
  
  int only_center = 0;

  //if(argc > 5)
    sscanf(argv[6], "%d", &only_center);
    
	int max_radius;
	sscanf(argv[7], "%d", &max_radius);

  dist_field = new float[X*Y*Z];
  marked_dist_field = new float[X*Y*Z];
  
  int num_elems=X*Y*Z;
  
  FILE* fin = fopen(argv[5], "rb");
  fread(dist_field, sizeof(float), num_elems, fin);
  fclose(fin);
  
  printf("start %s !\n", argv[1]);
  
  // read in the data
  std::vector<xyzi> line_pointset;
  std::map<int,std::vector<xyzi> > lines;
  
  char* lines_filename = argv[4];
  
  printf("reading file %s\n", lines_filename);

  std::string folder_path(lines_filename);


  folder_path = folder_path.substr(0, folder_path.find_last_of("/\\"));

  // read input into pointsets and lines
  vtkSmartPointer<vtkPolyData> polyData = ReadVTKLigaments(lines_filename, line_pointset, lines);
  
  vtkSmartPointer<vtkPointData> pointdata = polyData->GetPointData();
  vtkSmartPointer<vtkIntArray> ids = vtkIntArray::SafeDownCast(pointdata->GetArray("Id"));
  
#if 0
  // hack for testing specific ligament
  std::map<int,std::vector<xyzi> > lines_test;
  
  int check = 0;
  int lig_id = 23;
  //std::vector<xyzi> line;
  for(auto& l: lines){
    check += 1;
    if(l.first == lig_id){
      lines_test.insert( std::pair<int,std::vector<xyzi>>(l.first, l.second));
      break;
      //line = l.second;
    }
  }
#endif

  if(only_center){
    std::map<int,std::vector<xyzi> > lines_rnd;

    int n_lines=lines.size();
    std::random_device rd;
    std::mt19937 gen( rd());
    std::uniform_int_distribution<> dis(0, n_lines-1);

    for (int n=0; n<10; ++n){
      int idx=dis(gen);
      auto l=std::next(lines.begin(),dis(gen));
      lines_rnd[(*l).first]=(*l).second;
    }

    lines = lines_rnd;
    lines_rnd.clear();
    
  }

  cs_area_values->SetNumberOfComponents(1);
  cs_area_values->SetName("CrossArea");

  cs_minradius_values->SetNumberOfComponents(1);
  cs_minradius_values->SetName("MinRadius");


//cs_perimeter_values->SetNumberOfComponents(1);
//cs_perimeter_values->SetName("CrossPerimeter");
  
  //lig_curvature_values->SetNumberOfComponents(1);
  //lig_curvature_values->SetName("Curvature");
  
  // Iterate over ligaments filling cube with gradient values and taking slice vtkSlice
  // orthogonal to ligament l.

  int line_id = 0;
  int line_point_id=0;
  GInt::RegularGrid3D* datagrid = new GInt::RegularGrid3D({ X, Y, Z }, { 0,0,0 });
  GInt::RegularGridTrilinearFunction* tri =
	  new GInt::RegularGridTrilinearFunction(datagrid, dist_field);
   
  int dist = max_radius; // how far away to look
  int grid_points = dist * 2 + 1; // this will make 1 2d grid unit = 1 3d grid unit
// re-use some code to be able to do UF on 2d grid
  GInt::RegularGrid2D* small_grid = new GInt::RegularGrid2D({ grid_points, grid_points }, { 0,0 });
  float THRESHOLD = 0;



  printf("before loop lines %d\n", lines.size());
  for (auto& it_p : lines) {

	  std::vector<xyzi>& l = it_p.second;
          
	  // find lowest value - will be saddle position
	  int lowest_pos = 0; 
	  float low_val = tri->TriLinInterpValue(Vec3d(l[lowest_pos].v[0], l[lowest_pos].v[1], l[lowest_pos].v[2]));
	  for (int pos = 0; pos < l.size(); pos++) {
		  float otherval = tri->TriLinInterpValue(Vec3d(l[pos].v[0], l[pos].v[1], l[pos].v[2]));
		  if (otherval < low_val) {
			  lowest_pos = pos;
			  low_val = otherval;
		  }
	  }

	  //l.insert(l.begin(), it_p.second.begin(), it_p.second.end()); // copy list, so we don't destroy original geometry
	  //smooth(10000, l);
	  Vec3f origin(l[lowest_pos].v[0], l[lowest_pos].v[1], l[lowest_pos].v[2]);

	  // get sane points for normal computation - 1/3 the way to the ends
	  int num_to_start = lowest_pos;
	  int pos_to_start = lowest_pos - (num_to_start / 3);
	  Vec3f start_v(l[pos_to_start].v[0], l[pos_to_start].v[1], l[pos_to_start].v[2]);

	  int num_to_end = l.size() - lowest_pos;
	  int pos_to_end = lowest_pos + (num_to_end / 3);
	  Vec3f end_v(l[pos_to_end].v[0], l[pos_to_end].v[1], l[pos_to_end].v[2]);

	  Vec3f normal = end_v - start_v;
	  normal.Normalize();

	  GInt::PlaneGrid2D slice(origin, normal, grid_points, grid_points, dist, dist);

	  // now fill in grid
	  struct Quad {
		  vtkIdType id[4];
	  };
	  vector<Quad> quads;
	  vector<Vec3d> grid_vertices;
	  grid_vertices.reserve(grid_points*grid_points);
	  float* values = new float[grid_points*grid_points];
	
	  // y moves slowest
	  for (int jj = 0; jj < grid_points; jj++) {
		  // x moves fastest
		  for (int ii = 0; ii < grid_points; ii++) {
			  Vec3d sample_point = slice.Get3DPosition(ii, jj); // get the position of the ii,jj 2d grid sample
			  float val = tri->GetMinValue(); // use min value of dataset if sample point is outside grid
			  if (datagrid->Inside(sample_point)) {
				  val = tri->TriLinInterpValue(sample_point); // else sample the image
			  }
			  values[ii + jj*grid_points] = val;

			  // now record positions for the grid
			  grid_vertices.push_back(sample_point);
			  if (ii < grid_points - 1 && jj < grid_points - 1) {
				  // make a quad with ids
				  quads.push_back({
					  ii + jj*grid_points,
					  ii + 1 + jj*grid_points,
					  ii + 1 + (jj + 1)*grid_points,
					  ii + (jj + 1)*grid_points
				  });
			  }

		  }
	  }

	  // now do region-growing to count cross section
	  unordered_set<int> visited;
	  queue<Vec2l> frontier;
	  frontier.push({ dist, dist }); // start with the middle point since we KNOW it is in the right CC
	  while (!frontier.empty()) {
		  auto coords = frontier.front();
		  frontier.pop();
		  auto id = small_grid->Index2d(coords);
		  if (visited.count(id) != 0) continue; // skip already seen things
		  visited.insert(id);
		  Vec2l negs[8];
		  int num_neg = small_grid->GatherExistingNeighborsAll8(coords, negs);
		  for (int i = 0; i < num_neg; i++) {
			  auto nid = small_grid->Index2d(negs[i]);
			  if (values[nid] > THRESHOLD)
				  frontier.push(negs[i]);
		  }
	  }

	  // now visited holds the connected component from middle point
      double cc_size = visited.size();


	  // BEGIN VTK STUFF
	  // now we are ready to make vtk polydata for sanity checking
	  vtkNew<vtkPolyData> grid2d;
	  vtkNew<vtkPoints> points;
	  vtkNew<vtkCellArray> polys;
	  vtkNew<vtkFloatArray> scalars;
	  vtkNew<vtkUnsignedCharArray> cc_scalars;
	  scalars->SetName("dist_vals");
	  cc_scalars->SetName("cc_visited");
	  unsigned char* cc_array = new unsigned char[small_grid->NumElements()](); // *should zero initialize																				//memset(cc_array, 0, small_grid->NumElements()*sizeof(unsigned char))
	  for (auto id : visited) cc_array[id] = 1;

	  // Load the point, cell, and data attributes.
	  for (auto i = 0ul; i < grid_vertices.size(); ++i)
	  {
		  points->InsertPoint(i, grid_vertices[i].m_v);
		  scalars->InsertTuple1(i, values[i]);
		  cc_scalars->InsertTuple1(i, cc_array[i]);
	  }
	  for (auto& i : quads)
	  {
		  polys->InsertNextCell(vtkIdType(4), i.id);
	  }
	  // We now assign the pieces to the vtkPolyData.
	  grid2d->SetPoints(points);
	  grid2d->SetPolys(polys);
	  grid2d->GetPointData()->SetScalars(scalars);
	  grid2d->GetPointData()->AddArray(cc_scalars);
	  // now write to vtp

	  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      std::string l_filename(lines_filename);
      std::string f_name = l_filename +"_Slice_LIGG"+std::to_string(line_id)+".vtp";
          writer->SetInputData(grid2d);          
          writer->SetFileName(f_name.c_str());
	  printf("writing slice to %s\n", f_name.c_str());
	  writer->Write();


	  // END VTK STUFF

      // Add cross section to polydata
      for(int p=0 ; p < l.size(); p++){
		  cs_minradius_values->InsertComponent(line_point_id + p, 0, low_val);//

        cs_area_values->InsertComponent(line_point_id + p, 0, cc_size);//
      }



      line_point_id = line_point_id + l.size(); //increment point id over entire line
	  line_id++;
  }
  //all lines iterated save vtk polydata
  updateVTKLigaments(polyData, string(lines_filename) + "_cross_section.vtp");
  //else
  cout << "Poly data with cross section writen to " << string(lines_filename) << "_cross_section.vtp" <<endl;
#if 0
  int line_it = 0;
  unsigned long line_point_id=0;
  //for(auto& line: lines){
  //#pragma omp parallel for shared( line_it, line_point_id)//zip-zap-zoom
  for(int line_count = 0; line_count < lines.size(); line_count++){
    
    vtkImageData* cube = vtkImageData::New();
    vtkImageData* vtkSlice = vtkImageData::New();
    vtkMatrix4x4* SliceMatrix;
    vtkAbstractTransform *vtkTransform;

    std::map<int,std::vector<xyzi> >::iterator line = lines.begin();
    advance(line, line_count);

    line_it += 1;
    
    cout << "------------------------------------Processing Lig, ID: " << line->first << endl ;
    
    std::vector<xyzi> l = line->second;
    
    // smooth ligament sailing
    smooth(10000, l);
   
    // the following could be used to debug only a few 
    if(l.size()==1){
      for(int p = 0; p < l.size(); p++){ // we only compute stuff every "step" (see next line)                                                       
        cs_area_values->InsertComponent(line_point_id + p, 0, -1);  //#l[p_id]                                                                   
	cs_perimeter_values->InsertComponent(line_point_id + p, 0, -1);
      }
      line_point_id+=l.size();
      continue;
    }
    
    //used to mark points on ligament of interest and around point p
    float max_field_val = max_element(dist_field);
    
    //if we want to compute curvature of ligament via sum of cosine
    // similarity of normal vectors to slice
    //double curvature = 0;
    int lig_id = line->first;
    
    std::cout << "Lig Size: " << l.size() << std::endl;
    /*
     Iterate over points in middle of ligament. Currently 4 points
     centered at middle of ligament.
     */
    int p_ = 0;
    int step = l.size()/2; // step size will set the slice sampling rate (half the size will be one slice in the center)
    int stop = l.size();

    if(true){ //only_center==0){ //in case only cross section of center point is wanted
      p_ = 0;
      stop = l.size();
    }else{
      p_ = l.size()/2 - 2;
      stop = l.size()/2 + 2;
    }

    for(int p = 0; p < l.size(); p++){ // we only compute stuff every "step" (see next line)
      if(p % step != 0){ // set -1 for points on the ligament we skipped
        cs_area_values->InsertComponent(line_point_id + p, 0, -1);//#l[p_id]
        cs_perimeter_values->InsertComponent(line_point_id + p, 0, -1);
        continue;
      }

      // Ensure not stepping out of ligament for small ligaments
      if(p<l.size() && p+step<l.size()){
            // two points some distance along msc away from one another
            int p_store = p;
            xyzi p1 = l[p];
            xyzi p2 = l[p+step];
            vtkIdType pid = l[p].pid; //because ids remapped during transformations
            std::cout << "p " << p << " p_step " << p+step <<std::endl;
        
            // normal to two points on skeleton
        std::vector<double> nz = Normal(p1,p2);

        int steps = step;
        //std::cout << "STEP " << step << " STEPS " << steps <<std::endl;

        //Mark all points in line with fabricated highest gradient value to be identified
        //later through connected components.
        marked_dist_field = dist_field;
        int marker_base = WINDOW_SIZE;
        int marker_top = WINDOW_SIZE;
        // distance along line we can travel to set markers
        while( p - marker_base < 0)
          marker_base += 1;
        while( p + marker_top > l.size())
          marker_top -= 1;
        while( p + marker_top > X-1)
          marker_top -= 1;
        while( p + marker_top > Y-1)
          marker_top -= 1;
        while( p + marker_top > Z-1)
          marker_top -= 1;
        
        // mark volume at p
        marked_dist_field[LinearIndexFromCoordinate(l[p].v[0],l[p].v[1] ,l[p].v[2],X,Y)] = max_field_val + 0.0001;
        // mark volume at end of vector used for normal
        marked_dist_field[LinearIndexFromCoordinate(l[steps].v[0],l[steps].v[1] ,l[steps].v[2],X,Y)] = max_field_val + 0.0001;
        for(int pp = p - marker_base; pp < p + marker_top; pp++){
          if( pp>0 && pp < l.size() ){
            xyzi c = l[pp];
            if(marked_dist_field[LinearIndexFromCoordinate(c.v[0],c.v[1] ,c.v[2],X,Y)] > 0)
              marked_dist_field[LinearIndexFromCoordinate(c.v[0],c.v[1] ,c.v[2],X,Y)] = max_field_val + 0.0001; //+ 1.0;//(float)97.0; //if large value used false negative gradient values and false 0 gradient values placed inside ligament for some reason.
          }
        }
        marker_base = 5;
        marker_top = 5;
        // distance along line we can travel to set markers
        while( steps - marker_base < 0)
          marker_base += 1;
        while( steps + marker_top > l.size())
          marker_top -= 1;
        while( steps + marker_top > X-1)
          marker_top -= 1;
        while( steps + marker_top > Y-1)
          marker_top -= 1;
        while( steps + marker_top > Z-1)
          marker_top -= 1;
        for(int bb = steps - marker_base; bb < steps + marker_top; bb++){
          if( bb>0 && bb < l.size() ){
            xyzi c = l[bb];
            if(marked_dist_field[LinearIndexFromCoordinate(c.v[0],c.v[1] ,c.v[2],X,Y)] > 0)
              marked_dist_field[LinearIndexFromCoordinate(c.v[0],c.v[1] ,c.v[2],X,Y)] = max_field_val + 0.0001; //+ 1.0;//(float)97.0; //if large value used false negative gradient values and false 0 gradient values placed inside ligament for some reason.
          }
        }
        // compute cotangent and tangent of normal
        uint8_t m = fabs(nz[0]) < fabs(nz[1]) ? 0 : 1;
        m = fabs(nz[m]) < fabs(nz[2]) ? m : 2;
        
        std::vector<double> w = {0,0,0};
        w[m] = 1;
        
        //tan vector to normal
        std::vector<double> vx = Cross(w, nz);
        
        vx = Unit(vx);
        
        //cotan vector to normal
        std::vector<double> vy = Cross(nz, vx);
        
        vy = Unit(vy);
        
        //vtkMatrix4x4* SliceMatrix;

        //#Set inital state of matrix to identity

        SliceMatrix = vtkMatrix4x4::New();
        SliceMatrix->Identity();
        
        //#x axis column
        SliceMatrix->SetElement(0, 0, vx[0]);
        SliceMatrix->SetElement(1, 0, vx[1]);
        SliceMatrix->SetElement(2, 0, vx[2]);
        SliceMatrix->SetElement(3, 0, 0.0);
        
        //#y axis column
        SliceMatrix->SetElement(0, 1, vy[0]);
        SliceMatrix->SetElement(1, 1, vy[1]);
        SliceMatrix->SetElement(2, 1, vy[2]);
        SliceMatrix->SetElement(3, 1, 0.0);
        
        //#z axis column
        SliceMatrix->SetElement(0, 2, nz[0]);
        SliceMatrix->SetElement(1, 2, nz[1]);
        SliceMatrix->SetElement(2, 2, nz[2]);
        SliceMatrix->SetElement(3, 2, 0.0);

        
        xyzi unit_p = Unit(p1);
        int bounds_x = WINDOW_SIZE;
        int bounds_y = WINDOW_SIZE;
        int bounds_z = WINDOW_SIZE;
        std::vector<double> bbox = bounding_box(p1, bounds_x, bounds_y, bounds_z);

        xyzi scaled_p = unit_p;
        scaled_p.v[0] = bbox[0]; //(p1.v[0] - bbox[0]); //min of bounding box around vector p
        scaled_p.v[1] = bbox[2]; //(p1.v[1] - bbox[2]); // multiplied by unit of p to give box of size 'bounds'
        scaled_p.v[2] = bbox[4]; //(p1.v[2] - bbox[4]); // around p
        
        //#origin column // is this the correct one?
        SliceMatrix->SetElement(0, 3, scaled_p.v[0]);
        SliceMatrix->SetElement(1, 3, scaled_p.v[1]);
        SliceMatrix->SetElement(2, 3, scaled_p.v[2]);
        SliceMatrix->SetElement(3, 3, 1.0);
        
        //std::vector<double> rotate_p = RotateX({p1.v[0], p1.v[1], p1.v[2]});
        //double max_radius = PerimeterSearch(p1, toXYZI(nz,97), toXYZI(vx, 97), 1.0);
        //std::cout << "og_p: "<< p1.v[0] <<" "<<p1.v[1]<< " "<<p1.v[1]<<std::endl;
        //std::cout << "scaled+p: "<< scaled_p.v[0] <<" "<<scaled_p.v[1]<< " "<<scaled_p.v[1]<<std::endl;
        //std::cout << "bbox: " << bbox[0] <<" "<<bbox[1] <<" "<<bbox[2] <<" "<<bbox[3] <<" "<<bbox[4] <<" "<<bbox[5] <<std::endl;
        
        crosssection cs;
        
        //  if(line_it < 20) {
        
        create_vtkFunction(scaled_p, p1, cube, bbox);
        
        // Compute slice of cube given coordinated defining perpendicularly intersecting
        // plane to interior line

        compute_slice(cube, SliceMatrix, vtkTransform, vtkSlice);
        
        //bool save_image_slice = only_center==1;
        bool save_image_slice = 1;
        if(save_image_slice){
	  if(lig_id % 2 == 0 ){
          vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
          std::string f_name = folder_path+"/Slice_LIGG"+std::to_string(lig_id)+".vti";
          writer->SetFileName(f_name.c_str());
          writer->SetInputData(vtkSlice);
          writer->Write();
		}
	}

        //vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
        //writer->SetFileName("slice.vti");
        //writer->SetInputData(vtkSlice);
        //writer->Write();
        
        
        cs = cross_section(vtkSlice, scaled_p, lig_id, bbox, save_image_slice, folder_path);        
        cs_area_values->InsertComponent(line_point_id + p, 0, cs.area);//#l[p_id]
        cs_perimeter_values->InsertComponent(line_point_id + p, 0, cs.perimeter);

        if(p+step+step >= l.size()){//bc iteration stops when p+step > l.size()
            for( int tail = p+1; tail< l.size(); tail = tail + 1){
                cs_area_values->InsertComponent(line_point_id + tail, 0, cs.area);//#l[p_id]
                cs_perimeter_values->InsertComponent(line_point_id + tail, 0, cs.perimeter);
            }
        }

        
        if(cs.area == 0){
          cout << " .........................##############$$$$%%%%%%%%%%" << endl;
          cout << "norm: "<< nz[0] <<" "<< nz[1]<< " " << nz[2] << endl;
          cout <<" p "<< p << endl;
          cout << " .........................##############$$$$%%%%%%%%%%" << endl;
        }
        /*
         vtkSmartPointer<vtkXMLImageDataWriter> writer2 = vtkSmartPointer<vtkXMLImageDataWriter>::New();
         writer2->SetFileName("crop.vti");
         writer2->SetInputData(rawImage);
         writer2->Write();
         
         printf(" Done!!\n");
         */

      }//end conditional iterating over points

    }//end intereating over points in line

    line_point_id+=l.size(); //increment by line length to start indexing of next point on next lig

  }//end iterating over all lines

  printf("Last point in vtp file, line_point_id = %lld\n", line_point_id);
    //all lines iterated save vtk polydata
    updateVTKLigaments(polyData, string(lines_filename) + "_cross_section.vtp");
    //else
    cout << "Poly data with cross section writen to " << string(lines_filename) << "_cross_section.vtp" <<endl;
#endif
    return 1;
  }
  
