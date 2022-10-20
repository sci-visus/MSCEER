
#include <vector>
#include <map>
#include <string>
#include <mutex>
#include <functional>

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
#include "gi_bifiltration_pairing.h"
#include "gi_morse_smale_complex_basic.h"
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
 * currently iterates over 20 points centered at middle of ligament
 * to identify cross sectional slice of interest points in interior
 * of ligament marked in order to identify correct connected component
 * to sum interior points for area.
 * Boundary cases pose issue if marked points not found. results in large
 * incorrect area calculation.
 */
typedef GInt::RegularGrid3D GridType;

#define COMPUTE_CROSS 0
#define COMPUTE_CURV  1

using namespace GInt;

using namespace std;

int X, Y, Z;
int lig_boundary_x, lig_boundary_y, lig_boundary_z;

vtkSmartPointer<vtkDoubleArray> lig_curvature_values = vtkSmartPointer<vtkDoubleArray>::New();

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
  //polyData_out->GetPointData()->RemoveArray("CrossArea");
  //polyData_out->GetPointData()->RemoveArray("CrossPerimeter");
  polyData_out->GetPointData()->RemoveArray("Curvature");
  
  // Add updated arrays to the file
  //polyData_out->GetPointData()->AddArray(cs_area_values);
  //polyData_out->GetPointData()->AddArray(cs_perimeter_values);
  polyData_out->GetPointData()->AddArray(lig_curvature_values);
  
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(polyData_out);
  writer->SetFileName(filename_out.c_str());
  writer->Write();
  
  printf("Ligaments updated\n");
};


#endif



// filename, x,y,z , sx, sy, sz, ex, ey, ez, sx, sy, sz
int main(int argc, char** argv) {
  
  if(argc < 7){
    printf("usage: extractcurvature <lines.vtp>\n");
    return 1;
  }
  printf("start curvature\n");
  
  // read in the data
  std::vector<xyzi> line_pointset;
  std::map<int,std::vector<xyzi> > lines;
  std::map<int,std::vector<xyzi> > lines_test;
  
  char* lines_filename = argv[1];
  
  // read input into pointsets and lines
  vtkSmartPointer<vtkPolyData> polyData = ReadVTKLigaments(lines_filename, line_pointset, lines);
  
  vtkSmartPointer<vtkPointData> pointdata = polyData->GetPointData();
  vtkSmartPointer<vtkIntArray> ids = vtkIntArray::SafeDownCast(pointdata->GetArray("Id"));
  
  lig_curvature_values->SetNumberOfComponents(1);
  lig_curvature_values->SetName("Curvature");
  
  int line_it = 0;
  
  unsigned long line_point_id=0;

    // sum cosine similarity of normal vectors to perpendicularly intersecting
    // planes to ligament

  cout << "---- === Compute curvature === ----" << endl;


  for (auto& line : lines) {
	  cout << "-----New Lig, ID: " << line.first << endl;

	  std::vector<xyzi> l = line.second;

	  smooth(2000, l);

	  int lig_id = line.first;
	  int ten_p = int(l.size()*0.1);
	  int normal_step = 4;
	  double avg_curv = 0;
	  int pcount = 0;
	  int pstep = 10;


		  // PSEUDOCODE FOR CURVATURE MEASUREMENT HERE
		  // start loop: for (int p = pstep; p < l.size()-pstep; p+= pstep) {
		  // A = L[p-pstep], B = L[p], C = L[p+pstep]
		  // we want to accumulate angle between (B-A) and (C-B):
		  //	theta = arccos( Dot( (B-A).Normalize(), (C-B).Normalize() ))
		  // count the number of points actually seen
		  // after for loop divide by number of points seen. 
		  // done
	  double curvature=0;
	  for (int p = pstep; p < l.size()-pstep; p+= pstep) {
	    xyzi a = l[p-pstep];
	    xyzi b = l[p];
	    xyzi c = l[p+pstep];
	    
            double curve = acos(Norm(b-a)*Norm(c-b));
            if (::isnan(curve) == false)
	      curvature += curve;

		  pcount++;
		  cout << "Curvature: " << curvature << endl;
		  //}
	  }

	  avg_curv = curvature;///pcount;

	  // write average value on all points 
	  for (int p = 0; p < l.size(); p++) {
		  lig_curvature_values->InsertComponent(l[p].pid, 0, avg_curv);//line_point_id+p, 0, avg_curv);
	  }

	  line_point_id += l.size();

  }


  // boolean for saving results: save_results should we add a flag here?
  updateVTKLigaments(polyData, string(lines_filename) + "_curve.vtp");

  return 1;
}

