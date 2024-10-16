#include <stdio.h>
#include <vector>
#include <map>

#define VTK_ENABLED
#ifdef VTK_ENABLED

#include <vtkXMLPolyDataReader.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkXMLPolyDataWriter.h>
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
#include <list>
#include <algorithm>
#include <iterator>

/**
This executable is used to extract features from the vtp files 
generated by the extract_graph executable and save them into a CSV form
**/

using namespace std;

int X, Y, Z;

struct xyzi {
    double v[3];
    int i;
};
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

xyzi toXYZI(std::vector<double> p, int i){
  xyzi v;
  v.v[0] = p[0];
  v.v[1] = p[1];
  v.v[2] = p[2];
  v.i = i;
  return v;
}



struct lig_metrics{
  double length;
  double orientation[3];
  double cross_area;
  double cross_perimeter;
  double curvature;
  int valid;
};

std::map<int, lig_metrics> all_metrics;
bool has_area=false;
bool has_perimeter=false;
bool has_curvature=false;
bool has_valid=false;

#ifdef VTK_ENABLED

void ReadVTKLigaments(char* filename, std::vector<xyzi>& pointset,
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
  
  vtkSmartPointer<vtkIntArray> valids;
  has_valid=pointdata->GetArray("Valid");
  if(has_valid)
    valids = vtkIntArray::SafeDownCast(pointdata->GetArray("Valid"));
  else
    valids = vtkSmartPointer<vtkIntArray>::New();

  vtkSmartPointer<vtkDoubleArray> lengths;
  if(pointdata->GetArray("Length"))
    lengths = vtkDoubleArray::SafeDownCast(pointdata->GetArray("Length"));
  else
    lengths=vtkSmartPointer<vtkDoubleArray>::New();

  vtkSmartPointer<vtkDoubleArray> orientation[3];

  if(pointdata->GetArray("Orientation")){ // retro compatibility
    orientation[0] = vtkDoubleArray::SafeDownCast(pointdata->GetArray("Orientation"));
    orientation[1] = orientation[0];
    orientation[2] = orientation[0];
  }
  else if(pointdata->GetArray("Orientation-X")){
    orientation[0]  = vtkDoubleArray::SafeDownCast(pointdata->GetArray("Orientation-X"));
    orientation[1]  = vtkDoubleArray::SafeDownCast(pointdata->GetArray("Orientation-Y"));
    orientation[2]  = vtkDoubleArray::SafeDownCast(pointdata->GetArray("Orientation-Z"));
  } 
  else{
    orientation[0]=vtkSmartPointer<vtkDoubleArray>::New();
    orientation[1]=vtkSmartPointer<vtkDoubleArray>::New();
    orientation[2]=vtkSmartPointer<vtkDoubleArray>::New();

  }

  has_area = pointdata->GetArray("CrossArea");
  has_perimeter = pointdata->GetArray("CrossPerimeter");
  has_curvature = pointdata->GetArray("Curvature");

  vtkSmartPointer<vtkDoubleArray> cross_area;
  if(has_area) cross_area = vtkDoubleArray::SafeDownCast(pointdata->GetArray("CrossArea"));
  vtkSmartPointer<vtkDoubleArray> cross_perimeter; 
  if(has_perimeter) cross_perimeter = vtkDoubleArray::SafeDownCast(pointdata->GetArray("CrossPerimeter"));
  vtkSmartPointer<vtkDoubleArray> curvature;
  if(has_curvature) curvature = vtkDoubleArray::SafeDownCast(pointdata->GetArray("Curvature"));  

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

      if(pointId == 0) {
        int valid = 0;
        double length = 0;
        double orient[3];
          
        if(has_valid){
          valid = valids->GetValue(pid);
	}
        
        length = lengths->GetValue(pid);
        
        for (int oo = 0; oo < 3; oo++) orient[oo] = orientation[oo]->GetValue(pid);
	

        all_metrics[lid].length = length;
        all_metrics[lid].valid = valid;
        for (int oo = 0; oo < 3; oo++) all_metrics[lid].orientation[oo] = orient[oo];
	
	
        if(has_area){
          double area=cross_area->GetValue(pid);
          all_metrics[lid].cross_area=area;
        }
        
        if(has_perimeter) {
          double perimeter=cross_perimeter->GetValue(pid);
          all_metrics[lid].cross_perimeter=perimeter;
	}

        if(has_curvature) {
          double curvature_v=curvature->GetValue(pid);
          all_metrics[lid].curvature=curvature_v;
	}
      
      }

      xyzi t;
      for (int j = 0; j < 3; j++) t.v[j] = p[j];
      t.i = lid;
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

  //printf("Read complete correctly\n");
}
#endif

// filename, x,y,z , sx, sy, sz, ex, ey, ez, sx, sy, sz
int main(int argc, char** argv) {

  if(argc < 2){
    printf("usage: extraccentroids <lines.vtp>\n");
    return 1;
  }

  // read in the data
  std::vector<xyzi> line_pointset;
  std::map<int,std::vector<xyzi> > lines;

  // read input lines and populate metrics
  ReadVTKLigaments(argv[1], line_pointset, lines);

  std::string filename(argv[1]);
  ofstream outfile;
  outfile.open(filename+".ligcenters.csv");

  vtkSmartPointer<vtkPoints> points =
          vtkSmartPointer<vtkPoints>::New();

  printf("n ligs %d\n", all_metrics.size());
  
  outfile << "LID,"<<"centroid-X,"<<"centroid-Y,"<<"centroid-Z"<<std::endl;

  for(auto& lm : lines){
    auto& l = lm.second;
    xyzi centroid;
    for (int d = 0; d < 3; d++) centroid.v[d]=0;

    for(auto& p:l) {
      for (int d = 0; d < 3; d++)
        centroid.v[d] += p.v[d];
    }

    for (int d = 0; d < 3; d++) centroid.v[d]/=l.size();

    points->InsertNextPoint ( centroid.v[0],centroid.v[1],centroid.v[2]);

    outfile << lm.first << ","<< centroid.v[0]<<","<< centroid.v[1]<<","<< centroid.v[2]<< std::endl;

  }
  outfile.close();

  // Create a polydata object and add the points to it.
  vtkSmartPointer<vtkPolyData> polydata =
          vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);

  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
          vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName((filename+".ligcenters.vtp").c_str());

#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(polydata);
#else
  writer->SetInputData(polydata);
#endif
  writer->Write();

  printf("results saved in %s\n", std::string(std::string(argv[1])+".ligcenters.csv").c_str());
  return 0;
}