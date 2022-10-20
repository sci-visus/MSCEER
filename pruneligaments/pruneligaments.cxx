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
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
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

// For if connected components need to be computed from raw
// e.g. raw field passed and not connected comp file
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

typedef GInt::RegularGrid3D GridType;

using namespace GInt;

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
bool has_length=false;

//for retro compatability to
// avoid passing ligaments to write
// of prune lig
bool retro_orientation = false;

#ifdef VTK_ENABLED

void print(std::string exp){
    std::cout << exp <<std::endl;
}

void ReadVTKLigaments(char* filename, std::vector<xyzi>& pointset,
                      std::map<int,std::vector<xyzi> >& lines, std::map<int,int>* point_to_line_map=NULL) {

  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename);
  reader->Update();

  vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
    std::string f(filename);

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
    retro_orientation = true;
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
  has_length = pointdata->GetArray("Length");

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

//
// pruned observed boundary ligs
void PruneVTKLigaments(char* filename, std::vector<xyzi>& pointset,
                      std::map<int,std::vector<xyzi> >& lines,
                       std::vector<int> pruned_ligaments,
                       std::map<int,int> kept_lig_id_remap) { // std::map<int,int>* point_to_line_map=NULL,

    //instantiate polydata and smart pointers to keep
    // only ligaments not pruned
    vtkSmartPointer<vtkPolyData> polyData_pruned = vtkSmartPointer<vtkPolyData>::New();

    //area
    vtkSmartPointer<vtkDoubleArray> cs_area_values = vtkSmartPointer<vtkDoubleArray>::New();
        cs_area_values->SetNumberOfComponents(1);
        cs_area_values->SetName("CrossArea");
    //perimeter
    vtkSmartPointer<vtkDoubleArray> cs_perimeter_values = vtkSmartPointer<vtkDoubleArray>::New();
        cs_perimeter_values->SetNumberOfComponents(1);
        cs_perimeter_values->SetName("CrossPerimeter");
    //cells consist of polylines
    vtkSmartPointer<vtkCellArray> pruned_lines = vtkSmartPointer<vtkCellArray>::New();
    //points
    vtkSmartPointer<vtkPoints> points_mapped = vtkSmartPointer<vtkPoints>::New();
    std::vector<xyzi> prune_points;
    //ids
    vtkSmartPointer<vtkIntArray> ids_pruned = vtkSmartPointer<vtkIntArray>::New();
        ids_pruned->SetNumberOfComponents(1);
        ids_pruned->SetName("Id");
    //attributes
    vtkSmartPointer<vtkIntArray> valids_pruned =   vtkSmartPointer<vtkIntArray>::New();         //valid
        valids_pruned->SetNumberOfComponents(1);
        valids_pruned->SetName("Valid");
    vtkSmartPointer<vtkDoubleArray> length_pruned = vtkSmartPointer<vtkDoubleArray>::New();     //length
        length_pruned->SetNumberOfComponents(1);
        length_pruned->SetName("Length");
    vtkSmartPointer<vtkDoubleArray> curvature_pruned = vtkSmartPointer<vtkDoubleArray>::New();  //curvature
        curvature_pruned ->SetNumberOfComponents(1);
        curvature_pruned ->SetName("Curvature");
    vtkSmartPointer<vtkDoubleArray> orientation_pruned[3];                                      //orientation
        orientation_pruned[0]=vtkSmartPointer<vtkDoubleArray>::New();
        orientation_pruned[1]=vtkSmartPointer<vtkDoubleArray>::New();
        orientation_pruned[2]=vtkSmartPointer<vtkDoubleArray>::New();
        if(retro_orientation){
            orientation_pruned[0] ->SetNumberOfComponents(1);
            orientation_pruned[0]->SetName("Orientation");
            orientation_pruned[1]=orientation_pruned[0];
            orientation_pruned[2]=orientation_pruned[0];
        }
        else{
            orientation_pruned[0] ->SetNumberOfComponents(1);
            orientation_pruned[0]->SetName("Orientation-X");
            orientation_pruned[1] ->SetNumberOfComponents(1);
            orientation_pruned[1]->SetName("Orientation-Y");
            orientation_pruned[2] ->SetNumberOfComponents(1);
            orientation_pruned[2]->SetName("Orientation-Z");
        }

  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename);
  reader->Update();

  //original, non-pruned attributes
  vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

  //points
  vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
  //lines (cells)
  vtkSmartPointer<vtkCellArray> cells = polyData->GetLines();
  //point attributes
  vtkSmartPointer<vtkPointData> pointdata = polyData->GetPointData();
  //id
  vtkSmartPointer<vtkIntArray> ids = vtkIntArray::SafeDownCast(pointdata->GetArray("Id"));
  //valid
  vtkSmartPointer<vtkIntArray> valids;
  if(pointdata->GetArray("Valid"))
    valids = vtkIntArray::SafeDownCast(pointdata->GetArray("Valid"));
  else
    valids = vtkSmartPointer<vtkIntArray>::New();
  //length
  vtkSmartPointer<vtkDoubleArray> lengths;
  if(pointdata->GetArray("Length"))
    lengths = vtkDoubleArray::SafeDownCast(pointdata->GetArray("Length"));
  else
    lengths=vtkSmartPointer<vtkDoubleArray>::New();
  //orientation (with retro compatibility)
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
  // cross sec area
  vtkSmartPointer<vtkDoubleArray> cross_area;
  if(has_area) cross_area = vtkDoubleArray::SafeDownCast(pointdata->GetArray("CrossArea"));
  // cross sec perim
  vtkSmartPointer<vtkDoubleArray> cross_perimeter;
  if(has_perimeter) cross_perimeter = vtkDoubleArray::SafeDownCast(pointdata->GetArray("CrossPerimeter"));
  //curvature
  vtkSmartPointer<vtkDoubleArray> curvature;
  if(has_curvature) curvature = vtkDoubleArray::SafeDownCast(pointdata->GetArray("Curvature"));

  //begin iterating and updating
  unsigned long point_id=0;
  unsigned int point_count = 0;

  polyData->GetLines()->InitTraversal();
  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();

  //collect ligament attributes of only non-pruned ligs
  while(polyData->GetLines()->GetNextCell(idList)){

    int npoints = idList->GetNumberOfIds();
    std::cout << "Line has " << npoints << " points." << std::endl;

    if(npoints==0)
      continue;

    int id = ids->GetValue(idList->GetId(0)); //line id
    int mapped_id = kept_lig_id_remap[id]; //new id decremented by running some of pruned lig

    // Check if id was pruned
    std::vector<int>::iterator it = std::find(pruned_ligaments.begin(), pruned_ligaments.end(), id);

    //skip pruned add lig/points otherwise with new ids
    if(it == pruned_ligaments.end()){

        //set new lig
        vtkSmartPointer<vtkPolyLine> mapped_polyLine = vtkSmartPointer<vtkPolyLine>::New();
        mapped_polyLine->GetPointIds()->SetNumberOfIds(npoints);

        if(lines.count(id) > 0)
          continue;

        // add attributes to mapped polyline's points
        unsigned int lp = 0;
        for(vtkIdType pointId = 0; pointId < npoints; pointId++){

          // add point to lig
          vtkIdType pid = idList->GetId(pointId);
          double *p = points->GetPoint(pid);
          double point[3] = {p[0], p[1], p[2]};

          points_mapped->InsertNextPoint(point);

          mapped_polyLine->GetPointIds()->SetId(lp, point_count ); //mapped point in lig id: mapped_id + point_count
          ids_pruned->InsertComponent(point_count, 0, mapped_id); //lig id: mapped_id
                                                                              //point id: lp
          // add attributes to first point of lig
          //if(pointId == 0) {
            unsigned int valid = 0;
            double length = 0;
            double orient[3];

            if(pointdata->GetArray("Valid"))
              valid = valids->GetValue(pid);
            valids_pruned -> InsertComponent(point_count, 0, valid);
            if(pointdata->GetArray("Length")){
                length = lengths->GetValue(pid);
                std::cout << "LENGTH " << length << std::endl;
            }
            length_pruned->InsertComponent(point_count, 0, length);
            for (int oo = 0; oo < 3; oo++) orient[oo] = orientation[oo]->GetValue(pid);
            for (int oo = 0; oo < 3; oo++)
                orientation_pruned[oo]->InsertComponent(point_count, 0, orient[oo]);
            if(pointdata->GetArray("CrossArea")){
              double area=cross_area->GetValue(pid);
              cs_area_values ->InsertComponent(point_count, 0 , area);
            }
            if(pointdata->GetArray("CrossPerimeter")) {
              double perimeter=cross_perimeter->GetValue(pid);
              cs_perimeter_values -> InsertComponent(point_count, 0 , perimeter);
            }
            if(pointdata->GetArray("Curvature")) {
              double curvature_v=curvature->GetValue(pid);
              curvature_pruned -> InsertComponent(point_count , 0, curvature_v);
            }
          //}// end adding attributes to first point


          //add point to polyline

          lp++;
          point_id++;
          point_count++;
    } //end iteration over points in line

    pruned_lines->InsertNextCell(mapped_polyLine);

  }
}
    // Add the points to the dataset
    polyData_pruned->SetPoints(points_mapped);

    // Add the lines to the dataset
    polyData_pruned->SetLines(pruned_lines);

    //add ids
    polyData_pruned->GetPointData()->AddArray(ids_pruned);

    //add attr
    if(pointdata->GetArray("Valid"))
      polyData_pruned->GetPointData()->AddArray(valids_pruned);

    if(pointdata->GetArray("CrossArea"))
      polyData_pruned->GetPointData()->AddArray(cs_area_values);

    if(pointdata->GetArray("CrossPerimeter"))
      polyData_pruned->GetPointData()->AddArray(cs_perimeter_values);

    if(pointdata->GetArray("Curvature"))
      polyData_pruned->GetPointData()->AddArray(curvature_pruned);

    if(pointdata->GetArray("Length"))
        polyData_pruned->GetPointData()->AddArray(length_pruned);

    if(pointdata->GetArray("Orientation-Z")){
        polyData_pruned->GetPointData()->AddArray(orientation_pruned[0]);        //x
        polyData_pruned->GetPointData()->AddArray(orientation_pruned[1]);        //y
        polyData_pruned->GetPointData()->AddArray(orientation_pruned[2]);        //z
    }
    else
        polyData_pruned->GetPointData()->AddArray(orientation_pruned[0]);

    //file name
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(polyData_pruned);
#else
  writer->SetInputData(polyData_pruned);
#endif

    std::string filename_(filename);
    writer->SetFileName((filename_+".pruned.vtp").c_str());
    vtkSmartPointer<vtkPolyData> polyData_out = vtkSmartPointer<vtkPolyData>::New();
    //polyData_out->DeepCopy(polyData_pruned);

    writer->Write();

//    //testing:

//    // check write
//    vtkSmartPointer<vtkXMLPolyDataReader> reader_check = vtkSmartPointer<vtkXMLPolyDataReader>::New();
//    reader_check->SetFileName((filename_+".pruned.vtp").c_str());
//    reader_check->Update();

//    //original, non-pruned attributes
//    vtkSmartPointer<vtkPolyData> polyData_check = reader_check->GetOutput();

//    //points
//    vtkSmartPointer<vtkPoints> points_check = polyData_check->GetPoints();
//    //lines (cells)
//    vtkSmartPointer<vtkCellArray> cells_check = polyData_check->GetLines();

//    std::cout << "There were " << reader->GetNumberOfLines() << " lines." << std::endl;
//    std::cout << "There are now " << reader_check->GetNumberOfLines() << " lines." << std::endl;


//    polyData_pruned->GetLines()->InitTraversal();
//    vtkSmartPointer<vtkIdList> idList_pruned = vtkSmartPointer<vtkIdList>::New();
//    //collect ligament attributes of only non-pruned lig
//    vtkSmartPointer<vtkPointData> pointdata_check = polyData->GetPointData();
//    vtkSmartPointer<vtkDoubleArray> lengths_check = vtkDoubleArray::SafeDownCast(pointdata->GetArray("Length"));
//    while(polyData_pruned->GetLines()->GetNextCell(idList_pruned)){

//    int npoints_pruned = idList_pruned->GetNumberOfIds();

//    std::cout << "points in pruned lines " << npoints_pruned <<std::endl;
//    vtkIdType pid_check = idList->GetId(0);

//    int length_check = lengths_check->GetValue(pid_check);
//    std::cout << "length " << length_check <<std::endl;
//    }

  printf("Ligaments Pruned\n");
}

#endif

// filename, x,y,z , sx, sy, sz, ex, ey, ez, sx, sy, sz
int main(int argc, char** argv) {

  if(argc < 2){
    printf("usage: pruneligaments X Y Z <lines.vtp> <connectedcomponents[raworiginal/precomputed].raw> [threshold (if to compute connected components from original raw)]\\n");
    return 1;
  }



  // READ IN THE COMMAND LINE ARGUMENTS
  int X, Y, Z;
  std::string filename;
  if (argc < 6) { printf("Usage: X Y Z vtp_filename connectedcomponents[raworiginal/precomputed] [threshold (if to compute connected components from original raw)]\n"); return 0; }
  sscanf(argv[1], "%d", &X);
  sscanf(argv[2], "%d", &Y);
  sscanf(argv[3], "%d", &Z);
  filename = std::string(argv[4]);
  std::string concom_or_raw_filename = std::string(argv[5]);
  float threshold;
  sscanf(argv[6], "%f", &threshold);

  bool write_connected_comp;
  write_connected_comp = false;

  float* mask_field = NULL;

  mask_field = new float[X*Y*Z];

  //decide to use passed file or compute connected comp
  float* connected_comp;
  bool to_compute = true;
  if (concom_or_raw_filename.find("components.raw") != std::string::npos){
    connected_comp = new float[X*Y*Z];
    FILE* fin = fopen(argv[5], "rb");
    fread(connected_comp, sizeof(float), X*Y*Z, fin);
    fclose(fin);
    to_compute = false;
    std::cout << "Using component file" << std::endl;
  }
  // Need to compute connected comp
  if(to_compute){
    GridType* underlying_grid;
      // set up structures to navigate grid, and load the 3d image
    underlying_grid = new GridType(GInt::Vec3i(X, Y, Z), GInt::Vec3b(0, 0, 0));

    printf("loaded field function\n");

    VolumeConnectedComponents cc(underlying_grid);
    DenseLabeling<char> *maskvol = new DenseLabeling<char>(underlying_grid->NumElements());

    for(INDEX_TYPE i=0; i< underlying_grid->NumElements();i++)
      maskvol->SetLabel(i, mask_field[i] > threshold);

    cc.PerformUnionFind(maskvol);
    connected_comp = new float[underlying_grid->NumElements()];
    cc.mIDVol->ReMapIds(connected_comp);

    if(write_connected_comp){
        std::string outfilename = std::string(filename+".components.raw");
        printf("saving output to %s\n", outfilename.c_str());
        ofstream outFile;
        outFile.open(outfilename.c_str(), ios::out|ios::binary);
        outFile.write((char*)connected_comp, underlying_grid->NumElements()*sizeof(int));

    }
 }

  // read in the data
  std::vector<xyzi> line_pointset;
  std::map<int,std::vector<xyzi> > lines;

  // map old lig id to new lig id if ligament is kept
  // and based off of the number ligs so far removed.
  std::map<int, int> cleaned_lines;
  // lig ids of pruned lines
  std::vector<int > lines_pruned;

  //count removed to maintane running deduction of lig id numbers
  int removed_count = 0;

  // read input lines and populate metrics
  ReadVTKLigaments(argv[4], line_pointset, lines);                                     //read orginal


  for(auto& lm : lines){

      auto& l = lm.second;
      xyzi p_left = l.front();

      double x_left = p_left.v[0];
      double y_left = p_left.v[1];
      double z_left = p_left.v[2];

      xyzi p_right = l.back();
      double x_right = p_right.v[0];
      double y_right = p_right.v[1];
      double z_right = p_right.v[2];


      //mark lig as ending at boundry
      // could pose problem for ligs very near boundry
      // but still connected
      bool boundary_lig = false;
      for(int walk = 1; walk < 4; walk++){
          if(x_left - walk < 0 || x_left + walk > X)
              boundary_lig = true;
          if(y_left - walk < 0 || y_left + walk > Y)
              boundary_lig = true;
          if(z_left - walk < 0 || y_left + walk > Y)
              boundary_lig = true;
      }
      // if save to walk outside of lig, check at least one
      // direction stays connected.
      if( !boundary_lig){
          //walk around end and see if no connected component exists
          for (int x = 0; x < 4; x++){
            for (int y = 0; y < 4; y++){
              for (int z = 0; z < 4 ; z++){
                  if(connected_comp[  LinearIndexFromCoordinate(x_left+x, y_left, z_left, X, Y) ] < 0
                          || connected_comp[  LinearIndexFromCoordinate(x_right+x, y_right, z_right, X, Y) ] < 0)
                       boundary_lig = true;
              }
            }
           }
          }

      cleaned_lines.insert(std::pair<int,int>(lm.first , lm.first - removed_count));

      if(boundary_lig == true){
          removed_count += 1;
          lines_pruned.push_back(lm.first);
      }

  }

  // read in the data
  std::vector<xyzi> line_pointset2;
  std::map<int,std::vector<xyzi> > lines2;
  PruneVTKLigaments(argv[4], line_pointset2, lines2, lines_pruned,cleaned_lines);


  std::cout << "Removed " << lines_pruned.size() << " ligaments." << std::endl;
  printf("results saved in %s\n", std::string(std::string(argv[4])+".pruned.vtp").c_str());
  return 0;
}
