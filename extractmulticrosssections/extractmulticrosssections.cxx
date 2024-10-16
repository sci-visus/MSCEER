
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
static float* dist_field_slice;


// crosssection atribute list to fill


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




int main(int argc, char** argv) {
  
  if(argc < 9){
    printf("usage: extractmulticrosssections <X> <Y> <Z> <lines.vtp> <dist_field.raw> <float sample rate> <max slice radius> <output directory>\n");
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
  
  float sample_rate = 0;
    sscanf(argv[6], "%f", &sample_rate);
	int max_radius;
	sscanf(argv[7], "%d", &max_radius);

  dist_field = new float[X*Y*Z];

  
  int num_elems=X*Y*Z;
  
  FILE* fin = fopen(argv[5], "rb");
  fread(dist_field, sizeof(float), num_elems, fin);
  fclose(fin);
  
  printf("start %s !\n", argv[1]);
  
  // read in the data
  std::vector<xyzi> line_pointset;
  std::map<int,std::vector<xyzi> > lines;
  
  char* lines_filename = argv[4];
  char* dir_name = argv[8];

  printf("reading file %s\n", lines_filename);

  std::string folder_path(lines_filename);


  folder_path = folder_path.substr(0, folder_path.find_last_of("/\\"));

  // read input into pointsets and lines
  vtkSmartPointer<vtkPolyData> polyData = ReadVTKLigaments(lines_filename, line_pointset, lines);
  
  vtkSmartPointer<vtkPointData> pointdata = polyData->GetPointData();
  vtkSmartPointer<vtkIntArray> ids = vtkIntArray::SafeDownCast(pointdata->GetArray("Id"));
  
  
  // Iterate over ligaments filling cube with gradient values and taking slice vtkSlice
  // orthogonal to ligament l.

  //int line_id = 0;
 // int line_point_id=0;
  GInt::RegularGrid3D* datagrid = new GInt::RegularGrid3D({ X, Y, Z }, { 0,0,0 });
  GInt::RegularGridTrilinearFunction* tri =
	  new GInt::RegularGridTrilinearFunction(datagrid, dist_field);
   
  int dist = max_radius; // how far away to look
  int grid_points = dist * 2 + 1; // this will make 1 2d grid unit = 1 3d grid unit
// re-use some code to be able to do UF on 2d grid
  GInt::RegularGrid2D* small_grid = new GInt::RegularGrid2D({ grid_points, grid_points }, { 0,0 });
  float THRESHOLD = 0;

  printf("before loop lines %d\n", lines.size());
  vector<int> ordered_line_indices;
  for (auto& it_p : lines) {
	  ordered_line_indices.push_back(it_p.first);
  }

  int num_lines = ordered_line_indices.size();

  struct line_record {
	  float area;
	  float dist_from_lowest;
	  float dist_along_line;
	  float radius;
	  float radius_of_gyration;
	  float angle;
	  float chord_displacement;
  };
  vector<line_record>* record_array =
	  new vector<line_record>[num_lines];
 
  
  // make struct for dist squared lookup
  float* dist2_to_mid = new float[grid_points * grid_points];
  GInt::Vec2f mid(int(grid_points / 2), int(grid_points / 2));
  // y moves slowest
  for (int jj = 0; jj < grid_points; jj++) {
	  // x moves fastest
	  for (int ii = 0; ii < grid_points; ii++) {
		  GInt::Vec2f pt(ii, jj);
		  float dist = (pt - mid).MagSq();
		  printf("(%d,%d)-->%f\n", ii, jj, dist);
		  dist2_to_mid[ii + jj * grid_points] = dist;
	  }
  }

#pragma omp parallel for schedule(dynamic)
  for (int line_index = 0; line_index < num_lines; line_index++) {
	  auto& it_p = *lines.find(line_index);
	  //if (line_index > 3) continue; // hack

	  std::vector<xyzi>& l = it_p.second;
	  if (l.size() < 3) continue;

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
	 // printf("lig-%d: size=%d, lowest_pos=%d", it_p.first, l.size(), lowest_pos);
	  // pre-cache distances along the line
	  std::vector<float> dist_from_start;
	  dist_from_start.push_back(0);
	  //printf("    dists: ");

	  for (int pos = 1; pos < l.size(); pos++) {
		  Vec3d pcoods(l[pos - 1].v[0], l[pos - 1].v[1], l[pos - 1].v[2]);
		  Vec3d coords(l[pos].v[0], l[pos].v[1], l[pos].v[2]);
		  float last_dist = dist_from_start[dist_from_start.size() - 1];
		  dist_from_start.push_back(last_dist + (coords - pcoods).Mag());
		 // printf("%f, ", *dist_from_start.rbegin());
	  }
	 // printf("\n");
	  

	  // get the start and end positions
	  // use the lowest_pos as the anchor (0) for distance vs. cross section area
	  std::list<int> slice_position_list;
	  std::list<double> slice_dist_list;
	  slice_position_list.push_back(lowest_pos);// put midpoint on list
	  slice_dist_list.push_back(0);
	  float dist_left = sample_rate; // initialize

	  // count back from middle
	  double dist_so_far = 0;
	  for (int pos = lowest_pos-1; pos >= 0; pos--) {
		  float tdist = dist_from_start[pos + 1] - dist_from_start[pos];
		  dist_so_far += tdist;
		  if (dist_left <= tdist) {
			  // push back
			  slice_position_list.push_front(pos);
			  slice_dist_list.push_front(- dist_so_far);
			  dist_left += sample_rate - tdist;
		  }
		  else {
			  dist_left -= tdist;
		  }
	  }
	  if (slice_position_list.front() != 0) {
		  slice_position_list.push_front(0);
		  slice_dist_list.push_front(- dist_so_far);
	  }
	  // count forward from the middle
	  dist_so_far = 0;
	  dist_left = sample_rate;
	  for (int pos = lowest_pos + 1; pos < l.size(); pos++) {
		  float tdist = dist_from_start[pos] - dist_from_start[pos-1];
		  dist_so_far += tdist;
		  if (dist_left <= tdist) {
			  // push back
			  slice_position_list.push_back(pos);
			  slice_dist_list.push_back(dist_so_far);
			  dist_left += sample_rate - tdist;
		  }
		  else {
			  dist_left -= tdist;
		  }
	  }
	  if (slice_position_list.back() != l.size() - 1) {
		  slice_position_list.push_back(l.size() - 1);
		  slice_dist_list.push_back(dist_so_far);
	  }

	  printf("lig %d: doing %d slices\n", it_p.first, slice_position_list.size());
	  std::vector<int> slice_positions;
	  std::vector<double> slice_dists;

	  slice_positions.insert(slice_positions.begin(), slice_position_list.begin(), slice_position_list.end());
	  slice_dists.insert(slice_dists.begin(), slice_dist_list.begin(), slice_dist_list.end());

	  // BEGIN VTK STUFF
	  // now we are ready to make vtk polydata for sanity checking
	  vtkNew<vtkPolyData> grid2d;
	  vtkNew<vtkPoints> points;
	  vtkNew<vtkCellArray> polys;
	  vtkNew<vtkFloatArray> scalars;
	  vtkNew<vtkUnsignedCharArray> cc_scalars;
	  scalars->SetName("dist_vals");
	  cc_scalars->SetName("cc_visited");



	  int offset_positions = 0;
	  int last_index = l.size() - 1;
	  Vec3f start_line(l[0].v[0], l[0].v[1], l[0].v[2]);
	  Vec3f end_line(l[last_index].v[0], l[last_index].v[1], l[last_index].v[2]);
	  Vec3f line_ori = (end_line - start_line); line_ori.Normalize();

	  for (int pos = 0; pos < slice_positions.size(); pos++) {
		  
		  int prev_index = max(0, pos - 1);
		  int next_index = min((int) slice_positions.size() - 1, pos + 1);

		  int line_pos = slice_positions[pos];
		  int line_prev_pos = slice_positions[prev_index];
		  int line_next_pos = slice_positions[next_index];

		  Vec3f origin(l[line_pos].v[0], l[line_pos].v[1], l[line_pos].v[2]);
		  Vec3f start_v(l[line_prev_pos].v[0], l[line_prev_pos].v[1], l[line_prev_pos].v[2]);
		  Vec3f end_v(l[line_next_pos].v[0], l[line_next_pos].v[1], l[line_next_pos].v[2]);

		  Vec3f normal = end_v - start_v;
		  normal.Normalize();
		  GInt::PlaneGrid2D slice(origin, normal, grid_points, grid_points, dist, dist);

		  // measure angle
		  // measure chord disp;
		  float chord_displacement;
		  float angle;
		  if (pos == 0 || pos == slice_positions.size() - 1) {
			  angle = 0;
			  chord_displacement = 0;
		  }
		  else {
			  // compute angle
			  auto v1 = (origin - start_v); v1.Normalize();
			  auto v2 = (end_v - origin); v2.Normalize();
			  angle = fabs(acosf(v1.Dot(v2)));

			  // compute dist to line
			  auto moved_p = origin - start_line; // move start_line to zero, origin 
			  chord_displacement = ( (line_ori * line_ori.Dot(moved_p)) - moved_p).Mag();
		  }


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
						  offset_positions + ii + jj*grid_points,
						  offset_positions + ii + 1 + jj*grid_points,
						  offset_positions + ii + 1 + (jj + 1)*grid_points,
						  offset_positions + ii + (jj + 1)*grid_points
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

		  unsigned char* cc_array = new unsigned char[small_grid->NumElements()](); // *should zero initialize																				//memset(cc_array, 0, small_grid->NumElements()*sizeof(unsigned char))
		  // also accumulate the r^2 values
		  float rad2_accum = 0;
		  for (auto id : visited) {
			  cc_array[id] = 1;
			  rad2_accum += dist2_to_mid[id];
		  }
		  float radius_of_gyration = sqrtf(rad2_accum / visited.size());

		  // Load the point, cell, and data attributes.
		  for (auto i = 0ul; i < grid_vertices.size(); ++i)
		  {
			  points->InsertPoint(offset_positions + i, grid_vertices[i].m_v);
			  scalars->InsertTuple1(offset_positions + i, values[i]);
			  cc_scalars->InsertTuple1(offset_positions + i, cc_array[i]);
		  }
		  for (auto& qu : quads)
		  {
			  polys->InsertNextCell(vtkIdType(4), qu.id);
		  }

		  line_record lr;
		  lr.area = cc_size;
		  lr.dist_from_lowest = slice_dists[pos];
		  lr.dist_along_line = lr.dist_from_lowest - slice_dists[0];
		  lr.radius = tri->TriLinInterpValue(origin);
		  lr.radius_of_gyration = radius_of_gyration;
		  lr.angle = angle;
		  lr.chord_displacement = chord_displacement;

		  record_array[line_index].push_back(lr);

		  delete[] values;
		  delete[] cc_array;

		  offset_positions += grid_vertices.size();
		 
	  }

	  // We now assign the pieces to the vtkPolyData.
	  grid2d->SetPoints(points);
	  grid2d->SetPolys(polys);
	  grid2d->GetPointData()->SetScalars(scalars);
	  grid2d->GetPointData()->AddArray(cc_scalars);
	  // now write to vtp

	  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
      std::string l_filename(lines_filename);
	  std::string dirname(dir_name);
      std::string f_name = dirname +"Multi_Slice_LIGG_"+std::to_string(line_index)+".vtp";
          writer->SetInputData(grid2d);          
          writer->SetFileName(f_name.c_str());
	  printf("writing slice to %s\n", f_name.c_str());
	  writer->Write();


	  // END VTK STUFF


	  std::string rec_name = dirname + "Multi_Slice_LIGG_" + std::to_string(line_index) + "_rec.txt";

	  FILE* fout = fopen(rec_name.c_str(), "w");
	  for (int rc = 0; rc < record_array[line_index].size(); rc++) {
		  auto& lr = record_array[line_index][rc];
		  fprintf(fout, "%d, %f, %f, %f, %f, %f, %f, %f\n", 
			  slice_positions[rc], lr.area, lr.radius, 
			  lr.dist_along_line, lr.dist_from_lowest, 
			  lr.radius_of_gyration, lr.chord_displacement, 
			  lr.angle);
	  }
	  fclose(fout);

      //line_point_id = line_point_id + l.size(); //increment point id over entire line
	  //line_id++;
  } // end for loop over all lines
  //all lines iterated save vtk polydata
  
  //else 
  cout << "Poly data with cross section writen to " << string(lines_filename) << "_cross_section.vtp" <<endl;

    return 1;
  }
  
