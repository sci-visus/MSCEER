/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef EXTRACT_GRAPH_H
#define EXTRACT_GRAPH_H

#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <queue>
#include <time.h>
#include <fstream>

#ifdef VTK_ENABLED
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkNamedColors.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#endif

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
#include "gi_geodesic_density.h"

#include "nanoflann.h"

// using a regular grid
typedef GInt::RegularGrid3D GridType;
typedef GInt::TopologicalRegularGrid3D MeshType;
typedef GInt::RegularGridTrilinearFunction GridFuncType;
typedef GInt::DiscreteGradientLabeling<MeshType> GradType;

// for the topology we will use the lazy version of the topology-function, so no extra
// storage of computation
typedef GInt::LazyMaximumVertexLabeling<MeshType, GridFuncType> MaxVLType;
typedef GInt::TopologicalMaxVertexMeshFunction<MeshType, MaxVLType, GridFuncType, float> TopoFuncType;

// the Morse-Smale complex type based on these prior types
typedef GInt::MorseSmaleComplexBasic<float, MeshType, TopoFuncType, GradType> MscType;

struct ligament {
		INT_TYPE saddle_id;
		INT_TYPE arc1_id;
		INT_TYPE arc2_id;
    INDEX_TYPE max_id1;
    INDEX_TYPE max_id2;
    INDEX_TYPE remapped_max_id1;
    INDEX_TYPE remapped_max_id2;
};

typedef std::pair<int, INDEX_TYPE> ligid_cellid_pair;
typedef std::pair<GInt::Vec3d, ligid_cellid_pair> point_type;
class MyPointCloud : public std::vector<point_type> {
public:
	// Must return the number of data poins
	inline size_t kdtree_get_point_count() const { return size(); }
	// Must return the dim'th component of the idx'th point in the class:
	inline double kdtree_get_pt(const size_t idx, int dim) const { return this->operator[](idx).first[dim]; }
	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const
	{
		return false;
	}


};



#ifdef VTK_ENABLED

class SkeletonGraph {
  
public:
  
  SkeletonGraph(MscType* _msc, MeshType *_topological_grid, GridFuncType* _grid_function, int* _data_dims, int* _mask_field=NULL)
  : msc(_msc), topological_grid(_topological_grid), grid_function(_grid_function), data_dims(_data_dims), mask_field(_mask_field)
  {
    
    points = vtkSmartPointer<vtkPoints>::New();
    // Create a cell array to store the lines in and add the lines to it
    cells = vtkSmartPointer<vtkCellArray>::New();
    function_values = vtkSmartPointer<vtkDoubleArray>::New();
    function_values->SetNumberOfComponents(1);
    function_values->SetName("Function");
    
    id_values = vtkSmartPointer<vtkIntArray>::New();
    id_values->SetNumberOfComponents(1);
    id_values->SetName("Id");

	geodesicdensity_values = vtkSmartPointer<vtkDoubleArray>::New();
	geodesicdensity_values->SetNumberOfComponents(1);
	geodesicdensity_values->SetName("GeoDens");


    /*
    valid_values = vtkSmartPointer<vtkIntArray>::New();
    valid_values->SetNumberOfComponents(1);
    valid_values->SetName("Valid");
    */

    length_values = vtkSmartPointer<vtkDoubleArray>::New();
    length_values->SetNumberOfComponents(1);
    length_values->SetName("Length");
    
    orientation_x_values = vtkSmartPointer<vtkDoubleArray>::New();
    orientation_x_values->SetNumberOfComponents(1);
    orientation_x_values->SetName("Orientation-X");

		orientation_y_values = vtkSmartPointer<vtkDoubleArray>::New();
		orientation_y_values->SetNumberOfComponents(1);
		orientation_y_values->SetName("Orientation-Y");

		orientation_z_values = vtkSmartPointer<vtkDoubleArray>::New();
		orientation_z_values->SetNumberOfComponents(1);
		orientation_z_values->SetName("Orientation-Z");
    
    point_id=0;
  }
  
  double distance(GInt::Vec3d a, GInt::Vec3d b) {
    double total = 0;
    for (int i = 0; i < 3; i++){
      total += (a[i] - b[i])*(a[i] - b[i]);
    }
    return sqrt(total);
  }
  
  void connecting_vector(const GInt::Vec3d a, const GInt::Vec3d b, GInt::Vec3d& vec){
    for (int i = 0; i<3; i++){
      vec[i] = a[i] - b[i];
    }
  }
  
  double orientation(const GInt::Vec3d l, GInt::Vec3d a){
    double inner_prod = l.Dot(a);
    double norm = sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);
    // divisor = norm * norm(z);
    double orient = acos(inner_prod/norm)*57.2957;
    
    if(orient >= 90)
      orient = 180-orient;
    
    return orient;
  }
  
  void populatePolyline(unsigned long& point_id, unsigned int& ligament_id, vector<GInt::Vec3d>& linegeom, vtkSmartPointer<vtkPolyLine> polyLine){
    
    if(linegeom.size() == 0)
      return;
    
    int i=0;
    
    // smooth the input line (NOTE: this function edits the original line)
    smooth(250, linegeom);
    
    double running_length = 0;
    for(int i=1; i < linegeom.size(); i++){
      GInt::Vec3d p1 = linegeom[i];
      GInt::Vec3d p0 = linegeom[i-1];
      
      running_length += distance(p0, p1);
    }
    
    GInt::Vec3d conn_vector;
    connecting_vector(linegeom[0], linegeom.back(), conn_vector);
    
    double orientation_x = orientation(conn_vector,  GInt::Vec3d(1,0,0));
		double orientation_y = orientation(conn_vector,  GInt::Vec3d(0,1,0));
		double orientation_z = orientation(conn_vector,  GInt::Vec3d(0,0,1));
    

		printf("about to do populate lig geodesic, size=%d\n", linegeom.size());

		//id_values->Resize(linegeom.size());
    // crawling the ligament
    for(auto pl: linegeom){
      // position of point in ligament
      double p[3]={pl[0],pl[1],pl[2]};
           
      points->InsertNextPoint(p);
	  //points->Print(cout);
	  //--printf("doing "); pl.PrintFloat();

	  //-- printf("a %lld, %d\n", point_id, ligament_id);
      // setting id of the point inside the line
      polyLine->GetPointIds()->SetId(i,point_id);
	  //-- printf("b\n");

	  //--id_values->Print(cout);
      // setting ligament Id into array
	  //id_values->InsertComponent(point_id, 0, ligament_id);
	 //-- printf("pointid %lld, ligid %d\n", point_id, ligament_id);
	  double lig_id = ligament_id;
	  vtkIdType llpoint_id = point_id;
	  id_values->InsertNextValue(ligament_id);
	  //id_values->InsertComponent(llpoint_id, 0, lig_id);
	 //-- id_values->Print(cout);
	  //--printf("c\n");

      // interpolate value in the grid function
      double fvalue = grid_function->InterpolatedValue(pl);
	  //--function_values->Print(cout);
      //function_values->InsertComponent(llpoint_id, 0, fvalue);
	  function_values->InsertNextValue(fvalue);
	  
	  //-- function_values->Print(cout);
	  //--printf("d\n");

	  // do geodesic density here
	  GInt::Vec3i ipos(pl);
	  float dist_val = fmin(fvalue * 3, 20.0);
	  float density = 1;//geo_computer->GeodesicDensity(dist_val, ipos, m_data_box);
	  //geodesicdensity_values->InsertComponent(point_id, 0, density);
	  geodesicdensity_values->InsertNextValue(density);
	  //--printf("e\n");


      //length_values->InsertComponent(point_id, 0, running_length);
	  length_values->InsertNextValue(running_length);

	  //--printf("f\n");

      //orientation_x_values->InsertComponent(point_id, 0, orientation_x);
	//		orientation_y_values->InsertComponent(point_id, 0, orientation_y);
	//		orientation_z_values->InsertComponent(point_id, 0, orientation_z);
			orientation_x_values->InsertNextValue(orientation_x);
			orientation_y_values->InsertNextValue(orientation_y);
			orientation_z_values->InsertNextValue(orientation_z);
	  //--printf("g\n");

      i++;
      point_id++;
    }
	printf("done -- here\n");
  }
  
  bool isPointOnLigament(GInt::Vec3d vd){

    // check if point is masked (in a pancake region)                                                                      
    if(mask_field!=NULL){
      double p[3] = {vd[0],vd[1],vd[2]};
      bool invalid_point = (mask_field[linearIndexFromCoordinate(p, data_dims)] >= 0);
      if(invalid_point) return false;
    }
    return true;
  }

  GInt::GeodesicDensityComputer* geo_computer;
  GInt::GeodesicDensityComputer::DataBox* m_data_box;
  float m_max_dist;
  bool m_use_dist_radius;

  void SetUseDistRadius(bool value) { m_use_dist_radius = value; }

  void addNonOverlappingLigaments(vector<ligament>& ligaments, TopoFuncType* topological_grid_function) {
	  
	  // geodesic distance sampler
	  float* dist_field = this->grid_function->GetImage();
	  INDEX_TYPE num_elems = this->grid_function->GetGrid()->NumElements();
	  m_max_dist = dist_field[0];
	  for (INDEX_TYPE i = 1; i < num_elems; i++) {
		  if (m_max_dist < dist_field[i]) m_max_dist = dist_field[i];
	  }
	  if (m_max_dist > 20) m_max_dist = 20;
	  printf("using max_dist = %f\n", m_max_dist);

	  auto XYZ = this->grid_function->GetGrid()->XYZ();
	  geo_computer = new GInt::GeodesicDensityComputer(XYZ, dist_field);
	  geo_computer->initialize(m_max_dist, 0.0000, -0.1);
	  m_data_box = geo_computer->make_box();
	  
	  unsigned int ligament_id = 0;

	  map<INDEX_TYPE, vector<int>> conn_map;


	  printf("adding lig poitns to kd tree...\n");
	  MyPointCloud lig_points;
	  for (auto l : ligaments) {
		  // HERE WE FILL THE GEOMETRIC COORDINATES OF AN ARC
		  vector<INDEX_TYPE> linegeomids1;
		  vector<INDEX_TYPE> linegeomids2;
		  msc->fillArcGeometry(l.arc1_id, linegeomids1);
		  msc->fillArcGeometry(l.arc2_id, linegeomids2);
		 
		  for (auto id : linegeomids1) {
			  GInt::Vec3l vl;
			  topological_grid->cellid2Coords(id, vl);
			  GInt::Vec3d vd = vl * 0.5; // topological index space goes from 0 : 2X-1 , not 0 : X
			  lig_points.push_back({ vd, {l.saddle_id, id} });
		  }
		  for (auto id : linegeomids2) {
			  GInt::Vec3l vl;
			  topological_grid->cellid2Coords(id, vl);
			  GInt::Vec3d vd = vl * 0.5;
			  lig_points.push_back({ vd,{ l.saddle_id, id } });
		  }
	  }
	  printf("%d points added\n", lig_points.size());
	  printf("building kd tree...\n");

	  // construct a kd-tree index:
	  typedef nanoflann::KDTreeSingleIndexAdaptor<
		  nanoflann::L2_Simple_Adaptor<double, MyPointCloud >,
		  MyPointCloud,
		  3 /* dim */
	  > my_kd_tree_t;

	  my_kd_tree_t   index(3 /*dim*/, lig_points, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	  index.buildIndex();

	  printf("making ligaments...\n");
	  
	  // for results of knn search
	  double rad;
	  vector<pair<size_t, double>> results;
	  results.reserve(200);
	  nanoflann::SearchParams sp;
	  size_t num_trimmed = 0;
	  size_t num_accepted = 0;

	  unordered_set<INDEX_TYPE> super_ignore_set;
	  // first loop go through all arcs and mark every point 
	  for (auto l : ligaments) {
		  // HERE WE FILL THE GEOMETRIC COORDINATES OF AN ARC
		  vector<INDEX_TYPE> linegeomids1;
		  vector<INDEX_TYPE> linegeomids2;
		  msc->fillArcGeometry(l.arc1_id, linegeomids1);
		  msc->fillArcGeometry(l.arc2_id, linegeomids2);

		  std::vector<GInt::Vec3d> linegeom1;
		  std::vector<GInt::Vec3d> linegeom2;
		  std::vector<bool> valid1;
		  std::vector<bool> valid2;

		  for (auto id : linegeomids1) {
			  GInt::Vec3l vl;
			  topological_grid->cellid2Coords(id, vl);
			  GInt::Vec3d vd = vl * 0.5; // topological index space goes from 0 : 2X-1 , not 0 : X
			// do check here
			  if (m_use_dist_radius) {
				  rad = topological_grid_function->cellValue(id); // get radius along line
			  }
			  else {
				  rad = 0.2;
			  }
			  auto num_res = index.radiusSearch(vd.m_v, rad, results, sp);
			  bool has_different = false;
			  for (int i = 0; i < num_res; i++) {
				  if (lig_points[results[i].first].second.first != l.saddle_id) {
					  has_different = true;
					  super_ignore_set.insert(lig_points[results[i].first].second.second); // the cell id of the point to ignore
				  }
			  }
			  if (has_different) super_ignore_set.insert(id);
		  }
		  for (auto id : linegeomids2) {
			  GInt::Vec3l vl;
			  topological_grid->cellid2Coords(id, vl);
			  GInt::Vec3d vd = vl * 0.5;
			  // do check here
			  if (m_use_dist_radius) {
				  rad = topological_grid_function->cellValue(id); // get radius along line
			  }
			  else {
				  rad = 0.2;
			  }			  auto num_res = index.radiusSearch(vd.m_v, rad, results, sp);
			  bool has_different = false;
			  for (int i = 0; i < num_res; i++) {
				  if (lig_points[results[i].first].second.first != l.saddle_id) {
					  has_different = true;
					  super_ignore_set.insert(lig_points[results[i].first].second.second); // the cell id of the point to ignore
				  }
			  }
			  if (has_different) super_ignore_set.insert(id);
		  }
	  }

	  for (auto l : ligaments) {
		  // HERE WE FILL THE GEOMETRIC COORDINATES OF AN ARC
		  vector<INDEX_TYPE> linegeomids1;
		  vector<INDEX_TYPE> linegeomids2;
		  msc->fillArcGeometry(l.arc1_id, linegeomids1);
		  msc->fillArcGeometry(l.arc2_id, linegeomids2);

		  std::vector<GInt::Vec3d> linegeom1;
		  std::vector<GInt::Vec3d> linegeom2;
		  std::vector<bool> valid1;
		  std::vector<bool> valid2;
		  for (auto id : linegeomids1) {
			  GInt::Vec3l vl;
			  topological_grid->cellid2Coords(id, vl);
			  GInt::Vec3d vd = vl * 0.5; // topological index space goes from 0 : 2X-1 , not 0 : X
			  if (super_ignore_set.count(id) == 0) {
				  linegeom1.push_back(vd);
				  num_accepted++;
			  } else {
				  num_trimmed++;
			  }
		  }
		  for (auto id : linegeomids2) {
			  GInt::Vec3l vl;
			  topological_grid->cellid2Coords(id, vl);
			  GInt::Vec3d vd = vl * 0.5;
			  if (super_ignore_set.count(id) == 0) {
				  linegeom2.push_back(vd);
				  num_accepted++;
			  }
			  else {
				  num_trimmed++;
			  }
		  }
		  // now linegeom1 and linegeom2 have the 3d coordinates of the stair-stepping line

		  // reverse the second segment and merge together to make one ligament
		  std::reverse(linegeom2.begin(), linegeom2.end());

		  vector<GInt::Vec3d> merged_geom;
		  merged_geom.insert(merged_geom.begin(), linegeom1.begin(), linegeom1.end());
		  merged_geom.insert(merged_geom.end(), linegeom2.begin(), linegeom2.end());


		  if (merged_geom.size() == 0) continue;

		  vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();

		  polyLine->GetPointIds()->SetNumberOfIds(merged_geom.size());
		  printf("adding polyline %d, size=%d\n", ligament_id, merged_geom.size());
		  populatePolyline(point_id, ligament_id, merged_geom, polyLine);

		  cells->InsertNextCell(polyLine);

		  ligament_id++;
	  }
	  printf("num_trimmed = %lld, num_accepted = %lld\n", num_trimmed, num_accepted);

  }


  void addUnclippedLigamentsOutsideBlobs(vector<ligament>& ligaments, GInt::DenseLabeling<char> *maskvol = NULL) {
	  unsigned int ligament_id = 0;

	  if (mask_field)
		  printf("adding ligaments with field mask\n");
	  else
		  printf("adding ligaments...\n");

	  map<INDEX_TYPE, vector<int>> conn_map;

	  INDEX_TYPE max_nodes = msc->numNodes();


	  for (auto l : ligaments) {

		  //// if we have a blob mask, then skip ligaments that start and end in the same blob
		  //if (maskvol) {
			 // // get the blob ids
			 // auto n1 = msc->getNode(l.max_id1);
			 // INDEX_TYPE v1 = topological_grid->VertexNumberFromCellID(n1.cellindex);
			 // long long blob_id_1 = maskvol->GetLabel(v1);

			 // auto n2 = msc->getNode(l.max_id2);
			 // INDEX_TYPE v2 = topological_grid->VertexNumberFromCellID(n2.cellindex);
			 // long long blob_id_2 = maskvol->GetLabel(v2);

			 // // skip this ligament if we start and end int eh same blob
			 // if (blob_id_1 == blob_id_2) {
				//  printf("skipping lig %d, %d\n", blob_id_1, blob_id_2);
				//  continue;
			 // }
			 // // else record the connection
			 // if (blob_id_1 != -1) {
				//  INDEX_TYPE new_id = blob_id_1 + max_nodes;
				//  conn_map[new_id].push_back(ligament_id);
				//  l.remapped_max_id1 = new_id;

			 // }
			 // if (blob_id_2 != -1) {
				//  INDEX_TYPE new_id = blob_id_2 + max_nodes;
				//  conn_map[new_id].push_back(ligament_id);
				//  l.remapped_max_id2 = new_id;
			 // }
		  //}
		  // HERE WE FILL THE GEOMETRIC COORDINATES OF AN ARC
		  vector<INDEX_TYPE> linegeomids1;
		  vector<INDEX_TYPE> linegeomids2;
		  msc->fillArcGeometry(l.arc1_id, linegeomids1);
		  msc->fillArcGeometry(l.arc2_id, linegeomids2);

		  std::vector<GInt::Vec3d> linegeom1;
		  std::vector<GInt::Vec3d> linegeom2;
		  std::vector<bool> valid1;
		  std::vector<bool> valid2;
		  for (auto id : linegeomids1) {
			  GInt::Vec3l vl;
			  topological_grid->cellid2Coords(id, vl);
			  GInt::Vec3d vd = vl * 0.5; // topological index space goes from 0 : 2X-1 , not 0 : X
			  linegeom1.push_back(vd);
		  }
		  for (auto id : linegeomids2) {
			  GInt::Vec3l vl;
			  topological_grid->cellid2Coords(id, vl);
			  GInt::Vec3d vd = vl * 0.5;

			  linegeom2.push_back(vd);
		  }
		  // now linegeom1 and linegeom2 have the 3d coordinates of the stair-stepping line

		  // reverse the second segment and merge together to make one ligament
		  std::reverse(linegeom2.begin(), linegeom2.end());

		  vector<GInt::Vec3d> merged_geom;
		  merged_geom.insert(merged_geom.begin(), linegeom1.begin(), linegeom1.end());
		  merged_geom.insert(merged_geom.end(), linegeom2.begin(), linegeom2.end());


		  if (merged_geom.size() == 0) continue;

		  vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();

		  polyLine->GetPointIds()->SetNumberOfIds(merged_geom.size());
		  printf("adding polyline %d, size=%d\n", ligament_id, merged_geom.size());
		  populatePolyline(point_id, ligament_id, merged_geom, polyLine);

		  cells->InsertNextCell(polyLine);

		  ligament_id++;
	  }

  }

  void addLigaments(vector<ligament>& ligaments, GInt::DenseLabeling<char> *maskvol=NULL){
    unsigned int ligament_id=0;
    
    if(mask_field)
      printf("adding ligaments with field mask\n");
    else
      printf("adding ligaments...\n");
    
    map<INDEX_TYPE, vector<int>> conn_map;
    
    INDEX_TYPE max_nodes = msc->numNodes();
    
    for (auto l : ligaments) {
      
      {
        auto n1 = msc->getNode(l.max_id1);
        INDEX_TYPE v1 = topological_grid->VertexNumberFromCellID(n1.cellindex);
        
        if(maskvol){
        long long blob_id = maskvol->GetLabel(v1);
        
        if(blob_id != -1){
          INDEX_TYPE new_id = blob_id+max_nodes;
          conn_map[new_id].push_back(ligament_id);
          l.remapped_max_id1=new_id;
        }
        }
      }
      {
        auto n2 = msc->getNode(l.max_id2);
        INDEX_TYPE v2 = topological_grid->VertexNumberFromCellID(n2.cellindex);
        
	if(maskvol){
        long long blob_id = maskvol->GetLabel(v2);
        
        if(blob_id != -1){
          INDEX_TYPE new_id = blob_id+max_nodes;
          conn_map[new_id].push_back(ligament_id);
          l.remapped_max_id2=new_id;
        }
        }
      }
      
      // HERE WE FILL THE GEOMETRIC COORDINATES OF AN ARC
      vector<INDEX_TYPE> linegeomids1;
      vector<INDEX_TYPE> linegeomids2;
      msc->fillArcGeometry(l.arc1_id, linegeomids1);
      msc->fillArcGeometry(l.arc2_id, linegeomids2);

      std::vector<GInt::Vec3d> linegeom1;
      std::vector<GInt::Vec3d> linegeom2;
      std::vector<bool> valid1;
      std::vector<bool> valid2;
      for (auto id : linegeomids1) {
        GInt::Vec3l vl;
        topological_grid->cellid2Coords(id, vl);
        GInt::Vec3d vd = vl * 0.5; // topological index space goes from 0 : 2X-1 , not 0 : X
        
        if(isPointOnLigament(vd))
          linegeom1.push_back(vd);
        
      }
      for (auto id : linegeomids2) {
        GInt::Vec3l vl;
        topological_grid->cellid2Coords(id, vl);
        GInt::Vec3d vd = vl * 0.5;
        
	if(isPointOnLigament(vd))
           linegeom2.push_back(vd);
        
      }
      // now linegeom1 and linegeom2 have the 3d coordinates of the stair-stepping line
      
      // reverse the second segment and merge together to make one ligament
      std::reverse(linegeom2.begin(), linegeom2.end());
        
      vector<GInt::Vec3d> merged_geom;
      merged_geom.insert(merged_geom.begin(),linegeom1.begin(), linegeom1.end());
      merged_geom.insert(merged_geom.end(),linegeom2.begin(), linegeom2.end());
      
      
      if(merged_geom.size() == 0) continue;
      
      vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
      
      polyLine->GetPointIds()->SetNumberOfIds(merged_geom.size());
      
      populatePolyline(point_id, ligament_id, merged_geom,polyLine);
      
      cells->InsertNextCell(polyLine);
      
      ligament_id++;
    }
    
  }
  
  void save(std::string filename){
    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	printf("about to set points\n");
    // Add the points to the dataset
	points->Print(cout);
	//polyData->AllocatePointGhostArray();
    polyData->SetPoints(points);
    
	printf("b\n");
	// Add the lines to the dataset
    polyData->SetLines(cells);
	printf("c\n");

    polyData->GetPointData()->AddArray(id_values);
    polyData->GetPointData()->AddArray(function_values);
	polyData->GetPointData()->AddArray(geodesicdensity_values);
    //if(mask_field!=NULL)
    //  polyData->GetPointData()->AddArray(valid_values);
    
    polyData->GetPointData()->AddArray(length_values);
    polyData->GetPointData()->AddArray(orientation_x_values);
		polyData->GetPointData()->AddArray(orientation_y_values);
		polyData->GetPointData()->AddArray(orientation_z_values);
		printf("d\n");

    
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	printf("e\n");
	writer->SetInputData(polyData);
    writer->SetFileName(filename.c_str());
    writer->Write();
	printf("wrote and done!\n");
  };
  
private:
  MscType* msc;
  MeshType *topological_grid;
  GridFuncType* grid_function;
  
  vtkSmartPointer<vtkPoints> points;
  vtkSmartPointer<vtkCellArray> cells;
  
  vtkSmartPointer<vtkDoubleArray> function_values;
  vtkSmartPointer<vtkIntArray> id_values;
  //vtkSmartPointer<vtkIntArray> valid_values;
  vtkSmartPointer<vtkDoubleArray> length_values;
  vtkSmartPointer<vtkDoubleArray> orientation_x_values;
  vtkSmartPointer<vtkDoubleArray> orientation_y_values;
  vtkSmartPointer<vtkDoubleArray> orientation_z_values;
  vtkSmartPointer<vtkDoubleArray> geodesicdensity_values;

  int* data_dims;
  unsigned long point_id;
  int* mask_field;
  
  void smooth(int times, vector<GInt::Vec3d>& line) const{
	vector<GInt::Vec3d> oline = line;
    for (int k = 0; k < times; k++) {
      vector<GInt::Vec3d> copyline = line;
      for (int i = 1; i < line.size() - 1; i++) {
        copyline[i] = line[i - 1]* 0.25 + line[i]* 0.5 + line[i + 1] * 0.25;
		for (int cc = 0; cc < 3; cc++) {
			auto& v = copyline[i];
			auto& v2 = oline[i];
			if (v[cc] - v2[cc] > 1) v[cc] = v2[cc] + 1.0;
			if (v2[cc] - v[cc] > 1) v[cc] = v2[cc] - 1.0;
		}
      }
      line = copyline;
    }
    
  }
  
  inline INDEX_TYPE linearIndexFromCoordinate(const double* coords, int* dims) const {
    return ((INDEX_TYPE)coords[0])
    + ((INDEX_TYPE)coords[1]) * dims[0]
    + ((INDEX_TYPE)coords[2]) * dims[0] * dims[1];
  }
  
};

#endif // if VTK enabled

#endif
