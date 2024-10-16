#ifdef WIN32
#include <io.h>
#include <stdio.h>
#include <stdlib.h>
#define access _access
#define F_OK 0
#else
#include <unistd.h>
#endif
#include "extractgraph.h"
#include <fstream>
using namespace GInt;

#define COMPUTE_BLOBS 0
#define USE_BLUR_FIELD_AND_DIST 1

inline unsigned long linearIndexFromCoordinate(int i, int j, int k, unsigned long M, unsigned long N) {
  return i*(N*M) + j*M + k;
}

template<typename I>
void dumpArray(I* f, int* size, string filename){
  unsigned long outsize = size[0]*size[1]*size[2];

  I* outdata = new I[outsize];
  for (int i = 0; i < size[0]; ++i)
    for (int j = 0; j < size[1]; ++j) 
      for (int k = 0; k < size[2]; ++k) {

        unsigned long idx = linearIndexFromCoordinate(i,j,k,size[0],size[1]);
        outdata[idx] = f[idx];
      }

  ofstream outFile; 
  outFile.open(filename, ios::out|ios::binary|ios::ate);
  outFile.write((char*)outdata, outsize*sizeof(I));
  outFile.close();

  delete [] outdata;
  printf("data written\n");
}


// input is a set of 3d points
// output only valid if there are more than 3, for genericity
// -- center is the center of a box
// -- orientation is 3 normal vectors, with longes, 2nd longest, and 3rd longest axes
// -- spread is the lengths in each direction 
void ApproximateAspectRatio(vector<Vec3d> points, Vec3d& center, Vec3d* orientation, Vec3d& spread) {
	int num_points = points.size();
	if (num_points <= 3) return; // do nothing... degenerate condition - could be planar

	// Get the 2 farthest points
	int a, b;
	double longest_axis_distance = -1;
	for (int i = 0; i < num_points; i++) {
		for (int j = i + 1; j < num_points; j++) {
			Vec3d& i_p = points[i];
			Vec3d& j_p = points[j];
			double local_dist = (i_p - j_p).MagSq();

			if (local_dist > longest_axis_distance) {
				a = i; b = j; longest_axis_distance = local_dist;
			}
		}
	}
	
	// project all points onto plane orthogonal to a-b, midpoint (a+b)/2
	Vec3d& a_p = points[a];
	Vec3d& b_p = points[b];
	Vec3d midpoint = (a_p + b_p) * 0.5;
	Vec3d normal = (a_p - b_p); 
	normal.Normalize();
	orientation[0] = normal;

	vector<Vec3d> proj_points; 
	proj_points.reserve(num_points);
	for (int i = 0; i < num_points; i++) {
		Vec3d& p = points[i];
		Vec3d orig_to_point = p - midpoint;
		Vec3d proj_p = p - (normal * normal.Dot(orig_to_point));
		proj_points.push_back(proj_p);
	}

	// now find our longest axis again
	int a1, b1;
	double longest_axis_distance1 = -1;
	for (int i = 0; i < num_points; i++) {
		for (int j = i + 1; j < num_points; j++) {
			Vec3d& i_p = proj_points[i];
			Vec3d& j_p = proj_points[j];
			double local_dist = (i_p - j_p).MagSq();

			if (local_dist > longest_axis_distance1) {
				a1 = i; b1 = j; longest_axis_distance1 = local_dist;
			}
		}
	}

	// now project all proj_points to the 2ndary axis
	// project all points onto plane orthogonal to a-b, midpoint (a+b)/2
	Vec3d& a_p1 = proj_points[a1];
	Vec3d& b_p1 = proj_points[b1];
	Vec3d midpoint1 = (a_p1 + b_p1) * 0.5;
	Vec3d normal1 = (a_p1 - b_p1); normal1.Normalize();
	orientation[1] = normal1;

	vector<Vec3d> proj_points_1;
	proj_points_1.reserve(num_points);
	for (int i = 0; i < num_points; i++) {
		Vec3d& p = proj_points[i];
		Vec3d orig_to_point = p - midpoint1;
		Vec3d proj_p = p - (normal1 * normal1.Dot(orig_to_point));
		proj_points_1.push_back(proj_p);
	}

	// now find our last longest axis 
	int a2, b2;
	double longest_axis_distance2 = -1;
	for (int i = 0; i < num_points; i++) {
		for (int j = i + 1; j < num_points; j++) {
			Vec3d& i_p = proj_points_1[i];
			Vec3d& j_p = proj_points_1[j];
			double local_dist = (i_p - j_p).MagSq();

			if (local_dist > longest_axis_distance2) {
				a2 = i; b2 = j; longest_axis_distance2 = local_dist;
			}
		}
	}

	// now we know it all!
	Vec3d& a_p2 = proj_points_1[a2];
	Vec3d& b_p2 = proj_points_1[b2];
	Vec3d midpoint2 = (a_p2 + b_p2) * 0.5;
	Vec3d normal2 = (a_p2 - b_p2); normal2.Normalize();
	
	orientation[2] = normal2;
	center = midpoint2;
	spread = Vec3d(sqrt(longest_axis_distance), sqrt(longest_axis_distance1), sqrt(longest_axis_distance2));


}




void MakeVoids(MscType* msc, MeshType* mesh, std::string filename) {
        
         ofstream outfile;
	 outfile.open(filename+".voids.csv");
	 
	 vtkSmartPointer<vtkPoints> points = 
	   vtkSmartPointer<vtkPoints>::New();
  
	outfile << "NID,coord-X,coord-Y,coord-Z,center-X,center-Y,center-Z,spread-0,spread-1,spread-2,tori_0-X,tori_0-Y,tori_0-Z,tori_1-X,tori_1-Y,tori_1-Z,tori_2-X,tori_2-Y,tori_2-Z";
	
        std::vector<Vec3d> spreads;
        int counter=0;

	MscType::LivingNodesIterator nodes(msc);
	for (nodes.begin(); nodes.valid(); nodes.advance()) {
		auto nid = nodes.value();
		auto n = msc->getNode(nid);
		if (n.dim != 0) continue;			// skip all but minima
		if (n.boundary != 0) continue;		// don't want the boundary minima

		// NOW FIND ALL 1-saddles attached to it
		set<INDEX_TYPE> saddles1;
		MscType::SurroundingLivingArcsIterator arcsit(msc);
		for (arcsit.begin(nid); arcsit.valid(); arcsit.advance()) {
			auto aid = arcsit.value();
			auto a = msc->getArc(aid);
			saddles1.insert(a.upper);		// add the 1 saddle to the set
		}

		// NOW FIND ALL 2-SADDLES ATTACHED
		set<INDEX_TYPE> saddles2;
		for (auto sid : saddles1) {
			for (arcsit.begin(sid); arcsit.valid(); arcsit.advance()) {
				auto aid = arcsit.value();
				auto a = msc->getArc(aid);
				saddles2.insert(a.upper);		// add the 2 saddle to the set
			}
		}

		// NOW FIND ALL MAXIMA ATTACHED
		set<INDEX_TYPE> maxima;
		for (auto sid : saddles2) {
			for (arcsit.begin(sid); arcsit.valid(); arcsit.advance()) {
				auto aid = arcsit.value();
				auto a = msc->getArc(aid);
				maxima.insert(a.upper);		// add the 1 saddle to the set
			}
		}

		// NOW MAXIMA HAS THE MAX IDS AROUND OUR MINIMUM

		vector<Vec3d> max_coords;
		for (auto mid : maxima) {
			auto m = msc->getNode(mid);
			Vec3d max_c;
			mesh->centroid(m.cellindex, max_c);
			max_coords.push_back(max_c);
		}

		// skip small things
		if (max_coords.size() <= 3) continue;
		// now max_coords has coordinates of all the maxima around the minimum
		// do our aspect ratio...
		
		Vec3d center, spread;
		Vec3d orientation[3];

		ApproximateAspectRatio(max_coords, center, orientation, spread);

		/// now output them?
  
		//printf("node %d:\n", nid);
		Vec3d nc; mesh->centroid(n.cellindex, nc);
		//printf("\tcoords: ");	nc.PrintFloat();
		//printf("\tcenter: ");	center.PrintFloat();
		//printf("\tspread: ");	spread.PrintFloat();
		//printf("\tori_0: ");	orientation[0].PrintFloat();
		//printf("\tori_1: ");	orientation[1].PrintFloat();
		//printf("\tori_2: ");	orientation[2].PrintFloat();
		//printf("\t -->Sanitycheck: o1 . o2 = %f, o1 . o3 = %f, o2 . o3 = %f\n", orientation[0].Dot(orientation[1]), orientation[1].Dot(orientation[2]), orientation[1].Dot(orientation[2]));
		
		points->InsertNextPoint ( nc[0]*0.5,nc[1]*0.5,nc[2]*0.5 );
                
                spreads.push_back(spread);  
                counter++;
                outfile<<"\n"<<nid<<","<<nc[0]<<","<<nc[1]<<","<<nc[2]<<","<<center[0]<<","<<center[1]<<","<<center[2]<<","<<spread[0]<<","<<spread[1]<<","<<spread[2]<<","<<orientation[0][0]<<","<<orientation[0][1]<<","<<orientation[0][2]<<","<<orientation[1][0]<<","<<orientation[1][1]<<","<<orientation[1][2]<<","<<orientation[2][0]<<","<<orientation[2][1]<<","<<orientation[2][2];
	}
	outfile.close();

	 // Create a polydata object and add the points to it.
	 vtkSmartPointer<vtkPolyData> polydata = 
	   vtkSmartPointer<vtkPolyData>::New();
	 polydata->SetPoints(points);

	 vtkSmartPointer<vtkDoubleArray> pointNormalsArray =
	   vtkSmartPointer<vtkDoubleArray>::New();
	 pointNormalsArray->SetNumberOfComponents(3); //3d normals (ie x,y,z)
	 pointNormalsArray->SetNumberOfTuples(polydata->GetNumberOfPoints());

	 for(int i=0; i<polydata->GetNumberOfPoints();i++){
           double pN[3]={spreads[i][0],spreads[i][1],spreads[i][2]};
	   pointNormalsArray->SetTuple(i, pN) ;
	 }

	 // Add the normals to the points in the polydata
	 polydata->GetPointData()->SetNormals(pointNormalsArray);

	 // Write the file
	 vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
	   vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	 writer->SetFileName((filename+".centers.vtp").c_str());
#if VTK_MAJOR_VERSION <= 5
	 writer->SetInput(polydata);
#else
	 writer->SetInputData(polydata);
#endif
	 writer->Write();

}

#define USE_ABSOLUTE_PERSISTENCE


int main(int argc, char** argv) {

	// READ IN THE COMMAND LINE ARGUMENTS
	int X, Y, Z;
	std::string filename;
	if (argc < 6) { printf("Usage: X Y Z filename persistence [conected_comp_file]\n"); return 0; }
	sscanf(argv[1], "%d", &X);
	sscanf(argv[2], "%d", &Y);
	sscanf(argv[3], "%d", &Z);
	filename = std::string(argv[4]);
#ifndef USE_ABSOLUTE_PERSISTENCE
	int persistence;
	sscanf(argv[5], "%d", &persistence);
#else
	float persistence;
	sscanf(argv[5], "%f", &persistence);
#endif

	std::string mask_filename;

	int data_dims[3] = { X,Y,Z };

	int* mask_field = NULL;

	if (argc > 6) {
		mask_field = new int[X*Y*Z];
		FILE* fin = fopen(argv[6], "rb");
		fread(mask_field, sizeof(int), X*Y*Z, fin);
		fclose(fin);
	}
	else {
		char maskname[1024];
		string new_name(filename);

		// Use raw blurred connected components file instead of dist field
#if USE_BLUR_FIELD_AND_DIST
		new_name.replace(new_name.end() - 9, new_name.end(), ".raw");
#endif
		sprintf(maskname, "%s.pancakes.raw.components.raw", new_name.c_str());

		if (access(maskname, F_OK) != -1) {
			mask_field = new int[X*Y*Z];
			printf("found field mask here: %s\n", maskname);
			mask_field = new int[X*Y*Z];
			FILE* fin = fopen(maskname, "rb");
			fread(mask_field, sizeof(int), X*Y*Z, fin);
			fclose(fin);
		}
		else {
			new_name = filename;
			std::string toremove = ".distance.raw";
			size_t found = new_name.find(toremove);

			//if(found != std::string::npos)
			//  new_name.replace(found,toremove.length(), "");

			sprintf(maskname, "%s.pancakes.raw.components.raw", new_name.c_str());

			if (access(maskname, F_OK) != -1) {
				mask_field = new int[X*Y*Z];
				printf("found field mask here: %s\n", maskname);
				mask_field = new int[X*Y*Z];
				FILE* fin = fopen(maskname, "rb");
				fread(mask_field, sizeof(int), X*Y*Z, fin);
				fclose(fin);
			}
			else printf("mask file %s not found\n", maskname);
		}


	}

	GridType* underlying_grid;
	GridFuncType* grid_function;
	MeshType *topological_grid;
	TopoFuncType* topological_grid_function;
	GradType *discrete_gradient;

	// set up structures to navigate grid, and load the 3d image
	underlying_grid = new GridType(GInt::Vec3i(X, Y, Z), GInt::Vec3b(0, 0, 0));
	grid_function = new GridFuncType(underlying_grid);
	grid_function->LoadImageFromFile(filename.c_str());

	printf("loaded cont function\n");

	// now set up an indexing scheme to use in a topological interpretation of the 
	// regular grid
	topological_grid = new MeshType(underlying_grid);

	// we will use a lazy-evaluation max vertex mesh function - i.e. the value of a 
	// cell in the topological grid is the maximum value of its vertices in the input
	// image
	MaxVLType* maxv = new MaxVLType(topological_grid, grid_function);
	maxv->ComputeOutput();
	topological_grid_function = new TopoFuncType();
	topological_grid_function->setMeshAndFuncAndMaxVLabeling(topological_grid, grid_function, maxv);

	// read the discrete gradient from disk
	discrete_gradient = new GradType(topological_grid);
	printf("created dgrad struct\n");
	char gradname[2048];
	sprintf(gradname, "%s.grad", argv[4]);
	discrete_gradient->load_from_file(gradname);

	// now compute the Morse-Smale complex
	MscType* msc = new MscType(discrete_gradient, topological_grid, topological_grid_function);
	msc->SetBuildArcGeometry(Vec3b(false, false, true)); // we only need geometric realizations of 2saddle-max arcs
	msc->ComputeFromGrad();

	printf("going to simplify the persistance, %f\n", persistence);
	// simplify to persistence and dump record to file
	char persname[2048];
	sprintf(persname, "%s.2.pers", argv[4]);
	msc->set_output_cancellation_records(persname);

#ifndef USE_ABSOLUTE_PERSISTENCE
	// get persistence to get to the number of maxima requested
	msc->SetPersistanceByNumOfMaxima(persistence);
	printf("set persistence with %d n of maxima\n", persistence);
	float max_persistance = (grid_function->GetMaxValue() - grid_function->GetMinValue()) / 6.0f;
	printf("using max persistance %f\n", max_persistance);
	msc->ComputeHierarchy(max_persistance);
	printf("computed hierarchy\n");
#else
	// use absolute persistence value
	printf("using max persistance %f\n", persistence);
	msc->ComputeHierarchy(persistence);
	printf("computed hierarchy\n");
	printf("set persistence to %f\n", persistence);
	msc->SetSelectPersAbs(persistence);
#endif // !USE_ABSOLUTE_PERSISTENCE

	printf("got here\n");


	// set the persistence threshold
	//msc->SetSelectPersAbs(persistence);

	vector<ligament> ligaments;

	// now begin iteration - use the iterator for the msc class
	MscType::SurroundingLivingArcsIterator arcsit(msc);
	MscType::LivingNodesIterator nodes(msc);

	//sanity check arcs
	unordered_set<INT_TYPE> living_arcs;
	MscType::LivingArcsIterator liv_arcs_it(msc);
	for (liv_arcs_it.begin(); liv_arcs_it.valid(); liv_arcs_it.advance()) {
		living_arcs.insert(liv_arcs_it.value());
	}

	char pairsname[1024];
	sprintf(pairsname, "%s.pairs.txt", argv[4]);
	FILE* f_pairs = fopen(pairsname, "w");
	for (auto arcid : living_arcs) {
		auto& a = msc->getArc(arcid);
		if (a.dim != 2) continue;
		float low_val = msc->getNode(a.lower).value;
		float high_val = msc->getNode(a.upper).value;
		fprintf(f_pairs, "%f %f\n", low_val, high_val);
	}
	fclose(f_pairs);

	int count_maxima = 0;

	for (nodes.begin(); nodes.valid(); nodes.advance()) {
		INT_TYPE node_id = nodes.value();		//  get the id of a living node
		node<float>& n = msc->getNode(node_id); //  get the actual node
		// only look at 2-saddles
		if (n.dim == 3) count_maxima++;
		if (n.dim != 2) continue;
		if (n.value < 0.0) continue; // only want internal ones

		// look at surrounding arcs, only keep 2-saddles that have arcs going to 2 
		// distinct maxima
		int count_neighbors = 0;
		INT_TYPE neighbor_max_ids[2] = { -1, -1 };
		INT_TYPE neighbor_arc_ids[2] = { -1, -1 };
		//printf("tryign stuff arcs...\n");
		for (arcsit.begin(node_id); arcsit.valid(); arcsit.advance()) {
			INT_TYPE arc_id = arcsit.value();
			
			if (living_arcs.count(arc_id) == 0) {
				printf("WHOA insanity detected!\n");
			}

			//printf("nid = %d, aid = %d\n", node_id, arc_id);
			arc<float>& a = msc->getArc(arc_id);
			
			// skip any arcs where the lower node is not our 2-saddle, as
			// it is a 1-2 arc
			if (a.lower == node_id) {
				neighbor_arc_ids[count_neighbors] = arc_id;
				neighbor_max_ids[count_neighbors++] = a.upper;
			}
		}
		//printf("count = %d\t%d\t%d\n", count_neighbors, neighbor_max_ids[0], neighbor_max_ids[1]);
		// now check if this saddle goes between 2 distinct maxs
		if (count_neighbors == 2 && neighbor_max_ids[0] != neighbor_max_ids[1]) {
			// add this saddle to the ligaments list
			ligaments.push_back(ligament{ node_id, neighbor_arc_ids[0], neighbor_arc_ids[1], neighbor_max_ids[0], neighbor_max_ids[1],
      neighbor_max_ids[0], neighbor_max_ids[1]});
		}

	}
	printf("Found %d maxima\nFound %lu unique ligaments\n", count_maxima, ligaments.size());
	//DenseLabeling<char> *maskvol=NULL;
 //       
	//if (mask_field != NULL) {
	//	VolumeConnectedComponents cc(underlying_grid);
	//	maskvol = new DenseLabeling<char>(underlying_grid->NumElements());

	//	for (INDEX_TYPE i = 0; i < underlying_grid->NumElements(); i++)
	//		maskvol->SetLabel(i, mask_field[i] >= 0);
	//}
	//
	//// MEASUREING VOIDS HERE
	MakeVoids(msc, topological_grid, filename);


#ifdef VTK_ENABLED
  printf("writing lines...\n");  
  SkeletonGraph skeleton(msc, topological_grid, grid_function, data_dims, mask_field);
 
#ifndef USE_ABSOLUTE_PERSISTENCE
  skeleton.addLigaments(ligaments, maskvol);
#else
  //skeleton.addUnclippedLigamentsOutsideBlobs(ligaments, maskvol);
  skeleton.SetUseDistRadius(false);
  skeleton.addNonOverlappingLigaments(ligaments, topological_grid_function);
#endif // !USE_ABSOLUTE_PERSISTENCE
  
  skeleton.save(filename+".lig"+argv[5]+".vtp");
#endif

#if COMPUTE_BLOBS
        
  printf("Computing blobs...\n");

	// now gather all the maxima and give them new indices
	unordered_map<INT_TYPE, int> maximum_nodes;
	int nodeid = 0;
	for (auto l : ligaments) {
		arc<float>& a1 = msc->getArc(l.arc1_id);
		arc<float>& a2 = msc->getArc(l.arc2_id);
		if (maximum_nodes.count(a1.upper) == 0) {
			maximum_nodes[a1.upper] = nodeid++;
		}
		if (maximum_nodes.count(a2.upper) == 0) {
			maximum_nodes[a2.upper] = nodeid++;
		}
	}

	int* global_dm_ids = new int[underlying_grid->NumElements()];
	vector<pair<INT_TYPE, int> > vmaxnodes;
	for (auto p : maximum_nodes) vmaxnodes.push_back(p);
	int cmaxnodes = vmaxnodes.size();
  
  #pragma omp parallel for
	for (int i = 0; i < cmaxnodes; i++) {
		auto p = vmaxnodes[i];
		INT_TYPE maxid = p.first;

		set<INDEX_TYPE> geom;
		msc->fillGeometry(maxid, geom, false);

		for (auto id : geom) {
			INDEX_TYPE gridid = topological_grid->VertexNumberFromCellID(id);
			if (grid_function->SampleImage(gridid) > 0.0)
				global_dm_ids[gridid] = maxid;
			else
				global_dm_ids[gridid] = -1;
		}
	}

  int dims[3]={X,Y,Z};
	dumpArray(global_dm_ids, dims, filename+".blobs"+argv[5]+".raw");
#endif
  
#if 0 // the following is an example code to iterate through the ligaments
  
	// iterate through ligaments and build geometry, count incidence, 
	// get function values, etc.
	unordered_map<INT_TYPE, int> node_incidence;

	for (auto l : ligaments) {

		// HERE WE FILL THE GEOMETRIC COORDINATES OF AN ARC
		vector<INDEX_TYPE> linegeomids1;
		vector<INDEX_TYPE> linegeomids2;
		msc->fillArcGeometry(l.arc1_id, linegeomids1);
		msc->fillArcGeometry(l.arc2_id, linegeomids2);

		vector<Vec3d> linegeom1;
		vector<Vec3d> linegeom2;
		for (auto id : linegeomids1) {
			Vec3l vl;
			topological_grid->cellid2Coords(id, vl);
			Vec3d vd = vl * 0.5; // topological index space goes from 0 : 2X-1 , not 0 : X
			linegeom1.push_back(vd);
		}
		for (auto id : linegeomids2) {
			Vec3l vl;
			topological_grid->cellid2Coords(id, vl);
			Vec3d vd = vl * 0.5;
			linegeom2.push_back(vd);
		}
		// now linegeom1 and linegeom2 have the 3d coordinates of the stair-stepping line

		// HERE WE GET SOME NODE ATTRIBUTES

	}
  
#endif

	return 0;

}


