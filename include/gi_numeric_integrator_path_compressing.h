#ifndef NUMERIC_INTEGRATOR_PATH_COMPRESSING_H
#define NUMERIC_INTEGRATOR_PATH_COMPRESSING_H

#include <set>
#include <queue>
#include <stack>
#include <omp.h>

#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_labeling.h"
#include "gi_regular_grid.h"
#include "gi_regular_grid_trilinear_function.h"
#include "gi_adaptive_euler_advector.h"
#include "gi_timing.h"
//#include "gi_union_find_labeling.h"
//#include "gi_numeric_integrator_expanding_region_stop.h"
#include "gi_topological_regular_grid.h"
#include "gi_array_index_partition.h"



//#define OUTPUTINTERMEDIATE
namespace GInt {


    template< class Advector, class GridFuncType >
    class NumericIntegratorPathCompressingToTerminal {

    protected:

		
        DenseLabeling<DestType>* m_desttype;
        DenseLabeling<int>* m_dest_label;
        RegularGrid3D* m_grid;
		GridFuncType* m_func;
        int m_num_iterations_left;

        Vec3i m_xyz;
        Vec3b m_periodic;
        FLOATTYPE m_error_threshold;
        FLOATTYPE m_gradient_threshold;
		FLOATTYPE m_filterValue;

    public:

 
		NumericIntegratorPathCompressingToTerminal(GridFuncType* func, RegularGrid3D* grid, FLOATTYPE error_threshold, FLOATTYPE gradient_threshold, int interation_limit) :
            m_num_iterations_left(interation_limit), m_xyz(func->GetGrid()->XYZ()), m_periodic(func->GetGrid()->Periodic()), m_func(func), m_grid(grid),
            m_gradient_threshold(gradient_threshold), m_error_threshold(error_threshold) {          
			m_filterValue = std::numeric_limits<FLOATTYPE>::min(); 

        }

		~NumericIntegratorPathCompressingToTerminal() {
            delete m_desttype;
            delete m_dest_label;
        }

        DenseLabeling<int>* GetOutputLabels() { return m_dest_label; }
        RegularGrid3D* GetGrid() { return m_grid; }
		GridFuncType* GetFunction() { return m_func; }
 
        template <typename T>
        static void accumulate_and_clear_sets(std::vector<std::set<T>> & sets, std::set<T> &accSet) {

            accSet.clear();
            for(unsigned int i = 0; i < sets.size(); i++) {
                accSet.insert( sets[i].begin(), sets[i].end() );
                sets[i].clear();
            }
        }

#if 1
        void BeginIntegration(const std::unordered_map<INT_TYPE, std::vector<INDEX_TYPE> >& remap, DenseLabeling<int>* teminals, bool verbose = false) {

            ThreadedTimer gtimer(1);
            gtimer.StartGlobal();

            if(verbose)
                printf(" -- Performing numeric integration for volume assignment (%f)...\n", m_filterValue);

			m_dest_label = teminals;
            m_desttype = new DenseLabeling<DestType>(m_grid->NumElements());
            const INDEX_TYPE t_num_vertices = m_grid->NumElements();

            AdvectionChecker* inside_voxel_critical_advection_checker = new TerminateNearPathCompressedRegion(m_desttype, m_grid);

            ThreadedTimer ltimer0(1);
            ltimer0.StartGlobal();

            // ---------------------------------------------------------------------
            if (verbose){
                printf("   -- finding extrema to terminate integral lines...\n");
                fflush(stdout);
            }

#pragma omp parallel for
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {

                // Harsh added a new label for filtering
                //if (fabs(m_func->SampleImage(i)) <= m_filterValue) {
                //    m_desttype->SetLabel(i, DestType::BACKGROUND);
                //    m_dest_label->SetLabel(i, -2);
                //    continue;
                //}
				switch (m_dest_label->GetLabel(i)) {
				case -1:
					m_desttype->SetLabel(i, DestType::UNASSIGNED);
					break;
				default:
					m_desttype->SetLabel(i, DestType::CERTAIN_TERMINAL);
				}
            }

            ltimer0.EndGlobal();

			for (auto lp : remap) {
			
				for (auto id : lp.second) {
					m_func->SetGradExplicit(id, Vec3d(0, 0, 0));
				}
			}


#ifdef OUTPUTINTERMEDIATE
            m_dest_label->OutputToFile("certains.raw");
            m_desttype->OutputToIntFile("certains_type.raw");

            {
                m_dest_label->OutputToFile("certain_expansion.raw");
                m_dest_label->OutputToIntFile("dest_original.raw");
            }

            TopologicalRegularGridRestricted* ttgrid = new TopologicalRegularGridRestricted(m_grid);
            VertexLabelingToBoundaryLabeling<int>* tedge = new VertexLabelingToBoundaryLabeling<int>(m_dest_label, ttgrid);
            tedge->ComputeBoundary();
            tedge->OutputEdgesToFile("certain_edges.txt");
            FILE* fout = fopen("Linesout.txt", "w");
#endif

            if (verbose){
                //printf(" expansion done!\n", mExtrema.size());
                printf("   -- doing numerical integration first pass with path compression...");
                fflush(stdout);
            }

            // ---------------------------------------------------------------------

            ThreadedTimer ltimer1(1);
            ltimer1.StartGlobal();

            int t1, t2; t1 = t2 = 0;
#pragma omp parallel
            {
                Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_critical_advection_checker);
                std::vector<INDEX_TYPE> t_path;
                t_path.reserve(100);

                int num_threads = omp_get_num_threads();
                int thread_num = omp_get_thread_num();

                std::vector<INDEX_TYPE> partition;
                GInt::ArrayIndexPartitioner::EvenChunkSplit(t_num_vertices, num_threads, partition);

                INDEX_TYPE num_to_do = (partition[thread_num + 1] - partition[thread_num]);
                for (INDEX_TYPE kk = 1; kk <= num_to_do * 4; kk*=2) {

                    INDEX_TYPE startsize = num_to_do / kk;
                    INDEX_TYPE stepsize = (num_to_do*2) / kk;

                    if (stepsize == 0) continue;

                    for (INDEX_TYPE i = partition[thread_num] + startsize; i < partition[thread_num + 1]; i += stepsize) {

                        // early skip if this is already a maximum
                        if (m_desttype->GetLabel(i) != DestType::UNASSIGNED) {
                            continue;
                        }

                        t_path.clear();
                        t_path.push_back(i);

                        Vec3l t_coords = m_grid->Coords(i); // get the coordinates of the point
                        Vec3d t_current_point = t_coords;
                        int t_num_iterations_left = m_num_iterations_left;

                        bool t_continue = true;

#ifdef OUTPUTINTERMEDIATE
                        std::vector<Vec3d> line_soup;
                        line_soup.push_back(t_current_point);
#endif
                        while (t_continue) {

                            ADVECTION_EVENT t_return_code;
                            if (m_grid->DistToBoundary(t_coords) <= 1) {
                                t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
                                t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
                                t1++;
                            }
                            else {
                                t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
                                t_coords = (t_current_point + 0.5);
                                t2++;
                            }

                            INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
                            t_path.push_back(t_next_id);

                            if (t_return_code == ADVECTION_EVENT::OUT_OF_VOXEL) continue;
#ifdef OUTPUTINTERMEDIATE
                            line_soup.push_back(t_current_point);
#endif
                            // if we terminated or hit a critical point, then we are done
                            if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
                                t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
                                t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
                                t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {

                                int t_dest_label = m_dest_label->GetLabel(t_next_id);

                                //#pragma omp critical
                                {
                                for (int j = 0; j < t_path.size(); j++) {
                                    INDEX_TYPE jj = t_path[j];

                                    if (m_desttype->GetLabel(jj) == DestType::UNASSIGNED) {
                                        m_dest_label->SetLabel(jj, t_dest_label);
                                        m_desttype->SetLabel(jj, DestType::ASSIGNED);
                                    }
                                }
                                }
#ifdef OUTPUTINTERMEDIATE
#pragma omp critical
                            {
                                int tn = omp_get_thread_num();
                                fprintf(fout, "%d %d %d %d\n", i, tn, line_soup.size(), t_dest_label);
                                for (int j = 0; j < line_soup.size(); j++) {
                                    fprintf(fout, "%f %f %f\n", line_soup[j][0], line_soup[j][1], line_soup[j][2]);
                                }
                            }
#endif
                                t_continue = false;
                            }
                        }
                    }
                }
            }
#ifdef OUTPUTINTERMEDIATE
            fclose(fout);
            m_dest_label->OutputToIntFile("first_integration.raw");
            m_dest_label->OutputToIntFile("dests_after_first_integration.raw");

            m_dest_label->OutputToIntFile("first_integration.raw");
            m_desttype->OutputToIntFile("first_integration_type.raw");
#endif

            ltimer1.EndGlobal();
            // ---------------------------------------------------------------------

            if (verbose){
                printf(" done!");             ltimer1.PrintAll();
                printf("   -- checking unambiguous voxels...");
                fflush(stdout);
            }

            ThreadedTimer ltimer2(1);
            ltimer2.StartGlobal();

            // we will process vertices in iterations
            // this set will contain verts to be processed in the next iteration
            std::set<INDEX_TYPE> verts_2b_processed_set;

            // to avoid locking by threads, we will use local copies where threads will add verts
            // later, we will merge these into the main set defined above
            std::vector< std::set<INDEX_TYPE> > verts_thrds ( omp_get_max_threads() );

#pragma omp parallel for
            for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {

                if (m_desttype->GetLabel(i) == DestType::BACKGROUND)
                    continue;
				//if (m_desttype->GetLabel(i) == DestType::ASSIGNED) {
				//	verts_thrds[omp_get_thread_num()].insert(i);
				//}
                Vec3l t_coords = m_grid->Coords(i); // get the coordinates of the poitn
                Vec3l negs[6];

                int nn = m_grid->GatherExistingNeighborsAll6(t_coords, negs);

                for (int j = 0; j < nn; j++) {

                    INDEX_TYPE v2 = m_grid->Index3d(negs[j]);
                    if (v2 > i) continue;

                    if (m_dest_label->GetLabel(i) != m_dest_label->GetLabel(v2)) {

                        if (m_desttype->GetLabel(i) == DestType::ASSIGNED) {
                            verts_thrds[ omp_get_thread_num() ].insert(i);
                        }
                        if (m_desttype->GetLabel(v2) == DestType::ASSIGNED) {
                            verts_thrds[ omp_get_thread_num() ].insert(v2);
                        }
                    }
                }
            }

            // collect the results of all threads
            accumulate_and_clear_sets( verts_thrds, verts_2b_processed_set );
            ltimer2.EndGlobal();
            // ---------------------------------------------------------------------

            if (verbose){
                printf(" done!");   ltimer2.PrintAll();
                fflush(stdout);
                //printf("   -- found %d points needed correction...", verts_2b_processed_set.size());
                //fflush(stdout);
            }


            ThreadedTimer ltimer3(1);
            ltimer3.StartGlobal();

            AdvectionChecker* inside_voxel_nostop_advection_checker = new TerminateNearOriginalCertain(m_desttype, m_grid);

            size_t totalfixed = 0;

            // this loop will iterate until no more verts need to be processed
            for(unsigned int itern = 0; !verts_2b_processed_set.empty(); itern++) {

                // transfer from set to vector to start processing
                // use the set to store the verts needed in the next iteration
                std::vector<INDEX_TYPE> verts_2b_processed_vector ( verts_2b_processed_set.begin(), verts_2b_processed_set.end() );

                if (verbose){
                    printf("         -- iteration %d will process %d vertices\n", itern, verts_2b_processed_vector.size());
                }

                totalfixed += verts_2b_processed_vector.size();
#pragma omp parallel for
                for(int i = 0; i < verts_2b_processed_vector.size(); i++) {

                    Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_nostop_advection_checker);
                    INDEX_TYPE current_vertex = verts_2b_processed_vector[i];

                    int init_label = m_dest_label->GetLabel(current_vertex);

                    // INTEGRATE
                    // INTEGRATE
                    Vec3l t_coords = m_grid->Coords(current_vertex); // get the coordinates of the poitn
                    //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);

                    Vec3d t_current_point = t_coords;
                    int t_num_iterations_left = m_num_iterations_left;
                    bool t_continue = true;
                    int new_label;

                    while (t_continue) {

                        ADVECTION_EVENT t_return_code;
                        if (m_grid->DistToBoundary(t_coords) <= 1) {
                            t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
                            t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
                            t1++;
                        }
                        else {
                            t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
                            t_coords = (t_current_point + 0.5);
                            t2++;
                        }
                        INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
                        // if we terminated or hit a critical point, then we are done
                        if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
                            t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
                            t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
                            t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
                            new_label = m_dest_label->GetLabel(t_next_id);
                            if (m_desttype->GetLabel(t_next_id) != DestType::CERTAIN_TERMINAL){
                                //printf("whoatherenelly %d %d\n", m_desttype->GetLabel(t_next_id), t_return_code);
                            }
                            t_continue = false;
                        }
                    }
                    // INTEGRATE
                    // INTEGRATE

					//
                    m_dest_label->SetLabel(current_vertex, new_label);
					m_desttype->SetLabel(current_vertex, DestType::CERTAIN_NONTERMINAL);

                    if (new_label != init_label){

                        // ENQUEUE NEIGHBORS
                        Vec3l t_coords = m_grid->Coords(current_vertex); // get the coordinates of the poitn
                        Vec3l negs[6];
                        //INDEX_TYPE negids[6];
						int nn = m_grid->GatherExistingNeighborsAll6(t_coords, negs);

                        // for each neigbhor
                        for (int j = 0; j < nn; j++){

                            INDEX_TYPE negid = m_grid->Index3d(negs[j]);

                            // only if it has not yet been added to our update set
                            if (m_desttype->GetLabel(negid) == DestType::ASSIGNED &&
                                m_dest_label->GetLabel(negid) != new_label) {
                                    verts_thrds[ omp_get_thread_num() ].insert(negid);
                            }
                        }
                    }
                }

                accumulate_and_clear_sets( verts_thrds, verts_2b_processed_set );
            }
            

			ltimer3.EndGlobal();


#ifdef OUTPUTINTERMEDIATE
            m_desttype->OutputToIntFile("classes_type.raw");
#endif

            if (verbose){
                printf("   -- done! fixed a total of %d vertices!", totalfixed);              ltimer3.PrintAll();
            }

            gtimer.EndGlobal();
            if(verbose){
                printf(" -- done numerical integration!");  gtimer.PrintAll();
            }
        }
#else
        void BeginIntegration(bool verbose = false) {

                    verbose = true;

                    ThreadedTimer gtimer(1);
                    gtimer.StartGlobal();

                    if(verbose)
                        printf(" -- Performing numeric integration for volume assignment (%f)...\n", m_filterValue);

                    //m_func->ComputeGradFromImage(m_rkindex);

                    m_dest_label = new DenseLabeling<int>(m_grid->NumElements());
                    m_desttype = new DenseLabeling<DestType>(m_grid->NumElements());
                    const INDEX_TYPE t_num_vertices = m_grid->NumElements();

                    // THIS WILL NEED TO CHANGE
                    AdvectionChecker* inside_voxel_critical_advection_checker = new TerminateNearPathCompressedRegion(m_desttype, m_grid);
                    //AdvectionChecker* no_check = new NoTermination();//AdvectionChecker* inside_voxel_advection_checker = new TerminateNearAssigned(m_destinations, m_grid);


                    ThreadedTimer ltimer0(1);
                    ltimer0.StartGlobal();

                    if (verbose){
                        printf("   -- finding extrema to terminate integral lines...");
                        fflush(stdout);
                    }

                    // set all potential extrema, so we terminate near them
        //#pragma omp parallel for
                    for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {

                        //printf(" thread %d of %d max %d does vertx %d\n",
                        //       omp_get_thread_num(), omp_get_num_threads(), omp_get_max_threads(), i );

                        // Harsh added a new label for filtering
                        if (fabs(m_func->SampleImage(i)) <= m_filterValue) {
                            m_desttype->SetLabel(i, DestType::BACKGROUND);
                            m_dest_label->SetLabel(i, -2);
                            continue;
                        }

                        m_desttype->SetLabel(i, DestType::UNASSIGNED);
                        m_dest_label->SetLabel(i, -1);

                        if (IsExtremeVertexIn6Neighborhood(i)) {
        //#pragma omp critical
                            {
                                mExtrema.push_back(i);
                                //m_extrema.insert(i);
                                //m_dest_label->SetLabel(i, 0);
                            }
                        }

                    }


        #ifdef OUTPUTINTERMEDIATE
                    m_dest_label->OutputToIntFile("crits.raw");
        #endif
                    ltimer0.EndGlobal();
                    ltimer0.PrintAll();

                    if (verbose){
                        printf(" done! found %d extrema!\n", mExtrema.size());
                        printf("   -- expanding extrema certain regions...");
                        fflush(stdout);
                    }

                    int num_extrema = mExtrema.size();
        //#pragma omp parallel shared(mExtrema)
                    {
        //#pragma omp for schedule(dynamic)  nowait
                        for (int m = 0; m < num_extrema; m++) {
                            INDEX_TYPE maximum = mExtrema[m];
                            Expand_Lower_Neighborhood(maximum, m);
                            m_func->SetGradExplicit(maximum, Vec3d(0, 0, 0));
                        }
                    }

        #ifdef OUTPUTINTERMEDIATE
                    m_dest_label->OutputToFile("certains.raw");
                    m_desttype->OutputToIntFile("certains_type.raw");

                    {
                        m_dest_label->OutputToFile("certain_expansion.raw");
                        m_dest_label->OutputToIntFile("dest_original.raw");
                    }

                    TopologicalRegularGridRestricted* ttgrid = new TopologicalRegularGridRestricted(m_grid);
                    VertexLabelingToBoundaryLabeling<int>* tedge = new VertexLabelingToBoundaryLabeling<int>(m_dest_label, ttgrid);
                    tedge->ComputeBoundary();
                    tedge->OutputEdgesToFile("certain_edges.txt");
                    FILE* fout = fopen("Linesout.txt", "w");
        #endif
                    if (verbose){
                        printf(" expansion done!\n", mExtrema.size());
                        printf("   -- doing numerical integration first pass with path compression...");
                        fflush(stdout);
                    }

                    ThreadedTimer ltimer1(1);
                    ltimer1.StartGlobal();

                    int t1, t2; t1 = t2 = 0;
        //#pragma omp parallel
                    {
                        Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_critical_advection_checker);
                        std::vector<INDEX_TYPE> t_path;
                        t_path.reserve(100);

                        int num_threads = omp_get_num_threads();
                        int thread_num = omp_get_thread_num();

                        printf(" thread %d of %d\n", thread_num, num_threads);

                        std::vector<INDEX_TYPE> partition;
                        GInt::ArrayIndexPartitioner::EvenChunkSplit(t_num_vertices, num_threads, partition);
                        INDEX_TYPE num_to_do = (partition[thread_num + 1] - partition[thread_num]);
                        for (INDEX_TYPE kk = 1; kk <= num_to_do * 4; kk*=2) {
                            INDEX_TYPE startsize = num_to_do / kk;
                            INDEX_TYPE stepsize = (num_to_do*2) / kk;
                            if (stepsize == 0) continue;
                            //printf("%d %d %d\n", num_to_do, kk, stepsize);
                            for (INDEX_TYPE i = partition[thread_num] + startsize; i < partition[thread_num + 1]; i += stepsize) {


                                //#pragma omp for schedule(guided)  nowait
                                //				for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {

                                // early skip if this is already a maximum
                                if (m_desttype->GetLabel(i) != DestType::UNASSIGNED) {
                                    continue;
                                }

                                t_path.clear();
                                t_path.push_back(i);
                                Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
                                //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);
                                Vec3d t_current_point = t_coords;
                                int t_num_iterations_left = m_num_iterations_left;

                                bool t_continue = true;
        #ifdef OUTPUTINTERMEDIATE
                                std::vector<Vec3d> line_soup;
                                line_soup.push_back(t_current_point);
        #endif
                                while (t_continue) {
                                    Vec3d t_next_point;
                                    ADVECTION_EVENT t_return_code;
                                    if (m_grid->DistToBoundary(t_coords) <= 1) {
                                        t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
                                        t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
                                        t1++;
                                    }
                                    else {
                                        t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
                                        t_coords = (t_current_point + 0.5);
                                        t2++;
                                    }
                                    INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
                                    t_path.push_back(t_next_id);
                                    if (t_return_code == ADVECTION_EVENT::OUT_OF_VOXEL) continue;
        #ifdef OUTPUTINTERMEDIATE
                                    line_soup.push_back(t_current_point);
        #endif					// if we terminated or hit a critical point, then we are done
                                    if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
                                        t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
                                        t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
                                        t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
                                        int t_dest_label = m_dest_label->GetLabel(t_next_id);
                                        //DestType t_certain_label = m_certains->GetLabel(t_next_id);
                                        //#pragma omp critical
                                        //if (t_dest_label == -1) {
                                        //	printf("who there got here %d %d %d %d %d\n",
                                        //		m_desttype->GetLabel(t_next_id), t_next_id, t_return_code, t_num_iterations_left, t_path.size());
                                        //	m_grid->XYZ3d(i).PrintInt(); printf("->");  t_coords.PrintInt();
                                        //}
                                        //#pragma omp critical
                                {
                                    for (int j = 0; j < t_path.size(); j++) {
                                        INDEX_TYPE jj = t_path[j];
                                        if (m_desttype->GetLabel(jj) == DestType::UNASSIGNED) {
                                            m_dest_label->SetLabel(jj, t_dest_label);
                                            m_desttype->SetLabel(jj, DestType::ASSIGNED);
                                        }
                                        //m_certains->SetLabel(t_path[j], t_certain_label);
                                    }
        #ifdef OUTPUTINTERMEDIATE
        #pragma omp critical
                                    {
                                        int tn = omp_get_thread_num();
                                        fprintf(fout, "%d %d %d %d\n", i, tn, line_soup.size(), t_dest_label);
                                        for (int j = 0; j < line_soup.size(); j++) {
                                            fprintf(fout, "%f %f %f\n", line_soup[j][0], line_soup[j][1], line_soup[j][2]);
                                        }
                                    }
        #endif
                                }
                                        t_continue = false;
                                    }
                                }
                            }
                        }
                    }
        #ifdef OUTPUTINTERMEDIATE
                    fclose(fout);
                    m_dest_label->OutputToIntFile("first_integration.raw");
                    m_dest_label->OutputToIntFile("dests_after_first_integration.raw");

                    m_dest_label->OutputToIntFile("first_integration.raw");
                    m_desttype->OutputToIntFile("first_integration_type.raw");
        #endif

                    ltimer1.EndGlobal();
                    ltimer1.PrintAll();


                    if (verbose){
                        printf(" done!\n");
                        printf("   -- checking unambiguous voxels...");
                        fflush(stdout);
                    }

                    ThreadedTimer ltimer2(1);
                    ltimer2.StartGlobal();

                    std::unordered_set<INDEX_TYPE> added_vertices;
                    std::stack<INDEX_TYPE> vertex_stack;

        //#pragma omp parallel for
                    for (INDEX_TYPE i = 0; i < t_num_vertices; i++) {


                        Vec3l t_coords = m_grid->XYZ3d(i); // get the coordinates of the poitn
                        Vec3l negs[6];
                        //INDEX_TYPE negids[6];
                        int nn = m_grid->GatherExistingNeighborsSameBdry6(t_coords, negs);
                        for (int j = 0; j < nn; j++) {
                            INDEX_TYPE v2 = m_grid->Index3d(negs[j]);
                            if (v2 > i) continue;

                            //if (m_desttype->GetLabel(i) != DestType::ASSIGNED || m_desttype->GetLabel(v2) != DestType::ASSIGNED) continue;
                            if (m_dest_label->GetLabel(i) != m_dest_label->GetLabel(v2)) {

                                if (m_desttype->GetLabel(i) == DestType::ASSIGNED) {
        //#pragma omp critical
                                    {
                                        if (added_vertices.count(i) == 0) {
                                            added_vertices.insert(i);
                                            vertex_stack.push(i);
                                        }
                                    }
                                }
                                if (m_desttype->GetLabel(v2) == DestType::ASSIGNED) {
        //#pragma omp critical
                                    {
                                        if (added_vertices.count(v2) == 0) {
                                            added_vertices.insert(v2);
                                            vertex_stack.push(v2);
                                        }
                                    }
                                }

                            }
                        }

                    }




        //			TopologicalRegularGrid3D* tgrid = new TopologicalRegularGrid3D(m_grid);
        //			int num_cells = tgrid->numCells();
        //
        //
        //			/// in parallel gather the vertices that need to be updated
        //#pragma omp parallel
        //			{
        //
        //				int num_threads = omp_get_num_threads();
        //				int thread_num = omp_get_thread_num();
        //
        //				std::vector<INDEX_TYPE> partition;
        //				ArrayIndexPartitioner::EvenChunkSplit(tgrid->numCells(), num_threads, partition);
        //				TopologicalRegularGrid3D::DCellsIterator edges(tgrid, 1, partition[thread_num], partition[thread_num + 1]);
        //				for (edges.begin(); edges.valid(); edges.advance()) {
        //					TopologicalRegularGrid3D::FacetsIterator fit(tgrid);
        //					fit.begin(edges.value());
        //					INDEX_TYPE tv1 = fit.value();
        //					fit.advance();
        //					INDEX_TYPE tv2 = fit.value();
        //
        //					INDEX_TYPE v1 = tgrid->VertexNumberFromCellID(tv1);
        //					INDEX_TYPE v2 = tgrid->VertexNumberFromCellID(tv2);
        //
        //					if (m_desttype->GetLabel(v1) != DestType::ASSIGNED || m_desttype->GetLabel(v2) != DestType::ASSIGNED) continue;
        //					if (m_dest_label->GetLabel(v1) != m_dest_label->GetLabel(v2)) {
        //
        //						if (m_desttype->GetLabel(v1) == DestType::ASSIGNED) {
        //#pragma omp critical
        //							{
        //								if (added_vertices.count(v1) == 0) {
        //									added_vertices.insert(v1);
        //									vertex_stack.push(v1);
        //								}
        //							}
        //						}
        //						if (m_desttype->GetLabel(v2) == DestType::ASSIGNED) {
        //#pragma omp critical
        //							{
        //								if (added_vertices.count(v2) == 0) {
        //									added_vertices.insert(v2);
        //									vertex_stack.push(v2);
        //								}
        //							}
        //						}
        //					}
        //				}
        //			} // END PARALLEL SECTION
        //

                    //for (auto id : added_vertices) {
                    //	m_desttype->SetLabel(id, DestType::BACKGROUND);
                    //}
                    //m_desttype->OutputToIntFile("reintegrate.raw");

                    if (verbose){
                        printf(" done!");
                        ltimer2.EndGlobal();
                        ltimer2.PrintAll();

                        printf("   -- found %d points needed correction %d...", added_vertices.size(), vertex_stack.size());
                        fflush(stdout);
                    }

                    ThreadedTimer ltimer3(1);
                    ltimer3.StartGlobal();

                    // NOW FIX LABELS
                    AdvectionChecker* inside_voxel_nostop_advection_checker = new TerminateNearOriginalCertain(m_desttype, m_grid);

        //#pragma omp parallel
                    {
                        int cnt = 0;
                        Advector t_advector(m_grid, m_func, m_gradient_threshold, m_error_threshold, inside_voxel_nostop_advection_checker);
                        bool keep_going = true;
                        while (keep_going) {
                            INDEX_TYPE current_vertex;
        //#pragma omp critical
                            {
                                if (vertex_stack.size() > 0)
                                {
                                    current_vertex = vertex_stack.top();
                                    vertex_stack.pop();
                                }
                                else {
                                    keep_going = false;
                                }
                            }

                            if (cnt ++ < 10 ){
                                printf(" [%d] = %d\n", cnt, current_vertex);
                            }
                            if (keep_going) {
                                // INITIAL VALUE
                                int init_label = m_dest_label->GetLabel(current_vertex);

                                // INTEGRATE
                                // INTEGRATE
                                Vec3l t_coords = m_grid->XYZ3d(current_vertex); // get the coordinates of the poitn
                                //if (t_coords[0] == 0 && t_coords[1] == 0) printf("doing %d\n", t_coords[2]);

                                Vec3d t_current_point = t_coords;
                                int t_num_iterations_left = m_num_iterations_left;
                                bool t_continue = true;
                                int new_label;

                                while (t_continue) {
                                    Vec3d t_next_point;
                                    ADVECTION_EVENT t_return_code;
                                    if (m_grid->DistToBoundary(t_coords) <= 1) {
                                        t_return_code = t_advector.AdvectThroughVoxelNearBoundary(t_current_point, t_num_iterations_left);
                                        t_coords = m_grid->Inbounds(t_current_point + 0.5); // get nearest integer voxel
                                        t1++;
                                    }
                                    else {
                                        t_return_code = t_advector.AdvectThroughVoxelNoCheck(t_current_point, t_num_iterations_left);
                                        t_coords = (t_current_point + 0.5);
                                        t2++;
                                    }
                                    INDEX_TYPE t_next_id = m_grid->Index3d(t_coords);
                                    // if we terminated or hit a critical point, then we are done
                                    if (t_return_code == ADVECTION_EVENT::LOW_GRADIENT ||
                                        t_return_code == ADVECTION_EVENT::HIT_EXTREMUM ||
                                        t_return_code == ADVECTION_EVENT::HIT_PREASSIGNED ||
                                        t_return_code == ADVECTION_EVENT::OVER_MAX_ITERATIONS) {
                                        new_label = m_dest_label->GetLabel(t_next_id);
                                        if (m_desttype->GetLabel(t_next_id) != DestType::CERTAIN_TERMINAL) printf("whoatherenelly %d %d\n", m_desttype->GetLabel(t_next_id), t_return_code);
                                        t_continue = false;
                                    }
                                }
                                // INTEGRATE
                                // INTEGRATE

                                m_dest_label->SetLabel(current_vertex, new_label);
                                m_desttype->SetLabel(current_vertex, DestType::CERTAIN_NONTERMINAL);

                                if (new_label != init_label) {
                                    // ENQUEUE NEIGHBORS
                                    Vec3l t_coords = m_grid->XYZ3d(current_vertex); // get the coordinates of the poitn
                                    Vec3l negs[6];
                                    INDEX_TYPE negids[6];
                                    int nn = m_grid->GatherExistingNeighborsSameBdry6(t_coords, negs);
                                    for (int j = 0; j < nn; j++)  negids[j] = m_grid->Index3d(negs[j]);
        //#pragma omp critical
                                    {
                                        // for each neigbhor
                                        for (int j = 0; j < nn; j++) {
                                            INDEX_TYPE negid = negids[j];
                                            // only if it has not yet been added to our update set
                                            if (added_vertices.count(negid) == 0) {
                                                if (m_desttype->GetLabel(negid) == DestType::ASSIGNED &&
                                                    m_dest_label->GetLabel(negid) != new_label) {
                                                    added_vertices.insert(negid);
                                                    vertex_stack.push(negid);
                                                }
                                            }
                                        }
                                    }



                                }

                            }
                        } // END WHILE
                    } // END PARALLEL

                    ltimer3.EndGlobal();
                    ltimer3.PrintAll();


        #ifdef OUTPUTINTERMEDIATE
                    m_desttype->OutputToIntFile("classes_type.raw");
        #endif

                    if (verbose){
                        printf(" done! fixed a total of %d vertices\n", added_vertices.size());
                    }
                    gtimer.EndGlobal();
                    gtimer.PrintAll();

                }
#endif
            };


}
#endif
