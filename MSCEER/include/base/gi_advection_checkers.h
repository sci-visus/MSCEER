#ifndef ADVECTION_CHECKER_H
#define ADVECTION_CHECKER_H

#include <unordered_set>

#include "base/gi_basic_types.h"
#include "base/gi_vectors.h"
#include "base/gi_regular_grid_2d.h"
#include "base/gi_regular_grid_3d.h"
#include "base/gi_labeling.h"
#include "base/gi_advection_events.h"
#include "base/gi_topological_regular_grid_2d.h"
#include "base/gi_topological_regular_grid_3d.h"


namespace GInt {
    class AdvectionChecker {
    public:

        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary){ return NONE; }
    };

    // no early termination during advection interior to a voxel
    class NoTermination : public AdvectionChecker {
    public:
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            return NONE;
        }
    };


    // this variant
    class TerminateNearCrits : public AdvectionChecker {
    protected:
        DenseLabeling<INDEX_TYPE>* m_labeling;
        RegularGrid3D* m_grid;
    public:
        TerminateNearCrits(DenseLabeling<INDEX_TYPE>* labeling, RegularGrid3D*grid) : m_labeling(labeling), m_grid(grid) {}
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) == id) { // THIS IS WRONG!!!
                    return HIT_EXTREMUM; // reached a critical point
                }
            }
            else {
                Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) == id) { // THIS IS WRONG!!!
                    return HIT_EXTREMUM; // reached a critical point
                }

            }
            return NONE;
        }

    };

	// only works for non-periodic meshes right now
	template < class MeshType>
	class TerminateNearAssignedHex : public AdvectionChecker {
	protected:
		DenseLabeling<char>* m_labeling;
		MeshType* m_mesh;
		const int mTargetValue;
	public:
		TerminateNearAssignedHex(DenseLabeling<char>* labeling, MeshType* mesh, int target) : 
			m_labeling(labeling), m_mesh(mesh), mTargetValue(target) {}
		virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
			
			// closest_hex
			if (nearboundary) return HIT_PREASSIGNED;
			Vec3l closest_hex_coords = (new_point.IntFloor() * 2) + Vec3l(1, 1, 1);
			INDEX_TYPE hex_id = m_mesh->coords2Cellid(closest_hex_coords);
			if (m_labeling->GetLabel(hex_id) == mTargetValue) return HIT_PREASSIGNED;
			return NONE;
		}

	};

	// only works for non-periodic meshes right now
	template < class MeshType>
	class TerminateNearAssignedVert : public AdvectionChecker {
	protected:
		DenseLabeling<char>* m_labeling;
		MeshType* m_mesh;
		const int mTargetValue;
	public:
		TerminateNearAssignedVert(DenseLabeling<char>* labeling, MeshType* mesh, int target) :
			m_labeling(labeling), m_mesh(mesh), mTargetValue(target) {}
		virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {

			// closest_hex
			if (nearboundary) return HIT_PREASSIGNED;
			Vec3l closest_vert_coords = ((new_point + Vec3d(0.5, 0.5, 0.5)).IntFloor() * 2);
			INDEX_TYPE vert_id = m_mesh->coords2Cellid(closest_vert_coords);
			if (m_labeling->GetLabel(vert_id) == mTargetValue) return HIT_PREASSIGNED;
			return NONE;
		}

	};

	class AdvectionChecker2D {
	public:

		virtual ADVECTION_EVENT CheckAndDoStuff(const Vec2d& old_point, const Vec2d& new_point, bool nearboundary){ return NONE; }
	};
	// only works for non-periodic meshes right now
	template < class MeshType>
	class TerminateNearAssignedQuad2D : public AdvectionChecker2D {
	protected:
		DenseLabeling<char>* m_labeling;
		MeshType* m_mesh;
		const int mTargetValue;
	public:
		TerminateNearAssignedQuad2D(DenseLabeling<char>* labeling, MeshType* mesh, int target) :
			m_labeling(labeling), m_mesh(mesh), mTargetValue(target) {}
		virtual ADVECTION_EVENT CheckAndDoStuff(const Vec2d& old_point, const Vec2d& new_point, bool nearboundary) {

			// closest_hex
			if (nearboundary) return HIT_PREASSIGNED;
			Vec2l closest_hex_coords = (new_point.IntFloor() * 2) + Vec2l(1, 1);
			INDEX_TYPE hex_id = m_mesh->coords2Cellid(closest_hex_coords);
			if (m_labeling->GetLabel(hex_id) == mTargetValue) return HIT_PREASSIGNED;
			return NONE;
		}

	};

	// only works for non-periodic meshes right now
	template < class MeshType>
	class TerminateNearAssignedVert2D : public AdvectionChecker2D {
	protected:
		DenseLabeling<char>* m_labeling;
		MeshType* m_mesh;
		const int mTargetValue;
	public:
		TerminateNearAssignedVert2D(DenseLabeling<char>* labeling, MeshType* mesh, int target) :
			m_labeling(labeling), m_mesh(mesh), mTargetValue(target) {}
		virtual ADVECTION_EVENT CheckAndDoStuff(const Vec2d& old_point, const Vec2d& new_point, bool nearboundary) {

			// closest_hex
			if (nearboundary) return HIT_PREASSIGNED;
			Vec2l closest_vert_coords = ((new_point + Vec2d(0.5, 0.5)).IntFloor() * 2);
			INDEX_TYPE vert_id = m_mesh->coords2Cellid(closest_vert_coords);
			if (m_labeling->GetLabel(vert_id) == mTargetValue) return HIT_PREASSIGNED;
			return NONE;
		}

	};



    class TerminateNearAssigned : public AdvectionChecker {
    protected:
        DenseLabeling<INDEX_TYPE>* m_labeling;
        RegularGrid3D* m_grid;
    public:
        TerminateNearAssigned(DenseLabeling<INDEX_TYPE>* labeling, RegularGrid3D*grid) : m_labeling(labeling), m_grid(grid) {}
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) != -1) {
                    return HIT_PREASSIGNED; // reached a critical point
                }
            }
            else {
                Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) != -1) {
                    return HIT_PREASSIGNED; // reached a critical point
                }

            }
            return NONE;
        }

    };
    class TerminateNearCertain : public AdvectionChecker {
    protected:
        DenseLabeling<int>* m_labeling;
        RegularGrid3D* m_grid;
    public:
        TerminateNearCertain(DenseLabeling<int>* labeling, RegularGrid3D*grid) : m_labeling(labeling), m_grid(grid) {}
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) != -1) {
                    return HIT_PREASSIGNED; // reached a critical point
                }
            }
            else {
                Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_labeling->GetLabel(id) != -1) {
                    return HIT_PREASSIGNED; // reached a critical point
                }

            }
            return NONE;
        }

    };

	//labels are held on the VERTICES of the grid, and each region extends to 1/2 the
	// width of a voxel in each direction
	class TerminateNearPathCompressedRegion : public AdvectionChecker {
	protected:
		DenseLabeling<DestType>* m_labeling;
		RegularGrid3D* m_grid;
	public:
		TerminateNearPathCompressedRegion(DenseLabeling<DestType>* labeling, RegularGrid3D*grid) : m_labeling(labeling), m_grid(grid) {}
		virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
			if (!nearboundary) {
				Vec3l closest_vertex = new_point + 0.5;
				INDEX_TYPE id = m_grid->Index3d(closest_vertex);
				if (m_labeling->GetLabel(id) == DestType::ASSIGNED || m_labeling->GetLabel(id) == DestType::CERTAIN_TERMINAL) {
					return HIT_PREASSIGNED; // reached a path-compressed region
				}
			}
			else {
				Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
				INDEX_TYPE id = m_grid->Index3d(closest_vertex);
				if (m_labeling->GetLabel(id) == DestType::ASSIGNED || m_labeling->GetLabel(id) == DestType::CERTAIN_TERMINAL) {
					return HIT_PREASSIGNED; // reached a path-compressed region
				}

			}
			return NONE;
		}

	};

	class TerminateNearOriginalCertain : public AdvectionChecker {
	protected:
		DenseLabeling<DestType>* m_labeling;
		RegularGrid3D* m_grid;
	public:
		TerminateNearOriginalCertain(DenseLabeling<DestType>* labeling, RegularGrid3D*grid) : m_labeling(labeling), m_grid(grid) {}
		virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
			if (!nearboundary) {
				Vec3l closest_vertex = new_point + 0.5;
				INDEX_TYPE id = m_grid->Index3d(closest_vertex);
				if (m_labeling->GetLabel(id) == DestType::CERTAIN_TERMINAL) {
					return HIT_PREASSIGNED; 
				}
			}
			else {
				Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
				INDEX_TYPE id = m_grid->Index3d(closest_vertex);
				if (m_labeling->GetLabel(id) == DestType::CERTAIN_TERMINAL) {
					return HIT_PREASSIGNED; 
				}

			}
			return NONE;
		}

	};

    // THIS IS NOT FINISHED!!!
    class TerminateNearAssignedAndUpdate : public AdvectionChecker {
    protected:
        DenseLabeling<bool>* m_labeling;
        RegularGrid3D* m_grid;
    public:
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                if (m_labeling->GetLabel(m_grid->Index3d(closest_vertex))) {
                    return HIT_EXTREMUM; // reached a critical point
                }
            }
            return NONE;
        }

    };






    // no early termination during advection interior to a voxel
    class NoTermination2D : public AdvectionChecker2D {
    public:
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec2d& old_point, const Vec2d& new_point, bool nearboundary) {
            return NONE;
        }
    };



    class TerminateNearExtrema : public AdvectionChecker {
    protected:
        RegularGrid3D* m_grid;
    public:
        std::unordered_set<INDEX_TYPE> m_extrema;
        TerminateNearExtrema(const std::unordered_set<INDEX_TYPE>& extrema, RegularGrid3D* grid) : m_grid(grid), m_extrema(extrema) {}
        TerminateNearExtrema(const std::vector<INDEX_TYPE>& extrema, RegularGrid3D* grid) : m_grid(grid) {
            for(size_t i = 0; i < extrema.size(); i++)
                m_extrema.insert(extrema[i]);
        }
        void Printstuff() {
            printf("Termination minima are: \n");
            for (auto it = m_extrema.begin(); it != m_extrema.end(); it++) {
                Vec3i co = m_grid->XYZ3d(*it);
                printf("id %d = (%d, %d, %d)\n", *it, co[0], co[1], co[2]);
            }
        }

        INDEX_TYPE WhichExtremum(Vec3d point) {
            Vec3l closest_vertex = m_grid->Inbounds(point + 0.5); // dont know if our point is near boundary!
            //printf("(%f, %f, %f)->(%d, %d, %d)\n", point[0], point[1], point[2], closest_vertex[0], closest_vertex[1], closest_vertex[2]);
            return m_grid->Index3d(closest_vertex);
        }

        //void AddExtremum(INDEX_TYPE id) { m_extrema.insert(id); }
        virtual ADVECTION_EVENT CheckAndDoStuff(const Vec3d& old_point, const Vec3d& new_point, bool nearboundary) {
            if (!nearboundary) {
                Vec3l closest_vertex = new_point + 0.5;
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_extrema.count(id) != 0) {
                    return HIT_PREASSIGNED; // reached a critical point
                }
            }
            else {
                Vec3l closest_vertex = m_grid->Inbounds(new_point + 0.5);
                INDEX_TYPE id = m_grid->Index3d(closest_vertex);
                if (m_extrema.count(id) != 0) {
                    return HIT_PREASSIGNED; // reached a critical point
                }

            }
            return NONE;
        }

    };
}
#endif
