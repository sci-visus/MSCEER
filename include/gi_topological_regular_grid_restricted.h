/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef TOPOLOGICAL_REGULAR_GRID_RESTRICTED_H
#define TOPOLOGICAL_REGULAR_GRID_RESTRICTED_H


#include <map>
#include <set>
#include "gi_basic_types.h"
#include "gi_vectors.h"
#include "gi_regular_grid.h"
#include "gi_topological_regular_grid.h"
#include "gi_labeling.h"



namespace GInt {

    class TopologicalRegularGridRestricted : virtual public TopologicalRegularGrid3D {

    protected:
        const DenseLabeling<char> *m_restriction;

    public:
        TopologicalRegularGridRestricted(RegularGrid3D* base_grid) :
            TopologicalRegularGrid3D(base_grid), m_restriction(0) {

            printf(" -- Created TopologicalRegularGridRestricted\n");
        }

        void set_restriction(const DenseLabeling<char> *restriction) {
            //printf("TopologicalRegularGridRestricted::set_restriction(%p)\n", restriction);
            m_restriction = restriction;
        }
		char restrictionLabel(INDEX_TYPE cellid) const {
			return m_restriction->GetLabel(cellid);
		}

        BOUNDARY_TYPE boundaryValue(INDEX_TYPE cellid) const {

            // TopologicalRegularGrid3D::boundaryValue returns 0,1,2,3
            // m_restriction has labels 0 or 1.. so this function returns (0,1,2,3) or (4,5,6,7)

            //printf("TopologicalRegularGridRestricted::boundaryValue()\n");
            return (m_restriction->GetLabel (cellid) * (maxDim()+1)) + TopologicalRegularGrid3D::boundaryValue (cellid);
        }

		BOUNDARY_TYPE domainBoundaryValue(INDEX_TYPE cellid) const {

			// TopologicalRegularGrid3D::boundaryValue returns 0,1,2,3
			// m_restriction has labels 0 or 1.. so this function returns (0,1,2,3) or (4,5,6,7)

			//printf("TopologicalRegularGridRestricted::boundaryValue()\n");
			return TopologicalRegularGrid3D::boundaryValue(cellid);
		}
    };
}


#endif
