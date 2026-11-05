// Boundary types for the 4 sides of the domain.
#pragma once
#include "cfd/flow_fields.h"
#include "cfd/body_geometry.h"
namespace cfd {
    // Wall: no slip (cavity walls, the north wall can move at lid_speed)
    // Inlet: west side only, u = inlet_speed and v = 0
    // Outlet: east side only, zero gradient velocity, scaled so outflow = inflow
    // Slip: south/north only, v = 0 and no shear (far field "symmetry" boundary)
    enum class BoundaryType { Wall, Inlet, Outlet, Slip };
    struct BoundaryConditions {
        BoundaryType west=BoundaryType::Wall,east=BoundaryType::Wall,south=BoundaryType::Wall,north=BoundaryType::Wall;
        double lid_speed=1.0;
        double inlet_speed=0.0;
        static BoundaryConditions forCase(const Config& c);
        // set the known normal velocities on the domain edges, and zero velocity on body faces
        void enforceNormalVelocities(const StaggeredGrid& g,const BodyMask& body,FlowFields& f) const;
        // outlet faces copy the face just upstream, then get scaled so mass is conserved
        void updateOutflow(const StaggeredGrid& g,FlowFields& f) const;
        // known tangential wall values are taken care of in the momentum assembler (half a cell away)
    };
}
