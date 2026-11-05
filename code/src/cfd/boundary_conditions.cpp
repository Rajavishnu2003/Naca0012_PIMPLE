// Boundary conditions: which type each side is, and setting the known normal velocities.
#include "cfd/boundary_conditions.h"
#include <stdexcept>
namespace cfd {
    BoundaryConditions BoundaryConditions::forCase(const Config& c) {
        BoundaryConditions bc;
        // airfoil case: inflow on the left, outflow on the right, slip walls top and bottom
        // (the cavity keeps the default walls with a moving lid)
        if(c.case_type=="external") {
            bc.west=BoundaryType::Inlet;
            bc.east=BoundaryType::Outlet;
            bc.south=bc.north=BoundaryType::Slip;
            bc.lid_speed=0.0;
            bc.inlet_speed=1.0;
        }
        return bc;
    }
    void BoundaryConditions::enforceNormalVelocities(const StaggeredGrid& g,const BodyMask& body,FlowFields& f) const {
        // left and right edge u faces: inlet gets U = 1, walls get 0, the outlet is left alone
        for(int j=0;j<g.ny;++j) {
            f.u[g.uf(0,j)]=west==BoundaryType::Inlet?inlet_speed:0.0;
            if(east!=BoundaryType::Outlet) f.u[g.uf(g.nx,j)]=0;
        }
        // bottom and top v faces are always 0 (wall or slip)
        for(int i=0;i<g.nx;++i) {
            f.v[g.vf(i,0)]=0;
            f.v[g.vf(i,g.ny)]=0;
        }
        // faces inside the airfoil are set to zero velocity
        for(std::size_t k=0;k<f.u.size();++k) if(body.u[k]) f.u[k]=0;
        for(std::size_t k=0;k<f.v.size();++k) if(body.v[k]) f.v[k]=0;
    }
    void BoundaryConditions::updateOutflow(const StaggeredGrid& g,FlowFields& f) const {
        if(east!=BoundaryType::Outlet) return;
        double inflow=0,outflow=0;
        // zero gradient: copy u from the last inner face, and add up the flow in and out
        for(int j=0;j<g.ny;++j) {
            f.u[g.uf(g.nx,j)]=f.u[g.uf(g.nx-1,j)];
            inflow+=f.u[g.uf(0,j)]*g.dy;
            outflow+=f.u[g.uf(g.nx,j)]*g.dy;
        }
        // top and bottom are walls or slip (v = 0), so the only other flux is the inlet
        // scaling the outlet so outflow = inflow makes the pressure correction RHS sum to zero,
        // otherwise the Neumann pressure problem would have no solution
        if(!(outflow>0)) throw std::runtime_error("Nonpositive outflow; the flow has reversed at the outlet.");
        for(int j=0;j<g.ny;++j) f.u[g.uf(g.nx,j)]*=inflow/outflow;
    }
}
