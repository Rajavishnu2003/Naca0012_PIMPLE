// airfoil / cylinder shape and the solid mask on the grid
#pragma once
#include "cfd/staggered_grid.h"
#include "cfd/config.h"
#include <vector>
namespace cfd {
    // the body is made of grid cells (staircase / blocking off).
    // a pressure cell is solid if its centre is inside the body.
    // a velocity face is solid if the cell on either side is solid,
    // so the body surface always follows the cell faces, that is the staircase.
    struct BodyMask {
        std::vector<char> p,u,v;
        // same indexing as FlowFields: g.p(i,j), g.uf(i,j), g.vf(i,j)
        int solid_cells=0;
        explicit BodyMask(const StaggeredGrid& g);
        bool cell(const StaggeredGrid& g,int i,int j) const {
            return p[g.p(i,j)];
        }
        bool uFace(const StaggeredGrid& g,int i,int j) const {
            return u[g.uf(i,j)];
        }
        bool vFace(const StaggeredGrid& g,int i,int j) const {
            return v[g.vf(i,j)];
        }
    };
    // NACA 4 digit airfoil in its own frame, leading edge at (0,0) and trailing edge at (1,0)
    bool insideNaca4(const std::string& code,double x,double y);
    // is the point (x,y) inside the body, in domain coordinates
    bool insideBody(const Config& c,double x,double y);
    // builds the mask, it stays empty for the cavity or body = none
    BodyMask buildBodyMask(const StaggeredGrid& g,const Config& c);
}
