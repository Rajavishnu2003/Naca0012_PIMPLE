// grid spacing is just domain length / number of cells
#include "cfd/staggered_grid.h"
#include <stdexcept>
namespace cfd {
    StaggeredGrid::StaggeredGrid(int x,int y,double x0_,double y0_,double lx,double ly):nx(x),ny(y),x0(x0_),y0(y0_),dx(lx/x),dy(ly/y) {
        if(x<3||y<3||x>2048||y>2048) throw std::invalid_argument("Unsupported grid size.");
        if(!(lx>0)||!(ly>0)) throw std::invalid_argument("Domain lengths must be positive.");
    }
}
