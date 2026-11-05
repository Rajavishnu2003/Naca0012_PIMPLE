// the flow variables u, v, p stored on the staggered grid
#pragma once
#include "cfd/staggered_grid.h"
#include <vector>
namespace cfd {
    struct FlowFields {
        std::vector<double> u,v,p;
        explicit FlowFields(const StaggeredGrid& g):u((g.nx+1)*g.ny,0.0),v(g.nx*(g.ny+1),0.0),p(g.np(),0.0){}
    };
    // copy u or v between the full face arrays and the vector of matrix unknowns
    std::vector<double> packU(const StaggeredGrid& g,const FlowFields& f);
    std::vector<double> packV(const StaggeredGrid& g,const FlowFields& f);
    void unpackU(const StaggeredGrid& g,const std::vector<double>& u,FlowFields& f);
    void unpackV(const StaggeredGrid& g,const std::vector<double>& v,FlowFields& f);
}
