// pack and unpack: the matrix only has interior faces, the full arrays also keep the boundary faces
#include "cfd/flow_fields.h"
#include <stdexcept>
namespace cfd {
    std::vector<double> packU(const StaggeredGrid& g,const FlowFields& f) {
        std::vector<double> x(g.nu());
        for(int j=0;j<g.ny;++j) for(int i=1;i<g.nx;++i) x[g.u(i,j)]=f.u[g.uf(i,j)];
        return x;
    }
    std::vector<double> packV(const StaggeredGrid& g,const FlowFields& f) {
        std::vector<double> x(g.nv());
        for(int j=1;j<g.ny;++j) for(int i=0;i<g.nx;++i) x[g.v(i,j)]=f.v[g.vf(i,j)];
        return x;
    }
    void unpackU(const StaggeredGrid& g,const std::vector<double>& x,FlowFields& f) {
        if(static_cast<int>(x.size())!=g.nu()) throw std::invalid_argument("Wrong u vector size.");
        for(int j=0;j<g.ny;++j) for(int i=1;i<g.nx;++i) f.u[g.uf(i,j)]=x[g.u(i,j)];
    }
    void unpackV(const StaggeredGrid& g,const std::vector<double>& x,FlowFields& f) {
        if(static_cast<int>(x.size())!=g.nv()) throw std::invalid_argument("Wrong v vector size.");
        for(int j=1;j<g.ny;++j) for(int i=0;i<g.nx;++i) f.v[g.vf(i,j)]=x[g.v(i,j)];
    }
}
