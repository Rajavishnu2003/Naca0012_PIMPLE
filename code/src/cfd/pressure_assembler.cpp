// Pressure correction equation and the velocity/pressure update after solving it.
#include "cfd/pressure_assembler.h"
#include "cfd/csr_tools.h"
#include <algorithm>
#include <stdexcept>
namespace cfd {
    PressureAssembler::PressureAssembler(const StaggeredGrid& g,const BodyMask& body):A(fivePointPattern(g.nx,g.ny)),b(g.np(),0),g_(g),body_(body){}
    void PressureAssembler::assembleMatrix(const MomentumAssembler& m) {
        zeroValues(A);
        const auto& g=g_;
        // each link between two pressure cells is added once to both rows, so A stays symmetric
        auto edge=[&](int p,int q,double a) {
            if(a==0)return;  // face blocked by the body (d=0): no link between p and q
            if(!(a>0))throw std::runtime_error("Pressure conductance must be positive.");
            A.values[A.diagonalIndex[p]]+=a;
            A.values[A.diagonalIndex[q]]+=a;
            addCoefficient(A,p,q,-a);
            addCoefficient(A,q,p,-a);
        };
        // link coefficient = face area * d from the momentum equation (d = area / a_P)
        for(int j=0;j<g.ny;++j)for(int i=1;i<g.nx;++i)edge(g.p(i-1,j),g.p(i,j),g.dy*m.u.d[g.u(i,j)]);
        for(int j=1;j<g.ny;++j)for(int i=0;i<g.nx;++i)edge(g.p(i,j-1),g.p(i,j),g.dx*m.v.d[g.v(i,j)]);
        // solid cells have no links left, 1*p' = 0 keeps them out and A stays SPD (so CG works)
        for(int k=0;k<g.np();++k)if(body_.p[k])A.values[A.diagonalIndex[k]]=1.0;
        std::fill(b.begin(),b.end(),0);
        // fix p' at cell 0, otherwise the pressure system has no unique answer (only gradients matter)
        pinZeroReference(A,b,0);
    }
    // net volume flux out of each pressure cell, (u_e - u_w) dy + (v_n - v_s) dx
    std::vector<double> massImbalance(const StaggeredGrid& g,const FlowFields& f) {
        std::vector<double> b(g.np());
        for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i)
        b[g.p(i,j)]=(f.u[g.uf(i+1,j)]-f.u[g.uf(i,j)])*g.dy+(f.v[g.vf(i,j+1)]-f.v[g.vf(i,j)])*g.dx;
        return b;
    }
    // RHS of the p' equation is minus the mass imbalance of the predicted velocity
    void PressureAssembler::assembleRHS(const FlowFields& f) {
        b=massImbalance(g_,f);
        for(double& a:b)a=-a;
        b[0]=0;  // the pinned reference cell
    }
    void correctPressureVelocity(const StaggeredGrid& g,const MomentumAssembler& m,const std::vector<double>& pp,double alpha_p,FlowFields& f) {
        if(static_cast<int>(pp.size())!=g.np())throw std::invalid_argument("Wrong pressure correction size.");
        // u' = d * (p'_west - p'_east), the velocity gets the full correction so it becomes divergence free
        for(int j=0;j<g.ny;++j)for(int i=1;i<g.nx;++i)
        f.u[g.uf(i,j)]+=m.u.d[g.u(i,j)]*(pp[g.p(i-1,j)]-pp[g.p(i,j)]);
        for(int j=1;j<g.ny;++j)for(int i=0;i<g.nx;++i)
        f.v[g.vf(i,j)]+=m.v.d[g.v(i,j)]*(pp[g.p(i,j-1)]-pp[g.p(i,j)]);
        // pressure only gets alpha_p of p' (SIMPLE relaxation, alpha_p = 1 for PISO/PIMPLE)
        for(int k=0;k<g.np();++k)f.p[k]+=alpha_p*pp[k];
        // shift p so that cell 0 stays at zero
        const double ref=f.p[0];
        for(double& p:f.p)p-=ref;
    }
}
