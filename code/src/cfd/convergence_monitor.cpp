// Residuals that we log every step, and the check for "steady state reached".
#include "cfd/convergence_monitor.h"
#include "cfd/pressure_assembler.h"
#include "cfd/csr_tools.h"
#include <cmath>
#include <algorithm>
namespace cfd {
    FlowResiduals evaluateResiduals(const StaggeredGrid& g,const FlowFields& f,const FlowFields& old,const BoundaryConditions& bc,const Config& c,MomentumAssembler& m,double dt) {
        // build the steady momentum equations again without relaxation (dt = 0, alpha = 1),
        // not the relaxed matrix we just solved, so the residual means the real equation
        m.assemble(f,f,bc,c,0.0,1.0);
        FlowResiduals r;
        // max |b - A x| per unit volume for u and v
        r.u_momentum=normInf(computeResidual(m.u.A,packU(g,f),m.u.b))/g.volume();
        r.v_momentum=normInf(computeResidual(m.v.A,packV(g,f),m.v.b))/g.volume();
        // continuity: largest mass imbalance of any cell
        r.divergence=normInf(massImbalance(g,f))/g.volume();
        // biggest change in velocity since the last step
        for(std::size_t i=0;i<f.u.size();++i)r.velocity_change=std::max(r.velocity_change,std::abs(f.u[i]-old.u[i]));
        for(std::size_t i=0;i<f.v.size();++i)r.velocity_change=std::max(r.velocity_change,std::abs(f.v[i]-old.v[i]));
        // max Courant number dt*(|u|/dx + |v|/dy) using cell centre velocities
        if(dt>0)for(int j=0;j<g.ny;++j)for(int i=0;i<g.nx;++i) {
            const double co=dt*(0.5*(std::abs(f.u[g.uf(i,j)])+std::abs(f.u[g.uf(i+1,j)]))/g.dx+
            0.5*(std::abs(f.v[g.vf(i,j)])+std::abs(f.v[g.vf(i,j+1)]))/g.dy);
            r.courant=std::max(r.courant,co);
        }
        // how many faces used central and how many upwind (from this assembly)
        r.schemes=m.u.stats;
        r.schemes.merge(m.v.stats);
        return r;
    }
    // steady when momentum, continuity and the change per unit time are all under the tolerances
    bool steadyConverged(const FlowResiduals& r,const Config& c,double dt) {
        const double change=dt>0?r.velocity_change/dt:r.velocity_change;
        return std::max(r.u_momentum,r.v_momentum)<=c.momentum_tolerance&&r.divergence<=c.continuity_tolerance&&change<=c.change_tolerance;
    }
}
