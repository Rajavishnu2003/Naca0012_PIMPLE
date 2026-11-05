// Residuals printed and logged every step, and the steady state check.
#pragma once
#include "cfd/momentum_assembler.h"
namespace cfd {
    // max norms per unit volume, plus the largest velocity change and Courant number
    struct FlowResiduals {
        double u_momentum=0,v_momentum=0,divergence=0,velocity_change=0,courant=0;
        SchemeStatistics schemes;
    };
    FlowResiduals evaluateResiduals(const StaggeredGrid& g,const FlowFields& f,const FlowFields& previous,
    const BoundaryConditions& bc,const Config& c,MomentumAssembler& workspace,double dt);
    bool steadyConverged(const FlowResiduals& r,const Config& c,double dt);
}
