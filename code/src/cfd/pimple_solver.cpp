// PIMPLE time step = PISO inside a few outer (Picard) loops.
// This is what the airfoil cases run with.
#include "cfd/coupling_solver.h"
namespace cfd {
    void CouplingSolver::pimpleTimeStep(const FlowFields& old_time,double dt) {
        // each outer loop updates convection and pressure again, old_time (time level n) stays fixed
        // alpha_u under-relaxes the momentum, the airfoil cases use 2 outer loops with alpha_u = 0.7
        for(int outer=0;outer<c_.outer_correctors;++outer)
        transientOuterIteration(old_time,dt,c_.alpha_u);
    }
}
