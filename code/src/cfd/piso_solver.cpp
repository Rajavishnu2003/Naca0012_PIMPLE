// PISO time step, and the outer iteration that PIMPLE also reuses.
#include "cfd/coupling_solver.h"
namespace cfd {
    // one outer pass: momentum predictor, then a few pressure corrections with the same matrix
    void CouplingSolver::transientOuterIteration(const FlowFields& old_time,double dt,double alpha) {
        // keep the pressure that went into the momentum source, the later correctors need it
        const auto assembly_pressure=f_.p;
        momentumPredictor(old_time,dt,alpha);
        // pressure matrix is built once here and reused by every corrector below
        pressure_.assembleMatrix(momentum_);
        PressureLinearSolver solver(pressure_.A,c_,g_.nx,g_.ny);
        summary_.pressure_setup_seconds+=solver.setupSeconds();
        for(int k=0;k<c_.pressure_correctors;++k) {
            // from the 2nd corrector on, update the neighbour part H(u) with the corrected velocity first
            if(k>0)updateMomentumH(g_,momentum_,assembly_pressure,f_);
            // full pressure correction (no alpha_p) like in normal PISO
            pressureCorrect(solver,1.0);
        }
    }
    void CouplingSolver::pisoTimeStep(const FlowFields& old_time,double dt) {
        // one momentum assembly and solve per step, then 2 or more neighbour/pressure corrections
        // PISO does not use the SIMPLE relaxation, so alpha_u = alpha_p = 1 here
        transientOuterIteration(old_time,dt,1.0);
    }
}
