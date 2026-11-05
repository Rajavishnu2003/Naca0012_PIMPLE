// One SIMPLE iteration (steady, no time term).
// Used by the cavity cases; the airfoil cases use PIMPLE.
#include "cfd/coupling_solver.h"
namespace cfd {
    void CouplingSolver::simpleIteration(const FlowFields& previous) {
        // solve u and v momentum with the current pressure, dt = 0 means steady, alpha_u relaxes it
        momentumPredictor(previous,0.0,c_.alpha_u);
        // build the pressure correction matrix from the momentum diagonals (the d = area/a_P values)
        pressure_.assembleMatrix(momentum_);
        PressureLinearSolver solver(pressure_.A,c_,g_.nx,g_.ny);
        summary_.pressure_setup_seconds+=solver.setupSeconds();
        // one pressure correction, only alpha_p of p' is added to the pressure
        pressureCorrect(solver,c_.alpha_p);
    }
}
