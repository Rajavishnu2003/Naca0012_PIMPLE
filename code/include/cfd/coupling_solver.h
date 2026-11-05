// Pressure-velocity coupling solver (SIMPLE, PISO and PIMPLE share the same parts).
#pragma once
#include "cfd/convergence_monitor.h"
#include "cfd/pressure_assembler.h"
#include "cfd/linear_solver_adapter.h"
#include "cfd/body_forces.h"
#include <fstream>
#include <cstdint>
namespace cfd {
    // what we know at the end of a run, it also goes into summary.json
    struct RunSummary {
        bool converged=false,completed=false;
        int steps=0;
        double time=0,wall_seconds=0,momentum_solve_seconds=0,pressure_solve_seconds=0,pressure_setup_seconds=0;
        std::int64_t momentum_iterations=0,pressure_iterations=0,pressure_solves=0,momentum_solves=0;
        FlowResiduals residuals;
        BodyForces forces;
        // forces are only for the external flow (airfoil) case
    };
    // the discretisation and linear solvers are shared, the three algorithms only differ
    // in how they call the predictor and corrector (see the private functions below)
    class CouplingSolver {
        public:
        explicit CouplingSolver(Config config,bool write_output=true);
        RunSummary run();
        const FlowFields& fields() const {
            return f_;
        }
        const StaggeredGrid& grid() const {
            return g_;
        }
        const BodyMask& body() const {
            return body_;
        }
        private:
        Config c_;
        StaggeredGrid g_;
        BodyMask body_;
        BoundaryConditions bc_;
        FlowFields f_;
        // residual_workspace_ is a second assembler only used to compute the residuals
        MomentumAssembler momentum_,residual_workspace_;
        PressureAssembler pressure_;
        RunSummary summary_;
        bool output_;
        void momentumPredictor(const FlowFields& old_time,double dt,double alpha);
        void pressureCorrect(PressureLinearSolver& solver,double alpha_p);
        void simpleIteration(const FlowFields& old_iter);
        void pisoTimeStep(const FlowFields& old_time,double dt);
        void pimpleTimeStep(const FlowFields& old_time,double dt);
        void transientOuterIteration(const FlowFields& old_time,double dt,double alpha);
    };
}
