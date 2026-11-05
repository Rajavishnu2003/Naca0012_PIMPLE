// Pressure correction (p') equation, built from the momentum d coefficients.
#pragma once
#include "cfd/momentum_assembler.h"
namespace cfd {
    class PressureAssembler {
        public:
        CSRMatrix A;
        std::vector<double> b;
        PressureAssembler(const StaggeredGrid& g,const BodyMask& body);
        // matrix only changes when the momentum diagonals change
        void assembleMatrix(const MomentumAssembler& momentum);
        // RHS = minus the mass imbalance of the predicted velocity
        void assembleRHS(const FlowFields& predicted);
        private:
        const StaggeredGrid& g_;
        const BodyMask& body_;
    };
    // net volume flux out of each pressure cell (rho = 1)
    std::vector<double> massImbalance(const StaggeredGrid& g,const FlowFields& f);
    // add the p' correction to u, v and p
    void correctPressureVelocity(const StaggeredGrid& g,const MomentumAssembler& m,
    const std::vector<double>& pp,double alpha_p,FlowFields& f);
}
