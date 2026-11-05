// Momentum equations for u and v (matrix A, right hand side b, and d = area/a_P per face).
#pragma once
#include "cfd/flow_fields.h"
#include "cfd/boundary_conditions.h"
#include "cfd/config.h"
#include "cfd/convection_scheme.h"
#include "csr_matrix.h"
namespace cfd {
    // one linear system: A u = b, plus d that the pressure correction uses
    struct MomentumSystem {
        CSRMatrix A;
        std::vector<double> b,d;
        SchemeStatistics stats;
        MomentumSystem(int nx,int ny);
    };
    class MomentumAssembler {
        public:
        MomentumAssembler(const StaggeredGrid& g,const BodyMask& body);
        MomentumSystem u,v;
        // convecting is the latest (Picard) velocity used for the fluxes, old_time is the fixed time level n
        // dt = 0 gives the steady equations, alpha = 1 gives the matrix without relaxation
        void assemble(const FlowFields& convecting,const FlowFields& old_time,
        const BoundaryConditions& bc,const Config& c,double dt,double alpha);
        private:
        const StaggeredGrid& g_;
        const BodyMask& body_;
    };
    // body faces get an identity row (u = 0) and d = 0, so the pressure correction never moves them
    // updateMomentumH is the neighbour correction with a frozen matrix, used after the first PISO corrector
    // it recalculates H(u), it is not just the same continuity solve again
    void updateMomentumH(const StaggeredGrid& g,const MomentumAssembler& m,
    const std::vector<double>& assembly_pressure,FlowFields& f);
}
