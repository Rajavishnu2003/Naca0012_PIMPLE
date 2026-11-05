// connects the CFD code to the linear solvers (CG, BiCGSTAB, steepest descent, multigrid)
#pragma once
#include "cfd/config.h"
#include "cfd/pressure_multigrid.h"
#include "preconditioner.h"
#include <memory>
namespace cfd {
    struct LinearSolveResult {
        bool converged=false;
        int iterations=0;
        double initial_residual=0,final_residual=0,seconds=0;
    };
    // Jacobi preconditioner, just divide by the diagonal
    class DiagonalPreconditioner final:public Preconditioner {
        public:
        explicit DiagonalPreconditioner(const CSRMatrix& A);
        std::vector<double> applyInverse(const std::vector<double>& r) const override;
        private: std::vector<double> inverse_diagonal_;
    };
    // momentum matrix is not symmetric, so BiCGSTAB is used for it
    LinearSolveResult solveMomentum(const CSRMatrix& A,const std::vector<double>& b,std::vector<double>& x,const Config& c);
    // set up once for each new pressure matrix. in PISO only the rhs changes, so the same object is used again
    class PressureLinearSolver {
        public:
        PressureLinearSolver(const CSRMatrix& A,const Config& c,int nx,int ny);
        LinearSolveResult solve(const std::vector<double>& b,std::vector<double>& x) const;
        double setupSeconds() const {
            return setup_seconds_;
        }
        private:
        const CSRMatrix& A_;
        Config c_;
        std::unique_ptr<Preconditioner> pc_;
        std::unique_ptr<PressureMultigrid> mg_;
        double setup_seconds_=0;
    };
}
