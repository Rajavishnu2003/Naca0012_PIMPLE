// picks the linear solver and checks that it really converged
#include "cfd/linear_solver_adapter.h"
#include "cfd/csr_tools.h"
#include "bicgstab.h"
#include "conjugate_gradient.h"
#include "steepest_descent.h"
#include "ic0_preconditioner.h"
#include "linear_algebra.h"
#include <chrono>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace cfd {
    using Clock=std::chrono::steady_clock;
    namespace {
        double seconds(Clock::time_point a) {
            return std::chrono::duration<double>(Clock::now()-a).count();
        }
        // we want |r| <= max(atol, rtol*|b|)
        double target(const Config& c,const std::vector<double>& b) {
            return std::max(c.linear_atol,c.linear_rtol*vectorNorm(b));
        }
        // the solvers only take a relative tolerance, so turn our target into one
        double solverTolerance(const Config& c,const std::vector<double>& b) {
            double n=vectorNorm(b);
            return target(c,b)/(n>0?n:1.0);
        }
        // compute the true residual again after the solve, stop the run if it is not small enough
        void verify(LinearSolveResult& r,const CSRMatrix& A,const std::vector<double>& b,const std::vector<double>& x,const Config& c) {
            r.final_residual=vectorNorm(computeResidual(A,x,b));
            r.converged=std::isfinite(r.final_residual)&&r.final_residual<=target(c,b)*1.001;
            if(!r.converged)throw std::runtime_error("Linear solve reached limit without convergence. True residual="+std::to_string(r.final_residual));
        }
    }
    DiagonalPreconditioner::DiagonalPreconditioner(const CSRMatrix& A):inverse_diagonal_(A.rows) {
        for(int i=0;i<A.rows;++i) {
            double a=A.values.at(A.diagonalIndex.at(i));
            if(!std::isfinite(a)||a<=0)throw std::invalid_argument("Nonpositive diagonal.");
            inverse_diagonal_[i]=1/a;
        }
    }
    std::vector<double> DiagonalPreconditioner::applyInverse(const std::vector<double>& r) const {
        if(r.size()!=inverse_diagonal_.size())throw std::invalid_argument("Diagonal preconditioner vector mismatch.");
        auto z=r;
        for(std::size_t i=0;i<z.size();++i)z[i]*=inverse_diagonal_[i];
        return z;
    }
    LinearSolveResult solveMomentum(const CSRMatrix& A,const std::vector<double>& b,std::vector<double>& x,const Config& c) {
        auto start=Clock::now();
        LinearSolveResult r;
        r.initial_residual=vectorNorm(computeResidual(A,x,b));
        // skip the solve if the old guess is already good enough
        if(r.initial_residual>target(c,b)) {
            DiagonalPreconditioner pc(A);
            bicgstab_solver solver(solverTolerance(c,b),c.max_linear_iterations);
            solver.solve(A,b,x,pc);
            r.iterations=solver.iterations();
        }
        verify(r,A,b,x,c);
        r.seconds=seconds(start);
        return r;
    }
    PressureLinearSolver::PressureLinearSolver(const CSRMatrix& A,const Config& c,int nx,int ny):A_(A),c_(c) {
        auto start=Clock::now();
        // the pressure can use its own (looser) tolerance
        if(c_.pressure_rtol>0)c_.linear_rtol=c_.pressure_rtol;
        // the preconditioner is built here once, the airfoil cases use cg_ic0
        if(c.pressure_solver=="cg_ic0"||c.pressure_solver=="sd_ic0")pc_=std::make_unique<IC0Preconditioner>(A);
        else if(c.pressure_solver=="bicgstab")pc_=std::make_unique<DiagonalPreconditioner>(A);
        else pc_=std::make_unique<IdentityPreconditioner>();
        if(c.pressure_solver=="multigrid")mg_=std::make_unique<PressureMultigrid>(A,nx,ny,c.mg_pre,c.mg_post);
        setup_seconds_=seconds(start);
    }
    LinearSolveResult PressureLinearSolver::solve(const std::vector<double>& b,std::vector<double>& x) const {
        auto start=Clock::now();
        LinearSolveResult r;
        r.initial_residual=vectorNorm(computeResidual(A_,x,b));
        if(r.initial_residual>target(c_,b)) {
            const double tol=solverTolerance(c_,b);
            if(mg_) {
                if(c_.mg_fmg_start)mg_->applyFMG(x,b,c_.mg_cycle[0]);
                while(r.iterations<c_.max_linear_iterations&&vectorNorm(computeResidual(A_,x,b))>target(c_,b)) {
                    mg_->applyCycle(x,b,c_.mg_cycle[0]);
                    ++r.iterations;
                }
            }
            else if(c_.pressure_solver=="bicgstab") {
                bicgstab_solver s(tol,c_.max_linear_iterations);
                s.solve(A_,b,x,*pc_);
                r.iterations=s.iterations();
            }
            else if(c_.pressure_solver=="sd_ic0") {
                steepest_descent_solver s;
                s.tolerance=tol;
                s.maxIterations=c_.max_linear_iterations;
                s.verbose=false;
                s.solve(A_,b,x,*pc_);
                r.iterations=s.iterations();
            }
            else {
                // default: conjugate gradient, the pressure matrix is symmetric
                conjugate_gradient_solver s(tol,c_.max_linear_iterations);
                s.solve(A_,b,x,*pc_);
                r.iterations=s.iterations();
            }
        }
        verify(r,A_,b,x,c_);
        r.seconds=seconds(start);
        return r;
    }
}
