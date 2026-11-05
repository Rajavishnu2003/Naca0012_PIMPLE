#pragma once

#include "csr_matrix.h"
#include "linear_algebra.h"
#include "preconditioner.h"
#include <vector>

// (Preconditioned) steepest descent for SPD systems. Much slower than CG,
// it is here as a baseline to compare the solvers against.
class steepest_descent_solver{

    public : 

    // Convergence criterion: ||r|| / ||b|| <= tolerance
    double tolerance = 1e-8;

    int maxIterations = 1000000;

    // print the convergence message at the end, the comparison driver turns it off
    bool verbose = true;


    double computeAlpha(const CSRMatrix& A,
        const std::vector<double>& r);

    // preconditioned version: alpha = (r^T z) / (z^T A z) with z = M^-1 r
    double computeAlpha(const CSRMatrix& A,
        const std::vector<double>& r,
        const std::vector<double>& z);

    void updatePhi( std::vector<double>& x, 
         double alpha, const std::vector<double>& r);    

    void solve(
    const CSRMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x);

    void solve(
        const CSRMatrix& A,
        const std::vector<double>& b,
        std::vector<double>& x,
        const Preconditioner& preconditioner);

    int iterations() const { return iters_; }


private:
    int iters_ = 0;

};

