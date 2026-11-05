// small helper functions for CSR sparse matrices
#pragma once
#include "csr_matrix.h"
#include <string>
namespace cfd {
    // empty 5 point stencil pattern (P, E, W, N, S), values filled later
    CSRMatrix fivePointPattern(int nx,int ny);
    int coefficientPosition(const CSRMatrix& A,int row,int col);
    void addCoefficient(CSRMatrix& A,int row,int col,double value);
    void zeroValues(CSRMatrix& A);
    void validateCSR(const CSRMatrix& A);
    // fix p = 0 at the reference cell: clear its row and column and put 1 on the diagonal.
    // the zeros stay in the pattern, nothing is removed
    void pinZeroReference(CSRMatrix& A,std::vector<double>& b,int reference=0);
    // used to build the coarse multigrid matrices
    CSRMatrix transposeCSR(const CSRMatrix& A);
    CSRMatrix productCSR(const CSRMatrix& A,const CSRMatrix& B);
    void dumpCSR(const CSRMatrix& A,const std::vector<double>& b,const std::string& prefix);
    double normInf(const std::vector<double>& x);
}
