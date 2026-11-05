#pragma once

#include "csr_matrix.h"
#include "preconditioner.h"

#include <vector>


// IC0 = incomplete Cholesky, M = L*L^T with the same sparsity as A, no fill in.
// Used as the preconditioner for CG on the pressure correction equation.
//
// A has to be square, symmetric, with sorted unique columns per row and a
// positive diagonal stored in every row. The factor is built once in the
// constructor and A is not changed (its diagonalIndex is not used).
//
// L only keeps the positions that A already has in its lower triangle
// (stored zeros too). No shift or pivoting. If a pivot goes <= 0 we throw,
// but that alone does not mean A is not positive definite.
class IC0Preconditioner final : public Preconditioner
{
public:
    explicit IC0Preconditioner(const CSRMatrix& A);

    std::vector<double> applyInverse(
        const std::vector<double>& rhs
    ) const override;

    int size() const;

    // nonzeros stored in L (diagonal and kept zeros included), L^T not counted again
    int factorNonzeros() const;

private:
    CSRMatrix lowerFactor;
    CSRMatrix transposedLowerFactor;

    void buildTranspose();
};
