#pragma once

#include "csr_matrix.h"

#include <vector>

// Gauss-Seidel sweeps on A x = b, x is updated in place.
// Used as the smoother inside the pressure multigrid.
void gaussSeidel(
    const CSRMatrix& A,
    std::vector<double>& x,
    const std::vector<double>& b,
    int sweeps
);
