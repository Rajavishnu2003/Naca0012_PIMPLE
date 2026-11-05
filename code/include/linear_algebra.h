#pragma once

// Small vector helpers used by all the Krylov solvers (CG, BiCGSTAB, SD)

#include <vector>

// Dot product: a^T * b
double dotProduct(
    const std::vector<double>& a,
    const std::vector<double>& b
);

// Euclidean / L2 norm: ||a||_2
double vectorNorm(
    const std::vector<double>& a
);

// y = y + alpha * x
void axpy(
    std::vector<double>& y,
    double alpha,
    const std::vector<double>& x
);
