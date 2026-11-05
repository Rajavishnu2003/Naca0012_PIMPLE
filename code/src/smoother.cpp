// Gauss-Seidel smoother, works straight on the CSR rows
#include "smoother.h"

#include <cmath>
#include <stdexcept>

void gaussSeidel(
    const CSRMatrix& A,
    std::vector<double>& x,
    const std::vector<double>& b,
    int sweeps)
{
    // quick checks on sizes and the diagonal index first
    if (A.rows != A.cols)
    {
        throw std::invalid_argument(
            "Gauss-Seidel requires a square matrix."
        );
    }

    if (static_cast<int>(x.size()) != A.rows)
    {
        throw std::invalid_argument(
            "Solution vector size does not match matrix size."
        );
    }

    if (static_cast<int>(b.size()) != A.rows)
    {
        throw std::invalid_argument(
            "RHS vector size does not match matrix size."
        );
    }

    if (sweeps < 0)
    {
        throw std::invalid_argument(
            "Number of Gauss-Seidel sweeps cannot be negative."
        );
    }

    if (static_cast<int>(A.diagonalIndex.size()) != A.rows)
    {
        throw std::invalid_argument(
            "CSR matrix does not contain a valid diagonalIndex array."
        );
    }

    // the sweeps
    for (int sweep = 0; sweep < sweeps; ++sweep)
    {
        // go through the rows in order
        for (int i = 0; i < A.rows; ++i)
        {
            const int diagonalPosition =
                A.diagonalIndex[i];

            if (diagonalPosition < 0 ||
                diagonalPosition >= static_cast<int>(A.values.size()))
            {
                throw std::runtime_error(
                    "Missing diagonal entry in CSR matrix."
                );
            }

            const double diagonal =
                A.values[diagonalPosition];

            if (std::abs(diagonal) < 1.0e-14)
            {
                throw std::runtime_error(
                    "Zero or near-zero diagonal encountered in Gauss-Seidel."
                );
            }

            double offDiagonalSum = 0.0;

            // sum A(i,j) x_j over the nonzeros of row i
            for (int k = A.rowPtr[i];
                 k < A.rowPtr[i + 1];
                 ++k)
            {
                const int j = A.colIndex[k];

                // leave out A(i,i)
                if (j == i)
                {
                    continue;
                }

                offDiagonalSum +=
                    A.values[k] * x[j];
            }

            // x_i = (b_i - sum_{j != i} A_ij x_j) / A_ii
            // x is overwritten straight away, so for j < i we already use the
            // new values from this sweep and for j > i the old ones
            x[i] =
                (b[i] - offDiagonalSum)
                / diagonal;
        }
    }
}
