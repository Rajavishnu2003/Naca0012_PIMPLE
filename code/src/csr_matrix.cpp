// Basic CSR operations: y = A*x and the residual r = b - A*x
#include "csr_matrix.h"

#include <stdexcept>


std::vector<double> multiply(
    const CSRMatrix& A,
    const std::vector<double>& x)
{
    if (static_cast<int>(x.size()) != A.cols)
    {
        throw std::invalid_argument(
            "Matrix-vector size mismatch"
        );
    }

    std::vector<double> result(
        A.rows,
        0.0
    );


    // go row by row
    for (int i = 0; i < A.rows; ++i)
    {
        // nonzeros of row i sit at positions rowPtr[i] .. rowPtr[i+1]-1

        for (int k = A.rowPtr[i];
             k < A.rowPtr[i + 1];
             ++k)
        {
            const int j = A.colIndex[k];

            result[i] +=
                A.values[k] * x[j];
        }
    }

    return result;
}


std::vector<double> computeResidual(
    const CSRMatrix& A,
    const std::vector<double>& x,
    const std::vector<double>& b)
{
    if (static_cast<int>(b.size()) != A.rows)
    {
        throw std::invalid_argument(
            "RHS size does not match matrix"
        );
    }

    // First calculate A*x
    std::vector<double> Ax =
        multiply(A, x);


    std::vector<double> residual(
        A.rows,
        0.0
    );


    // r = b - A*x
    for (int i = 0; i < A.rows; ++i)
    {
        residual[i] =
            b[i] - Ax[i];
    }

    return residual;
}