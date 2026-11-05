#pragma once

// Sparse matrix in CSR format. All the big systems (momentum and pressure)
// are stored like this, so we never build a full nx*ny by nx*ny matrix.

#include <vector>

// CSR keeps only the nonzeros: the values, their column index, and where each row starts
struct CSRMatrix
{
    int rows = 0;
    int cols = 0;

    // All non-zero matrix values
    std::vector<double> values;

    // Column number corresponding to each value
    std::vector<int> colIndex;

    // rowPtr[i] = starting location of row i
    // Size = rows + 1
    std::vector<int> rowPtr;

    // position of A[i][i] inside values, -1 if row i has no diagonal
    // handy for Gauss-Seidel and for changing the diagonal directly
    std::vector<int> diagonalIndex;


    // Default constructor
    CSRMatrix() = default;


    // Constructor
    CSRMatrix(int numberOfRows, int numberOfCols)
        : rows(numberOfRows),
          cols(numberOfCols),
          rowPtr(numberOfRows + 1, 0),
          diagonalIndex(numberOfRows, -1)
    {
    }


    // Number of non-zero elements
    int nnz() const
    {
        return static_cast<int>(values.size());
    }
};


// Matrix-vector multiplication
// result = A * x
std::vector<double> multiply(
    const CSRMatrix& A,
    const std::vector<double>& x
);


// Residual calculation
// r = b - A*x
std::vector<double> computeResidual(
    const CSRMatrix& A,
    const std::vector<double>& x,
    const std::vector<double>& b
);