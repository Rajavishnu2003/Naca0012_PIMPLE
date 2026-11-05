// IC0 preconditioner: checks A, builds L row by row, then applies
// z = (L L^T)^-1 r with one forward and one backward solve.
#include "ic0_preconditioner.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>


namespace
{
    constexpr double symmetryTolerance = 1e-12;


    // check the CSR arrays and find the diagonals ourselves,
    // since A.diagonalIndex may be missing or out of date
    std::vector<int> validateMatrix(const CSRMatrix& A)
    {
        if (A.rows <= 0 || A.rows != A.cols)
        {
            throw std::invalid_argument(
                "IC(0) requires a nonempty square matrix."
            );
        }

        if (A.rows == std::numeric_limits<int>::max()
            || A.values.size()
                > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::invalid_argument(
                "IC(0) matrix exceeds the supported CSR integer index range."
            );
        }

        const int n = A.rows;
        const int nnz = static_cast<int>(A.values.size());

        if (A.values.size() != A.colIndex.size()
            || A.rowPtr.size() != static_cast<std::size_t>(n) + 1)
        {
            throw std::invalid_argument(
                "IC(0) received inconsistent CSR array sizes."
            );
        }

        if (A.rowPtr.front() != 0 || A.rowPtr.back() != nnz)
        {
            throw std::invalid_argument(
                "IC(0) requires rowPtr[0] = 0 and rowPtr[n] = nnz."
            );
        }

        std::vector<int> diagonalLocations(n, -1);

        for (int i = 0; i < n; ++i)
        {
            const int start = A.rowPtr[i];
            const int end = A.rowPtr[i + 1];

            if (start < 0 || end < start || end > nnz)
            {
                throw std::invalid_argument(
                    "IC(0) received invalid row pointers at row "
                    + std::to_string(i) + "."
                );
            }

            int previousColumn = -1;

            for (int position = start; position < end; ++position)
            {
                const int column = A.colIndex[position];

                if (column < 0 || column >= n || column <= previousColumn)
                {
                    throw std::invalid_argument(
                        "IC(0) requires valid, strictly increasing, unique "
                        "column indices in every CSR row; check row "
                        + std::to_string(i) + "."
                    );
                }

                if (!std::isfinite(A.values[position]))
                {
                    throw std::invalid_argument(
                        "IC(0) requires finite matrix coefficients."
                    );
                }

                if (column == i)
                {
                    diagonalLocations[i] = position;
                }

                previousColumn = column;
            }

            if (diagonalLocations[i] < 0
                || A.values[diagonalLocations[i]] <= 0.0)
            {
                throw std::invalid_argument(
                    "IC(0) requires an explicit positive diagonal at row "
                    + std::to_string(i) + "."
                );
            }
        }

        // symmetry check: compare A(i,j) with A(j,i), missing entry counts as zero.
        // both are divided by the bigger one first, so the test is relative
        // and big values can't overflow
        for (int i = 0; i < n; ++i)
        {
            for (int position = A.rowPtr[i];
                 position < A.rowPtr[i + 1]; ++position)
            {
                const int j = A.colIndex[position];

                if (j == i)
                {
                    continue;
                }

                const auto rowStart = A.colIndex.begin() + A.rowPtr[j];
                const auto rowEnd = A.colIndex.begin() + A.rowPtr[j + 1];
                const auto counterpart = std::lower_bound(rowStart, rowEnd, i);

                double aji = 0.0;

                if (counterpart != rowEnd && *counterpart == i)
                {
                    const std::size_t index =
                        static_cast<std::size_t>(counterpart - A.colIndex.begin());
                    aji = A.values[index];
                }

                const double aij = A.values[position];
                const double scale = std::max(std::abs(aij), std::abs(aji));

                if (scale > 0.0
                    && std::abs(aij / scale - aji / scale) > symmetryTolerance)
                {
                    throw std::invalid_argument(
                        "IC(0) requires numerical symmetry; mismatch at ("
                        + std::to_string(i) + ", " + std::to_string(j) + ")."
                    );
                }
            }
        }

        return diagonalLocations;
    }


    // turn the long double result back to double, stop if it blew up
    double checkedSolveValue(long double value, int row, const char* stage)
    {
        const long double largestDouble = std::numeric_limits<double>::max();

        if (!std::isfinite(value) || std::abs(value) > largestDouble)
        {
            throw std::runtime_error(
                std::string("IC(0) produced a nonfinite value during ")
                + stage + " at row " + std::to_string(row) + "."
            );
        }

        return static_cast<double>(value);
    }
}


IC0Preconditioner::IC0Preconditioner(const CSRMatrix& A)
{
    const std::vector<int> inputDiagonal = validateMatrix(A);
    const int n = A.rows;

    lowerFactor = CSRMatrix(n, n);

    // count how many entries the lower triangle has, so we can reserve once
    std::size_t lowerPositions = 0;

    for (int i = 0; i < n; ++i)
    {
        lowerPositions +=
            static_cast<std::size_t>(inputDiagonal[i] - A.rowPtr[i] + 1);
    }

    lowerFactor.values.reserve(lowerPositions);
    lowerFactor.colIndex.reserve(lowerPositions);

    // copy the lower triangle of A (up to and including the diagonal) into L
    for (int i = 0; i < n; ++i)
    {
        lowerFactor.rowPtr[i] =
            static_cast<int>(lowerFactor.values.size());

        for (int position = A.rowPtr[i];
             position <= inputDiagonal[i]; ++position)
        {
            lowerFactor.colIndex.push_back(A.colIndex[position]);
            lowerFactor.values.push_back(A.values[position]);
        }

        lowerFactor.diagonalIndex[i] =
            static_cast<int>(lowerFactor.values.size()) - 1;
    }

    lowerFactor.rowPtr[n] =
        static_cast<int>(lowerFactor.values.size());

    // the factorisation itself, row by row:
    // L(i,j) = (A(i,j) - sum_{k<j} L(i,k)*L(j,k)) / L(j,j)
    // the sum walks rows i and j together and only multiplies matching columns
    for (int i = 0; i < n; ++i)
    {
        const int diagonal = lowerFactor.diagonalIndex[i];

        for (int position = lowerFactor.rowPtr[i];
             position < diagonal; ++position)
        {
            const int j = lowerFactor.colIndex[position];
            int left = lowerFactor.rowPtr[i];
            int right = lowerFactor.rowPtr[j];
            const int jDiagonal = lowerFactor.diagonalIndex[j];
            long double sum = 0.0L;

            while (left < position && right < jDiagonal)
            {
                const int leftColumn = lowerFactor.colIndex[left];
                const int rightColumn = lowerFactor.colIndex[right];

                if (leftColumn == rightColumn)
                {
                    sum += static_cast<long double>(lowerFactor.values[left])
                        * lowerFactor.values[right];
                    ++left;
                    ++right;
                }
                else if (leftColumn < rightColumn)
                {
                    ++left;
                }
                else
                {
                    ++right;
                }
            }

            const long double value =
                (static_cast<long double>(lowerFactor.values[position]) - sum)
                / lowerFactor.values[jDiagonal];

            lowerFactor.values[position] =
                checkedSolveValue(value, i, "factor setup");
        }

        // diagonal: L(i,i) = sqrt(A(i,i) - sum_k L(i,k)^2)
        long double squareSum = 0.0L;

        for (int position = lowerFactor.rowPtr[i];
             position < diagonal; ++position)
        {
            const long double value = lowerFactor.values[position];
            squareSum += value * value;
        }

        const long double pivot =
            static_cast<long double>(lowerFactor.values[diagonal]) - squareSum;

        if (!std::isfinite(pivot) || pivot <= 0.0L)
        {
            throw std::runtime_error(
                "IC(0) breakdown: nonpositive or nonfinite incomplete pivot "
                "at row " + std::to_string(i)
                + ". Incomplete-factorization breakdown does not by itself "
                  "prove that A is indefinite."
            );
        }

        lowerFactor.values[diagonal] =
            checkedSolveValue(std::sqrt(pivot), i, "factor setup");

        if (lowerFactor.values[diagonal] <= 0.0)
        {
            throw std::runtime_error(
                "IC(0) breakdown: the factor diagonal underflowed at row "
                + std::to_string(i) + "."
            );
        }
    }

    buildTranspose();
}


// store L^T as its own CSR matrix so the backward solve can go row by row
void IC0Preconditioner::buildTranspose()
{
    const int n = lowerFactor.rows;
    const int nnz = lowerFactor.nnz();

    transposedLowerFactor = CSRMatrix(n, n);
    transposedLowerFactor.values.resize(nnz);
    transposedLowerFactor.colIndex.resize(nnz);

    // count entries per column of L, then a running sum gives rowPtr of L^T
    for (int column : lowerFactor.colIndex)
    {
        ++transposedLowerFactor.rowPtr[column + 1];
    }

    for (int i = 0; i < n; ++i)
    {
        transposedLowerFactor.rowPtr[i + 1] +=
            transposedLowerFactor.rowPtr[i];
    }

    std::vector<int> nextPosition = transposedLowerFactor.rowPtr;

    // going through the rows of L in order means each row of L^T comes out sorted,
    // with the diagonal first and then the upper entries
    for (int i = 0; i < n; ++i)
    {
        for (int position = lowerFactor.rowPtr[i];
             position < lowerFactor.rowPtr[i + 1]; ++position)
        {
            const int column = lowerFactor.colIndex[position];
            const int destination = nextPosition[column]++;

            transposedLowerFactor.colIndex[destination] = i;
            transposedLowerFactor.values[destination] =
                lowerFactor.values[position];

            if (column == i)
            {
                transposedLowerFactor.diagonalIndex[column] = destination;
            }
        }
    }
}


std::vector<double> IC0Preconditioner::applyInverse(
    const std::vector<double>& rhs) const
{
    const int n = size();

    if (rhs.size() != static_cast<std::size_t>(n))
    {
        throw std::invalid_argument(
            "IC(0) right-hand-side size does not match its factors."
        );
    }

    for (double value : rhs)
    {
        if (!std::isfinite(value))
        {
            throw std::invalid_argument(
                "IC(0) requires a finite right-hand-side vector."
            );
        }
    }

    std::vector<double> y(n, 0.0);
    std::vector<double> z(n, 0.0);

    // forward substitution, L*y = rhs
    for (int i = 0; i < n; ++i)
    {
        const int diagonal = lowerFactor.diagonalIndex[i];
        long double value = rhs[i];

        for (int position = lowerFactor.rowPtr[i];
             position < diagonal; ++position)
        {
            value -= static_cast<long double>(lowerFactor.values[position])
                * y[lowerFactor.colIndex[position]];
        }

        y[i] = checkedSolveValue(
            value / lowerFactor.values[diagonal], i, "forward substitution"
        );
    }

    // backward substitution, L^T*z = y
    // so z = (L L^T)^-1 rhs, and we never form an inverse
    for (int i = n - 1; i >= 0; --i)
    {
        const int diagonal = transposedLowerFactor.diagonalIndex[i];
        long double value = y[i];

        for (int position = diagonal + 1;
             position < transposedLowerFactor.rowPtr[i + 1]; ++position)
        {
            value -=
                static_cast<long double>(transposedLowerFactor.values[position])
                * z[transposedLowerFactor.colIndex[position]];
        }

        z[i] = checkedSolveValue(
            value / transposedLowerFactor.values[diagonal],
            i,
            "backward substitution"
        );
    }

    return z;
}


int IC0Preconditioner::size() const
{
    return lowerFactor.rows;
}


int IC0Preconditioner::factorNonzeros() const
{
    return lowerFactor.nnz();
}
