// Dot product, 2-norm and axpy. Sizes are checked so a wrong vector length
// fails right away instead of reading past the end.
#include "linear_algebra.h"

#include <cmath>
#include <stdexcept>

double dotProduct(
    const std::vector<double>& a,
    const std::vector<double>& b)
{
    if (a.size() != b.size())
    {
        throw std::invalid_argument(
            "dotProduct: vector sizes do not match"
        );
    }

    double result = 0.0;

    for (std::size_t i = 0; i < a.size(); ++i)
    {
        result += a[i] * b[i];
    }

    return result;
}


double vectorNorm(
    const std::vector<double>& a)
{
    return std::sqrt(dotProduct(a, a));
}

void axpy(
    std::vector<double>& y,
    double alpha,
    const std::vector<double>& x)
{
    if (y.size() != x.size())
    {
        throw std::invalid_argument(
            "axpy: vector sizes do not match"
        );
    }

    for (std::size_t i = 0; i < y.size(); ++i)
    {
        y[i] += alpha * x[i];
    }
}