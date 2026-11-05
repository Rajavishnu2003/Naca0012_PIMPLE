// Preconditioned CG. Without a preconditioner it just uses M = I.
#include "conjugate_gradient.h"

#include <cmath>
#include <stdexcept>

namespace
{
// check sizes, CSR arrays and that nothing is inf/nan before we start
void validateSolve(
    const CSRMatrix& A,
    const std::vector<double>& b,
    const std::vector<double>& x,
    double tolerance,
    int maxIterations)
{
    if (!std::isfinite(tolerance) || tolerance <= 0.0 || maxIterations < 0)
    {
        throw std::invalid_argument(
            "conjugate_gradient_solver: tolerance must be finite and positive; max iterations must be nonnegative");
    }
    const std::size_t n = b.size();
    if (A.rows < 0 || A.rows != A.cols ||
        static_cast<std::size_t>(A.rows) != n || x.size() != n)
    {
        throw std::invalid_argument(
            "conjugate_gradient_solver: A must be square and match the sizes of b and x");
    }
    if (A.rowPtr.size() != n + 1 || A.values.size() != A.colIndex.size() ||
        A.rowPtr.front() != 0 || A.rowPtr.back() < 0 ||
        static_cast<std::size_t>(A.rowPtr.back()) != A.values.size())
    {
        throw std::invalid_argument("conjugate_gradient_solver: invalid CSR storage dimensions");
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        if (A.rowPtr[i] < 0 || A.rowPtr[i] > A.rowPtr[i + 1])
        {
            throw std::invalid_argument("conjugate_gradient_solver: invalid CSR row pointers");
        }
        if (!std::isfinite(b[i]) || !std::isfinite(x[i]))
        {
            throw std::invalid_argument("conjugate_gradient_solver: b and x must be finite");
        }
    }
    for (std::size_t k = 0; k < A.values.size(); ++k)
    {
        if (A.colIndex[k] < 0 || A.colIndex[k] >= A.cols ||
            !std::isfinite(A.values[k]))
        {
            throw std::invalid_argument("conjugate_gradient_solver: invalid CSR column or matrix value");
        }
    }
}

double checkedNorm(const std::vector<double>& v)
{
    const double result = vectorNorm(v);
    if (!std::isfinite(result))
    {
        throw std::runtime_error("conjugate_gradient_solver: nonfinite vector norm");
    }
    return result;
}

std::vector<double> applyChecked(
    const Preconditioner& preconditioner,
    const std::vector<double>& rhs)
{
    std::vector<double> result = preconditioner.applyInverse(rhs);
    if (result.size() != rhs.size())
    {
        throw std::invalid_argument("conjugate_gradient_solver: preconditioner returned the wrong vector size");
    }
    for (double value : result)
    {
        if (!std::isfinite(value))
        {
            throw std::runtime_error("conjugate_gradient_solver: preconditioner returned a nonfinite value");
        }
    }
    return result;
}

// rho = (r,z) has to stay > 0 for CG with an SPD preconditioner
double positiveRho(
    const std::vector<double>& r,
    const std::vector<double>& z)
{
    const double rz = dotProduct(r, z);
    if (!std::isfinite(rz) || rz <= 0.0)
    {
        throw std::runtime_error(
            "conjugate_gradient_solver: r.z must be positive; use an SPD preconditioner");
    }
    return rz;
}

} // namespace

double conjugate_gradient_solver::computeBeta(double rzNew, double rzOld) const
{
    if (!std::isfinite(rzNew) || !std::isfinite(rzOld) || rzNew <= 0.0 || rzOld <= 0.0)
    {
        throw std::runtime_error("conjugate_gradient_solver: invalid r.z in beta");
    }
    // beta = (r,z)_new / (r,z)_old
    const double beta = rzNew / rzOld;
    if (!std::isfinite(beta))
    {
        throw std::runtime_error("conjugate_gradient_solver: nonfinite beta");
    }
    return beta;
}

void conjugate_gradient_solver::updateDirection(
    const std::vector<double>& z,
    double beta,
    std::vector<double>& d) const
{
    for (std::size_t i = 0; i < d.size(); ++i)
    {
        d[i] = z[i] + beta * d[i]; // PCG: d_(k+1) = z_(k+1) + beta d_k.
    }
}

double conjugate_gradient_solver::computeAlpha(
    double rzOld,
    const std::vector<double>& d,
    const std::vector<double>& Ad) const
{
    const double dAd = dotProduct(d, Ad);
    if (!std::isfinite(dAd) || dAd <= 0.0 || !std::isfinite(rzOld) || rzOld <= 0.0)
    {
        throw std::runtime_error(
            "conjugate_gradient_solver: nonpositive or nonfinite curvature; A and M must be SPD");
    }
    // step length, alpha = (r,z)/(d,Ad). d.Ad <= 0 would mean A is not SPD
    const double alpha = rzOld / dAd;
    if (!std::isfinite(alpha) || alpha <= 0.0)
    {
        throw std::runtime_error("conjugate_gradient_solver: invalid step size");
    }
    return alpha;
}

void conjugate_gradient_solver::updatePhi(
    std::vector<double>& x,
    double alpha,
    const std::vector<double>& d) const
{
    axpy(x, alpha, d);
}

void conjugate_gradient_solver::updateResidual(
    std::vector<double>& r,
    double alpha,
    const std::vector<double>& Ad) const
{
    axpy(r, -alpha, Ad);
}

void conjugate_gradient_solver::solve(
    const CSRMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x)
{
    const IdentityPreconditioner identity;
    solve(A, b, x, identity);
}

void conjugate_gradient_solver::solve(
    const CSRMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    const Preconditioner& preconditioner)
{
    iters_ = 0;
    validateSolve(A, b, x, tol_, maxIter_);

    // relative tolerance on ||b||, if b = 0 we just use the absolute residual
    double bNorm = checkedNorm(b);
    if (bNorm == 0.0) bNorm = 1.0;
    std::vector<double> r = computeResidual(A, x, b);
    if (checkedNorm(r) / bNorm <= tol_ || maxIter_ == 0) return;

    std::vector<double> z;
    std::vector<double> d;
    std::vector<double> Ad;
    double rzOld = 0.0;

    // first search direction d = z = M^-1 r. Also used again if we restart
    // after replacing r with the true residual
    const auto resetDirection = [&]()
    {
        z = applyChecked(preconditioner, r); // M z = r.
        d = z;
        rzOld = positiveRho(r, z);
    };
    resetDirection();

    // main CG loop: x += alpha d, r -= alpha Ad, then new z and new direction
    while (iters_ < maxIter_)
    {
        Ad = multiply(A, d);
        const double alpha = computeAlpha(rzOld, d, Ad);
        updatePhi(x, alpha, d);
        updateResidual(r, alpha, Ad);
        ++iters_;

        // stop test uses the plain ||r||, not the preconditioned norm
        if (checkedNorm(r) / bNorm <= tol_)
        {
            r = computeResidual(A, x, b);
            if (checkedNorm(r) / bNorm <= tol_) return;

            // the updated r drifted away from b - A*x (round off), so restart from the true one
            if (iters_ < maxIter_) resetDirection();
            continue;
        }
        if (iters_ == maxIter_) break;

        z = applyChecked(preconditioner, r);
        const double rzNew = positiveRho(r, z);
        const double beta = computeBeta(rzNew, rzOld);
        updateDirection(z, beta, d);
        rzOld = rzNew;
    }
    // hitting maxIter does not mean it converged, the caller can check
    // computeResidual(A, x, b) itself (the comparison driver does that)
}
