// BiCGSTAB (van der Vorst) with right preconditioning. Each iteration has two
// half steps: a BiCG step with alpha, then a minimum residual step with omega.
#include "bicgstab.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace
{
// same input checks as in CG: sizes, CSR arrays, no inf/nan
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
            "bicgstab_solver: tolerance must be finite and positive; max iterations must be nonnegative");
    }
    const std::size_t n = b.size();
    if (A.rows < 0 || A.rows != A.cols ||
        static_cast<std::size_t>(A.rows) != n || x.size() != n)
    {
        throw std::invalid_argument(
            "bicgstab_solver: A must be square and match the sizes of b and x");
    }
    if (A.rowPtr.size() != n + 1 || A.values.size() != A.colIndex.size() ||
        A.rowPtr.front() != 0 || A.rowPtr.back() < 0 ||
        static_cast<std::size_t>(A.rowPtr.back()) != A.values.size())
    {
        throw std::invalid_argument("bicgstab_solver: invalid CSR storage dimensions");
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        if (A.rowPtr[i] < 0 || A.rowPtr[i] > A.rowPtr[i + 1])
        {
            throw std::invalid_argument("bicgstab_solver: invalid CSR row pointers");
        }
        if (!std::isfinite(b[i]) || !std::isfinite(x[i]))
        {
            throw std::invalid_argument("bicgstab_solver: b and x must be finite");
        }
    }
    for (std::size_t k = 0; k < A.values.size(); ++k)
    {
        if (A.colIndex[k] < 0 || A.colIndex[k] >= A.cols ||
            !std::isfinite(A.values[k]))
        {
            throw std::invalid_argument("bicgstab_solver: invalid CSR column or matrix value");
        }
    }
}

double checkedNorm(const std::vector<double>& v)
{
    const double result = vectorNorm(v);
    if (!std::isfinite(result))
    {
        throw std::runtime_error("bicgstab_solver: nonfinite vector norm");
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
        throw std::invalid_argument("bicgstab_solver: preconditioner returned the wrong vector size");
    }
    for (double value : result)
    {
        if (!std::isfinite(value))
        {
            throw std::runtime_error("bicgstab_solver: preconditioner returned a nonfinite value");
        }
    }
    return result;
}

// says if a dot product is basically zero compared to the two vectors.
// that is the BiCGSTAB breakdown case. Divide by the bigger norm first so
// |a|*|b| can't overflow
bool nearlyOrthogonal(
    double product,
    const std::vector<double>& a,
    const std::vector<double>& b)
{
    if (!std::isfinite(product) || product == 0.0) return true;
    const double aNorm = checkedNorm(a);
    const double bNorm = checkedNorm(b);
    if (aNorm == 0.0 || bNorm == 0.0) return true;
    const double large = std::max(aNorm, bNorm);
    const double small = std::min(aNorm, bNorm);
    return std::abs(product) / large <= std::numeric_limits<double>::epsilon() * small;
}

} // namespace

double bicgstab_solver::computeBeta(
    double rho,
    double rhoOld,
    double alpha,
    double omega) const
{
    if (!std::isfinite(rho) || !std::isfinite(rhoOld) ||
        !std::isfinite(alpha) || !std::isfinite(omega) ||
        rho == 0.0 || rhoOld == 0.0 || omega == 0.0)
    {
        throw std::runtime_error("bicgstab_solver: invalid scalar in beta");
    }
    const double beta = (rho / rhoOld) * (alpha / omega);
    if (!std::isfinite(beta))
    {
        throw std::runtime_error("bicgstab_solver: nonfinite beta");
    }
    return beta;
}

void bicgstab_solver::updateDirection(
    const std::vector<double>& r,
    double beta,
    double omega,
    const std::vector<double>& v,
    std::vector<double>& p) const
{
    for (std::size_t i = 0; i < p.size(); ++i)
    {
        p[i] = r[i] + beta * (p[i] - omega * v[i]);
    }
}

double bicgstab_solver::computeAlpha(
    double rho,
    const std::vector<double>& rHat,
    const std::vector<double>& v) const
{
    const double rHatV = dotProduct(rHat, v);
    if (nearlyOrthogonal(rHatV, rHat, v))
    {
        throw std::runtime_error(
            "bicgstab_solver: rHat.v is numerically zero; alpha breakdown");
    }
    const double alpha = rho / rHatV;
    if (!std::isfinite(alpha) || alpha == 0.0)
    {
        throw std::runtime_error("bicgstab_solver: invalid alpha");
    }
    return alpha;
}

void bicgstab_solver::computeIntermediateResidual(
    const std::vector<double>& r,
    double alpha,
    const std::vector<double>& v,
    std::vector<double>& s) const
{
    s = r;
    axpy(s, -alpha, v);
}

double bicgstab_solver::computeOmega(
    const std::vector<double>& t,
    const std::vector<double>& s) const
{
    const double tt = dotProduct(t, t);
    if (!std::isfinite(tt) || tt <= 0.0)
    {
        throw std::runtime_error("bicgstab_solver: invalid t.t; omega breakdown");
    }
    const double ts = dotProduct(t, s);
    if (nearlyOrthogonal(ts, t, s))
    {
        throw std::runtime_error(
            "bicgstab_solver: t.s is numerically zero; omega breakdown");
    }
    const double omega = ts / tt;
    if (!std::isfinite(omega) || omega == 0.0)
    {
        throw std::runtime_error("bicgstab_solver: invalid omega");
    }
    return omega;
}

void bicgstab_solver::updatePhi(
    std::vector<double>& x,
    double alpha,
    const std::vector<double>& pPre,
    double omega,
    const std::vector<double>& sPre) const
{
    axpy(x, alpha, pPre);
    axpy(x, omega, sPre);
}

void bicgstab_solver::updateResidual(
    std::vector<double>& r,
    const std::vector<double>& s,
    double omega,
    const std::vector<double>& t) const
{
    r = s;
    axpy(r, -omega, t);
}

void bicgstab_solver::solve(
    const CSRMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x)
{
    const IdentityPreconditioner identity;
    solve(A, b, x, identity);
}

void bicgstab_solver::solve(
    const CSRMatrix& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    const Preconditioner& preconditioner)
{
    iters_ = 0;
    validateSolve(A, b, x, tol_, maxIter_);

    double bNorm = checkedNorm(b);
    if (bNorm == 0.0) bNorm = 1.0;
    std::vector<double> r = computeResidual(A, x, b);
    if (checkedNorm(r) / bNorm <= tol_ || maxIter_ == 0) return;

    const std::size_t n = b.size();
    std::vector<double> rHat = r; // shadow residual, stays fixed unless we restart
    std::vector<double> p(n, 0.0), v(n, 0.0), s(n, 0.0), t;
    double rhoOld = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    bool firstStep = true;

    // start over from the current x (new rHat, p and v back to zero)
    const auto restart = [&]()
    {
        rHat = r;
        std::fill(p.begin(), p.end(), 0.0);
        std::fill(v.begin(), v.end(), 0.0);
        rhoOld = alpha = omega = 1.0;
        firstStep = true;
    };

    while (iters_ < maxIter_)
    {
        // rho = (rHat, r), if it goes to zero the method breaks down
        const double rho = dotProduct(rHat, r);
        if (nearlyOrthogonal(rho, rHat, r))
        {
            throw std::runtime_error(
                "bicgstab_solver: rHat.r is numerically zero; rho breakdown");
        }
        // new search direction, first time p = r, after that p = r + beta (p - omega v)
        if (firstStep)
        {
            p = r;
        }
        else
        {
            const double beta = computeBeta(rho, rhoOld, alpha, omega);
            updateDirection(r, beta, omega, v, p);
        }

        // first half step: v = A M^-1 p, alpha, and s = r - alpha v
        // (M^-1 p is a preconditioner solve, no inverse is formed)
        const std::vector<double> pPre = applyChecked(preconditioner, p);
        v = multiply(A, pPre);
        alpha = computeAlpha(rho, rHat, v);
        computeIntermediateResidual(r, alpha, v, s);

        if (checkedNorm(s) / bNorm <= tol_)
        {
            // s is already small enough, so take only the alpha part of the update and stop
            axpy(x, alpha, pPre);
            ++iters_;
            r = computeResidual(A, x, b);
            if (checkedNorm(r) / bNorm <= tol_) return;
            restart(); // the true residual was not small after all, so start again
            continue;
        }

        // second half step: t = A M^-1 s, omega minimises ||s - omega t||
        const std::vector<double> sPre = applyChecked(preconditioner, s);
        t = multiply(A, sPre);
        omega = computeOmega(t, s);
        updatePhi(x, alpha, pPre, omega, sPre);
        updateResidual(r, s, omega, t); // r = s - omega t, should still equal b - A*x
        ++iters_;

        // converged? check again with the true residual b - A*x before leaving
        if (checkedNorm(r) / bNorm <= tol_)
        {
            r = computeResidual(A, x, b);
            if (checkedNorm(r) / bNorm <= tol_) return;
            restart();
            continue;
        }
        rhoOld = rho;
        firstStep = false;
    }
    // ran out of iterations, that is not converged, the caller can check the true residual
}
