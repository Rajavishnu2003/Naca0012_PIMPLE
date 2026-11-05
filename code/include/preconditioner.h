#pragma once

#include <vector>

// Base class for preconditioners. The solver only needs z = M^-1 r,
// so that is the one function here. Building M is done separately.
class Preconditioner
{
public:
    virtual ~Preconditioner() = default;

    virtual std::vector<double> applyInverse(
        const std::vector<double>& rhs
    ) const = 0;
};


// M = I just returns r, so CG/SD/BiCGSTAB run without preconditioning
class IdentityPreconditioner final : public Preconditioner
{
public:
    std::vector<double> applyInverse(
        const std::vector<double>& rhs
    ) const override
    {
        return rhs;
    }
};
