// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "OptimumProblem.hpp"

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>

namespace Reaktor {

OptimumProblem::OptimumProblem(unsigned n, unsigned m)
: n(n), m(m)
{}

auto OptimumProblem::setObjective(const ObjectiveFunction& f) -> void
{
    this->f = f;
}

auto OptimumProblem::setObjectiveGrad(const ObjectiveGradFunction& g) -> void
{
    this->g = g;
}

auto OptimumProblem::setObjectiveHessian(const ObjectiveHessianFunction& H) -> void
{
    this->H = H;
    hessian_scheme = HessianScheme::Regular;
}

auto OptimumProblem::setObjectiveDiagonalHessian(const ObjectiveDiagonalHessianFunction& diagH) -> void
{
    this->diagH = diagH;
    hessian_scheme = HessianScheme::Diagonal;
}

auto OptimumProblem::setObjectiveInverseHessian(const ObjectiveInverseHessianFunction& invH) -> void
{
    this->invH = invH;
    hessian_scheme = HessianScheme::Inverse;
}

auto OptimumProblem::setConstraint(const ConstraintFunction& h) -> void
{
    this->h = h;
}

auto OptimumProblem::setConstraintGrad(const ConstraintGradFunction& A) -> void
{
    this->A = A;
}

auto OptimumProblem::setLowerBounds(const Vector& lower) -> void
{
    l = lower;
}

auto OptimumProblem::setLowerBounds(double lower) -> void
{
    l = lower * ones(n);
}

auto OptimumProblem::setUpperBounds(const Vector& upper) -> void
{
    u = upper;
}

auto OptimumProblem::setUpperBounds(double upper) -> void
{
    u = upper * ones(n);
}

auto OptimumProblem::numVariables() const -> unsigned
{
    return n;
}

auto OptimumProblem::numConstraints() const -> unsigned
{
    return m;
}

auto OptimumProblem::lowerBounds() const -> const Vector&
{
    return l;
}

auto OptimumProblem::upperBounds() const -> const Vector&
{
    return u;
}

auto OptimumProblem::objective(const Vector& x) const -> double
{
    if(not f)
        error("Cannot evaluate OptimumProblem::objective.",
              "Have you called OptimumProblem::setObjective before?");
    return f(x);
}

auto OptimumProblem::objectiveGrad(const Vector& x) const -> Vector
{
    if(not g)
        error("Cannot evaluate OptimumProblem::objectiveGrad.",
              "Have you called OptimumProblem::setObjectiveGrad before?");
    return g(x);
}

auto OptimumProblem::objectiveHessian(const Vector& x) const -> Matrix
{
    if(not H)
        error("Cannot evaluate OptimumProblem::objectiveHessian.",
              "Have you called OptimumProblem::setObjectiveHessian before?");
    return H(x);
}

auto OptimumProblem::objectiveDiagonalHessian(const Vector& x) const -> Vector
{
    if(not diagH)
        error("Cannot evaluate OptimumProblem::objectiveDiagonalHessian.",
              "Have you called OptimumProblem::setObjectiveDiagonalHessian before?");
    return diagH(x);
}

auto OptimumProblem::objectiveInverseHessian(const Vector& x, const Vector& g) const -> Matrix
{
    if(not invH)
        error("Cannot evaluate OptimumProblem::objectiveInverseHessian.",
              "Have you called OptimumProblem::setObjectiveInverseHessian before?");
    return invH(x, g);
}

auto OptimumProblem::constraint(const Vector& x) const -> Vector
{
    if(not h)
        error("Cannot evaluate OptimumProblem::constraint.",
              "Have you called OptimumProblem::setConstraint before?");
    return h(x);
}

auto OptimumProblem::constraintGrad(const Vector& x) const -> Matrix
{
    if(not A)
        error("Cannot evaluate OptimumProblem::constraintGrad.",
              "Have you called OptimumProblem::setConstraintGrad before?");
    return A(x);
}

auto OptimumProblem::hessianScheme() const -> HessianScheme
{
    return hessian_scheme;
}

auto bfgs() -> ObjectiveInverseHessianFunction
{
    Vector x0;
    Vector g0;
    Matrix B0;
    Vector dx;
    Vector dg;

    ObjectiveInverseHessianFunction f = [=](const Vector& x, const Vector& g) mutable
    {
        if(x0.size() == 0)
        {
            x0.noalias() = x;
            g0.noalias() = g;
            B0 = diag(x);
            return B0;
        }

        dx.noalias() = x - x0;
        dg.noalias() = g - g0;
        x0.noalias() = x;
        g0.noalias() = g;

        const unsigned n = x.size();
        const double a = dx.dot(dg);
        const auto I = identity(n, n);
        Matrix B = (I - dx*tr(dg)/a)*B0*(I - dg*tr(dx)/a) + dx*tr(dx)/a;

        return B;
    };

    return f;
}

} // namespace Reaktor
