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
#include <Reaktor/Optimization/Hessian.hpp>

namespace Reaktor {

OptimumProblem::OptimumProblem()
: n(), m()
{}

auto OptimumProblem::setNumVariables(unsigned n) -> void
{
    this->n = n;
}

auto OptimumProblem::setNumConstraints(unsigned m) -> void
{
    this->m = m;
}

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

auto OptimumProblem::objectiveHessian(const Vector& x, const Vector& g) const -> Hessian
{
    if(not H)
        error("Cannot evaluate OptimumProblem::objectiveHessian.",
              "Have you called OptimumProblem::setObjectiveHessian before?");
    return H(x, g);
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

auto bfgs() -> ObjectiveHessianFunction
{
    Vector x0;
    Vector g0;
    Matrix B0;
    Vector dx;
    Vector dg;

    ObjectiveHessianFunction f = [=](const Vector& x, const Vector& g) mutable
    {
        Hessian H;
        H.mode = Hessian::Inverse;

        if(x0.size() == 0)
        {
            x0.noalias() = x;
            g0.noalias() = g;
            B0 = diag(x);
            H.inverse.noalias() = B0;
            return H;
        }

        dx.noalias() = x - x0;
        dg.noalias() = g - g0;
        x0.noalias() = x;
        g0.noalias() = g;

        const unsigned n = x.size();
        const double a = dx.dot(dg);
        const auto I = identity(n, n);

        H.inverse = (I - dx*tr(dg)/a)*B0*(I - dg*tr(dx)/a) + dx*tr(dx)/a;

        return H;
    };

    return f;
}

} // namespace Reaktor
