/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "Derivatives.hpp"

namespace Reaktor {
namespace internal {

// The square root of the machine epsilon (to be used in 1st-order finite difference schemes)
const double eps1storder = 1.0e-8;

// The cubic root of the machine epsilon (to be used in 2nd-order finite difference schemes)
const double eps2ndorder = 1.0e-6;

} /* namespace internal */

using namespace internal;

auto derivativeForward(const ScalarFunction& f, const Vector& x) -> Vector
{
    const double fx = f(x);

    const unsigned nrows = x.rows();

    Vector dfdx = zeros(nrows);
    Vector xh(nrows);
    for(unsigned i = 0; i < nrows; ++i)
    {
        const double h = eps1storder * std::max(std::abs(x[i]), 1.0);

        xh.noalias() = x + h * Vector::Unit(nrows, i);

        dfdx[i] = (f(xh) - fx)/h;
    }

    return dfdx;
}

auto derivativeBackward(const ScalarFunction& f, const Vector& x) -> Vector
{
    const double fx = f(x);

    const unsigned nrows = x.rows();

    Vector dfdx = zeros(nrows);
    Vector xh(nrows);
    for(unsigned i = 0; i < nrows; ++i)
    {
        const double h = eps1storder * std::max(std::abs(x[i]), 1.0);

        xh.noalias() = x - h * Vector::Unit(nrows, i);

        dfdx[i] = (fx - f(xh))/h;
    }

    return dfdx;
}

auto derivativeCentral(const ScalarFunction& f, const Vector& x) -> Vector
{
    const unsigned nrows = x.rows();

    Vector dfdx = zeros(nrows);
    Vector xh1(nrows), xh2(nrows);
    for(unsigned i = 0; i < nrows; ++i)
    {
        const double h = eps2ndorder * std::max(std::abs(x[i]), 1.0);

        xh1.noalias() = x + h * Vector::Unit(nrows, i);
        xh2.noalias() = x - h * Vector::Unit(nrows, i);

        dfdx[i] = (f(xh1) - f(xh2))/(2*h);
    }

    return dfdx;
}

auto derivativeForward(const VectorFunction& f, const Vector& x) -> Matrix
{
    Vector fx = f(x);

    const unsigned nrows = fx.rows();
    const unsigned ncols = x.rows();

    Matrix dfdx = zeros(nrows, ncols);
    Vector xh(ncols);
    for(unsigned i = 0; i < ncols; ++i)
    {
        const double h = eps1storder * std::max(std::abs(x[i]), 1.0);

        xh.noalias() = x + h * Vector::Unit(ncols, i);

        dfdx.col(i) = (f(xh) - fx)/h;
    }

    return dfdx;
}

auto derivativeBackward(const VectorFunction& f, const Vector& x) -> Matrix
{
    Vector fx = f(x);

    const unsigned nrows = fx.rows();
    const unsigned ncols = x.rows();

    Matrix dfdx = zeros(nrows, ncols);
    Vector xh(ncols);
    for(unsigned i = 0; i < ncols; ++i)
    {
        const double h = eps1storder * std::max(std::abs(x[i]), 1.0);

        xh.noalias() = x - h * Vector::Unit(ncols, i);

        dfdx.col(i) = (fx - f(xh))/h;
    }

    return dfdx;
}

auto derivativeCentral(const VectorFunction& f, const Vector& x) -> Matrix
{
    Vector fx = f(x);

    const unsigned nrows = fx.rows();
    const unsigned ncols = x.rows();

    Matrix dfdx = zeros(nrows, ncols);
    Vector xh1(ncols), xh2(ncols);
    for(unsigned i = 0; i < ncols; ++i)
    {
        const double h = eps2ndorder * std::max(std::abs(x[i]), 1.0);

        xh1 = x + h * Vector::Unit(ncols, i);
        xh2 = x - h * Vector::Unit(ncols, i);

        dfdx.col(i) = (f(xh1) - f(xh2))/(2*h);
    }

    return dfdx;
}

} /* namespace Reaktor */
