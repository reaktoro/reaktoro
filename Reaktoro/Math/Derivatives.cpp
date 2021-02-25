// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "Derivatives.hpp"

namespace Reaktoro {
namespace {

// The square root of the machine epsilon (to be used in 1st-order finite difference schemes)
const double eps1storder = 1.0e-8;

// The cubic root of the machine epsilon (to be used in 2nd-order finite difference schemes)
const double eps2ndorder = 1.0e-6;

} // namespace

auto derivativeForward(const ScalarFunction& f, VectorXdConstRef x) -> VectorXd
{
    const double fx = f(x);
    const unsigned nrows = x.rows();
    VectorXd dfdx = zeros(nrows);
    VectorXd xh(nrows);
    for(unsigned i = 0; i < nrows; ++i)
    {
        const double h = eps1storder * (x[i] ? std::abs(x[i]) : 1.0);
        xh = x;
        xh[i] += h;
        dfdx[i] = (f(xh) - fx)/h;
    }
    return dfdx;
}

auto derivativeBackward(const ScalarFunction& f, VectorXdConstRef x) -> VectorXd
{
    const double fx = f(x);
    const unsigned nrows = x.rows();
    VectorXd dfdx = zeros(nrows);
    VectorXd xh(nrows);
    for(unsigned i = 0; i < nrows; ++i)
    {
        const double h = eps1storder * (x[i] ? std::abs(x[i]) : 1.0);
        xh = x;
        xh[i] -= h;
        dfdx[i] = (fx - f(xh))/h;
    }
    return dfdx;
}

auto derivativeCentral(const ScalarFunction& f, VectorXdConstRef x) -> VectorXd
{
    const unsigned nrows = x.rows();
    VectorXd dfdx = zeros(nrows);
    VectorXd xh1(nrows), xh2(nrows);
    for(unsigned i = 0; i < nrows; ++i)
    {
        const double h = eps2ndorder * (x[i] ? std::abs(x[i]) : 1.0);
        xh1 = x;
        xh2 = x;
        xh1[i] += h;
        xh2[i] -= h;
        dfdx[i] = (f(xh1) - f(xh2))/(2*h);
    }
    return dfdx;
}

auto derivativeForward(const VectorFunction& f, VectorXdConstRef x) -> MatrixXd
{
    const VectorXd fx = f(x);
    const unsigned nrows = fx.rows();
    const unsigned ncols = x.rows();
    MatrixXd dfdx = zeros(nrows, ncols);
    VectorXd xh(ncols);
    for(unsigned i = 0; i < ncols; ++i)
    {
        const double h = eps1storder * (x[i] ? std::abs(x[i]) : 1.0);
        xh = x;
        xh[i] += h;
        dfdx.col(i) = (f(xh) - fx)/h;
    }
    return dfdx;
}

auto derivativeBackward(const VectorFunction& f, VectorXdConstRef x) -> MatrixXd
{
    const VectorXd fx = f(x);
    const unsigned nrows = fx.rows();
    const unsigned ncols = x.rows();
    MatrixXd dfdx = zeros(nrows, ncols);
    VectorXd xh(ncols);
    for(unsigned i = 0; i < ncols; ++i)
    {
        const double h = eps1storder * (x[i] ? std::abs(x[i]) : 1.0);
        xh = x;
        xh[i] -= h;
        dfdx.col(i) = (fx - f(xh))/h;
    }
    return dfdx;
}

auto derivativeCentral(const VectorFunction& f, VectorXdConstRef x) -> MatrixXd
{
    const VectorXd fx = f(x);
    const unsigned nrows = fx.rows();
    const unsigned ncols = x.rows();
    MatrixXd dfdx = zeros(nrows, ncols);
    VectorXd xh1(ncols), xh2(ncols);
    for(unsigned i = 0; i < ncols; ++i)
    {
        const double h = eps2ndorder * (x[i] ? std::abs(x[i]) : 1.0);
        xh1 = x;
        xh2 = x;
        xh1[i] += h;
        xh2[i] -= h;
        dfdx.col(i) = (f(xh1) - f(xh2))/(2*h);
    }
    return dfdx;
}

} // namespace Reaktoro
