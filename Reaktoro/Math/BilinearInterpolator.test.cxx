// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Math/BilinearInterpolator.hpp>
using namespace Reaktoro;

TEST_CASE("Testing BilinearInterpolator", "[BilinearInterpolator]")
{
    auto nx = GENERATE(10, 5, 1);
    auto ny = GENERATE(10, 5, 1);

    INFO("nx = " << nx);
    INFO("ny = " << ny);

    Vec<double> x(nx);
    Vec<double> y(ny);
    Vec<double> z(nx * ny);

    for(auto i = 0; i < nx; ++i)
        x[i] = i;

    for(auto j = 0; j < ny; ++j)
        y[j] = j*j;

    for(auto i = 0; i < nx; ++i) for(auto j = 0; j < ny; ++j)
    {
        const auto k = i + nx*j;
        z[k] = 1 + 3*x[i] + 5*y[j] + 7*x[i]*y[j];
    }

    BilinearInterpolator f(x, y, z);

    for(auto i = 0; i < nx; ++i) for(auto j = 0; j < ny; ++j)
        CHECK( f(x[i], y[j]) == Approx(z[i + nx*j]) );

    for(auto i = 0; i < 2*nx; ++i) for(auto j = 0; j < 2*ny; ++j)
    {
        const auto xi = x.front() + i*(x.back() - x.front())/(2*nx);
        const auto yj = y.front() + j*(y.back() - y.front())/(2*ny);

        CHECK( f(xi, yj) == Approx(1 + 3*xi + 5*yj + 7*xi*yj) );
    }
}
