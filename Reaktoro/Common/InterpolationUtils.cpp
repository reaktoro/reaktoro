// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

#include "InterpolationUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Math/BilinearInterpolator.hpp>

namespace Reaktoro {

auto interpolate(
    const Vec<double>& temperatures,
    const Vec<double>& pressures,
    const Vec<double>& scalars) -> Fn<real(real, real)>
{
    return BilinearInterpolator(temperatures, pressures, scalars);
}

auto interpolate(
    const Vec<double>& temperatures,
    const Vec<double>& pressures,
    const Fn<double(double, double)>& func) -> Fn<real(real, real)>
{
    return BilinearInterpolator(temperatures, pressures, func);
}

auto interpolate(
    const Vec<double>& temperatures,
    const Vec<double>& pressures,
    const Vec<Fn<double(double, double)>>& fs) -> Fn<ArrayXr(real, real)>
{
    const auto size = fs.size();

    Vec<BilinearInterpolator> interps(size);
    for(unsigned i = 0; i < size; ++i)
        interps[i] = BilinearInterpolator(temperatures, pressures, fs[i]);

    ArrayXr res(size);

    auto func = [=](real T, real P) mutable
    {
        for(unsigned i = 0; i < size; ++i)
            res[i] = interps[i](T, P);
        return res;
    };

    return func;
}

} // namespace Reaktoro
