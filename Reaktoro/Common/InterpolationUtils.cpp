// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
    const std::vector<real>& temperatures,
    const std::vector<real>& pressures,
    const std::vector<real>& scalars) -> std::function<real(real, real)>
{
    return BilinearInterpolator(temperatures, pressures, scalars);
}

auto interpolate(
    const std::vector<real>& temperatures,
    const std::vector<real>& pressures,
    const std::function<real(real, real)>& func) -> std::function<real(real, real)>
{
    return BilinearInterpolator(temperatures, pressures, func);
}

auto interpolate(
    const std::vector<real>& temperatures,
    const std::vector<real>& pressures,
    const std::vector<std::function<real(real, real)>>& fs) -> std::function<VectorXr(real, real)>
{
    const auto size = fs.size();

    std::vector<BilinearInterpolator> interps(size);
    for(unsigned i = 0; i < size; ++i)
        interps[i] = BilinearInterpolator(temperatures, pressures, fs[i]);

    VectorXr res(size);

    auto func = [=](real T, real P) mutable
    {
        for(unsigned i = 0; i < size; ++i)
            res[i] = interps[i](T, P);
        return res;
    };

    return func;
}

} // namespace Reaktoro
