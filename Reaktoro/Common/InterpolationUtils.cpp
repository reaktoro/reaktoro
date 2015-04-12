// Reaktoro is a C++ library for computational reaction modelling.
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

#include "InterpolationUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Math/BilinearInterpolator.hpp>

namespace Reaktoro {

auto interpolate(
    const std::vector<double>& temperatures,
    const std::vector<double>& pressures,
    const std::vector<ThermoScalar>& scalars) -> ThermoScalarFunction
{
    std::vector<double> vals, ddts, ddps;
    vals.reserve(scalars.size());
    ddts.reserve(scalars.size());
    ddps.reserve(scalars.size());
    for(const ThermoScalar& scalar : scalars)
    {
        vals.push_back(scalar.val);
        ddts.push_back(scalar.ddt);
        ddps.push_back(scalar.ddp);
    }
    BilinearInterpolator val(temperatures, pressures, vals);
    BilinearInterpolator ddt(temperatures, pressures, ddts);
    BilinearInterpolator ddp(temperatures, pressures, ddps);

    auto func = [=](double T, double P)
    {
        return ThermoScalar(val(T, P), ddt(T, P), ddp(T, P));
    };

    return func;
}

auto interpolate(
    const std::vector<double>& temperatures,
    const std::vector<double>& pressures,
    const ThermoScalarFunction& f) -> ThermoScalarFunction
{
    auto val_func = [=](double T, double P) { return f(T, P).val; };
    auto ddt_func = [=](double T, double P) { return f(T, P).ddt; };
    auto ddp_func = [=](double T, double P) { return f(T, P).ddp; };

    BilinearInterpolator val(temperatures, pressures, val_func);
    BilinearInterpolator ddt(temperatures, pressures, ddt_func);
    BilinearInterpolator ddp(temperatures, pressures, ddp_func);

    auto func = [=](double T, double P)
    {
        return ThermoScalar(val(T, P), ddt(T, P), ddp(T, P));
    };

    return func;
}

} // namespace Reaktoro
