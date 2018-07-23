// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>

namespace Reaktoro {

auto interpolate(
    const std::vector<double>& temperatures,
    const std::vector<double>& pressures,
    const std::vector<ThermoScalar>& scalars) -> ThermoScalarFunction
{
    std::vector<double> vals, ddTs, ddPs;
    vals.reserve(scalars.size());
    ddTs.reserve(scalars.size());
    ddPs.reserve(scalars.size());
    for(const ThermoScalar& scalar : scalars)
    {
        vals.push_back(scalar.val);
        ddTs.push_back(scalar.ddT);
        ddPs.push_back(scalar.ddP);
    }
    BilinearInterpolator val(temperatures, pressures, vals);
    BilinearInterpolator ddT(temperatures, pressures, ddTs);
    BilinearInterpolator ddP(temperatures, pressures, ddPs);

    auto func = [=](double T, double P)
    {
        return ThermoScalar(val(T, P), ddT(T, P), ddP(T, P));
    };

    return func;
}

auto interpolate(
    const std::vector<double>& temperatures,
    const std::vector<double>& pressures,
    const ThermoScalarFunction& f) -> ThermoScalarFunction
{
    auto val_func = [=](double T, double P) { return f(T, P).val; };
    auto ddT_func = [=](double T, double P) { return f(T, P).ddT; };
    auto ddP_func = [=](double T, double P) { return f(T, P).ddP; };

    BilinearInterpolator val(temperatures, pressures, val_func);
    BilinearInterpolator ddT(temperatures, pressures, ddT_func);
    BilinearInterpolator ddP(temperatures, pressures, ddP_func);

    auto func = [=](double T, double P)
    {
        return ThermoScalar(val(T, P), ddT(T, P), ddP(T, P));
    };

    return func;
}

auto interpolate(
    const std::vector<double>& temperatures,
    const std::vector<double>& pressures,
    const std::vector<ThermoScalarFunction>& fs) -> ThermoVectorFunction
{
    const unsigned size = fs.size();

    std::vector<BilinearInterpolator> val(size), ddT(size), ddP(size);

    for(unsigned i = 0; i < size; ++i)
    {
        auto val_func = [=](double T, double P) { return fs[i](T, P).val; };
        auto ddT_func = [=](double T, double P) { return fs[i](T, P).ddT; };
        auto ddP_func = [=](double T, double P) { return fs[i](T, P).ddP; };

        val[i] = BilinearInterpolator(temperatures, pressures, val_func);
        ddT[i] = BilinearInterpolator(temperatures, pressures, ddT_func);
        ddP[i] = BilinearInterpolator(temperatures, pressures, ddP_func);
    }

    ThermoVector res(size);

    auto func = [=](double T, double P) mutable
    {
        for(unsigned i = 0; i < size; ++i)
        {
            res.val[i] = val[i](T, P);
            res.ddT[i] = ddT[i](T, P);
            res.ddP[i] = ddP[i](T, P);
        }
        return res;
    };

    return func;
}

} // namespace Reaktoro
