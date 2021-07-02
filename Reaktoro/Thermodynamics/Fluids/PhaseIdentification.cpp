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

#include "PhaseIdentification.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Math/Roots.hpp>

// Eigen includes
#include <Eigen/Dense>
#include <unsupported/Eigen/Polynomials>

namespace Reaktoro {

using std::log;

auto identifyPhaseUsingVolume(
    const real& temperature,
    const real& pressure,
    const real& Z,
    const real& b) -> CubicEOSFluidType
{
    auto volume = Z * universalGasConstant * temperature / pressure;
    return volume / b > 1.75 ? CubicEOSFluidType::Vapor : CubicEOSFluidType::Liquid;
}

auto identifyPhaseUsingIsothermalCompressibility(
    const real& temperature, const real& pressure, const real& Z) -> CubicEOSFluidType
{
    error(true, "identifyPhaseUsingIsothermalCompressibility is currently not accessible.");
    // auto volume = Z * universalGasConstant * temperature / pressure;
    // auto dkdt = (1.0 / (volume * volume)) * volume.ddP * volume.ddT;

    // return (dkdt <= 0.0) ? CubicEOSFluidType::Vapor : CubicEOSFluidType::Liquid;
    return {};
}

auto pressureComparison(
    double Pressure,
    double Temperature,
    double amix,
    double bmix,
    double A,
    double B,
    double C,
    double epsilon,
    double sigma) -> CubicEOSFluidType
{
    // WARNING: The use of Eigen::PolynomialSolver with real seems to cause
    // errors in macOS:
    //
    // error: no matching function for call to
    // '__libcpp_isfinite_or_builtin' else if ((__libcpp_isinf_or_builtin(__a)
    // || __libcpp_isinf_or_builtin(__b)) && __libcpp_isfinite_or_builtin(__c)
    // && __libcpp_isfinite_or_builtin(__d))
    //
    // Using double instead.

    auto p = [&](double V) -> double
    {
        return ((universalGasConstant*Temperature) / (V - bmix)) - (amix / ((V + epsilon * bmix) * (V + sigma * bmix)));
    };

    auto k1 = epsilon * bmix;
    auto k2 = sigma * bmix;

    // Computing parameters AP, BP, CP, DP and EP of equation AP*P^4 + BP*P^3 + CP*P^2 + DP*P + EP = 0,
    // which gives the values of P where the EoS changes slope (Local max and min)
    const auto R = universalGasConstant;
    const auto T = Temperature;
    const auto AP = R * T;
    const auto BP = 2 * R * T * (k2 + k1) - 2 * amix;
    const auto CP = R * T * (k2 * k2 + 4.0 * k1 * k2 + k1 * k1) - amix * (k1 + k2 - 4 * bmix);
    const auto DP = 2 * R * T * (k1 * k2 * k2 + k1 * k1 * k2) - 2 * amix * (bmix * bmix - k2 * bmix - k1 * bmix);
    const auto EP = R * T * k1 * k1 * k2 * k2 - amix * (k1 + k2) * bmix * bmix;

    auto polynomial_solver = Eigen::PolynomialSolver<double, 4>(
        Eigen::Matrix<double, 5, 1>{EP, DP, CP, BP, AP});

    constexpr auto abs_imaginary_threshold = 1e-15;
    auto real_roots = std::vector<double>();
    real_roots.reserve(5);
    polynomial_solver.realRoots(real_roots, abs_imaginary_threshold);

    // removing roots lower than bmix, no physical meaning
    auto new_end = std::remove_if(real_roots.begin(), real_roots.end(), [&](double r){
        return r < bmix;
    });
    real_roots.resize(new_end - real_roots.begin());

    if (real_roots.size() == 0)
    {
        return CubicEOSFluidType::Vapor;
    }

    std::vector<double> pressures;
    for (const auto& volume : real_roots)
        pressures.push_back(p(volume));

    auto Pmin = std::min_element(pressures.begin(), pressures.end());
    auto Pmax = std::max_element(pressures.begin(), pressures.end());

    if (Pressure < *Pmin)
        return CubicEOSFluidType::Vapor;
    if (Pressure > *Pmax)
        return CubicEOSFluidType::Liquid;

    Exception exception;
    exception.error << "Could not define phase type.";
    exception.reason << "gibbsEnergyAndEquationOfStateMethod has received one Z but the pressure is between Pmin and Pmax.";
    RaiseError(exception);
}


auto gibbsResidualEnergyComparison(
    const real& pressure,
    const real& temperature,
    const real& amix,
    const real& bmix,
    const real& A,
    const real& B,
    const real& Z_min,
    const real& Z_max,
    const real epsilon,
    const real sigma) -> CubicEOSFluidType
{
    auto constexpr R = universalGasConstant;
    auto const& T = temperature;

    // Computing the values of residual Gibbs energy for all Zs
    std::vector<real> Gs;
    for (const auto Z : {Z_max, Z_min})
    {
        const real factor = -1.0 / (3 * Z*Z + 2 * A*Z + B);
        const real beta = pressure * bmix / (R * T);
        const real q = amix / (bmix * R * T);

        // Calculate the integration factor I and its temperature derivative IT
        real I = {};
        if (epsilon != sigma)
            I = log((Z + sigma * beta) / (Z + epsilon * beta)) / (sigma - epsilon);
        else
            I = beta / (Z + epsilon * beta);

        Gs.push_back(R * temperature*(Z - 1 - log(Z - beta) - q * I));
    }

    return (Gs[0] < Gs[1]) ? CubicEOSFluidType::Vapor : CubicEOSFluidType::Liquid;
}


auto identifyPhaseUsingGibbsEnergyAndEos(
    const real& pressure,
    const real& temperature,
    const real& amix,
    const real& bmix,
    const real& A,
    const real& B,
    const real& C,
    std::vector<real> Zs,
    const real epsilon,
    const real sigma) -> CubicEOSFluidType
{
    if (Zs.size() == 1)
    {
        return pressureComparison(pressure, temperature, amix, bmix, A, B, C, epsilon, sigma);
    }

    if (Zs.size() != 2) {
        Exception exception;
        exception.error << "identifyPhaseUsingGibbsEnergyAndEos received invalid input";
        exception.reason << "Zs should have size 1 or 2 in identifyPhaseUsingGibbsEnergyAndEos, "
            << "but has a size of " << Zs.size();
        RaiseError(exception);
    }

    const real& Z_min = Zs[0] < Zs[1] ? Zs[0] : Zs[1];
    const real& Z_max = Zs[0] < Zs[1] ? Zs[1] : Zs[0];

    return gibbsResidualEnergyComparison(pressure, temperature, amix, bmix, A, B, Z_min, Z_max, epsilon, sigma);
}

}
