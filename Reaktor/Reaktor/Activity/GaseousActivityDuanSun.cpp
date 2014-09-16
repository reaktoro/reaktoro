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

#include "GaseousActivityDuanSun.hpp"

// C++ includes
#include <functional>
using namespace std::placeholders;

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/ConvertUtils.hpp>
#include <Reaktor/Common/Exception.hpp>

namespace Reaktor {
namespace internal {

const double coeffs1[]  = {1.0, -7.1734882e-1, -6.5129019e-2, 5.0383896, -16.063152, -1.5693490e-1};
const double coeffs2[]  = {4.7586835e-3, 1.5985379e-4, -2.1429977e-4, -4.4257744e-3, -2.7057990e-3, 4.4621407e-4};
const double coeffs3[]  = {-3.3569963e-6, -4.9286471e-7, -1.1444930e-6, 0.0, 0.0, -9.1080591e-7};
const double coeffs4[]  = {0.0, 0.0, 0.0, 1.9572733, 1.4119239e-1, 0.0};
const double coeffs5[]  = {-1.3179396, 0.0, 0.0, 0.0, 0.0, 0.0};
const double coeffs6[]  = {-3.8389101e-6, -2.7855285e-7, -1.1558081e-7, 2.4223436e-6, 8.1132965e-7, 1.0647399e-7};
const double coeffs7[]  = {0.0, 1.1877015e-9, 1.1952370e-9, 0.0, 0.0, 2.4273357e-10};
const double coeffs8[]  = {2.2815104e-3, 0.0, 0.0, -9.3796135e-4, -1.1453082e-4, 0.0};
const double coeffs9[]  = {0.0, 0.0, 0.0, -1.5026030, 2.3895671, 3.5874255e-1};
const double coeffs10[] = {0.0, 0.0, 0.0, 3.0272240e-3, 5.0527457e-4, 6.3319710e-5};
const double coeffs11[] = {0.0, 0.0, 0.0, -31.377342, -17.763460, -249.89661};
const double coeffs12[] = {0.0, -96.539512, -221.34306, -12.847063, 985.92232, 0.0};
const double coeffs13[] = {0.0, 4.4774938e-1, 0.0, 0.0, 0.0, 0.0};
const double coeffs14[] = {0.0, 101.81078, 71.820393, 0.0, 0.0, 888.768};
const double coeffs15[] = {0.0, 5.3783879e-6, 6.6089246e-6, -1.5056648e-5, -5.4965256e-7, -6.6348003e-7};

/// Calculate the saturated pressure of CO2 using Poling et al. (2001) model
auto saturatedPressureCO2(double T) -> double
{
    const double TcrCO2 = 304.20;
    const double PcrCO2 = 73.83;

    const double x   = 1 - T/TcrCO2;
    const double x15 = std::pow(x, 1.5);
    const double x3  = x15*x15;
    const double x6  = x3*x3;

    return PcrCO2 * std::exp((-6.95626*x + 1.19695*x15 -
        3.12614*x3 + 2.99448*x6)/(1 - x));
}

/// Determine the range index at which the provided temperature and pressure belongs to according to Table 1 of Duan et al. (2006)
auto regionIndex(double T, double Pbar) -> Index
{
    const double P1 = (T < 305) ? saturatedPressureCO2(T) :
        (305 <= T and T <= 405) ? 75 + (T - 305)*1.25 : 200;

    if(237 <= T and T <= 573 and 0 < Pbar and Pbar <= P1)
        return 0;
    if(237 <= T and T <= 340 and P1 <= Pbar and Pbar <= 1000)
        return 1;
    if(273 <= T and T <= 340 and Pbar > 1000)
        return 2;
    if(340 <= T and T <= 435 and P1 <= Pbar and Pbar <= 1000)
        return 3;
    if(340 <= T and T <= 435 and Pbar > 1000)
        return 4;
    if(435 <= T and T <= 533 and P1 <= Pbar and Pbar <= 2000)
        return 5;

    Exception exception;
    exception.error << "Cannot determine the temperature and pressure region for the calculation of "
        "the fugacity coefficient of CO2(g) using the Duan et al. (2006) model";
    exception.reason << "The temperature and/or pressure is out of range. " <<
        "(T = " << T << " K and P = " << Pbar*1e5 << " Pa).";

    raise(exception);

    return unsigned(-1);
}

auto gaseousActivityDuanSunCO2(const GaseousMixtureState& params, Index iCO2) -> ThermoScalar
{
    // The temperature (in units of K) and pressure (in units of bar)
    const double T  = params.T;
    const double Pb = convert<Pa,bar>(params.P);

    // Determine the region index for the use of the Duan et al. (2006) coefficients
    const Index idx = regionIndex(T, Pb);

    const double c1    = coeffs1[idx];
    const double c2    = coeffs2[idx];
    const double c3    = coeffs3[idx];
    const double c4    = coeffs4[idx];
    const double c5    = coeffs5[idx];
    const double c6    = coeffs6[idx];
    const double c7    = coeffs7[idx];
    const double c8    = coeffs8[idx];
    const double c9    = coeffs9[idx];
    const double c10   = coeffs10[idx];
    const double c11   = coeffs11[idx];
    const double c12   = coeffs12[idx];
    const double c13   = coeffs13[idx];
    const double c14   = coeffs14[idx];
    const double c15   = coeffs15[idx];

    // Calculate the fugacity coefficient of species CO2(g)
    const double phi = c1 + (c2 + c3*T + c4/T + c5/(T - 150))*Pb +
        (c6 + c7*T + c8/T)*Pb*Pb + (c9 + c10*T + c11/T)*std::log(Pb) +
        (c12 + c13*T)/Pb + c14/T + c15*T*T;

    // The molar fractions of the gaseous species and its molar derivatives
    const auto& x = params.x;

    // The molar fraction of CO2(g) and its molar derivatives
    ThermoScalar xCO2 = x.row(iCO2);

    // Calculate the activity of CO2(g) and its molar derivatives
    ThermoScalar aCO2;
    aCO2.val = phi * Pb * xCO2.val;
    aCO2.ddn = phi * Pb * xCO2.ddn;

    return aCO2;
}

} /* namespace internal */

auto gaseousActivityDuanSunCO2(const GaseousMixture& mixture) -> GaseousActivity
{
    const Index iCO2 = indexSpecies(mixture, "CO2(g)");

    return std::bind(internal::gaseousActivityDuanSunCO2, _1, iCO2);
}

} // namespace Reaktor
