// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Models/StandardThermoModelNasa.hpp>
using namespace Reaktoro;

/// Manufacture coefficients a1, a2, ..., a7, b1, b2 in NasaThermoData object.
/// @param Tr The reference temperature at which Cp0(Tr) = Cp0x, H0(Tr) = H0x, and S0(Tr) = S0x
/// @param Tmin The minimum temperature in the interval
/// @param Tmax The maximum temperature in the interval
/// @param Cp0x The expected value of Cp0 at Tr
/// @param H0x The expected value of H0 at Tr
/// @param S0x The expected value of S0 at Tr
/// @return NasaThermoData
auto manufactureNasaPolynomial(real Tr, real Tmin, real Tmax, real Cp0x, real H0x, real S0x) -> StandardThermoModelParamsNasa::Polynomial
{
    const auto Tr2 = Tr*Tr;
    const auto Tr3 = Tr*Tr2;
    const auto Tr4 = Tr*Tr3;
    const auto lnTr = log(Tr);
    const auto R = universalGasConstant;
    const auto a1 = Tr2;
    const auto a2 = Tr;
    const auto a4 = 1/Tr;
    const auto a5 = 1/Tr2;
    const auto a6 = 1/Tr3;
    const auto a7 = 1/Tr4;
    const auto a3 = -(a1/Tr2 + a2/Tr + a4*Tr + a5*Tr2 + a6*Tr3 + a7*Tr4) + Cp0x/R;
    const auto b1 = -Tr*(-a1/Tr2 + a2/Tr*lnTr + a3 + a4*Tr/2 + a5*Tr2/3 + a6*Tr3/4 + a7*Tr4/5) + H0x/R;
    const auto b2 = -(-a1/Tr2/2 - a2/Tr + a3*lnTr + a4*Tr + a5*Tr2/2 + a6*Tr3/3 + a7*Tr4/4) + S0x/R;

    StandardThermoModelParamsNasa::Polynomial polynomial;
    polynomial.Tmin = Tmin;
    polynomial.Tmax = Tmax;
    polynomial.a1 = a1;
    polynomial.a2 = a2;
    polynomial.a3 = a3;
    polynomial.a4 = a4;
    polynomial.a5 = a5;
    polynomial.a6 = a6;
    polynomial.a7 = a7;
    polynomial.b1 = b1;
    polynomial.b2 = b2;
    return polynomial;
}

TEST_CASE("Testing StandardThermoModelNasa module", "[StandardThermoModelNasa]")
{
    StandardThermoModelParamsNasa params;

    const auto Cp0x = 345.6; // the expected Cp0 at reference temperatures for manufactured thermo params
    const auto H0x  = 123.4; // the expected  H0 at reference temperatures for manufactured thermo params
    const auto S0x  = 234.5; // the expected  S0 at reference temperatures for manufactured thermo params

    const auto Tr0 =  500.0; // the reference temperature at temperature interval 0
    const auto Tr1 = 2000.0; // the reference temperature at temperature interval 1
    const auto Tr2 = 8000.0; // the reference temperature at temperature interval 2

    params.polynomials.resize(3);
    params.polynomials[0] = manufactureNasaPolynomial(Tr0,  200.0, 1000.0, Cp0x, H0x, S0x);
    params.polynomials[1] = manufactureNasaPolynomial(Tr1, 1000.0, 6000.0, Cp0x, H0x, S0x);
    params.polynomials[2] = manufactureNasaPolynomial(Tr2, 6000.0, 9000.0, Cp0x, H0x, S0x);

    //======================================================================
    // Testing method detail::getNasaThermoDataForGivenTemperature
    //======================================================================

    CHECK( detail::indexTemperatureInterval(params.polynomials,  200.0) == 0 );
    CHECK( detail::indexTemperatureInterval(params.polynomials,  500.0) == 0 );
    CHECK( detail::indexTemperatureInterval(params.polynomials, 1000.0) == 0 );
    CHECK( detail::indexTemperatureInterval(params.polynomials, 2000.0) == 1 );
    CHECK( detail::indexTemperatureInterval(params.polynomials, 6000.0) == 1 );
    CHECK( detail::indexTemperatureInterval(params.polynomials, 8000.0) == 2 );
    CHECK( detail::indexTemperatureInterval(params.polynomials, 9000.0) == 2 );

    //======================================================================
    // Testing method detail::computeStandardThermoProps(polynomial, T)
    //======================================================================

    StandardThermoProps props0 = detail::computeStandardThermoProps(params.polynomials[0], Tr0);
    StandardThermoProps props1 = detail::computeStandardThermoProps(params.polynomials[1], Tr1);
    StandardThermoProps props2 = detail::computeStandardThermoProps(params.polynomials[2], Tr2);

    CHECK( props0.Cp0 == Approx(Cp0x)          );
    CHECK( props0.H0  == Approx(H0x)           );
    CHECK( props0.G0  == Approx(H0x - Tr0*S0x) );
    CHECK( props0.V0  == 0.0                   );
    CHECK( props0.Cv0 == 0.0                   );

    CHECK( props1.Cp0 == Approx(Cp0x)          );
    CHECK( props1.H0  == Approx(H0x)           );
    CHECK( props1.G0  == Approx(H0x - Tr1*S0x) );
    CHECK( props1.V0  == 0.0                   );
    CHECK( props1.Cv0 == 0.0                   );

    CHECK( props2.Cp0 == Approx(Cp0x)          );
    CHECK( props2.H0  == Approx(H0x)           );
    CHECK( props2.G0  == Approx(H0x - Tr2*S0x) );
    CHECK( props2.V0  == 0.0                   );
    CHECK( props2.Cv0 == 0.0                   );

    //======================================================================
    // Testing method detail::computeStandardThermoProps(params, T)
    //======================================================================

    CHECK( detail::computeStandardThermoProps(params, 650.0).G0 ==
           detail::computeStandardThermoProps(params.polynomials[0], 650.0).G0 );

    CHECK( detail::computeStandardThermoProps(params, 1650.0).G0 ==
           detail::computeStandardThermoProps(params.polynomials[1], 1650.0).G0 );

    CHECK( detail::computeStandardThermoProps(params, 8650.0).G0 ==
           detail::computeStandardThermoProps(params.polynomials[2], 8650.0).G0 );
}
