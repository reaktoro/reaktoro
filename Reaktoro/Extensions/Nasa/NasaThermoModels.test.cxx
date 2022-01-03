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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Extensions/Nasa/NasaThermoModels.hpp>
using namespace Reaktoro;

/// Manufacture coefficients a1, a2, ..., a7, b1, b2 in NasaThermoData object.
/// @param Tr The reference temperature at which Cp0(Tr) = Cp0x, H0(Tr) = H0x, and S0(Tr) = S0x
/// @param Tmin The minimum temperature in the interval
/// @param Tmax The maximum temperature in the interval
/// @param Cp0x The expected value of Cp0 at Tr
/// @param H0x The expected value of H0 at Tr
/// @param S0x The expected value of S0 at Tr
/// @return NasaThermoData
auto manufactureNasaThermoData(real Tr, real Tmin, real Tmax, real Cp0x, real H0x, real S0x) -> NasaThermoData
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

    NasaThermoData params;
    params.Tmin = Tmin;
    params.Tmax = Tmax;
    params.a1 = a1;
    params.a2 = a2;
    params.a3 = a3;
    params.a4 = a4;
    params.a5 = a5;
    params.a6 = a6;
    params.a7 = a7;
    params.b1 = b1;
    params.b2 = b2;
    return params;
}

TEST_CASE("Testing NasaThermoModels module", "[NasaThermoModels]")
{
    NasaSpecies species;
    species.name = "Abc(cr)";
    species.thermodata.resize(3);

    const auto Cp0x = 345.6; // the expected Cp0 at reference temperatures for manufactured thermo params
    const auto H0x  = 123.4; // the expected  H0 at reference temperatures for manufactured thermo params
    const auto S0x  = 234.5; // the expected  S0 at reference temperatures for manufactured thermo params

    const auto Tr0 =  500.0; // the reference temperature at temperature interval 0
    const auto Tr1 = 2000.0; // the reference temperature at temperature interval 1
    const auto Tr2 = 8000.0; // the reference temperature at temperature interval 2

    species.thermodata[0] = manufactureNasaThermoData(Tr0,  200.0, 1000.0, Cp0x, H0x, S0x);
    species.thermodata[1] = manufactureNasaThermoData(Tr1, 1000.0, 6000.0, Cp0x, H0x, S0x);
    species.thermodata[2] = manufactureNasaThermoData(Tr2, 6000.0, 9000.0, Cp0x, H0x, S0x);

    //======================================================================
    // Testing method NasaUtils::getNasaThermoDataForGivenTemperature
    //======================================================================

    CHECK( NasaUtils::getNasaThermoDataForGivenTemperature(species.thermodata,  200.0).Tmin == 200.0 );
    CHECK( NasaUtils::getNasaThermoDataForGivenTemperature(species.thermodata,  500.0).Tmin == 200.0 );
    CHECK( NasaUtils::getNasaThermoDataForGivenTemperature(species.thermodata, 1000.0).Tmin == 1000.0 );
    CHECK( NasaUtils::getNasaThermoDataForGivenTemperature(species.thermodata, 2000.0).Tmin == 1000.0 );
    CHECK( NasaUtils::getNasaThermoDataForGivenTemperature(species.thermodata, 6000.0).Tmin == 6000.0 );
    CHECK( NasaUtils::getNasaThermoDataForGivenTemperature(species.thermodata, 8000.0).Tmin == 6000.0 );
    CHECK( NasaUtils::getNasaThermoDataForGivenTemperature(species.thermodata, 9000.0).Tmin == 6000.0 );

    //======================================================================
    // Testing method NasaUtils::computeStandardThermoProps(params, T)
    //======================================================================

    StandardThermoProps props0 = NasaUtils::computeStandardThermoProps(species.thermodata[0], Tr0);
    StandardThermoProps props1 = NasaUtils::computeStandardThermoProps(species.thermodata[1], Tr1);
    StandardThermoProps props2 = NasaUtils::computeStandardThermoProps(species.thermodata[2], Tr2);

    CHECK( props0.Cp0 == Approx(Cp0x)          );
    CHECK( props0.H0  == Approx(H0x)           );
    CHECK( props0.G0  == Approx(H0x - Tr0*S0x) );
    CHECK( props0.V0  == 0.0                   );

    CHECK( props1.Cp0 == Approx(Cp0x)          );
    CHECK( props1.H0  == Approx(H0x)           );
    CHECK( props1.G0  == Approx(H0x - Tr1*S0x) );
    CHECK( props1.V0  == 0.0                   );

    CHECK( props2.Cp0 == Approx(Cp0x)          );
    CHECK( props2.H0  == Approx(H0x)           );
    CHECK( props2.G0  == Approx(H0x - Tr2*S0x) );
    CHECK( props2.V0  == 0.0                   );

    //======================================================================
    // Testing method NasaUtils::computeStandardThermoProps(species, T)
    //======================================================================

    CHECK( NasaUtils::computeStandardThermoProps(species, 650.0).G0 ==
           NasaUtils::computeStandardThermoProps(species.thermodata[0], 650.0).G0 );

    CHECK( NasaUtils::computeStandardThermoProps(species, 1650.0).G0 ==
           NasaUtils::computeStandardThermoProps(species.thermodata[1], 1650.0).G0 );

    CHECK( NasaUtils::computeStandardThermoProps(species, 8650.0).G0 ==
           NasaUtils::computeStandardThermoProps(species.thermodata[2], 8650.0).G0 );
}
