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
#include <Reaktoro/Core/PhaseChemicalProps.hpp>
using namespace Reaktoro;

auto createStandardThermoPropsFn(double param)
{
    StandardThermoPropsFn fn = [=](real T, real P)
    {
        StandardThermoProps props;
        props.G0  = 0.1 * param * T * P;
        props.H0  = 0.2 * param * T * P;
        props.V0  = 0.3 * param * T * P;
        props.Cp0 = 0.4 * param * T * P;
        props.Cv0 = 0.5 * param * T * P;
        return props;
    };
    return fn;
}

TEST_CASE("Testing PhaseChemicalProps class", "[PhaseChemicalProps]")
{
    const auto R = universalGasConstant;

    ActivityPropsFn activity_props_fn = [](ActivityPropsRef props, ActivityArgs args)
    {
        const auto [T, P, x, extra] = args;
        props.Vex  = 1.0 * (T * P);
        props.VexT = 2.0 * (T * P);
        props.VexP = 3.0 * (T * P);
        props.Gex  = 4.0 * (T * P);
        props.Hex  = 5.0 * (T * P);
        props.Cpex = 6.0 * (T * P);
        props.Cvex = 7.0 * (T * P);
        props.ln_g = 8.0 * x;
        props.ln_a = 9.0 * x;
    };

    Phase phase;
    phase = phase.withName("SomeGas");
    phase = phase.withActivityPropsFn(activity_props_fn);
    phase = phase.withStateOfMatter(StateOfMatter::Gas);
    phase = phase.withSpecies({
        Species("H2O(g)").withStandardThermoPropsFn(createStandardThermoPropsFn(10.0)), // param = 10.0
        Species("CO2(g)").withStandardThermoPropsFn(createStandardThermoPropsFn(20.0)), // param = 20.0
        Species("CH4(g)").withStandardThermoPropsFn(createStandardThermoPropsFn(30.0)), // param = 30.0
        Species("H2S(g)").withStandardThermoPropsFn(createStandardThermoPropsFn(40.0))  // param = 40.0
    });

    PhaseChemicalProps props(phase);

    SECTION("testing when species have non-zero amounts")
    {
        const real T = 300.0;
        const real P = 123.0e5;

        const ArrayXr n = ArrayXr{{ 1.0, 2.0, 3.0, 4.0 }};
        const ArrayXr x = ArrayXr{{ 0.1, 0.2, 0.3, 0.4 }};

        const ArrayXr G0  = 0.1 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr H0  = 0.2 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr V0  = 0.3 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr Cp0 = 0.4 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr Cv0 = 0.5 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr S0  = (H0 - G0)/T;
        const ArrayXr U0  = H0 - P*V0;
        const ArrayXr A0  = G0 - P*V0;

        const real Vex  = 1.0 * (T * P);
        const real VexT = 2.0 * (T * P);
        const real VexP = 3.0 * (T * P);
        const real Gex  = 4.0 * (T * P);
        const real Hex  = 5.0 * (T * P);
        const real Cpex = 6.0 * (T * P);
        const real Cvex = 7.0 * (T * P);

        const ArrayXr ln_g = 8.0 * x;
        const ArrayXr ln_a = 9.0 * x;
        const ArrayXr u    = G0 + R*T*ln_a;

        const real G  = (G0 * x).sum() + Gex;
        const real H  = (H0 * x).sum() + Hex;
        const real V  = (V0 * x).sum() + Vex;
        const real Cp = (Cp0 * x).sum() + Cpex;
        const real Cv = (Cv0 * x).sum() + Cvex;
        const real S  = (H - G)/T;
        const real U  = H - P*V;
        const real A  = G - P*V;

        const real amount = n.sum();

        real mass = 0.0;
        for(auto i = 0; i < phase.species().size(); ++i)
            mass += n[i] * phase.species(i).molarMass();

        REQUIRE_NOTHROW( props.update(T, P, n) );

        REQUIRE( props.temperature() == T );
        REQUIRE( props.pressure() == P );

        REQUIRE( props.speciesAmounts()        .isApprox(n)    );
        REQUIRE( props.moleFractions()         .isApprox(x)    );
        REQUIRE( props.lnActivityCoefficients().isApprox(ln_g) );
        REQUIRE( props.lnActivities()          .isApprox(ln_a) );
        REQUIRE( props.chemicalPotentials()    .isApprox(u)    );

        REQUIRE( (props.standardGibbsEnergies()        == G0  ).all() );
        REQUIRE( (props.standardEnthalpies()           == H0  ).all() );
        REQUIRE( (props.standardVolumes()              == V0  ).all() );
        REQUIRE( (props.standardEntropies()            == S0  ).all() );
        REQUIRE( (props.standardInternalEnergies()     == U0  ).all() );
        REQUIRE( (props.standardHelmholtzEnergies()    == A0  ).all() );
        REQUIRE( (props.standardHeatCapacitiesConstP() == Cp0 ).all() );
        REQUIRE( (props.standardHeatCapacitiesConstV() == Cv0 ).all() );

        REQUIRE( props.molarGibbsEnergy()        == Approx(G)          );
        REQUIRE( props.molarEnthalpy()           == Approx(H)          );
        REQUIRE( props.molarVolume()             == Approx(V)          );
        REQUIRE( props.molarEntropy()            == Approx(S)          );
        REQUIRE( props.molarInternalEnergy()     == Approx(U)          );
        REQUIRE( props.molarHelmholtzEnergy()    == Approx(A)          );
        REQUIRE( props.molarHeatCapacityConstP() == Approx(Cp)         );
        REQUIRE( props.molarHeatCapacityConstV() == Approx(Cv)         );
        REQUIRE( props.molarDensity()            == Approx(1.0 / V)    );
        REQUIRE( props.amount()                  == Approx(amount)     );
        REQUIRE( props.mass()                    == Approx(mass)       );
        REQUIRE( props.volume()                  == Approx(V * amount) );
    }

    SECTION("testing when species have zero amounts")
    {
        const real T = 300.0;
        const real P = 123.0e5;

        const ArrayXr n = ArrayXr{{ 0.0, 0.0, 0.0, 0.0 }};
        const ArrayXr x = ArrayXr{{ 0.0, 0.0, 0.0, 0.0 }};

        const ArrayXr G0  = 0.1 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr H0  = 0.2 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr V0  = 0.3 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr Cp0 = 0.4 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr Cv0 = 0.5 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr S0  = (H0 - G0)/T;
        const ArrayXr U0  = H0 - P*V0;
        const ArrayXr A0  = G0 - P*V0;

        const real Vex  = 0.0;
        const real VexT = 0.0;
        const real VexP = 0.0;
        const real Gex  = 0.0;
        const real Hex  = 0.0;
        const real Cpex = 0.0;
        const real Cvex = 0.0;

        const ArrayXr ln_g = ArrayXr{{ 0.0, 0.0, 0.0, 0.0 }};
        const ArrayXr ln_a = ArrayXr{{ 0.0, 0.0, 0.0, 0.0 }};
        const ArrayXr u    = G0;

        const real G  = 0.0;
        const real H  = 0.0;
        const real V  = 0.0;
        const real Cp = 0.0;
        const real Cv = 0.0;
        const real S  = 0.0;
        const real U  = 0.0;
        const real A  = 0.0;

        REQUIRE_NOTHROW( props.update(T, P, n) );

        REQUIRE( props.temperature() == T );
        REQUIRE( props.pressure()    == P );

        REQUIRE( props.speciesAmounts()        .isApprox(n)    );
        REQUIRE( props.moleFractions()         .isApprox(x)    );
        REQUIRE( props.lnActivityCoefficients().isApprox(ln_g) );
        REQUIRE( props.lnActivities()          .isApprox(ln_a) );
        REQUIRE( props.chemicalPotentials()    .isApprox(u)    );

        REQUIRE( (props.standardGibbsEnergies()        == G0  ).all() );
        REQUIRE( (props.standardEnthalpies()           == H0  ).all() );
        REQUIRE( (props.standardVolumes()              == V0  ).all() );
        REQUIRE( (props.standardEntropies()            == S0  ).all() );
        REQUIRE( (props.standardInternalEnergies()     == U0  ).all() );
        REQUIRE( (props.standardHelmholtzEnergies()    == A0  ).all() );
        REQUIRE( (props.standardHeatCapacitiesConstP() == Cp0 ).all() );
        REQUIRE( (props.standardHeatCapacitiesConstV() == Cv0 ).all() );

        REQUIRE( props.molarGibbsEnergy()        == Approx(0.0) );
        REQUIRE( props.molarEnthalpy()           == Approx(0.0) );
        REQUIRE( props.molarVolume()             == Approx(0.0) );
        REQUIRE( props.molarEntropy()            == Approx(0.0) );
        REQUIRE( props.molarInternalEnergy()     == Approx(0.0) );
        REQUIRE( props.molarHelmholtzEnergy()    == Approx(0.0) );
        REQUIRE( props.molarHeatCapacityConstP() == Approx(0.0) );
        REQUIRE( props.molarHeatCapacityConstV() == Approx(0.0) );
        REQUIRE( props.molarDensity()            == Approx(0.0) );
        REQUIRE( props.amount()                  == Approx(0.0) );
        REQUIRE( props.mass()                    == Approx(0.0) );
        REQUIRE( props.volume()                  == Approx(0.0) );
    }
}
