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
#include <Reaktoro/Core/ChemicalProps.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ChemicalProps class", "[ChemicalProps]")
{
    const auto R = universalGasConstant;

    StandardThermoPropsFn standard_thermo_props_fn_gas = [](real T, real P)
    {
        StandardThermoProps props;
        props.G0  = 0.1 * (T * P);
        props.H0  = 0.2 * (T * P);
        props.V0  = 0.3 * (T * P);
        props.Cp0 = 0.4 * (T * P);
        props.Cv0 = 0.5 * (T * P);
        return props;
    };

    StandardThermoPropsFn standard_thermo_props_fn_solid = [](real T, real P)
    {
        StandardThermoProps props;
        props.G0  = 1.1 * (T * P);
        props.H0  = 1.2 * (T * P);
        props.V0  = 1.3 * (T * P);
        props.Cp0 = 1.4 * (T * P);
        props.Cv0 = 1.5 * (T * P);
        return props;
    };

    ActivityPropsFn activity_props_fn_gas = [](ActivityProps props, real T, real P, ArrayXrConstRef x)
    {
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

    ActivityPropsFn activity_props_fn_solid = [](ActivityProps props, real T, real P, ArrayXrConstRef x)
    {
        props.Vex  = 1.1 * (T * P);
        props.VexT = 2.1 * (T * P);
        props.VexP = 3.1 * (T * P);
        props.Gex  = 4.1 * (T * P);
        props.Hex  = 5.1 * (T * P);
        props.Cpex = 6.1 * (T * P);
        props.Cvex = 7.1 * (T * P);
        props.ln_g = 8.1 * x;
        props.ln_a = 9.1 * x;
    };

    Vec<Phase> phases
    {
        Phase()
            .withName("SomeGas")
            .withActivityPropsFn(activity_props_fn_gas)
            .withStateOfMatter(StateOfMatter::Gas)
            .withSpecies({
                Species("H2O(g)").withStandardThermoPropsFn(standard_thermo_props_fn_gas),
                Species("CO2(g)").withStandardThermoPropsFn(standard_thermo_props_fn_gas)}),

        Phase()
            .withName("SomeSolid")
            .withActivityPropsFn(activity_props_fn_solid)
            .withStateOfMatter(StateOfMatter::Solid)
            .withSpecies({
                Species("CaCO3(s)").withStandardThermoPropsFn(standard_thermo_props_fn_solid)})
    };

    ChemicalSystem system(phases);

    ChemicalProps props(system);

    SECTION("testing when species have non-zero amounts")
    {
        const real T = 300.0;
        const real P = 123.0e5;

        const ArrayXr n = ArrayXr{{ 4.0, 6.0, 5.0 }};
        const ArrayXr x = ArrayXr{{ 0.4, 0.6, 1.0 }};

        const ArrayXr G0  = ArrayXr{{ 0.1, 0.1, 1.1 }} * T * P;
        const ArrayXr H0  = ArrayXr{{ 0.2, 0.2, 1.2 }} * T * P;
        const ArrayXr V0  = ArrayXr{{ 0.3, 0.3, 1.3 }} * T * P;
        const ArrayXr Cp0 = ArrayXr{{ 0.4, 0.4, 1.4 }} * T * P;
        const ArrayXr Cv0 = ArrayXr{{ 0.5, 0.5, 1.5 }} * T * P;
        const ArrayXr S0  = (H0 - G0)/T;
        const ArrayXr U0  = H0 - P*V0;
        const ArrayXr A0  = G0 - P*V0;

        const ArrayXr Vex  = ArrayXr{{ 1.0, 1.1 }} * T * P;
        const ArrayXr VexT = ArrayXr{{ 2.0, 2.1 }} * T * P;
        const ArrayXr VexP = ArrayXr{{ 3.0, 3.1 }} * T * P;
        const ArrayXr Gex  = ArrayXr{{ 4.0, 4.1 }} * T * P;
        const ArrayXr Hex  = ArrayXr{{ 5.0, 5.1 }} * T * P;
        const ArrayXr Cpex = ArrayXr{{ 6.0, 6.1 }} * T * P;
        const ArrayXr Cvex = ArrayXr{{ 7.0, 7.1 }} * T * P;

        const ArrayXr ln_g = ArrayXr{{ 8.0*x[0], 8.0*x[1], 8.1*x[2] }};
        const ArrayXr ln_a = ArrayXr{{ 9.0*x[0], 9.0*x[1], 9.1*x[2] }};
        const ArrayXr u    = G0 + R*T*ln_a;

        REQUIRE_NOTHROW( props.update(T, P, n) );

        REQUIRE( props.temperature() == T );
        REQUIRE( props.pressure() == P );

        REQUIRE( props.standardGibbsEnergies()       .isApprox(G0)  );
        REQUIRE( props.standardEnthalpies()          .isApprox(H0)  );
        REQUIRE( props.standardVolumes()             .isApprox(V0)  );
        REQUIRE( props.standardEntropies()           .isApprox(S0)  );
        REQUIRE( props.standardInternalEnergies()    .isApprox(U0)  );
        REQUIRE( props.standardHelmholtzEnergies()   .isApprox(A0)  );
        REQUIRE( props.standardHeatCapacitiesConstP().isApprox(Cp0) );
        REQUIRE( props.standardHeatCapacitiesConstV().isApprox(Cv0) );

        REQUIRE( props.speciesAmounts()        .isApprox(n)    );
        REQUIRE( props.moleFractions()         .isApprox(x)    );
        REQUIRE( props.lnActivityCoefficients().isApprox(ln_g) );
        REQUIRE( props.lnActivities()          .isApprox(ln_a) );
        REQUIRE( props.chemicalPotentials()    .isApprox(u)    );

        real volume = 0.0;
        for(auto i = 0; i < phases.size(); ++i)
            volume += props.phaseProps(i).volume();

        REQUIRE( props.volume() == Approx(volume) );
    }

    SECTION("testing when species have zero amounts")
    {
        const real T = 300.0;
        const real P = 123.0e5;

        const ArrayXr n = ArrayXr{{ 0.0, 0.0, 0.0 }};
        const ArrayXr x = ArrayXr{{ 0.0, 0.0, 1.0 }}; // the solid phase has only one species, so xi = 1 for that species

        const ArrayXr G0  = ArrayXr{{ 0.1, 0.1, 1.1 }} * T * P;
        const ArrayXr H0  = ArrayXr{{ 0.2, 0.2, 1.2 }} * T * P;
        const ArrayXr V0  = ArrayXr{{ 0.3, 0.3, 1.3 }} * T * P;
        const ArrayXr Cp0 = ArrayXr{{ 0.4, 0.4, 1.4 }} * T * P;
        const ArrayXr Cv0 = ArrayXr{{ 0.5, 0.5, 1.5 }} * T * P;
        const ArrayXr S0  = (H0 - G0)/T;
        const ArrayXr U0  = H0 - P*V0;
        const ArrayXr A0  = G0 - P*V0;

        const ArrayXr Vex  = ArrayXr{{ 0.0, 0.0 }};
        const ArrayXr VexT = ArrayXr{{ 0.0, 0.0 }};
        const ArrayXr VexP = ArrayXr{{ 0.0, 0.0 }};
        const ArrayXr Gex  = ArrayXr{{ 0.0, 0.0 }};
        const ArrayXr Hex  = ArrayXr{{ 0.0, 0.0 }};
        const ArrayXr Cpex = ArrayXr{{ 0.0, 0.0 }};
        const ArrayXr Cvex = ArrayXr{{ 0.0, 0.0 }};

        const ArrayXr ln_g = ArrayXr{{ 0.0, 0.0, 0.0 }};
        const ArrayXr ln_a = ArrayXr{{ 0.0, 0.0, 0.0 }};
        const ArrayXr u    = G0 + R*T*ln_a;

        REQUIRE_NOTHROW( props.update(T, P, n) );

        REQUIRE( props.temperature() == T );
        REQUIRE( props.pressure() == P );

        REQUIRE( props.standardGibbsEnergies()       .isApprox(G0)  );
        REQUIRE( props.standardEnthalpies()          .isApprox(H0)  );
        REQUIRE( props.standardVolumes()             .isApprox(V0)  );
        REQUIRE( props.standardEntropies()           .isApprox(S0)  );
        REQUIRE( props.standardInternalEnergies()    .isApprox(U0)  );
        REQUIRE( props.standardHelmholtzEnergies()   .isApprox(A0)  );
        REQUIRE( props.standardHeatCapacitiesConstP().isApprox(Cp0) );
        REQUIRE( props.standardHeatCapacitiesConstV().isApprox(Cv0) );

        REQUIRE( props.speciesAmounts()        .isApprox(n)    );
        REQUIRE( props.moleFractions()         .isApprox(x)    );
        REQUIRE( props.lnActivityCoefficients().isApprox(ln_g) );
        REQUIRE( props.lnActivities()          .isApprox(ln_a) );
        REQUIRE( props.chemicalPotentials()    .isApprox(u)    );

        real volume = 0.0;
        for(auto i = 0; i < phases.size(); ++i)
            volume += props.phaseProps(i).volume();

        REQUIRE( props.volume() == Approx(volume) );
    }
}
