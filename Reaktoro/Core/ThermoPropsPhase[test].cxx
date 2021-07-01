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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ThermoPropsPhase.hpp>
using namespace Reaktoro;

inline auto createStandardThermoModel(double param)
{
    StandardThermoModel model = [=](real T, real P)
    {
        StandardThermoProps props;
        props.G0  = 0.1 * param * T * P;
        props.H0  = 0.2 * param * T * P;
        props.V0  = 0.3 * param * T * P;
        props.Cp0 = 0.4 * param * T * P;
        props.Cv0 = 0.5 * param * T * P;
        return props;
    };
    return model;
}

TEST_CASE("Testing ThermoPropsPhase class", "[ThermoPropsPhase]")
{
    ActivityPropsFn activity_props_fn = [](ActivityPropsRef props, ActivityArgs args) {};

    Phase phase;
    phase = phase.withName("SomeGas");
    phase = phase.withActivityPropsFn(activity_props_fn);
    phase = phase.withStateOfMatter(StateOfMatter::Gas);
    phase = phase.withSpecies({
        Species("H2O(g)").withStandardThermoModel(createStandardThermoModel(10.0)), // param = 10.0
        Species("CO2(g)").withStandardThermoModel(createStandardThermoModel(20.0)), // param = 20.0
        Species("CH4(g)").withStandardThermoModel(createStandardThermoModel(30.0)), // param = 30.0
        Species("H2S(g)").withStandardThermoModel(createStandardThermoModel(40.0))  // param = 40.0
    });

    ThermoPropsPhase props(phase);

    real T = 300.0;
    real P = 123.0e5;

    const ArrayXr G0  = 0.1 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr H0  = 0.2 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr V0  = 0.3 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr Cp0 = 0.4 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr Cv0 = 0.5 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr S0  = (H0 - G0)/T;
    const ArrayXr U0  = H0 - P*V0;
    const ArrayXr A0  = G0 - P*V0;

    REQUIRE_NOTHROW( props.update(T, P) );

    REQUIRE( props.temperature() == T );
    REQUIRE( props.pressure()    == P );

    REQUIRE( props.standardGibbsEnergies()        .isApprox(G0)   );
    REQUIRE( props.standardEnthalpies()           .isApprox(H0)   );
    REQUIRE( props.standardVolumes()              .isApprox(V0)   );
    REQUIRE( props.standardEntropies()            .isApprox(S0)   );
    REQUIRE( props.standardInternalEnergies()     .isApprox(U0)   );
    REQUIRE( props.standardHelmholtzEnergies()    .isApprox(A0)   );
    REQUIRE( props.standardHeatCapacitiesConstP() .isApprox(Cp0)  );
    REQUIRE( props.standardHeatCapacitiesConstV() .isApprox(Cv0)  );

    //---------------------------------------------------------------------
    // Testing temperature derivatives of the properties
    //---------------------------------------------------------------------
    const ArrayXd  G0_T = 0.1 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
    const ArrayXd  H0_T = 0.2 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
    const ArrayXd  V0_T = 0.3 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
    const ArrayXd Cp0_T = 0.4 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
    const ArrayXd Cv0_T = 0.5 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
    const ArrayXd  S0_T = (H0_T - G0_T)/T - (H0 - G0)/(T*T);
    const ArrayXd  U0_T = H0_T - P*V0_T;
    const ArrayXd  A0_T = G0_T - P*V0_T;

    REQUIRE_NOTHROW( props.update(T, P, wrt(T)) );

    REQUIRE( grad(props.temperature()) == 1.0 );
    REQUIRE( grad(props.pressure())    == 0.0 );

    REQUIRE( grad(props.standardGibbsEnergies())        .isApprox(G0_T)   );
    REQUIRE( grad(props.standardEnthalpies())           .isApprox(H0_T)   );
    REQUIRE( grad(props.standardVolumes())              .isApprox(V0_T)   );
    REQUIRE( grad(props.standardEntropies())            .isApprox(S0_T)   );
    REQUIRE( grad(props.standardInternalEnergies())     .isApprox(U0_T)   );
    REQUIRE( grad(props.standardHelmholtzEnergies())    .isApprox(A0_T)   );
    REQUIRE( grad(props.standardHeatCapacitiesConstP()) .isApprox(Cp0_T)  );
    REQUIRE( grad(props.standardHeatCapacitiesConstV()) .isApprox(Cv0_T)  );

    //---------------------------------------------------------------------
    // Testing pressure derivatives of the properties
    //---------------------------------------------------------------------
    const ArrayXd  G0_P = 0.1 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
    const ArrayXd  H0_P = 0.2 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
    const ArrayXd  V0_P = 0.3 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
    const ArrayXd Cp0_P = 0.4 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
    const ArrayXd Cv0_P = 0.5 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
    const ArrayXd  S0_P = (H0_P - G0_P)/T;
    const ArrayXd  U0_P = H0_P - V0 - P*V0_P;
    const ArrayXd  A0_P = G0_P - V0 - P*V0_P;

    REQUIRE_NOTHROW( props.update(T, P, wrt(P)) );

    REQUIRE( grad(props.temperature()) == 0.0 );
    REQUIRE( grad(props.pressure())    == 1.0 );

    REQUIRE( grad(props.standardGibbsEnergies())        .isApprox(G0_P)   );
    REQUIRE( grad(props.standardEnthalpies())           .isApprox(H0_P)   );
    REQUIRE( grad(props.standardVolumes())              .isApprox(V0_P)   );
    REQUIRE( grad(props.standardEntropies())            .isApprox(S0_P)   );
    REQUIRE( grad(props.standardInternalEnergies())     .isApprox(U0_P)   );
    REQUIRE( grad(props.standardHelmholtzEnergies())    .isApprox(A0_P)   );
    REQUIRE( grad(props.standardHeatCapacitiesConstP()) .isApprox(Cp0_P)  );
    REQUIRE( grad(props.standardHeatCapacitiesConstV()) .isApprox(Cv0_P)  );
}
