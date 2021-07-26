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
#include <Reaktoro/Core/ThermoProps.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ThermoProps class", "[ThermoProps]")
{
    StandardThermoModel standard_thermo_model_gas = [](real T, real P)
    {
        StandardThermoProps props;
        props.G0  = 0.1 * (T*P)*(T*P);
        props.H0  = 0.2 * (T*P)*(T*P);
        props.V0  = 0.3 * (T*P)*(T*P);
        props.Cp0 = 0.4 * (T*P)*(T*P);
        props.Cv0 = 0.5 * (T*P)*(T*P);
        return props;
    };

    StandardThermoModel standard_thermo_model_solid = [](real T, real P)
    {
        StandardThermoProps props;
        props.G0  = 1.1 * (T*P)*(T*P);
        props.H0  = 1.2 * (T*P)*(T*P);
        props.V0  = 1.3 * (T*P)*(T*P);
        props.Cp0 = 1.4 * (T*P)*(T*P);
        props.Cv0 = 1.5 * (T*P)*(T*P);
        return props;
    };

    ActivityModel activity_props_fn_gas = [](ActivityPropsRef props, ActivityArgs args) {};

    ActivityModel activity_props_fn_solid = [](ActivityPropsRef props, ActivityArgs args) {};

    Database db;

    db.addSpecies( Species("H2O(g)").withStandardThermoModel(standard_thermo_model_gas) );
    db.addSpecies( Species("CO2(g)").withStandardThermoModel(standard_thermo_model_gas) );
    db.addSpecies( Species("CaCO3(s)").withStandardThermoModel(standard_thermo_model_solid) );

    Vec<Phase> phases
    {
        Phase()
            .withName("SomeGas")
            .withActivityPropsFn(activity_props_fn_gas)
            .withStateOfMatter(StateOfMatter::Gas)
            .withSpecies({
                db.species().get("H2O(g)"),
                db.species().get("CO2(g)")}),

        Phase()
            .withName("SomeSolid")
            .withActivityPropsFn(activity_props_fn_solid)
            .withStateOfMatter(StateOfMatter::Solid)
            .withSpecies({
                db.species().get("CaCO3(s)") })
    };

    ChemicalSystem system(db, phases);

    ThermoProps props(system);

    real T = 5.0;
    real P = 7.0;

    const ArrayXr G0  = ArrayXr{{ 0.1, 0.1, 1.1 }} * (T*P)*(T*P);
    const ArrayXr H0  = ArrayXr{{ 0.2, 0.2, 1.2 }} * (T*P)*(T*P);
    const ArrayXr V0  = ArrayXr{{ 0.3, 0.3, 1.3 }} * (T*P)*(T*P);
    const ArrayXr Cp0 = ArrayXr{{ 0.4, 0.4, 1.4 }} * (T*P)*(T*P);
    const ArrayXr Cv0 = ArrayXr{{ 0.5, 0.5, 1.5 }} * (T*P)*(T*P);
    const ArrayXr S0  = (H0 - G0)/T;
    const ArrayXr U0  = H0 - P*V0;
    const ArrayXr A0  = G0 - P*V0;

    CHECK_NOTHROW( props.update(T, P) );

    CHECK( props.temperature() == T );
    CHECK( props.pressure()    == P );

    CHECK( props.standardGibbsEnergies()       .isApprox(G0)   );
    CHECK( props.standardEnthalpies()          .isApprox(H0)   );
    CHECK( props.standardVolumes()             .isApprox(V0)   );
    CHECK( props.standardEntropies()           .isApprox(S0)   );
    CHECK( props.standardInternalEnergies()    .isApprox(U0)   );
    CHECK( props.standardHelmholtzEnergies()   .isApprox(A0)   );
    CHECK( props.standardHeatCapacitiesConstP().isApprox(Cp0)  );
    CHECK( props.standardHeatCapacitiesConstV().isApprox(Cv0)  );

    //---------------------------------------------------------------------
    // Testing temperature derivatives of the properties
    //---------------------------------------------------------------------
    const ArrayXd  G0_T = ArrayXd{{ 0.1, 0.1, 1.1 }} * 2*P*(T*P);
    const ArrayXd  H0_T = ArrayXd{{ 0.2, 0.2, 1.2 }} * 2*P*(T*P);
    const ArrayXd  V0_T = ArrayXd{{ 0.3, 0.3, 1.3 }} * 2*P*(T*P);
    const ArrayXd Cp0_T = ArrayXd{{ 0.4, 0.4, 1.4 }} * 2*P*(T*P);
    const ArrayXd Cv0_T = ArrayXd{{ 0.5, 0.5, 1.5 }} * 2*P*(T*P);
    const ArrayXd  S0_T = (H0_T - G0_T)/T - (H0 - G0)/(T*T);
    const ArrayXd  U0_T = H0_T - P*V0_T;
    const ArrayXd  A0_T = G0_T - P*V0_T;

    CHECK_NOTHROW( props.update(T, P, wrt(T)) );

    CHECK( grad(props.temperature()) == 1.0 );
    CHECK( grad(props.pressure())    == 0.0 );

    CHECK( grad(props.standardGibbsEnergies())        .isApprox(G0_T)   );
    CHECK( grad(props.standardEnthalpies())           .isApprox(H0_T)   );
    CHECK( grad(props.standardVolumes())              .isApprox(V0_T)   );
    CHECK( grad(props.standardEntropies())            .isApprox(S0_T)   );
    CHECK( grad(props.standardInternalEnergies())     .isApprox(U0_T)   );
    CHECK( grad(props.standardHelmholtzEnergies())    .isApprox(A0_T)   );
    CHECK( grad(props.standardHeatCapacitiesConstP()) .isApprox(Cp0_T)  );
    CHECK( grad(props.standardHeatCapacitiesConstV()) .isApprox(Cv0_T)  );

    //---------------------------------------------------------------------
    // Testing pressure derivatives of the properties
    //---------------------------------------------------------------------
    const ArrayXd  G0_P = ArrayXd{{ 0.1, 0.1, 1.1 }} * 2*T*(T*P);
    const ArrayXd  H0_P = ArrayXd{{ 0.2, 0.2, 1.2 }} * 2*T*(T*P);
    const ArrayXd  V0_P = ArrayXd{{ 0.3, 0.3, 1.3 }} * 2*T*(T*P);
    const ArrayXd Cp0_P = ArrayXd{{ 0.4, 0.4, 1.4 }} * 2*T*(T*P);
    const ArrayXd Cv0_P = ArrayXd{{ 0.5, 0.5, 1.5 }} * 2*T*(T*P);
    const ArrayXd  S0_P = (H0_P - G0_P)/T;
    const ArrayXd  U0_P = H0_P - V0 - P*V0_P;
    const ArrayXd  A0_P = G0_P - V0 - P*V0_P;

    CHECK_NOTHROW( props.update(T, P, wrt(P)) );

    CHECK( grad(props.temperature()) == 0.0 );
    CHECK( grad(props.pressure())    == 1.0 );

    CHECK( grad(props.standardGibbsEnergies())        .isApprox(G0_P)   );
    CHECK( grad(props.standardEnthalpies())           .isApprox(H0_P)   );
    CHECK( grad(props.standardVolumes())              .isApprox(V0_P)   );
    CHECK( grad(props.standardEntropies())            .isApprox(S0_P)   );
    CHECK( grad(props.standardInternalEnergies())     .isApprox(U0_P)   );
    CHECK( grad(props.standardHelmholtzEnergies())    .isApprox(A0_P)   );
    CHECK( grad(props.standardHeatCapacitiesConstP()) .isApprox(Cp0_P)  );
    CHECK( grad(props.standardHeatCapacitiesConstV()) .isApprox(Cv0_P)  );
}

