// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
        props.VT0 = 0.4 * (T*P)*(T*P);
        props.VP0 = 0.5 * (T*P)*(T*P);
        props.Cp0 = 0.6 * (T*P)*(T*P);
        return props;
    };

    StandardThermoModel standard_thermo_model_solid = [](real T, real P)
    {
        StandardThermoProps props;
        props.G0  = 1.1 * (T*P)*(T*P);
        props.H0  = 1.2 * (T*P)*(T*P);
        props.V0  = 1.3 * (T*P)*(T*P);
        props.VT0 = 1.4 * (T*P)*(T*P);
        props.VP0 = 1.5 * (T*P)*(T*P);
        props.Cp0 = 1.6 * (T*P)*(T*P);
        return props;
    };

    ActivityModel activity_model_gas = [](ActivityPropsRef props, ActivityModelArgs args) {};

    ActivityModel activity_model_solid = [](ActivityPropsRef props, ActivityModelArgs args) {};

    Database db;

    db.addSpecies( Species("H2O(g)").withStandardThermoModel(standard_thermo_model_gas) );
    db.addSpecies( Species("CO2(g)").withStandardThermoModel(standard_thermo_model_gas) );
    db.addSpecies( Species("CaCO3(s)").withStandardThermoModel(standard_thermo_model_solid) );

    Vec<Phase> phases
    {
        Phase()
            .withName("SomeGas")
            .withActivityModel(activity_model_gas)
            .withStateOfMatter(StateOfMatter::Gas)
            .withSpecies({
                db.species().get("H2O(g)"),
                db.species().get("CO2(g)")}),

        Phase()
            .withName("SomeSolid")
            .withActivityModel(activity_model_solid)
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
    const ArrayXr VT0 = ArrayXr{{ 0.4, 0.4, 1.4 }} * (T*P)*(T*P);
    const ArrayXr VP0 = ArrayXr{{ 0.5, 0.5, 1.5 }} * (T*P)*(T*P);
    const ArrayXr Cp0 = ArrayXr{{ 0.6, 0.6, 1.6 }} * (T*P)*(T*P);
    const ArrayXr Cv0 = Cp0 + T*VT0*VT0/VP0;
    const ArrayXr S0  = (H0 - G0)/T;
    const ArrayXr U0  = H0 - P*V0;
    const ArrayXr A0  = G0 - P*V0;

    CHECK_NOTHROW( props.update(T, P) );

    CHECK( props.temperature() == T );
    CHECK( props.pressure()    == P );

    CHECK( props.speciesStandardGibbsEnergies()       .isApprox(G0)   );
    CHECK( props.speciesStandardEnthalpies()          .isApprox(H0)   );
    CHECK( props.speciesStandardVolumes()             .isApprox(V0)   );
    CHECK( props.speciesStandardVolumesT()            .isApprox(VT0)  );
    CHECK( props.speciesStandardVolumesP()            .isApprox(VP0)  );
    CHECK( props.speciesStandardEntropies()           .isApprox(S0)   );
    CHECK( props.speciesStandardInternalEnergies()    .isApprox(U0)   );
    CHECK( props.speciesStandardHelmholtzEnergies()   .isApprox(A0)   );
    CHECK( props.speciesStandardHeatCapacitiesConstP().isApprox(Cp0)  );
    CHECK( props.speciesStandardHeatCapacitiesConstV().isApprox(Cv0)  );

    //---------------------------------------------------------------------
    // Testing temperature derivatives of the properties
    //---------------------------------------------------------------------
    const ArrayXd  G0_T = ArrayXd{{ 0.1, 0.1, 1.1 }} * 2*P*(T*P);
    const ArrayXd  H0_T = ArrayXd{{ 0.2, 0.2, 1.2 }} * 2*P*(T*P);
    const ArrayXd  V0_T = ArrayXd{{ 0.3, 0.3, 1.3 }} * 2*P*(T*P);
    const ArrayXd VT0_T = ArrayXd{{ 0.4, 0.4, 1.4 }} * 2*P*(T*P);
    const ArrayXd VP0_T = ArrayXd{{ 0.5, 0.5, 1.5 }} * 2*P*(T*P);
    const ArrayXd Cp0_T = ArrayXd{{ 0.6, 0.6, 1.6 }} * 2*P*(T*P);
    const ArrayXd Cv0_T = Cp0_T + VT0*VT0/VP0 + 2*T*VT0*VT0_T/VP0 - T*VT0*VT0/VP0/VP0*VP0_T;
    const ArrayXd  S0_T = (H0_T - G0_T)/T - (H0 - G0)/(T*T);
    const ArrayXd  U0_T = H0_T - P*V0_T;
    const ArrayXd  A0_T = G0_T - P*V0_T;

    CHECK_NOTHROW( props.update(T, P, wrt(T)) );

    CHECK( grad(props.temperature()) == 1.0 );
    CHECK( grad(props.pressure())    == 0.0 );

    CHECK( grad(props.speciesStandardGibbsEnergies())        .isApprox(G0_T)   );
    CHECK( grad(props.speciesStandardEnthalpies())           .isApprox(H0_T)   );
    CHECK( grad(props.speciesStandardVolumes())              .isApprox(V0_T)   );
    CHECK( grad(props.speciesStandardVolumesT())             .isApprox(VT0_T)  );
    CHECK( grad(props.speciesStandardVolumesP())             .isApprox(VP0_T)  );
    CHECK( grad(props.speciesStandardEntropies())            .isApprox(S0_T)   );
    CHECK( grad(props.speciesStandardInternalEnergies())     .isApprox(U0_T)   );
    CHECK( grad(props.speciesStandardHelmholtzEnergies())    .isApprox(A0_T)   );
    CHECK( grad(props.speciesStandardHeatCapacitiesConstP()) .isApprox(Cp0_T)  );
    CHECK( grad(props.speciesStandardHeatCapacitiesConstV()) .isApprox(Cv0_T)  );

    //---------------------------------------------------------------------
    // Testing pressure derivatives of the properties
    //---------------------------------------------------------------------
    const ArrayXd  G0_P = ArrayXd{{ 0.1, 0.1, 1.1 }} * 2*T*(T*P);
    const ArrayXd  H0_P = ArrayXd{{ 0.2, 0.2, 1.2 }} * 2*T*(T*P);
    const ArrayXd  V0_P = ArrayXd{{ 0.3, 0.3, 1.3 }} * 2*T*(T*P);
    const ArrayXd VT0_P = ArrayXd{{ 0.4, 0.4, 1.4 }} * 2*T*(T*P);
    const ArrayXd VP0_P = ArrayXd{{ 0.5, 0.5, 1.5 }} * 2*T*(T*P);
    const ArrayXd Cp0_P = ArrayXd{{ 0.6, 0.6, 1.6 }} * 2*T*(T*P);
    const ArrayXd Cv0_P = Cp0_P + 2*T*VT0*VT0_P/VP0 - T*VT0*VT0/VP0/VP0*VP0_P;
    const ArrayXd  S0_P = (H0_P - G0_P)/T;
    const ArrayXd  U0_P = H0_P - V0 - P*V0_P;
    const ArrayXd  A0_P = G0_P - V0 - P*V0_P;

    CHECK_NOTHROW( props.update(T, P, wrt(P)) );

    CHECK( grad(props.temperature()) == 0.0 );
    CHECK( grad(props.pressure())    == 1.0 );

    CHECK( grad(props.speciesStandardGibbsEnergies())        .isApprox(G0_P)   );
    CHECK( grad(props.speciesStandardEnthalpies())           .isApprox(H0_P)   );
    CHECK( grad(props.speciesStandardVolumes())              .isApprox(V0_P)   );
    CHECK( grad(props.speciesStandardVolumesT())             .isApprox(VT0_P)  );
    CHECK( grad(props.speciesStandardVolumesP())             .isApprox(VP0_P)  );
    CHECK( grad(props.speciesStandardEntropies())            .isApprox(S0_P)   );
    CHECK( grad(props.speciesStandardInternalEnergies())     .isApprox(U0_P)   );
    CHECK( grad(props.speciesStandardHelmholtzEnergies())    .isApprox(A0_P)   );
    CHECK( grad(props.speciesStandardHeatCapacitiesConstP()) .isApprox(Cp0_P)  );
    CHECK( grad(props.speciesStandardHeatCapacitiesConstV()) .isApprox(Cv0_P)  );
}

