// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

namespace {

auto createStandardThermoModel(double param)
{
    StandardThermoModel model = [=](real T, real P)
    {
        StandardThermoProps props;
        props.G0  = 0.1 * param * T * P;
        props.H0  = 0.2 * param * T * P;
        props.V0  = 0.3 * param * T * P;
        props.VT0 = 0.4 * param * T * P;
        props.VP0 = 0.5 * param * T * P;
        props.Cp0 = 0.6 * param * T * P;
        return props;
    };
    return model;
}

} // namespace

TEST_CASE("Testing ThermoPropsPhase class", "[ThermoPropsPhase]")
{
    ActivityModel activity_model = [](ActivityPropsRef props, ActivityModelArgs args) {};

    Phase phase;
    phase = phase.withName("SomeGas");
    phase = phase.withActivityModel(activity_model);
    phase = phase.withStateOfMatter(StateOfMatter::Gas);
    phase = phase.withSpecies({
        Species("H2O(g)").withStandardThermoModel(createStandardThermoModel(10.0)), // param = 10.0
        Species("CO2(g)").withStandardThermoModel(createStandardThermoModel(20.0)), // param = 20.0
        Species("CH4(g)").withStandardThermoModel(createStandardThermoModel(30.0)), // param = 30.0
        Species("H2S(g)").withStandardThermoModel(createStandardThermoModel(40.0))  // param = 40.0
    });

    ThermoPropsPhase props(phase);

    real T = 5.0;
    real P = 7.0;

    const ArrayXr G0  = 0.1 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr H0  = 0.2 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr V0  = 0.3 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr VT0 = 0.4 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr VP0 = 0.5 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr Cp0 = 0.6 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
    const ArrayXr Cv0 = Cp0 + T*VT0*VT0/VP0;
    const ArrayXr S0  = (H0 - G0)/T;
    const ArrayXr U0  = H0 - P*V0;
    const ArrayXr A0  = G0 - P*V0;

    CHECK_NOTHROW( props.update(T, P) );

    CHECK( props.temperature() == T );
    CHECK( props.pressure()    == P );

    CHECK( props.speciesStandardGibbsEnergies()        .isApprox(G0)   );
    CHECK( props.speciesStandardEnthalpies()           .isApprox(H0)   );
    CHECK( props.speciesStandardVolumes()              .isApprox(V0)   );
    CHECK( props.speciesStandardEntropies()            .isApprox(S0)   );
    CHECK( props.speciesStandardInternalEnergies()     .isApprox(U0)   );
    CHECK( props.speciesStandardHelmholtzEnergies()    .isApprox(A0)   );
    CHECK( props.speciesStandardHeatCapacitiesConstP() .isApprox(Cp0)  );
    CHECK( props.speciesStandardHeatCapacitiesConstV() .isApprox(Cv0)  );

    //---------------------------------------------------------------------
    // Testing temperature derivatives of the properties
    //---------------------------------------------------------------------
    const ArrayXd  G0_T = 0.1 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
    const ArrayXd  H0_T = 0.2 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
    const ArrayXd  V0_T = 0.3 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
    const ArrayXd VT0_T = 0.4 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
    const ArrayXd VP0_T = 0.5 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
    const ArrayXd Cp0_T = 0.6 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
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
    const ArrayXd  G0_P = 0.1 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
    const ArrayXd  H0_P = 0.2 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
    const ArrayXd  V0_P = 0.3 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
    const ArrayXd VT0_P = 0.4 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
    const ArrayXd VP0_P = 0.5 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
    const ArrayXd Cp0_P = 0.6 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
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
