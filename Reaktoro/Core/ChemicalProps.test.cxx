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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
using namespace Reaktoro;

namespace test { extern auto createDatabase() -> Database; }

TEST_CASE("Testing ChemicalProps class", "[ChemicalProps]")
{
    const auto R = universalGasConstant;

    StandardThermoPropsFn standard_thermo_props_fn_gas = [](real T, real P)
    {
        StandardThermoProps props;
        props.G0  = 0.1 * (T*P)*(T*P);
        props.H0  = 0.2 * (T*P)*(T*P);
        props.V0  = 0.3 * (T*P)*(T*P);
        props.Cp0 = 0.4 * (T*P)*(T*P);
        props.Cv0 = 0.5 * (T*P)*(T*P);
        return props;
    };

    StandardThermoPropsFn standard_thermo_props_fn_solid = [](real T, real P)
    {
        StandardThermoProps props;
        props.G0  = 1.1 * (T*P)*(T*P);
        props.H0  = 1.2 * (T*P)*(T*P);
        props.V0  = 1.3 * (T*P)*(T*P);
        props.Cp0 = 1.4 * (T*P)*(T*P);
        props.Cv0 = 1.5 * (T*P)*(T*P);
        return props;
    };

    ActivityPropsFn activity_props_fn_gas = [](ActivityPropsRef props, ActivityArgs args)
    {
        const auto [T, P, x, extra] = args;
        props.Vex  = 1.0 * (T*P)*(T*P);
        props.VexT = 2.0 * (T*P)*(T*P);
        props.VexP = 3.0 * (T*P)*(T*P);
        props.Gex  = 4.0 * (T*P)*(T*P);
        props.Hex  = 5.0 * (T*P)*(T*P);
        props.Cpex = 6.0 * (T*P)*(T*P);
        props.Cvex = 7.0 * (T*P)*(T*P);
        props.ln_g = 8.0 * x;
        props.ln_a = 9.0 * x;
    };

    ActivityPropsFn activity_props_fn_solid = [](ActivityPropsRef props, ActivityArgs args)
    {
        const auto [T, P, x, extra] = args;
        props.Vex  = 1.1 * (T*P)*(T*P);
        props.VexT = 2.1 * (T*P)*(T*P);
        props.VexP = 3.1 * (T*P)*(T*P);
        props.Gex  = 4.1 * (T*P)*(T*P);
        props.Hex  = 5.1 * (T*P)*(T*P);
        props.Cpex = 6.1 * (T*P)*(T*P);
        props.Cvex = 7.1 * (T*P)*(T*P);
        props.ln_g = 8.1 * x;
        props.ln_a = 9.1 * x;
    };

    Database db;

    db.addSpecies( Species("H2O(g)").withStandardThermoPropsFn(standard_thermo_props_fn_gas) );
    db.addSpecies( Species("CO2(g)").withStandardThermoPropsFn(standard_thermo_props_fn_gas) );
    db.addSpecies( Species("CaCO3(s)").withStandardThermoPropsFn(standard_thermo_props_fn_solid) );

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

    ChemicalProps props(system);

    SECTION("testing when species have non-zero amounts")
    {
        real T = 300.0;
        real P = 123.0e5;
        ArrayXr n = ArrayXr{{ 4.0, 6.0, 5.0 }};

        const ArrayXr x = ArrayXr{{ 0.4, 0.6, 1.0 }};

        const ArrayXr nsumphases = ArrayXr{{ 10.0, 5.0 }}; // the total amount of each phase

        const ArrayXr G0  = ArrayXr{{ 0.1, 0.1, 1.1 }} * (T*P)*(T*P);
        const ArrayXr H0  = ArrayXr{{ 0.2, 0.2, 1.2 }} * (T*P)*(T*P);
        const ArrayXr V0  = ArrayXr{{ 0.3, 0.3, 1.3 }} * (T*P)*(T*P);
        const ArrayXr Cp0 = ArrayXr{{ 0.4, 0.4, 1.4 }} * (T*P)*(T*P);
        const ArrayXr Cv0 = ArrayXr{{ 0.5, 0.5, 1.5 }} * (T*P)*(T*P);
        const ArrayXr S0  = (H0 - G0)/T;
        const ArrayXr U0  = H0 - P*V0;
        const ArrayXr A0  = G0 - P*V0;

        const ArrayXr Vex  = ArrayXr{{ 1.0, 1.1 }} * (T*P)*(T*P);
        const ArrayXr VexT = ArrayXr{{ 2.0, 2.1 }} * (T*P)*(T*P);
        const ArrayXr VexP = ArrayXr{{ 3.0, 3.1 }} * (T*P)*(T*P);
        const ArrayXr Gex  = ArrayXr{{ 4.0, 4.1 }} * (T*P)*(T*P);
        const ArrayXr Hex  = ArrayXr{{ 5.0, 5.1 }} * (T*P)*(T*P);
        const ArrayXr Cpex = ArrayXr{{ 6.0, 6.1 }} * (T*P)*(T*P);
        const ArrayXr Cvex = ArrayXr{{ 7.0, 7.1 }} * (T*P)*(T*P);

        const ArrayXr ln_g = ArrayXr{{ 8.0*x[0], 8.0*x[1], 8.1*x[2] }};
        const ArrayXr ln_a = ArrayXr{{ 9.0*x[0], 9.0*x[1], 9.1*x[2] }};
        const ArrayXr u    = G0 + R*T*ln_a;

        const real Gtot = (G0 * n).sum() + (nsumphases * Gex).sum();
        const real Htot = (H0 * n).sum() + (nsumphases * Hex).sum();
        const real Vtot = (V0 * n).sum() + (nsumphases * Vex).sum();
        const real Stot = (Htot - Gtot)/T;
        const real Utot = Htot - P*Vtot;
        const real Atot = Gtot - P*Vtot;

        REQUIRE_NOTHROW( props.update(T, P, n) );

        REQUIRE( props.temperature() == T );
        REQUIRE( props.pressure()    == P );

        REQUIRE( props.speciesAmounts()              .isApprox(n)    );
        REQUIRE( props.moleFractions()               .isApprox(x)    );
        REQUIRE( props.lnActivityCoefficients()      .isApprox(ln_g) );
        REQUIRE( props.lnActivities()                .isApprox(ln_a) );
        REQUIRE( props.chemicalPotentials()          .isApprox(u)    );
        REQUIRE( props.standardGibbsEnergies()       .isApprox(G0)   );
        REQUIRE( props.standardEnthalpies()          .isApprox(H0)   );
        REQUIRE( props.standardVolumes()             .isApprox(V0)   );
        REQUIRE( props.standardEntropies()           .isApprox(S0)   );
        REQUIRE( props.standardInternalEnergies()    .isApprox(U0)   );
        REQUIRE( props.standardHelmholtzEnergies()   .isApprox(A0)   );
        REQUIRE( props.standardHeatCapacitiesConstP().isApprox(Cp0)  );
        REQUIRE( props.standardHeatCapacitiesConstV().isApprox(Cv0)  );

        REQUIRE( props.gibbsEnergy()     == Approx(Gtot) );
        REQUIRE( props.enthalpy()        == Approx(Htot) );
        REQUIRE( props.volume()          == Approx(Vtot) );
        REQUIRE( props.entropy()         == Approx(Stot) );
        REQUIRE( props.internalEnergy()  == Approx(Utot) );
        REQUIRE( props.helmholtzEnergy() == Approx(Atot) );

        //---------------------------------------------------------------------
        // Testing temperature derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXd n_T = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd x_T = ArrayXd{{ 0.0, 0.0, 0.0 }};

        const ArrayXd nsumphases_T = ArrayXd{{ 0.0, 0.0 }};

        const ArrayXd  G0_T = ArrayXd{{ 0.1, 0.1, 1.1 }} * 2*P*(T*P);
        const ArrayXd  H0_T = ArrayXd{{ 0.2, 0.2, 1.2 }} * 2*P*(T*P);
        const ArrayXd  V0_T = ArrayXd{{ 0.3, 0.3, 1.3 }} * 2*P*(T*P);
        const ArrayXd Cp0_T = ArrayXd{{ 0.4, 0.4, 1.4 }} * 2*P*(T*P);
        const ArrayXd Cv0_T = ArrayXd{{ 0.5, 0.5, 1.5 }} * 2*P*(T*P);
        const ArrayXd  S0_T = (H0_T - G0_T)/T - (H0 - G0)/(T*T);
        const ArrayXd  U0_T = H0_T - P*V0_T;
        const ArrayXd  A0_T = G0_T - P*V0_T;

        const ArrayXd  Vex_T = ArrayXd{{ 1.0, 1.1 }} * 2*P*(T*P);
        const ArrayXd VexT_T = ArrayXd{{ 2.0, 2.1 }} * 2*P*(T*P);
        const ArrayXd VexP_T = ArrayXd{{ 3.0, 3.1 }} * 2*P*(T*P);
        const ArrayXd  Gex_T = ArrayXd{{ 4.0, 4.1 }} * 2*P*(T*P);
        const ArrayXd  Hex_T = ArrayXd{{ 5.0, 5.1 }} * 2*P*(T*P);
        const ArrayXd Cpex_T = ArrayXd{{ 6.0, 6.1 }} * 2*P*(T*P);
        const ArrayXd Cvex_T = ArrayXd{{ 7.0, 7.1 }} * 2*P*(T*P);

        const ArrayXd ln_g_T = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd ln_a_T = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd    u_T = G0_T + R*ln_a;

        const double Gtot_T = (G0_T * n).sum() + (nsumphases * Gex_T).sum();
        const double Htot_T = (H0_T * n).sum() + (nsumphases * Hex_T).sum();
        const double Vtot_T = (V0_T * n).sum() + (nsumphases * Vex_T).sum();
        const double Stot_T = (Htot_T - Gtot_T)/T - (Htot - Gtot)/(T*T);
        const double Utot_T = Htot_T - P*Vtot_T;
        const double Atot_T = Gtot_T - P*Vtot_T;

        REQUIRE_NOTHROW( props.update(T, P, n, wrt(T)) );

        REQUIRE( grad(props.temperature()) == 1.0 );
        REQUIRE( grad(props.pressure())    == 0.0 );

        REQUIRE( grad(props.speciesAmounts())               .isApprox(n_T)    );
        REQUIRE( grad(props.moleFractions())                .isApprox(x_T)    );
        REQUIRE( grad(props.lnActivityCoefficients())       .isApprox(ln_g_T) );
        REQUIRE( grad(props.lnActivities())                 .isApprox(ln_a_T) );
        REQUIRE( grad(props.chemicalPotentials())           .isApprox(u_T)    );
        REQUIRE( grad(props.standardGibbsEnergies())        .isApprox(G0_T)   );
        REQUIRE( grad(props.standardEnthalpies())           .isApprox(H0_T)   );
        REQUIRE( grad(props.standardVolumes())              .isApprox(V0_T)   );
        REQUIRE( grad(props.standardEntropies())            .isApprox(S0_T)   );
        REQUIRE( grad(props.standardInternalEnergies())     .isApprox(U0_T)   );
        REQUIRE( grad(props.standardHelmholtzEnergies())    .isApprox(A0_T)   );
        REQUIRE( grad(props.standardHeatCapacitiesConstP()) .isApprox(Cp0_T)  );
        REQUIRE( grad(props.standardHeatCapacitiesConstV()) .isApprox(Cv0_T)  );

        REQUIRE( grad(props.gibbsEnergy())     == Approx(Gtot_T) );
        REQUIRE( grad(props.enthalpy())        == Approx(Htot_T) );
        REQUIRE( grad(props.volume())          == Approx(Vtot_T) );
        REQUIRE( grad(props.entropy())         == Approx(Stot_T) );
        REQUIRE( grad(props.internalEnergy())  == Approx(Utot_T) );
        REQUIRE( grad(props.helmholtzEnergy()) == Approx(Atot_T) );

        //---------------------------------------------------------------------
        // Testing pressure derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXd n_P = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd x_P = ArrayXd{{ 0.0, 0.0, 0.0 }};

        const ArrayXd nsumphases_P = ArrayXd{{ 0.0, 0.0 }};

        const ArrayXd  G0_P = ArrayXd{{ 0.1, 0.1, 1.1 }} * 2*T*(T*P);
        const ArrayXd  H0_P = ArrayXd{{ 0.2, 0.2, 1.2 }} * 2*T*(T*P);
        const ArrayXd  V0_P = ArrayXd{{ 0.3, 0.3, 1.3 }} * 2*T*(T*P);
        const ArrayXd Cp0_P = ArrayXd{{ 0.4, 0.4, 1.4 }} * 2*T*(T*P);
        const ArrayXd Cv0_P = ArrayXd{{ 0.5, 0.5, 1.5 }} * 2*T*(T*P);
        const ArrayXd  S0_P = (H0_P - G0_P)/T;
        const ArrayXd  U0_P = H0_P - V0 - P*V0_P;
        const ArrayXd  A0_P = G0_P - V0 - P*V0_P;

        const ArrayXd  Vex_P = ArrayXd{{ 1.0, 1.1 }} * 2*T*(T*P);
        const ArrayXd VexT_P = ArrayXd{{ 2.0, 2.1 }} * 2*T*(T*P);
        const ArrayXd VexP_P = ArrayXd{{ 3.0, 3.1 }} * 2*T*(T*P);
        const ArrayXd  Gex_P = ArrayXd{{ 4.0, 4.1 }} * 2*T*(T*P);
        const ArrayXd  Hex_P = ArrayXd{{ 5.0, 5.1 }} * 2*T*(T*P);
        const ArrayXd Cpex_P = ArrayXd{{ 6.0, 6.1 }} * 2*T*(T*P);
        const ArrayXd Cvex_P = ArrayXd{{ 7.0, 7.1 }} * 2*T*(T*P);

        const ArrayXd ln_g_P = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd ln_a_P = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd    u_P = G0_P;

        const double Gtot_P = (G0_P * n).sum() + (nsumphases * Gex_P).sum();
        const double Htot_P = (H0_P * n).sum() + (nsumphases * Hex_P).sum();
        const double Vtot_P = (V0_P * n).sum() + (nsumphases * Vex_P).sum();
        const double Stot_P = (Htot_P - Gtot_P)/T;
        const double Utot_P = Htot_P - Vtot - P*Vtot_P;
        const double Atot_P = Gtot_P - Vtot - P*Vtot_P;

        REQUIRE_NOTHROW( props.update(T, P, n, wrt(P)) );

        REQUIRE( grad(props.temperature()) == 0.0 );
        REQUIRE( grad(props.pressure())    == 1.0 );

        REQUIRE( grad(props.speciesAmounts())               .isApprox(n_P)    );
        REQUIRE( grad(props.moleFractions())                .isApprox(x_P)    );
        REQUIRE( grad(props.lnActivityCoefficients())       .isApprox(ln_g_P) );
        REQUIRE( grad(props.lnActivities())                 .isApprox(ln_a_P) );
        REQUIRE( grad(props.chemicalPotentials())           .isApprox(u_P)    );
        REQUIRE( grad(props.standardGibbsEnergies())        .isApprox(G0_P)   );
        REQUIRE( grad(props.standardEnthalpies())           .isApprox(H0_P)   );
        REQUIRE( grad(props.standardVolumes())              .isApprox(V0_P)   );
        REQUIRE( grad(props.standardEntropies())            .isApprox(S0_P)   );
        REQUIRE( grad(props.standardInternalEnergies())     .isApprox(U0_P)   );
        REQUIRE( grad(props.standardHelmholtzEnergies())    .isApprox(A0_P)   );
        REQUIRE( grad(props.standardHeatCapacitiesConstP()) .isApprox(Cp0_P)  );
        REQUIRE( grad(props.standardHeatCapacitiesConstV()) .isApprox(Cv0_P)  );

        REQUIRE( grad(props.gibbsEnergy())     == Approx(Gtot_P) );
        REQUIRE( grad(props.enthalpy())        == Approx(Htot_P) );
        REQUIRE( grad(props.volume())          == Approx(Vtot_P) );
        REQUIRE( grad(props.entropy())         == Approx(Stot_P) );
        REQUIRE( grad(props.internalEnergy())  == Approx(Utot_P) );
        REQUIRE( grad(props.helmholtzEnergy()) == Approx(Atot_P) );

        //---------------------------------------------------------------------
        // Testing compositional derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXXd n_n = MatrixXd::Identity(3, 3);
        const ArrayXXd x_n = ArrayXXd{{
            // x = [xA, xB] = [ [0.4, 0.6], [1.0] ] for phases A and B
            { 0.06, -0.04, 0.00},
            {-0.06,  0.04, 0.00},
            { 0.00,  0.00, 0.00},
        }};

        const ArrayXXd nsumphases_n = ArrayXXd{{
            {1.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
        }};

        const ArrayXXd  G0_n = ArrayXXd::Zero(3, 3);
        const ArrayXXd  H0_n = ArrayXXd::Zero(3, 3);
        const ArrayXXd  V0_n = ArrayXXd::Zero(3, 3);
        const ArrayXXd Cp0_n = ArrayXXd::Zero(3, 3);
        const ArrayXXd Cv0_n = ArrayXXd::Zero(3, 3);
        const ArrayXXd  S0_n = ArrayXXd::Zero(3, 3);
        const ArrayXXd  U0_n = ArrayXXd::Zero(3, 3);
        const ArrayXXd  A0_n = ArrayXXd::Zero(3, 3);

        const ArrayXXd  Vex_n = ArrayXXd::Zero(2, 3);
        const ArrayXXd VexT_n = ArrayXXd::Zero(2, 3);
        const ArrayXXd VexP_n = ArrayXXd::Zero(2, 3);
        const ArrayXXd  Gex_n = ArrayXXd::Zero(2, 3);
        const ArrayXXd  Hex_n = ArrayXXd::Zero(2, 3);
        const ArrayXXd Cpex_n = ArrayXXd::Zero(2, 3);
        const ArrayXXd Cvex_n = ArrayXXd::Zero(2, 3);

        const ArrayXXd ln_g_n = ArrayXXd{{
            // ln_g = { 8.0*x[0], 8.0*x[1], 8.1*x[2] }
            { 8.0*x_n(0, 0), 8.0*x_n(0, 1), 8.1*x_n(0, 2) },
            { 8.0*x_n(1, 0), 8.0*x_n(1, 1), 8.1*x_n(1, 2) },
            { 8.0*x_n(2, 0), 8.0*x_n(2, 1), 8.1*x_n(2, 2) },
        }};

        const ArrayXXd ln_a_n = ArrayXXd{{
            // ln_a = { 9.0*x[0], 9.0*x[1], 9.1*x[2] }
            { 9.0*x_n(0, 0), 9.0*x_n(0, 1), 9.1*x_n(0, 2) },
            { 9.0*x_n(1, 0), 9.0*x_n(1, 1), 9.1*x_n(1, 2) },
            { 9.0*x_n(2, 0), 9.0*x_n(2, 1), 9.1*x_n(2, 2) },
        }};

        const ArrayXXd u_n = R*T*ln_a_n;

        auto dot = [](auto A, auto x)
        {
            return (A.matrix().transpose() * x.matrix()).array();
        };

        const ArrayXd Gtot_n = dot(n_n, G0) + dot(nsumphases_n, Gex);
        const ArrayXd Htot_n = dot(n_n, H0) + dot(nsumphases_n, Hex);
        const ArrayXd Vtot_n = dot(n_n, V0) + dot(nsumphases_n, Vex);
        const ArrayXd Stot_n = (Htot_n - Gtot_n)/T;
        const ArrayXd Utot_n = Htot_n - P*Vtot_n;
        const ArrayXd Atot_n = Gtot_n - P*Vtot_n;

        for(auto i = 0; i < 3; ++i)
        {
            REQUIRE_NOTHROW( props.update(T, P, n, wrt(n[i])) );

            REQUIRE( grad(props.temperature()) == 0.0 );
            REQUIRE( grad(props.pressure())    == 0.0 );

            REQUIRE( grad(props.speciesAmounts())               .isApprox(    n_n.col(i)) );
            REQUIRE( grad(props.moleFractions())                .isApprox(    x_n.col(i)) );
            REQUIRE( grad(props.lnActivityCoefficients())       .isApprox( ln_g_n.col(i)) );
            REQUIRE( grad(props.lnActivities())                 .isApprox( ln_a_n.col(i)) );
            REQUIRE( grad(props.chemicalPotentials())           .isApprox(    u_n.col(i)) );
            REQUIRE( grad(props.standardGibbsEnergies())        .isApprox(   G0_n.col(i)) );
            REQUIRE( grad(props.standardEnthalpies())           .isApprox(   H0_n.col(i)) );
            REQUIRE( grad(props.standardVolumes())              .isApprox(   V0_n.col(i)) );
            REQUIRE( grad(props.standardEntropies())            .isApprox(   S0_n.col(i)) );
            REQUIRE( grad(props.standardInternalEnergies())     .isApprox(   U0_n.col(i)) );
            REQUIRE( grad(props.standardHelmholtzEnergies())    .isApprox(   A0_n.col(i)) );
            REQUIRE( grad(props.standardHeatCapacitiesConstP()) .isApprox(  Cp0_n.col(i)) );
            REQUIRE( grad(props.standardHeatCapacitiesConstV()) .isApprox(  Cv0_n.col(i)) );

            REQUIRE( grad(props.gibbsEnergy())     == Approx(Gtot_n[i]) );
            REQUIRE( grad(props.enthalpy())        == Approx(Htot_n[i]) );
            REQUIRE( grad(props.volume())          == Approx(Vtot_n[i]) );
            REQUIRE( grad(props.entropy())         == Approx(Stot_n[i]) );
            REQUIRE( grad(props.internalEnergy())  == Approx(Utot_n[i]) );
            REQUIRE( grad(props.helmholtzEnergy()) == Approx(Atot_n[i]) );
        }
    }

    SECTION("testing when species have zero mole fractions")
    {
        const real T = 300.0;
        const real P = 123.0e5;

        // The gas phase has last species with zero mole fraction!
        const ArrayXr n1 = ArrayXr{{ 1.0, 0.0, 1.0 }};

        REQUIRE_THROWS( props.update(T, P, n1) );

        // The single-species solid phase has unit mole fraction!
        const ArrayXr n2 = ArrayXr{{ 1.0, 1.0, 0.0 }};

        REQUIRE_NOTHROW( props.update(T, P, n2) );
    }
}
