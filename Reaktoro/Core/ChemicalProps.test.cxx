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
#include <Reaktoro/Common/AutoDiff.hpp>
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
using namespace Reaktoro;

TEST_CASE("Testing ChemicalProps class", "[ChemicalProps]")
{
    const auto R = universalGasConstant;

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

    ActivityModel activity_model_gas = [](ActivityPropsRef props, ActivityArgs args)
    {
        const auto [T, P, x] = args;
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

    ActivityModel activity_model_solid = [](ActivityPropsRef props, ActivityArgs args)
    {
        const auto [T, P, x] = args;
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

    ChemicalProps props(system);

    const auto molar_masses = ArrayXd{{
        system.species(0).molarMass(),
        system.species(1).molarMass(),
        system.species(2).molarMass(),
    }};

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

        const real Ntot = n.sum();
        const real Mtot = (n * molar_masses).sum();
        const real Gtot = (G0 * n).sum() + (nsumphases * Gex).sum();
        const real Htot = (H0 * n).sum() + (nsumphases * Hex).sum();
        const real Vtot = (V0 * n).sum() + (nsumphases * Vex).sum();
        const real Stot = (Htot - Gtot)/T;
        const real Utot = Htot - P*Vtot;
        const real Atot = Gtot - P*Vtot;

        CHECK_NOTHROW( props.update(T, P, n) );

        CHECK( props.temperature() == T );
        CHECK( props.pressure()    == P );

        CHECK( props.speciesAmounts()              .isApprox(n)    );
        CHECK( props.moleFractions()               .isApprox(x)    );
        CHECK( props.lnActivityCoefficients()      .isApprox(ln_g) );
        CHECK( props.lnActivities()                .isApprox(ln_a) );
        CHECK( props.chemicalPotentials()          .isApprox(u)    );
        CHECK( props.standardVolumes()             .isApprox(V0)   );
        CHECK( props.standardGibbsEnergies()       .isApprox(G0)   );
        CHECK( props.standardEnthalpies()          .isApprox(H0)   );
        CHECK( props.standardEntropies()           .isApprox(S0)   );
        CHECK( props.standardInternalEnergies()    .isApprox(U0)   );
        CHECK( props.standardHelmholtzEnergies()   .isApprox(A0)   );
        CHECK( props.standardHeatCapacitiesConstP().isApprox(Cp0)  );
        CHECK( props.standardHeatCapacitiesConstV().isApprox(Cv0)  );

        CHECK( props.amount()          == Approx(Ntot) );
        CHECK( props.mass()            == Approx(Mtot) );
        CHECK( props.volume()          == Approx(Vtot) );
        CHECK( props.gibbsEnergy()     == Approx(Gtot) );
        CHECK( props.enthalpy()        == Approx(Htot) );
        CHECK( props.entropy()         == Approx(Stot) );
        CHECK( props.internalEnergy()  == Approx(Utot) );
        CHECK( props.helmholtzEnergy() == Approx(Atot) );

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

        const double Ntot_T = 0.0;
        const double Mtot_T = 0.0;
        const double Gtot_T = (G0_T * n).sum() + (nsumphases * Gex_T).sum();
        const double Htot_T = (H0_T * n).sum() + (nsumphases * Hex_T).sum();
        const double Vtot_T = (V0_T * n).sum() + (nsumphases * Vex_T).sum();
        const double Stot_T = (Htot_T - Gtot_T)/T - (Htot - Gtot)/(T*T);
        const double Utot_T = Htot_T - P*Vtot_T;
        const double Atot_T = Gtot_T - P*Vtot_T;

        autodiff::seed(T);
        props.update(T, P, n);
        autodiff::unseed(T);

        CHECK( grad(props.temperature()) == 1.0 );
        CHECK( grad(props.pressure())    == 0.0 );

        CHECK( grad(props.speciesAmounts())               .isApprox(n_T)    );
        CHECK( grad(props.moleFractions())                .isApprox(x_T)    );
        CHECK( grad(props.lnActivityCoefficients())       .isApprox(ln_g_T) );
        CHECK( grad(props.lnActivities())                 .isApprox(ln_a_T) );
        CHECK( grad(props.chemicalPotentials())           .isApprox(u_T)    );
        CHECK( grad(props.standardVolumes())              .isApprox(V0_T)   );
        CHECK( grad(props.standardGibbsEnergies())        .isApprox(G0_T)   );
        CHECK( grad(props.standardEnthalpies())           .isApprox(H0_T)   );
        CHECK( grad(props.standardEntropies())            .isApprox(S0_T)   );
        CHECK( grad(props.standardInternalEnergies())     .isApprox(U0_T)   );
        CHECK( grad(props.standardHelmholtzEnergies())    .isApprox(A0_T)   );
        CHECK( grad(props.standardHeatCapacitiesConstP()) .isApprox(Cp0_T)  );
        CHECK( grad(props.standardHeatCapacitiesConstV()) .isApprox(Cv0_T)  );

        CHECK( grad(props.amount())          == Approx(Ntot_T) );
        CHECK( grad(props.mass())            == Approx(Mtot_T) );
        CHECK( grad(props.volume())          == Approx(Vtot_T) );
        CHECK( grad(props.gibbsEnergy())     == Approx(Gtot_T) );
        CHECK( grad(props.enthalpy())        == Approx(Htot_T) );
        CHECK( grad(props.entropy())         == Approx(Stot_T) );
        CHECK( grad(props.internalEnergy())  == Approx(Utot_T) );
        CHECK( grad(props.helmholtzEnergy()) == Approx(Atot_T) );

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

        const double Ntot_P = 0.0;
        const double Mtot_P = 0.0;
        const double Gtot_P = (G0_P * n).sum() + (nsumphases * Gex_P).sum();
        const double Htot_P = (H0_P * n).sum() + (nsumphases * Hex_P).sum();
        const double Vtot_P = (V0_P * n).sum() + (nsumphases * Vex_P).sum();
        const double Stot_P = (Htot_P - Gtot_P)/T;
        const double Utot_P = Htot_P - Vtot - P*Vtot_P;
        const double Atot_P = Gtot_P - Vtot - P*Vtot_P;

        autodiff::seed(P);
        props.update(T, P, n);
        autodiff::unseed(P);

        CHECK( grad(props.temperature()) == 0.0 );
        CHECK( grad(props.pressure())    == 1.0 );

        CHECK( grad(props.speciesAmounts())               .isApprox(n_P)    );
        CHECK( grad(props.moleFractions())                .isApprox(x_P)    );
        CHECK( grad(props.lnActivityCoefficients())       .isApprox(ln_g_P) );
        CHECK( grad(props.lnActivities())                 .isApprox(ln_a_P) );
        CHECK( grad(props.standardVolumes())              .isApprox(V0_P)   );
        CHECK( grad(props.chemicalPotentials())           .isApprox(u_P)    );
        CHECK( grad(props.standardGibbsEnergies())        .isApprox(G0_P)   );
        CHECK( grad(props.standardEnthalpies())           .isApprox(H0_P)   );
        CHECK( grad(props.standardEntropies())            .isApprox(S0_P)   );
        CHECK( grad(props.standardInternalEnergies())     .isApprox(U0_P)   );
        CHECK( grad(props.standardHelmholtzEnergies())    .isApprox(A0_P)   );
        CHECK( grad(props.standardHeatCapacitiesConstP()) .isApprox(Cp0_P)  );
        CHECK( grad(props.standardHeatCapacitiesConstV()) .isApprox(Cv0_P)  );

        CHECK( grad(props.amount())          == Approx(Ntot_P) );
        CHECK( grad(props.mass())            == Approx(Mtot_P) );
        CHECK( grad(props.volume())          == Approx(Vtot_P) );
        CHECK( grad(props.gibbsEnergy())     == Approx(Gtot_P) );
        CHECK( grad(props.enthalpy())        == Approx(Htot_P) );
        CHECK( grad(props.entropy())         == Approx(Stot_P) );
        CHECK( grad(props.internalEnergy())  == Approx(Utot_P) );
        CHECK( grad(props.helmholtzEnergy()) == Approx(Atot_P) );

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
            return ((A.matrix().transpose() * x.matrix()).array()).eval();
        };

        const ArrayXd Ntot_n = ArrayXd::Ones(3);
        const ArrayXd Mtot_n = molar_masses;
        const ArrayXd Gtot_n = dot(n_n, G0) + dot(nsumphases_n, Gex);
        const ArrayXd Htot_n = dot(n_n, H0) + dot(nsumphases_n, Hex);
        const ArrayXd Vtot_n = dot(n_n, V0) + dot(nsumphases_n, Vex);
        const ArrayXd Stot_n = (Htot_n - Gtot_n)/T;
        const ArrayXd Utot_n = Htot_n - P*Vtot_n;
        const ArrayXd Atot_n = Gtot_n - P*Vtot_n;

        for(auto i = 0; i < 3; ++i)
        {
            autodiff::seed(n[i]);
            props.update(T, P, n);
            autodiff::unseed(n[i]);

            CHECK( grad(props.temperature()) == 0.0 );
            CHECK( grad(props.pressure())    == 0.0 );

            CHECK( grad(props.speciesAmounts())               .isApprox(    n_n.col(i)) );
            CHECK( grad(props.moleFractions())                .isApprox(    x_n.col(i)) );
            CHECK( grad(props.lnActivityCoefficients())       .isApprox( ln_g_n.col(i)) );
            CHECK( grad(props.lnActivities())                 .isApprox( ln_a_n.col(i)) );
            CHECK( grad(props.chemicalPotentials())           .isApprox(    u_n.col(i)) );
            CHECK( grad(props.standardVolumes())              .isApprox(   V0_n.col(i)) );
            CHECK( grad(props.standardGibbsEnergies())        .isApprox(   G0_n.col(i)) );
            CHECK( grad(props.standardEnthalpies())           .isApprox(   H0_n.col(i)) );
            CHECK( grad(props.standardEntropies())            .isApprox(   S0_n.col(i)) );
            CHECK( grad(props.standardInternalEnergies())     .isApprox(   U0_n.col(i)) );
            CHECK( grad(props.standardHelmholtzEnergies())    .isApprox(   A0_n.col(i)) );
            CHECK( grad(props.standardHeatCapacitiesConstP()) .isApprox(  Cp0_n.col(i)) );
            CHECK( grad(props.standardHeatCapacitiesConstV()) .isApprox(  Cv0_n.col(i)) );

            CHECK( grad(props.amount())          == Approx(Ntot_n[i]) );
            CHECK( grad(props.mass())            == Approx(Mtot_n[i]) );
            CHECK( grad(props.gibbsEnergy())     == Approx(Gtot_n[i]) );
            CHECK( grad(props.enthalpy())        == Approx(Htot_n[i]) );
            CHECK( grad(props.volume())          == Approx(Vtot_n[i]) );
            CHECK( grad(props.entropy())         == Approx(Stot_n[i]) );
            CHECK( grad(props.internalEnergy())  == Approx(Utot_n[i]) );
            CHECK( grad(props.helmholtzEnergy()) == Approx(Atot_n[i]) );
        }
    }

    SECTION("testing when species have zero mole fractions")
    {
        const real T = 11.0;
        const real P = 13.0;

        // The gas phase has last species with zero mole fraction!
        const ArrayXr n1 = ArrayXr{{ 1.0, 0.0, 1.0 }};

        CHECK_THROWS( props.update(T, P, n1) );

        // The single-species solid phase has unit mole fraction!
        const ArrayXr n2 = ArrayXr{{ 1.0, 1.0, 0.0 }};

        CHECK_NOTHROW( props.update(T, P, n2) );
    }

    SECTION("testing convenience methods")
    {
        const real T = 11.0;
        const real P = 13.0;
        const ArrayXr n = ArrayXr{{ 1.0, 2.0, 3.0 }};

        props.update(T, P, n);

        const auto A = system.formulaMatrixElements();

        const auto idxElement = [&](auto symbol) { return system.elements().index(symbol); };
        const auto idxSpecies = [&](auto name) { return system.species().index(name); };

        const VectorXr b = A * n.matrix();
        const VectorXr b0 = A.leftCols(2) * n.head(2).matrix();  // element amounts in phase 0: H2O(g) CO2(g)
        const VectorXr b1 = A.rightCols(1) * n.tail(1).matrix(); // element amounts in phase 1: CaCO3(s)

        CHECK( props.elementAmounts().matrix().isApprox(b) );
        CHECK( props.elementAmountsInPhase(0).matrix().isApprox(b0) );
        CHECK( props.elementAmountsInPhase(1).matrix().isApprox(b1) );
        CHECK( props.elementAmountsInPhase("SomeGas").matrix().isApprox(b0) );
        CHECK( props.elementAmountsInPhase("SomeSolid").matrix().isApprox(b1) );
        CHECK( props.elementAmountsAmongSpecies(ArrayXl{{0, 1}}).matrix().isApprox(b0) );
        CHECK( props.elementAmountsAmongSpecies(ArrayXl{{2}}).matrix().isApprox(b1) );

        CHECK( props.elementAmount("H")  == Approx(b[idxElement("H")]) );
        CHECK( props.elementAmount("O")  == Approx(b[idxElement("O")]) );
        CHECK( props.elementAmount("C")  == Approx(b[idxElement("C")]) );
        CHECK( props.elementAmount("Ca") == Approx(b[idxElement("Ca")]) );

        CHECK( props.elementAmount(0) == Approx(b[0]) );
        CHECK( props.elementAmount(1) == Approx(b[1]) );
        CHECK( props.elementAmount(2) == Approx(b[2]) );
        CHECK( props.elementAmount(3) == Approx(b[3]) );

        CHECK( props.elementAmountInPhase("H", "SomeGas")  == Approx(b0[idxElement("H")]) );
        CHECK( props.elementAmountInPhase("O", "SomeGas")  == Approx(b0[idxElement("O")]) );
        CHECK( props.elementAmountInPhase("C", "SomeGas")  == Approx(b0[idxElement("C")]) );
        CHECK( props.elementAmountInPhase("Ca", "SomeGas") == Approx(b0[idxElement("Ca")]) );

        CHECK( props.elementAmountInPhase("H", "SomeSolid")  == Approx(b1[idxElement("H")]) );
        CHECK( props.elementAmountInPhase("O", "SomeSolid")  == Approx(b1[idxElement("O")]) );
        CHECK( props.elementAmountInPhase("C", "SomeSolid")  == Approx(b1[idxElement("C")]) );
        CHECK( props.elementAmountInPhase("Ca", "SomeSolid") == Approx(b1[idxElement("Ca")]) );

        CHECK( props.elementAmountInPhase(0, 0) == Approx(b0[0]) );
        CHECK( props.elementAmountInPhase(1, 0) == Approx(b0[1]) );
        CHECK( props.elementAmountInPhase(2, 0) == Approx(b0[2]) );
        CHECK( props.elementAmountInPhase(3, 0) == Approx(b0[3]) );

        CHECK( props.elementAmountInPhase(0, 1) == Approx(b1[0]) );
        CHECK( props.elementAmountInPhase(1, 1) == Approx(b1[1]) );
        CHECK( props.elementAmountInPhase(2, 1) == Approx(b1[2]) );
        CHECK( props.elementAmountInPhase(3, 1) == Approx(b1[3]) );

        auto i = 0;
        for(const auto& species : system.species())
        {
            const auto name = species.name();
            const auto idx = system.species().index(name);
            CHECK( props.speciesAmount(name)              == Approx(props.speciesAmounts()[idx])               );
            CHECK( props.speciesMass(name)                == Approx(props.speciesMasses()[idx])                );
            CHECK( props.moleFraction(name)               == Approx(props.moleFractions()[idx])                );
            CHECK( props.activityCoefficient(name)        == Approx(exp(props.lnActivityCoefficients()[idx]))  );
            CHECK( props.lgActivityCoefficient(name)      == Approx(props.lnActivityCoefficients()[idx]/ln10)  );
            CHECK( props.lnActivityCoefficient(name)      == Approx(props.lnActivityCoefficients()[idx])       );
            CHECK( props.activity(name)                   == Approx(exp(props.lnActivities()[idx]))            );
            CHECK( props.lgActivity(name)                 == Approx(props.lnActivities()[idx]/ln10)            );
            CHECK( props.lnActivity(name)                 == Approx(props.lnActivities()[idx])                 );
            CHECK( props.chemicalPotential(name)          == Approx(props.chemicalPotentials()[idx])           );
            CHECK( props.standardVolume(name)             == Approx(props.standardVolumes()[idx])              );
            CHECK( props.standardGibbsEnergy(name)        == Approx(props.standardGibbsEnergies()[idx])        );
            CHECK( props.standardEnthalpy(name)           == Approx(props.standardEnthalpies()[idx])           );
            CHECK( props.standardEntropy(name)            == Approx(props.standardEntropies()[idx])            );
            CHECK( props.standardInternalEnergy(name)     == Approx(props.standardInternalEnergies()[idx])     );
            CHECK( props.standardHelmholtzEnergy(name)    == Approx(props.standardHelmholtzEnergies()[idx])    );
            CHECK( props.standardHeatCapacityConstP(name) == Approx(props.standardHeatCapacitiesConstP()[idx]) );
            CHECK( props.standardHeatCapacityConstV(name) == Approx(props.standardHeatCapacitiesConstV()[idx]) );

            CHECK( props.speciesAmount(i)              == Approx(props.speciesAmounts()[idx])               );
            CHECK( props.speciesMass(i)                == Approx(props.speciesMasses()[idx])                );
            CHECK( props.moleFraction(i)               == Approx(props.moleFractions()[idx])                );
            CHECK( props.activityCoefficient(i)        == Approx(exp(props.lnActivityCoefficients()[idx]))  );
            CHECK( props.lgActivityCoefficient(i)      == Approx(props.lnActivityCoefficients()[idx]/ln10)  );
            CHECK( props.lnActivityCoefficient(i)      == Approx(props.lnActivityCoefficients()[idx])       );
            CHECK( props.activity(i)                   == Approx(exp(props.lnActivities()[idx]))            );
            CHECK( props.lgActivity(i)                 == Approx(props.lnActivities()[idx]/ln10)            );
            CHECK( props.lnActivity(i)                 == Approx(props.lnActivities()[idx])                 );
            CHECK( props.chemicalPotential(i)          == Approx(props.chemicalPotentials()[idx])           );
            CHECK( props.standardVolume(i)             == Approx(props.standardVolumes()[idx])              );
            CHECK( props.standardGibbsEnergy(i)        == Approx(props.standardGibbsEnergies()[idx])        );
            CHECK( props.standardEnthalpy(i)           == Approx(props.standardEnthalpies()[idx])           );
            CHECK( props.standardEntropy(i)            == Approx(props.standardEntropies()[idx])            );
            CHECK( props.standardInternalEnergy(i)     == Approx(props.standardInternalEnergies()[idx])     );
            CHECK( props.standardHelmholtzEnergy(i)    == Approx(props.standardHelmholtzEnergies()[idx])    );
            CHECK( props.standardHeatCapacityConstP(i) == Approx(props.standardHeatCapacitiesConstP()[idx]) );
            CHECK( props.standardHeatCapacityConstV(i) == Approx(props.standardHeatCapacitiesConstV()[idx]) );

            i += 1;
        }
    }
}
