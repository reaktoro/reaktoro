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
#include <Reaktoro/Common/MoleFractionUtils.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
using namespace Reaktoro;

inline auto createStandardThermoPropsFn(double param)
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

TEST_CASE("Testing ChemicalPropsPhase class", "[ChemicalPropsPhase]")
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

    const auto molar_masses = ArrayXd{{
        phase.species(0).molarMass(),
        phase.species(1).molarMass(),
        phase.species(2).molarMass(),
        phase.species(3).molarMass(),
    }};

    ChemicalPropsPhase props(phase);

    SECTION("Testing when species have non-zero amounts")
    {
        real T = 300.0;
        real P = 123.0e5;
        ArrayXr n = ArrayXr{{ 1.0, 2.0, 3.0, 4.0 }};

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

        const real nsum = n.sum();
        const real mass = (n * molar_masses).sum();

        const real rho = 1.0 / V;

        const real Gtot = G * nsum;
        const real Htot = H * nsum;
        const real Vtot = V * nsum;
        const real Stot = S * nsum;
        const real Utot = U * nsum;
        const real Atot = A * nsum;

        REQUIRE_NOTHROW( props.update(T, P, n) );

        REQUIRE( props.temperature() == T );
        REQUIRE( props.pressure()    == P );

        REQUIRE( props.speciesAmounts()               .isApprox(n)    );
        REQUIRE( props.moleFractions()                .isApprox(x)    );
        REQUIRE( props.lnActivityCoefficients()       .isApprox(ln_g) );
        REQUIRE( props.lnActivities()                 .isApprox(ln_a) );
        REQUIRE( props.chemicalPotentials()           .isApprox(u)    );
        REQUIRE( props.standardGibbsEnergies()        .isApprox(G0)   );
        REQUIRE( props.standardEnthalpies()           .isApprox(H0)   );
        REQUIRE( props.standardVolumes()              .isApprox(V0)   );
        REQUIRE( props.standardEntropies()            .isApprox(S0)   );
        REQUIRE( props.standardInternalEnergies()     .isApprox(U0)   );
        REQUIRE( props.standardHelmholtzEnergies()    .isApprox(A0)   );
        REQUIRE( props.standardHeatCapacitiesConstP() .isApprox(Cp0)  );
        REQUIRE( props.standardHeatCapacitiesConstV() .isApprox(Cv0)  );

        REQUIRE( props.molarGibbsEnergy()        == Approx(G)    );
        REQUIRE( props.molarEnthalpy()           == Approx(H)    );
        REQUIRE( props.molarVolume()             == Approx(V)    );
        REQUIRE( props.molarEntropy()            == Approx(S)    );
        REQUIRE( props.molarInternalEnergy()     == Approx(U)    );
        REQUIRE( props.molarHelmholtzEnergy()    == Approx(A)    );
        REQUIRE( props.molarHeatCapacityConstP() == Approx(Cp)   );
        REQUIRE( props.molarHeatCapacityConstV() == Approx(Cv)   );
        REQUIRE( props.molarDensity()            == Approx(rho)  );
        REQUIRE( props.amount()                  == Approx(nsum) );
        REQUIRE( props.mass()                    == Approx(mass) );
        REQUIRE( props.gibbsEnergy()             == Approx(Gtot) );
        REQUIRE( props.enthalpy()                == Approx(Htot) );
        REQUIRE( props.volume()                  == Approx(Vtot) );
        REQUIRE( props.entropy()                 == Approx(Stot) );
        REQUIRE( props.internalEnergy()          == Approx(Utot) );
        REQUIRE( props.helmholtzEnergy()         == Approx(Atot) );

        //---------------------------------------------------------------------
        // Testing temperature derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXd n_T = ArrayXd{{ 0.0, 0.0, 0.0, 0.0 }};
        const ArrayXd x_T = ArrayXd{{ 0.0, 0.0, 0.0, 0.0 }};

        const ArrayXd  G0_T = 0.1 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd  H0_T = 0.2 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd  V0_T = 0.3 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd Cp0_T = 0.4 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd Cv0_T = 0.5 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd  S0_T = (H0_T - G0_T)/T - (H0 - G0)/(T*T);
        const ArrayXd  U0_T = H0_T - P*V0_T;
        const ArrayXd  A0_T = G0_T - P*V0_T;

        const double  Vex_T = 1.0 * P;
        const double VexT_T = 2.0 * P;
        const double VexP_T = 3.0 * P;
        const double  Gex_T = 4.0 * P;
        const double  Hex_T = 5.0 * P;
        const double Cpex_T = 6.0 * P;
        const double Cvex_T = 7.0 * P;

        const ArrayXd ln_g_T = ArrayXd{{0.0, 0.0, 0.0, 0.0}};
        const ArrayXd ln_a_T = ArrayXd{{0.0, 0.0, 0.0, 0.0}};
        const ArrayXd    u_T = G0_T + R*ln_a;

        const double  G_T = ( G0_T * x).sum() +  Gex_T;
        const double  H_T = ( H0_T * x).sum() +  Hex_T;
        const double  V_T = ( V0_T * x).sum() +  Vex_T;
        const double Cp_T = (Cp0_T * x).sum() + Cpex_T;
        const double Cv_T = (Cv0_T * x).sum() + Cvex_T;
        const double  S_T = (H_T - G_T)/T - (H - G)/(T*T);
        const double  U_T = H_T - P*V_T;
        const double  A_T = G_T - P*V_T;

        const double rho_T = -V_T/(V*V);

        const double Gtot_T = G_T * nsum;
        const double Htot_T = H_T * nsum;
        const double Vtot_T = V_T * nsum;
        const double Stot_T = S_T * nsum;
        const double Utot_T = U_T * nsum;
        const double Atot_T = A_T * nsum;

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

        REQUIRE( grad(props.molarGibbsEnergy())        == Approx(G_T)    );
        REQUIRE( grad(props.molarEnthalpy())           == Approx(H_T)    );
        REQUIRE( grad(props.molarVolume())             == Approx(V_T)    );
        REQUIRE( grad(props.molarEntropy())            == Approx(S_T)    );
        REQUIRE( grad(props.molarInternalEnergy())     == Approx(U_T)    );
        REQUIRE( grad(props.molarHelmholtzEnergy())    == Approx(A_T)    );
        REQUIRE( grad(props.molarHeatCapacityConstP()) == Approx(Cp_T)   );
        REQUIRE( grad(props.molarHeatCapacityConstV()) == Approx(Cv_T)   );
        REQUIRE( grad(props.molarDensity())            == Approx(rho_T)  );
        REQUIRE( grad(props.amount())                  == Approx(0.0)    );
        REQUIRE( grad(props.mass())                    == Approx(0.0)    );
        REQUIRE( grad(props.gibbsEnergy())             == Approx(Gtot_T) );
        REQUIRE( grad(props.enthalpy())                == Approx(Htot_T) );
        REQUIRE( grad(props.volume())                  == Approx(Vtot_T) );
        REQUIRE( grad(props.entropy())                 == Approx(Stot_T) );
        REQUIRE( grad(props.internalEnergy())          == Approx(Utot_T) );
        REQUIRE( grad(props.helmholtzEnergy())         == Approx(Atot_T) );

        //---------------------------------------------------------------------
        // Testing pressure derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXd n_P = ArrayXd{{ 0.0, 0.0, 0.0, 0.0 }};
        const ArrayXd x_P = ArrayXd{{ 0.0, 0.0, 0.0, 0.0 }};

        const ArrayXd  G0_P = 0.1 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd  H0_P = 0.2 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd  V0_P = 0.3 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd Cp0_P = 0.4 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd Cv0_P = 0.5 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd  S0_P = (H0_P - G0_P)/T;
        const ArrayXd  U0_P = H0_P - V0 - P*V0_P;
        const ArrayXd  A0_P = G0_P - V0 - P*V0_P;

        const double  Vex_P = 1.0 * T;
        const double VexT_P = 2.0 * T;
        const double VexP_P = 3.0 * T;
        const double  Gex_P = 4.0 * T;
        const double  Hex_P = 5.0 * T;
        const double Cpex_P = 6.0 * T;
        const double Cvex_P = 7.0 * T;

        const ArrayXd ln_g_P = ArrayXd{{0.0, 0.0, 0.0, 0.0}};
        const ArrayXd ln_a_P = ArrayXd{{0.0, 0.0, 0.0, 0.0}};
        const ArrayXd    u_P = G0_P;

        const double  G_P = ( G0_P * x).sum() +  Gex_P;
        const double  H_P = ( H0_P * x).sum() +  Hex_P;
        const double  V_P = ( V0_P * x).sum() +  Vex_P;
        const double Cp_P = (Cp0_P * x).sum() + Cpex_P;
        const double Cv_P = (Cv0_P * x).sum() + Cvex_P;
        const double  S_P = (H_P - G_P)/T;
        const double  U_P = H_P - V - P*V_P;
        const double  A_P = G_P - V - P*V_P;

        const double rho_P = -V_P/(V*V);

        const double Gtot_P = G_P * nsum;
        const double Htot_P = H_P * nsum;
        const double Vtot_P = V_P * nsum;
        const double Stot_P = S_P * nsum;
        const double Utot_P = U_P * nsum;
        const double Atot_P = A_P * nsum;

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

        REQUIRE( grad(props.molarGibbsEnergy())        == Approx(G_P)    );
        REQUIRE( grad(props.molarEnthalpy())           == Approx(H_P)    );
        REQUIRE( grad(props.molarVolume())             == Approx(V_P)    );
        REQUIRE( grad(props.molarEntropy())            == Approx(S_P)    );
        REQUIRE( grad(props.molarInternalEnergy())     == Approx(U_P)    );
        REQUIRE( grad(props.molarHelmholtzEnergy())    == Approx(A_P)    );
        REQUIRE( grad(props.molarHeatCapacityConstP()) == Approx(Cp_P)   );
        REQUIRE( grad(props.molarHeatCapacityConstV()) == Approx(Cv_P)   );
        REQUIRE( grad(props.molarDensity())            == Approx(rho_P)  );
        REQUIRE( grad(props.amount())                  == Approx(0.0)    );
        REQUIRE( grad(props.mass())                    == Approx(0.0)    );
        REQUIRE( grad(props.gibbsEnergy())             == Approx(Gtot_P) );
        REQUIRE( grad(props.enthalpy())                == Approx(Htot_P) );
        REQUIRE( grad(props.volume())                  == Approx(Vtot_P) );
        REQUIRE( grad(props.entropy())                 == Approx(Stot_P) );
        REQUIRE( grad(props.internalEnergy())          == Approx(Utot_P) );
        REQUIRE( grad(props.helmholtzEnergy())         == Approx(Atot_P) );

        //---------------------------------------------------------------------
        // Testing compositional derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXXd n_n = MatrixXd::Identity(4, 4);
        const ArrayXXd x_n = moleFractionsJacobian(n);

        const ArrayXXd  G0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd  H0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd  V0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd Cp0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd Cv0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd  S0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd  U0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd  A0_n = ArrayXXd::Zero(4, 4);

        const ArrayXd  Vex_n = ArrayXd::Zero(4);
        const ArrayXd VexT_n = ArrayXd::Zero(4);
        const ArrayXd VexP_n = ArrayXd::Zero(4);
        const ArrayXd  Gex_n = ArrayXd::Zero(4);
        const ArrayXd  Hex_n = ArrayXd::Zero(4);
        const ArrayXd Cpex_n = ArrayXd::Zero(4);
        const ArrayXd Cvex_n = ArrayXd::Zero(4);

        const ArrayXXd ln_g_n = 8.0 * x_n;
        const ArrayXXd ln_a_n = 9.0 * x_n;
        const ArrayXXd u_n    = R*T*ln_a_n;

        auto dot = [](auto A, auto x)
        {
            return (A.matrix().transpose() * x.matrix()).array();
        };

        const ArrayXd  G_n = dot(x_n,  G0) +  Gex_n;
        const ArrayXd  H_n = dot(x_n,  H0) +  Hex_n;
        const ArrayXd  V_n = dot(x_n,  V0) +  Vex_n;
        const ArrayXd Cp_n = dot(x_n, Cp0) + Cpex_n;
        const ArrayXd Cv_n = dot(x_n, Cv0) + Cvex_n;
        const ArrayXd  S_n = (H_n - G_n)/T;
        const ArrayXd  U_n = H_n - P*V_n;
        const ArrayXd  A_n = G_n - P*V_n;

        const ArrayXd nsum_n = ArrayXd::Ones(4);
        const ArrayXd mass_n = molar_masses;

        const ArrayXd rho_n = -V_n / (V*V);

        const ArrayXd Gtot_n = G_n * nsum + G * nsum_n;
        const ArrayXd Htot_n = H_n * nsum + H * nsum_n;
        const ArrayXd Vtot_n = V_n * nsum + V * nsum_n;
        const ArrayXd Stot_n = S_n * nsum + S * nsum_n;
        const ArrayXd Utot_n = U_n * nsum + U * nsum_n;
        const ArrayXd Atot_n = A_n * nsum + A * nsum_n;

        for(auto i = 0; i < 4; ++i)
        {
            REQUIRE_NOTHROW( props.update(T, P, n, wrt(n[i])) );

            REQUIRE( grad(props.temperature()) == 0.0 );
            REQUIRE( grad(props.pressure())    == 0.0 );

            REQUIRE( grad(props.speciesAmounts())               .isApprox(   n_n.col(i)) );
            REQUIRE( grad(props.moleFractions())                .isApprox(   x_n.col(i)) );
            REQUIRE( grad(props.lnActivityCoefficients())       .isApprox(ln_g_n.col(i)) );
            REQUIRE( grad(props.lnActivities())                 .isApprox(ln_a_n.col(i)) );
            REQUIRE( grad(props.chemicalPotentials())           .isApprox(   u_n.col(i)) );
            REQUIRE( grad(props.standardGibbsEnergies())        .isApprox(  G0_n.col(i)) );
            REQUIRE( grad(props.standardEnthalpies())           .isApprox(  H0_n.col(i)) );
            REQUIRE( grad(props.standardVolumes())              .isApprox(  V0_n.col(i)) );
            REQUIRE( grad(props.standardEntropies())            .isApprox(  S0_n.col(i)) );
            REQUIRE( grad(props.standardInternalEnergies())     .isApprox(  U0_n.col(i)) );
            REQUIRE( grad(props.standardHelmholtzEnergies())    .isApprox(  A0_n.col(i)) );
            REQUIRE( grad(props.standardHeatCapacitiesConstP()) .isApprox( Cp0_n.col(i)) );
            REQUIRE( grad(props.standardHeatCapacitiesConstV()) .isApprox( Cv0_n.col(i)) );

            REQUIRE( grad(props.molarGibbsEnergy())        == Approx(   G_n[i]) );
            REQUIRE( grad(props.molarEnthalpy())           == Approx(   H_n[i]) );
            REQUIRE( grad(props.molarVolume())             == Approx(   V_n[i]) );
            REQUIRE( grad(props.molarEntropy())            == Approx(   S_n[i]) );
            REQUIRE( grad(props.molarInternalEnergy())     == Approx(   U_n[i]) );
            REQUIRE( grad(props.molarHelmholtzEnergy())    == Approx(   A_n[i]) );
            REQUIRE( grad(props.molarHeatCapacityConstP()) == Approx(  Cp_n[i]) );
            REQUIRE( grad(props.molarHeatCapacityConstV()) == Approx(  Cv_n[i]) );
            REQUIRE( grad(props.molarDensity())            == Approx( rho_n[i]) );
            REQUIRE( grad(props.amount())                  == Approx(nsum_n[i]) );
            REQUIRE( grad(props.mass())                    == Approx(mass_n[i]) );
            REQUIRE( grad(props.gibbsEnergy())             == Approx(Gtot_n[i]) );
            REQUIRE( grad(props.enthalpy())                == Approx(Htot_n[i]) );
            REQUIRE( grad(props.volume())                  == Approx(Vtot_n[i]) );
            REQUIRE( grad(props.entropy())                 == Approx(Stot_n[i]) );
            REQUIRE( grad(props.internalEnergy())          == Approx(Utot_n[i]) );
            REQUIRE( grad(props.helmholtzEnergy())         == Approx(Atot_n[i]) );
        }
    }

    SECTION("Testing when species have zero amounts")
    {
        const real T = 300.0;
        const real P = 123.0e5;

        const ArrayXr n = ArrayXr{{ 0.0, 0.0, 0.0, 0.0 }};

        REQUIRE_THROWS( props.update(T, P, n) );
    }
}
