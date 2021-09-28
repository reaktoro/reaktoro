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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/MoleFractionUtils.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
using namespace Reaktoro;

template<typename T>
auto approx(const T& val)
{
    return Approx(val).scale(1.0);
}

inline auto createStandardThermoModel(double param)
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
        props.Cv0 = 0.7 * param * T * P;
        return props;
    };
    return model;
}

TEST_CASE("Testing ChemicalPropsPhase class", "[ChemicalPropsPhase]")
{
    const auto R = universalGasConstant;

    ActivityModel activity_model = [](ActivityPropsRef props, ActivityArgs args)
    {
        const auto [T, P, x] = args;
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
    phase = phase.withActivityModel(activity_model);
    phase = phase.withStateOfMatter(StateOfMatter::Gas);
    phase = phase.withSpecies({
        Species("H2O(g)").withStandardThermoModel(createStandardThermoModel(10.0)), // param = 10.0
        Species("CO2(g)").withStandardThermoModel(createStandardThermoModel(20.0)), // param = 20.0
        Species("CH4(g)").withStandardThermoModel(createStandardThermoModel(30.0)), // param = 30.0
        Species("H2S(g)").withStandardThermoModel(createStandardThermoModel(40.0))  // param = 40.0
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
        real T = 5.0;
        real P = 7.0;
        ArrayXr n = ArrayXr{{ 1.0, 2.0, 3.0, 4.0 }};

        const ArrayXr x = ArrayXr{{ 0.1, 0.2, 0.3, 0.4 }};

        const ArrayXr G0  = 0.1 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr H0  = 0.2 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr V0  = 0.3 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr VT0 = 0.4 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr VP0 = 0.5 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr Cp0 = 0.6 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
        const ArrayXr Cv0 = 0.7 * ArrayXr{{ 10.0, 20.0, 30.0, 40.0 }} * T * P;
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
        const real VT = (VT0 * x).sum() + VexT;
        const real VP = (VP0 * x).sum() + VexP;
        const real Cp = (Cp0 * x).sum() + Cpex;
        const real Cv = (Cv0 * x).sum() + Cvex;
        const real S  = (H - G)/T;
        const real U  = H - P*V;
        const real A  = G - P*V;

        const real nsum = n.sum();
        const real mass = (n * molar_masses).sum();

        const real rho = 1.0 / V;

        const real Gtot  = nsum * G;
        const real Htot  = nsum * H;
        const real Vtot  = nsum * V;
        const real VtotT = nsum * VT;
        const real VtotP = nsum * VP;
        const real Stot  = nsum * S;
        const real Utot  = nsum * U;
        const real Atot  = nsum * A;

        Map<String, Any> extra;

        CHECK_NOTHROW( props.update(T, P, n, extra) );

        CHECK( props.temperature() == T );
        CHECK( props.pressure()    == P );

        CHECK( props.speciesAmounts()               .isApprox(n)    );
        CHECK( props.moleFractions()                .isApprox(x)    );
        CHECK( props.lnActivityCoefficients()       .isApprox(ln_g) );
        CHECK( props.lnActivities()                 .isApprox(ln_a) );
        CHECK( props.chemicalPotentials()           .isApprox(u)    );
        CHECK( props.standardGibbsEnergies()        .isApprox(G0)   );
        CHECK( props.standardEnthalpies()           .isApprox(H0)   );
        CHECK( props.standardVolumes()              .isApprox(V0)   );
        CHECK( props.standardVolumesT()             .isApprox(VT0)  );
        CHECK( props.standardVolumesP()             .isApprox(VP0)  );
        CHECK( props.standardEntropies()            .isApprox(S0)   );
        CHECK( props.standardInternalEnergies()     .isApprox(U0)   );
        CHECK( props.standardHelmholtzEnergies()    .isApprox(A0)   );
        CHECK( props.standardHeatCapacitiesConstP() .isApprox(Cp0)  );
        CHECK( props.standardHeatCapacitiesConstV() .isApprox(Cv0)  );

        CHECK( props.molarGibbsEnergy()        == approx(G)     );
        CHECK( props.molarEnthalpy()           == approx(H)     );
        CHECK( props.molarVolume()             == approx(V)     );
        CHECK( props.molarVolumeT()            == approx(VT)    );
        CHECK( props.molarVolumeP()            == approx(VP)    );
        CHECK( props.molarEntropy()            == approx(S)     );
        CHECK( props.molarInternalEnergy()     == approx(U)     );
        CHECK( props.molarHelmholtzEnergy()    == approx(A)     );
        CHECK( props.molarHeatCapacityConstP() == approx(Cp)    );
        CHECK( props.molarHeatCapacityConstV() == approx(Cv)    );
        CHECK( props.molarDensity()            == approx(rho)   );
        CHECK( props.amount()                  == approx(nsum)  );
        CHECK( props.mass()                    == approx(mass)  );
        CHECK( props.gibbsEnergy()             == approx(Gtot)  );
        CHECK( props.enthalpy()                == approx(Htot)  );
        CHECK( props.volume()                  == approx(Vtot)  );
        CHECK( props.volumeT()                 == approx(VtotT) );
        CHECK( props.volumeP()                 == approx(VtotP) );
        CHECK( props.entropy()                 == approx(Stot)  );
        CHECK( props.internalEnergy()          == approx(Utot)  );
        CHECK( props.helmholtzEnergy()         == approx(Atot)  );

        //---------------------------------------------------------------------
        // Testing temperature derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXd n_T = ArrayXd{{ 0.0, 0.0, 0.0, 0.0 }};
        const ArrayXd x_T = ArrayXd{{ 0.0, 0.0, 0.0, 0.0 }};

        const ArrayXd  G0_T = 0.1 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd  H0_T = 0.2 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd  V0_T = 0.3 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd VT0_T = 0.4 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd VP0_T = 0.5 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd Cp0_T = 0.6 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
        const ArrayXd Cv0_T = 0.7 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * P;
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

        const double  G_T = ( G0_T * x).sum() + Gex_T;
        const double  H_T = ( H0_T * x).sum() + Hex_T;
        const double  V_T = ( V0_T * x).sum() + Vex_T;
        const double VT_T = (VT0_T * x).sum() + VexT_T;
        const double VP_T = (VP0_T * x).sum() + VexP_T;
        const double Cp_T = (Cp0_T * x).sum() + Cpex_T;
        const double Cv_T = (Cv0_T * x).sum() + Cvex_T;
        const double  S_T = (H_T - G_T)/T - (H - G)/(T*T);
        const double  U_T = H_T - P*V_T;
        const double  A_T = G_T - P*V_T;

        const double rho_T = -V_T/(V*V);

        const double Gtot_T  = nsum * G_T;
        const double Htot_T  = nsum * H_T;
        const double Vtot_T  = nsum * V_T;
        const double VtotT_T = nsum * VT_T;
        const double VtotP_T = nsum * VP_T;
        const double Stot_T  = nsum * S_T;
        const double Utot_T  = nsum * U_T;
        const double Atot_T  = nsum * A_T;

        autodiff::seed(T);
        props.update(T, P, n, extra);
        autodiff::unseed(T);

        CHECK( grad(props.temperature()) == 1.0 );
        CHECK( grad(props.pressure())    == 0.0 );

        CHECK( grad(props.speciesAmounts())               .isApprox(n_T)    );
        CHECK( grad(props.moleFractions())                .isApprox(x_T)    );
        CHECK( grad(props.lnActivityCoefficients())       .isApprox(ln_g_T) );
        CHECK( grad(props.lnActivities())                 .isApprox(ln_a_T) );
        CHECK( grad(props.chemicalPotentials())           .isApprox(u_T)    );
        CHECK( grad(props.standardGibbsEnergies())        .isApprox(G0_T)   );
        CHECK( grad(props.standardEnthalpies())           .isApprox(H0_T)   );
        CHECK( grad(props.standardVolumes())              .isApprox(V0_T)   );
        CHECK( grad(props.standardVolumesT())             .isApprox(VT0_T)  );
        CHECK( grad(props.standardVolumesP())             .isApprox(VP0_T)  );
        CHECK( grad(props.standardEntropies())            .isApprox(S0_T)   );
        CHECK( grad(props.standardInternalEnergies())     .isApprox(U0_T)   );
        CHECK( grad(props.standardHelmholtzEnergies())    .isApprox(A0_T)   );
        CHECK( grad(props.standardHeatCapacitiesConstP()) .isApprox(Cp0_T)  );
        CHECK( grad(props.standardHeatCapacitiesConstV()) .isApprox(Cv0_T)  );

        CHECK( grad(props.molarGibbsEnergy())        == approx(G_T)     );
        CHECK( grad(props.molarEnthalpy())           == approx(H_T)     );
        CHECK( grad(props.molarVolume())             == approx(V_T)     );
        CHECK( grad(props.molarEntropy())            == approx(S_T)     );
        CHECK( grad(props.molarInternalEnergy())     == approx(U_T)     );
        CHECK( grad(props.molarHelmholtzEnergy())    == approx(A_T)     );
        CHECK( grad(props.molarHeatCapacityConstP()) == approx(Cp_T)    );
        CHECK( grad(props.molarHeatCapacityConstV()) == approx(Cv_T)    );
        CHECK( grad(props.molarDensity())            == approx(rho_T)   );
        CHECK( grad(props.amount())                  == approx(0.0)     );
        CHECK( grad(props.mass())                    == approx(0.0)     );
        CHECK( grad(props.gibbsEnergy())             == approx(Gtot_T)  );
        CHECK( grad(props.enthalpy())                == approx(Htot_T)  );
        CHECK( grad(props.volume())                  == approx(Vtot_T)  );
        CHECK( grad(props.volumeT())                 == approx(VtotT_T) );
        CHECK( grad(props.volumeP())                 == approx(VtotP_T) );
        CHECK( grad(props.entropy())                 == approx(Stot_T)  );
        CHECK( grad(props.internalEnergy())          == approx(Utot_T)  );
        CHECK( grad(props.helmholtzEnergy())         == approx(Atot_T)  );

        //---------------------------------------------------------------------
        // Testing pressure derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXd n_P = ArrayXd{{ 0.0, 0.0, 0.0, 0.0 }};
        const ArrayXd x_P = ArrayXd{{ 0.0, 0.0, 0.0, 0.0 }};

        const ArrayXd  G0_P = 0.1 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd  H0_P = 0.2 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd  V0_P = 0.3 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd VT0_P = 0.4 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd VP0_P = 0.5 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd Cp0_P = 0.6 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
        const ArrayXd Cv0_P = 0.7 * ArrayXd{{ 10.0, 20.0, 30.0, 40.0 }} * T;
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
        const double VT_P = (VT0_P * x).sum() + VexT_P;
        const double VP_P = (VP0_P * x).sum() + VexP_P;
        const double Cp_P = (Cp0_P * x).sum() + Cpex_P;
        const double Cv_P = (Cv0_P * x).sum() + Cvex_P;
        const double  S_P = (H_P - G_P)/T;
        const double  U_P = H_P - V - P*V_P;
        const double  A_P = G_P - V - P*V_P;

        const double rho_P = -V_P/(V*V);

        const double  Gtot_P = nsum * G_P;
        const double  Htot_P = nsum * H_P;
        const double  Vtot_P = nsum * V_P;
        const double VtotT_P = nsum * VT_P;
        const double VtotP_P = nsum * VP_P;
        const double  Stot_P = nsum * S_P;
        const double  Utot_P = nsum * U_P;
        const double  Atot_P = nsum * A_P;

        autodiff::seed(P);
        props.update(T, P, n, extra);
        autodiff::unseed(P);

        CHECK( grad(props.temperature()) == 0.0 );
        CHECK( grad(props.pressure())    == 1.0 );

        CHECK( grad(props.speciesAmounts())               .isApprox(n_P)    );
        CHECK( grad(props.moleFractions())                .isApprox(x_P)    );
        CHECK( grad(props.lnActivityCoefficients())       .isApprox(ln_g_P) );
        CHECK( grad(props.lnActivities())                 .isApprox(ln_a_P) );
        CHECK( grad(props.chemicalPotentials())           .isApprox(u_P)    );
        CHECK( grad(props.standardGibbsEnergies())        .isApprox(G0_P)   );
        CHECK( grad(props.standardEnthalpies())           .isApprox(H0_P)   );
        CHECK( grad(props.standardVolumes())              .isApprox(V0_P)   );
        CHECK( grad(props.standardVolumesT())             .isApprox(VT0_P)  );
        CHECK( grad(props.standardVolumesP())             .isApprox(VP0_P)  );
        CHECK( grad(props.standardEntropies())            .isApprox(S0_P)   );
        CHECK( grad(props.standardInternalEnergies())     .isApprox(U0_P)   );
        CHECK( grad(props.standardHelmholtzEnergies())    .isApprox(A0_P)   );
        CHECK( grad(props.standardHeatCapacitiesConstP()) .isApprox(Cp0_P)  );
        CHECK( grad(props.standardHeatCapacitiesConstV()) .isApprox(Cv0_P)  );

        CHECK( grad(props.molarGibbsEnergy())        == approx(G_P)     );
        CHECK( grad(props.molarEnthalpy())           == approx(H_P)     );
        CHECK( grad(props.molarVolume())             == approx(V_P)     );
        CHECK( grad(props.molarEntropy())            == approx(S_P)     );
        CHECK( grad(props.molarInternalEnergy())     == approx(U_P)     );
        CHECK( grad(props.molarHelmholtzEnergy())    == approx(A_P)     );
        CHECK( grad(props.molarHeatCapacityConstP()) == approx(Cp_P)    );
        CHECK( grad(props.molarHeatCapacityConstV()) == approx(Cv_P)    );
        CHECK( grad(props.molarDensity())            == approx(rho_P)   );
        CHECK( grad(props.amount())                  == approx(0.0)     );
        CHECK( grad(props.mass())                    == approx(0.0)     );
        CHECK( grad(props.gibbsEnergy())             == approx(Gtot_P)  );
        CHECK( grad(props.enthalpy())                == approx(Htot_P)  );
        CHECK( grad(props.volume())                  == approx(Vtot_P)  );
        CHECK( grad(props.volumeT())                 == approx(VtotT_P) );
        CHECK( grad(props.volumeP())                 == approx(VtotP_P) );
        CHECK( grad(props.entropy())                 == approx(Stot_P)  );
        CHECK( grad(props.internalEnergy())          == approx(Utot_P)  );
        CHECK( grad(props.helmholtzEnergy())         == approx(Atot_P)  );

        //---------------------------------------------------------------------
        // Testing compositional derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXXd n_n = MatrixXd::Identity(4, 4);
        const ArrayXXd x_n = moleFractionsJacobian(n);

        const ArrayXXd  G0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd  H0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd  V0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd VT0_n = ArrayXXd::Zero(4, 4);
        const ArrayXXd VP0_n = ArrayXXd::Zero(4, 4);
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
            return ((A.matrix().transpose() * x.matrix()).array()).eval();
        };

        const ArrayXd  G_n = dot(x_n,  G0) + Gex_n;
        const ArrayXd  H_n = dot(x_n,  H0) + Hex_n;
        const ArrayXd  V_n = dot(x_n,  V0) + Vex_n;
        const ArrayXd VT_n = dot(x_n, VT0) + VexT_n;
        const ArrayXd VP_n = dot(x_n, VP0) + VexP_n;
        const ArrayXd Cp_n = dot(x_n, Cp0) + Cpex_n;
        const ArrayXd Cv_n = dot(x_n, Cv0) + Cvex_n;
        const ArrayXd  S_n = (H_n - G_n)/T;
        const ArrayXd  U_n = H_n - P*V_n;
        const ArrayXd  A_n = G_n - P*V_n;

        const ArrayXd nsum_n = ArrayXd::Ones(4);
        const ArrayXd mass_n = molar_masses;

        const ArrayXd rho_n = -V_n / (V*V);

        const ArrayXd  Gtot_n = nsum * G_n  + nsum_n * G;
        const ArrayXd  Htot_n = nsum * H_n  + nsum_n * H;
        const ArrayXd  Vtot_n = nsum * V_n  + nsum_n * V;
        const ArrayXd VtotT_n = nsum * VT_n + nsum_n * VT;
        const ArrayXd VtotP_n = nsum * VP_n + nsum_n * VP;
        const ArrayXd  Stot_n = nsum * S_n  + nsum_n * S;
        const ArrayXd  Utot_n = nsum * U_n  + nsum_n * U;
        const ArrayXd  Atot_n = nsum * A_n  + nsum_n * A;

        for(auto i = 0; i < 4; ++i)
        {
            INFO("i = " << i);
            autodiff::seed(n[i]);
            props.update(T, P, n, extra);
            autodiff::unseed(n[i]);

            CHECK( grad(props.temperature()) == 0.0 );
            CHECK( grad(props.pressure())    == 0.0 );

            CHECK( grad(props.speciesAmounts())               .isApprox(   n_n.col(i)) );
            CHECK( grad(props.moleFractions())                .isApprox(   x_n.col(i)) );
            CHECK( grad(props.lnActivityCoefficients())       .isApprox(ln_g_n.col(i)) );
            CHECK( grad(props.lnActivities())                 .isApprox(ln_a_n.col(i)) );
            CHECK( grad(props.chemicalPotentials())           .isApprox(   u_n.col(i)) );
            CHECK( grad(props.standardGibbsEnergies())        .isApprox(  G0_n.col(i)) );
            CHECK( grad(props.standardEnthalpies())           .isApprox(  H0_n.col(i)) );
            CHECK( grad(props.standardVolumes())              .isApprox(  V0_n.col(i)) );
            CHECK( grad(props.standardVolumesT())             .isApprox( VT0_n.col(i)) );
            CHECK( grad(props.standardVolumesP())             .isApprox( VP0_n.col(i)) );
            CHECK( grad(props.standardEntropies())            .isApprox(  S0_n.col(i)) );
            CHECK( grad(props.standardInternalEnergies())     .isApprox(  U0_n.col(i)) );
            CHECK( grad(props.standardHelmholtzEnergies())    .isApprox(  A0_n.col(i)) );
            CHECK( grad(props.standardHeatCapacitiesConstP()) .isApprox( Cp0_n.col(i)) );
            CHECK( grad(props.standardHeatCapacitiesConstV()) .isApprox( Cv0_n.col(i)) );

            CHECK( grad(props.molarGibbsEnergy())        == approx(   G_n[i])  );
            CHECK( grad(props.molarEnthalpy())           == approx(   H_n[i])  );
            CHECK( grad(props.molarVolume())             == approx(   V_n[i])  );
            CHECK( grad(props.molarEntropy())            == approx(   S_n[i])  );
            CHECK( grad(props.molarInternalEnergy())     == approx(   U_n[i])  );
            CHECK( grad(props.molarHelmholtzEnergy())    == approx(   A_n[i])  );
            CHECK( grad(props.molarHeatCapacityConstP()) == approx(  Cp_n[i])  );
            CHECK( grad(props.molarHeatCapacityConstV()) == approx(  Cv_n[i])  );
            CHECK( grad(props.molarDensity())            == approx( rho_n[i])  );
            CHECK( grad(props.amount())                  == approx(nsum_n[i])  );
            CHECK( grad(props.mass())                    == approx(mass_n[i])  );
            CHECK( grad(props.gibbsEnergy())             == approx(Gtot_n[i])  );
            CHECK( grad(props.enthalpy())                == approx(Htot_n[i])  );
            CHECK( grad(props.volume())                  == approx(Vtot_n[i])  );
            CHECK( grad(props.volumeT())                 == approx(VtotT_n[i]) );
            CHECK( grad(props.volumeP())                 == approx(VtotP_n[i]) );
            CHECK( grad(props.entropy())                 == approx(Stot_n[i])  );
            CHECK( grad(props.internalEnergy())          == approx(Utot_n[i])  );
            CHECK( grad(props.helmholtzEnergy())         == approx(Atot_n[i])  );
        }
    }

    SECTION("Testing when species have zero amounts")
    {
        const real T = 300.0;
        const real P = 123.0e5;

        const ArrayXr n = ArrayXr{{ 0.0, 0.0, 0.0, 0.0 }};

        Map<String, Any> extra;

        CHECK_THROWS( props.update(T, P, n, extra) );
    }
}
