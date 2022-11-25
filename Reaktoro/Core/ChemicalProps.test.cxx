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

    ActivityModel activity_model_gas = [](ActivityPropsRef props, ActivityModelArgs args)
    {
        const auto [T, P, x] = args;
        props.Vx  = 1.0 * (T*P)*(T*P);
        props.VxT = 2.0 * (T*P)*(T*P);
        props.VxP = 3.0 * (T*P)*(T*P);
        props.Gx  = 4.0 * (T*P)*(T*P);
        props.Hx  = 5.0 * (T*P)*(T*P);
        props.Cpx = 6.0 * (T*P)*(T*P);
        props.ln_g = 8.0 * x;
        props.ln_a = 9.0 * x;
        props.som = StateOfMatter::Gas;
    };

    ActivityModel activity_model_solid = [](ActivityPropsRef props, ActivityModelArgs args)
    {
        const auto [T, P, x] = args;
        props.Vx  = 1.1 * (T*P)*(T*P);
        props.VxT = 2.1 * (T*P)*(T*P);
        props.VxP = 3.1 * (T*P)*(T*P);
        props.Gx  = 4.1 * (T*P)*(T*P);
        props.Hx  = 5.1 * (T*P)*(T*P);
        props.Cpx = 6.1 * (T*P)*(T*P);
        props.ln_g = 8.1 * x;
        props.ln_a = 9.1 * x;
        props.som = StateOfMatter::Solid;
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
            .withIdealActivityModel(activity_model_gas)
            .withStateOfMatter(StateOfMatter::Gas)
            .withSpecies({
                db.species().get("H2O(g)"),
                db.species().get("CO2(g)")}),

        Phase()
            .withName("SomeSolid")
            .withActivityModel(activity_model_solid)
            .withIdealActivityModel(activity_model_solid)
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

    SECTION("Testing when species have non-zero amounts")
    {
        real T = 3.0;
        real P = 5.0;
        ArrayXr n = ArrayXr{{ 4.0, 6.0, 5.0 }};

        const ArrayXr x = ArrayXr{{ 0.4, 0.6, 1.0 }};

        const ArrayXr nsumphases = ArrayXr{{ 10.0, 5.0 }}; // the total amount of each phase

        const ArrayXr  G0 = ArrayXr{{ 0.1, 0.1, 1.1 }} * (T*P)*(T*P);
        const ArrayXr  H0 = ArrayXr{{ 0.2, 0.2, 1.2 }} * (T*P)*(T*P);
        const ArrayXr  V0 = ArrayXr{{ 0.3, 0.3, 1.3 }} * (T*P)*(T*P);
        const ArrayXr VT0 = ArrayXr{{ 0.4, 0.4, 1.4 }} * (T*P)*(T*P);
        const ArrayXr VP0 = ArrayXr{{ 0.5, 0.5, 1.5 }} * (T*P)*(T*P);
        const ArrayXr Cp0 = ArrayXr{{ 0.6, 0.6, 1.6 }} * (T*P)*(T*P);
        const ArrayXr Cv0 = Cp0 + T*VT0*VT0/VP0;
        const ArrayXr  S0 = (H0 - G0)/T;
        const ArrayXr  U0 = H0 - P*V0;
        const ArrayXr  A0 = G0 - P*V0;

        const ArrayXr Vx  = ArrayXr{{ 1.0, 1.1 }} * (T*P)*(T*P);
        const ArrayXr VxT = ArrayXr{{ 2.0, 2.1 }} * (T*P)*(T*P);
        const ArrayXr VxP = ArrayXr{{ 3.0, 3.1 }} * (T*P)*(T*P);
        const ArrayXr Gx  = ArrayXr{{ 4.0, 4.1 }} * (T*P)*(T*P);
        const ArrayXr Hx  = ArrayXr{{ 5.0, 5.1 }} * (T*P)*(T*P);
        const ArrayXr Cpx = ArrayXr{{ 6.0, 6.1 }} * (T*P)*(T*P);

        const ArrayXr ln_g = ArrayXr{{ 8.0*x[0], 8.0*x[1], 8.1*x[2] }};
        const ArrayXr ln_a = ArrayXr{{ 9.0*x[0], 9.0*x[1], 9.1*x[2] }};
        const ArrayXr u    = G0 + R*T*ln_a;

        const real Ntot  = n.sum();
        const real Mtot  = (n * molar_masses).sum();
        const real Vtot  = (V0 * n).sum() + (nsumphases * Vx).sum();
        const real VTtot = (VT0 * n).sum() + (nsumphases * VxT).sum();
        const real VPtot = (VP0 * n).sum() + (nsumphases * VxP).sum();
        const real Gtot  = (G0 * n).sum() + (nsumphases * Gx).sum();
        const real Htot  = (H0 * n).sum() + (nsumphases * Hx).sum();
        const real Cptot = (Cp0 * n).sum() + (nsumphases * Cpx).sum();
        const real Cvtot = Cptot + T*VTtot*VTtot/VPtot;
        const real Stot  = (Htot - Gtot)/T;
        const real Utot  = Htot - P*Vtot;
        const real Atot  = Gtot - P*Vtot;

        CHECK_NOTHROW( props.update(T, P, n) );

        CHECK( props.temperature() == T );
        CHECK( props.pressure()    == P );

        CHECK( props.speciesAmounts()                     .isApprox(n)    );
        CHECK( props.speciesMoleFractions()               .isApprox(x)    );
        CHECK( props.speciesActivityCoefficientsLn()      .isApprox(ln_g) );
        CHECK( props.speciesActivitiesLn()                .isApprox(ln_a) );
        CHECK( props.speciesChemicalPotentials()          .isApprox(u)    );
        CHECK( props.speciesStandardVolumes()             .isApprox(V0)   );
        CHECK( props.speciesStandardVolumesT()            .isApprox(VT0)  );
        CHECK( props.speciesStandardVolumesP()            .isApprox(VP0)  );
        CHECK( props.speciesStandardGibbsEnergies()       .isApprox(G0)   );
        CHECK( props.speciesStandardEnthalpies()          .isApprox(H0)   );
        CHECK( props.speciesStandardEntropies()           .isApprox(S0)   );
        CHECK( props.speciesStandardInternalEnergies()    .isApprox(U0)   );
        CHECK( props.speciesStandardHelmholtzEnergies()   .isApprox(A0)   );
        CHECK( props.speciesStandardHeatCapacitiesConstP().isApprox(Cp0)  );
        CHECK( props.speciesStandardHeatCapacitiesConstV().isApprox(Cv0)  );

        CHECK( props.molarVolume()                == Approx(Vtot/Ntot)  );
        CHECK( props.molarVolumeT()               == Approx(VTtot/Ntot) );
        CHECK( props.molarVolumeP()               == Approx(VPtot/Ntot) );
        CHECK( props.molarGibbsEnergy()           == Approx(Gtot/Ntot)  );
        CHECK( props.molarEnthalpy()              == Approx(Htot/Ntot)  );
        CHECK( props.molarEntropy()               == Approx(Stot/Ntot)  );
        CHECK( props.molarInternalEnergy()        == Approx(Utot/Ntot)  );
        CHECK( props.molarHelmholtzEnergy()       == Approx(Atot/Ntot)  );
        CHECK( props.molarHeatCapacityConstP()    == Approx(Cptot/Ntot) );
        CHECK( props.molarHeatCapacityConstV()    == Approx(Cvtot/Ntot) );
        CHECK( props.specificVolume()             == Approx(Vtot/Mtot)  );
        CHECK( props.specificVolumeT()            == Approx(VTtot/Mtot) );
        CHECK( props.specificVolumeP()            == Approx(VPtot/Mtot) );
        CHECK( props.specificGibbsEnergy()        == Approx(Gtot/Mtot)  );
        CHECK( props.specificEnthalpy()           == Approx(Htot/Mtot)  );
        CHECK( props.specificEntropy()            == Approx(Stot/Mtot)  );
        CHECK( props.specificInternalEnergy()     == Approx(Utot/Mtot)  );
        CHECK( props.specificHelmholtzEnergy()    == Approx(Atot/Mtot)  );
        CHECK( props.specificHeatCapacityConstP() == Approx(Cptot/Mtot) );
        CHECK( props.specificHeatCapacityConstV() == Approx(Cvtot/Mtot) );
        CHECK( props.density()                    == Approx(Mtot/Vtot)  );

        CHECK( props.amount()             == Approx(Ntot)  );
        CHECK( props.mass()               == Approx(Mtot)  );
        CHECK( props.volume()             == Approx(Vtot)  );
        CHECK( props.gibbsEnergy()        == Approx(Gtot)  );
        CHECK( props.enthalpy()           == Approx(Htot)  );
        CHECK( props.entropy()            == Approx(Stot)  );
        CHECK( props.internalEnergy()     == Approx(Utot)  );
        CHECK( props.helmholtzEnergy()    == Approx(Atot)  );
        CHECK( props.heatCapacityConstP() == Approx(Cptot) );
        CHECK( props.heatCapacityConstV() == Approx(Cvtot) );

        //---------------------------------------------------------------------
        // Testing temperature derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXd n_T = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd x_T = ArrayXd{{ 0.0, 0.0, 0.0 }};

        const ArrayXd nsumphases_T = ArrayXd{{ 0.0, 0.0 }};

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

        const ArrayXd  Vex_T = ArrayXd{{ 1.0, 1.1 }} * 2*P*(T*P);
        const ArrayXd VexT_T = ArrayXd{{ 2.0, 2.1 }} * 2*P*(T*P);
        const ArrayXd VexP_T = ArrayXd{{ 3.0, 3.1 }} * 2*P*(T*P);
        const ArrayXd  Gex_T = ArrayXd{{ 4.0, 4.1 }} * 2*P*(T*P);
        const ArrayXd  Hex_T = ArrayXd{{ 5.0, 5.1 }} * 2*P*(T*P);
        const ArrayXd Cpex_T = ArrayXd{{ 6.0, 6.1 }} * 2*P*(T*P);

        const ArrayXd ln_g_T = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd ln_a_T = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd    u_T = G0_T + R*ln_a;

        const double Ntot_T  = 0.0;
        const double Mtot_T  = 0.0;
        const double Gtot_T  = (G0_T * n).sum() + (nsumphases * Gex_T).sum();
        const double Htot_T  = (H0_T * n).sum() + (nsumphases * Hex_T).sum();
        const double Vtot_T  = (V0_T * n).sum() + (nsumphases * Vex_T).sum();
        const double VTtot_T = (VT0_T * n).sum() + (nsumphases * VexT_T).sum();
        const double VPtot_T = (VP0_T * n).sum() + (nsumphases * VexP_T).sum();
        const double Cptot_T = (Cp0_T * n).sum() + (nsumphases * Cpex_T).sum();
        const double Cvtot_T = Cptot_T + VTtot*VTtot/VPtot + 2*T*VTtot*VTtot_T/VPtot - T*VTtot*VTtot/VPtot/VPtot*VPtot_T;
        const double Stot_T  = (Htot_T - Gtot_T)/T - (Htot - Gtot)/(T*T);
        const double Utot_T  = Htot_T - P*Vtot_T;
        const double Atot_T  = Gtot_T - P*Vtot_T;

        autodiff::seed(T);
        props.update(T, P, n);
        autodiff::unseed(T);

        CHECK( grad(props.temperature()) == 1.0 );
        CHECK( grad(props.pressure())    == 0.0 );

        CHECK( grad(props.speciesAmounts())                      .isApprox(n_T)    );
        CHECK( grad(props.speciesMoleFractions())                .isApprox(x_T)    );
        CHECK( grad(props.speciesActivityCoefficientsLn())       .isApprox(ln_g_T) );
        CHECK( grad(props.speciesActivitiesLn())                 .isApprox(ln_a_T) );
        CHECK( grad(props.speciesChemicalPotentials())           .isApprox(u_T)    );
        CHECK( grad(props.speciesStandardVolumes())              .isApprox(V0_T)   );
        CHECK( grad(props.speciesStandardVolumesT())             .isApprox(VT0_T)  );
        CHECK( grad(props.speciesStandardVolumesP())             .isApprox(VP0_T)  );
        CHECK( grad(props.speciesStandardGibbsEnergies())        .isApprox(G0_T)   );
        CHECK( grad(props.speciesStandardEnthalpies())           .isApprox(H0_T)   );
        CHECK( grad(props.speciesStandardEntropies())            .isApprox(S0_T)   );
        CHECK( grad(props.speciesStandardInternalEnergies())     .isApprox(U0_T)   );
        CHECK( grad(props.speciesStandardHelmholtzEnergies())    .isApprox(A0_T)   );
        CHECK( grad(props.speciesStandardHeatCapacitiesConstP()) .isApprox(Cp0_T)  );
        CHECK( grad(props.speciesStandardHeatCapacitiesConstV()) .isApprox(Cv0_T)  );

        CHECK( grad(props.amount())             == Approx(Ntot_T)  );
        CHECK( grad(props.mass())               == Approx(Mtot_T)  );
        CHECK( grad(props.volume())             == Approx(Vtot_T)  );
        CHECK( grad(props.volumeT())            == Approx(VTtot_T) );
        CHECK( grad(props.volumeP())            == Approx(VPtot_T) );
        CHECK( grad(props.gibbsEnergy())        == Approx(Gtot_T)  );
        CHECK( grad(props.enthalpy())           == Approx(Htot_T)  );
        CHECK( grad(props.entropy())            == Approx(Stot_T)  );
        CHECK( grad(props.internalEnergy())     == Approx(Utot_T)  );
        CHECK( grad(props.helmholtzEnergy())    == Approx(Atot_T)  );
        CHECK( grad(props.heatCapacityConstP()) == Approx(Cptot_T) );
        CHECK( grad(props.heatCapacityConstV()) == Approx(Cvtot_T) );

        //---------------------------------------------------------------------
        // Testing pressure derivatives of the properties
        //---------------------------------------------------------------------
        const ArrayXd n_P = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd x_P = ArrayXd{{ 0.0, 0.0, 0.0 }};

        const ArrayXd nsumphases_P = ArrayXd{{ 0.0, 0.0 }};

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

        const ArrayXd  Vex_P = ArrayXd{{ 1.0, 1.1 }} * 2*T*(T*P);
        const ArrayXd VexT_P = ArrayXd{{ 2.0, 2.1 }} * 2*T*(T*P);
        const ArrayXd VexP_P = ArrayXd{{ 3.0, 3.1 }} * 2*T*(T*P);
        const ArrayXd  Gex_P = ArrayXd{{ 4.0, 4.1 }} * 2*T*(T*P);
        const ArrayXd  Hex_P = ArrayXd{{ 5.0, 5.1 }} * 2*T*(T*P);
        const ArrayXd Cpex_P = ArrayXd{{ 6.0, 6.1 }} * 2*T*(T*P);

        const ArrayXd ln_g_P = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd ln_a_P = ArrayXd{{ 0.0, 0.0, 0.0 }};
        const ArrayXd    u_P = G0_P;

        const double Ntot_P  = 0.0;
        const double Mtot_P  = 0.0;
        const double Gtot_P  = (G0_P * n).sum() + (nsumphases * Gex_P).sum();
        const double Htot_P  = (H0_P * n).sum() + (nsumphases * Hex_P).sum();
        const double Vtot_P  = (V0_P * n).sum() + (nsumphases * Vex_P).sum();
        const double VTtot_P = (VT0_P * n).sum() + (nsumphases * VexT_P).sum();
        const double VPtot_P = (VP0_P * n).sum() + (nsumphases * VexP_P).sum();
        const double Cptot_P = (Cp0_P * n).sum() + (nsumphases * Cpex_P).sum();
        const double Cvtot_P = Cptot_P + 2*T*VTtot*VTtot_P/VPtot - T*VTtot*VTtot/VPtot/VPtot*VPtot_P;
        const double Stot_P  = (Htot_P - Gtot_P)/T;
        const double Utot_P  = Htot_P - Vtot - P*Vtot_P;
        const double Atot_P  = Gtot_P - Vtot - P*Vtot_P;

        autodiff::seed(P);
        props.update(T, P, n);
        autodiff::unseed(P);

        CHECK( grad(props.temperature()) == 0.0 );
        CHECK( grad(props.pressure())    == 1.0 );

        CHECK( grad(props.speciesAmounts())                      .isApprox(n_P)    );
        CHECK( grad(props.speciesMoleFractions())                .isApprox(x_P)    );
        CHECK( grad(props.speciesActivityCoefficientsLn())       .isApprox(ln_g_P) );
        CHECK( grad(props.speciesActivitiesLn())                 .isApprox(ln_a_P) );
        CHECK( grad(props.speciesStandardVolumes())              .isApprox(V0_P)   );
        CHECK( grad(props.speciesStandardVolumesT())             .isApprox(VT0_P)  );
        CHECK( grad(props.speciesStandardVolumesP())             .isApprox(VP0_P)  );
        CHECK( grad(props.speciesChemicalPotentials())           .isApprox(u_P)    );
        CHECK( grad(props.speciesStandardGibbsEnergies())        .isApprox(G0_P)   );
        CHECK( grad(props.speciesStandardEnthalpies())           .isApprox(H0_P)   );
        CHECK( grad(props.speciesStandardEntropies())            .isApprox(S0_P)   );
        CHECK( grad(props.speciesStandardInternalEnergies())     .isApprox(U0_P)   );
        CHECK( grad(props.speciesStandardHelmholtzEnergies())    .isApprox(A0_P)   );
        CHECK( grad(props.speciesStandardHeatCapacitiesConstP()) .isApprox(Cp0_P)  );
        CHECK( grad(props.speciesStandardHeatCapacitiesConstV()) .isApprox(Cv0_P)  );

        CHECK( grad(props.amount())             == Approx(Ntot_P)  );
        CHECK( grad(props.mass())               == Approx(Mtot_P)  );
        CHECK( grad(props.volume())             == Approx(Vtot_P)  );
        CHECK( grad(props.volumeT())            == Approx(VTtot_P) );
        CHECK( grad(props.volumeP())            == Approx(VPtot_P) );
        CHECK( grad(props.gibbsEnergy())        == Approx(Gtot_P)  );
        CHECK( grad(props.enthalpy())           == Approx(Htot_P)  );
        CHECK( grad(props.entropy())            == Approx(Stot_P)  );
        CHECK( grad(props.internalEnergy())     == Approx(Utot_P)  );
        CHECK( grad(props.helmholtzEnergy())    == Approx(Atot_P)  );
        CHECK( grad(props.heatCapacityConstP()) == Approx(Cptot_P) );
        CHECK( grad(props.heatCapacityConstV()) == Approx(Cvtot_P) );

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
        const ArrayXXd VT0_n = ArrayXXd::Zero(3, 3);
        const ArrayXXd VP0_n = ArrayXXd::Zero(3, 3);
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

        const ArrayXd Ntot_n  = ArrayXd::Ones(3);
        const ArrayXd Mtot_n  = molar_masses;
        const ArrayXd Gtot_n  = dot(n_n, G0) + dot(nsumphases_n, Gx);
        const ArrayXd Htot_n  = dot(n_n, H0) + dot(nsumphases_n, Hx);
        const ArrayXd Vtot_n  = dot(n_n, V0) + dot(nsumphases_n, Vx);
        const ArrayXd VTtot_n = dot(n_n, VT0) + dot(nsumphases_n, VxT);
        const ArrayXd VPtot_n = dot(n_n, VP0) + dot(nsumphases_n, VxP);
        const ArrayXd Cptot_n = dot(n_n, Cp0) + dot(nsumphases_n, Cpx);
        const ArrayXd Cvtot_n = Cptot_n + 2*T*VTtot*VTtot_n/VPtot - T*VTtot*VTtot/VPtot/VPtot * VPtot_n;
        const ArrayXd Stot_n  = (Htot_n - Gtot_n)/T;
        const ArrayXd Utot_n  = Htot_n - P*Vtot_n;
        const ArrayXd Atot_n  = Gtot_n - P*Vtot_n;

        for(auto i = 0; i < 3; ++i)
        {
            autodiff::seed(n[i]);
            props.update(T, P, n);
            autodiff::unseed(n[i]);

            CHECK( grad(props.temperature()) == 0.0 );
            CHECK( grad(props.pressure())    == 0.0 );

            CHECK( grad(props.speciesAmounts())                      .isApprox(    n_n.col(i)) );
            CHECK( grad(props.speciesMoleFractions())                .isApprox(    x_n.col(i)) );
            CHECK( grad(props.speciesActivityCoefficientsLn())       .isApprox( ln_g_n.col(i)) );
            CHECK( grad(props.speciesActivitiesLn())                 .isApprox( ln_a_n.col(i)) );
            CHECK( grad(props.speciesChemicalPotentials())           .isApprox(    u_n.col(i)) );
            CHECK( grad(props.speciesStandardVolumes())              .isApprox(   V0_n.col(i)) );
            CHECK( grad(props.speciesStandardVolumesT())             .isApprox(  VT0_n.col(i)) );
            CHECK( grad(props.speciesStandardVolumesP())             .isApprox(  VP0_n.col(i)) );
            CHECK( grad(props.speciesStandardGibbsEnergies())        .isApprox(   G0_n.col(i)) );
            CHECK( grad(props.speciesStandardEnthalpies())           .isApprox(   H0_n.col(i)) );
            CHECK( grad(props.speciesStandardEntropies())            .isApprox(   S0_n.col(i)) );
            CHECK( grad(props.speciesStandardInternalEnergies())     .isApprox(   U0_n.col(i)) );
            CHECK( grad(props.speciesStandardHelmholtzEnergies())    .isApprox(   A0_n.col(i)) );
            CHECK( grad(props.speciesStandardHeatCapacitiesConstP()) .isApprox(  Cp0_n.col(i)) );
            CHECK( grad(props.speciesStandardHeatCapacitiesConstV()) .isApprox(  Cv0_n.col(i)) );

            CHECK( grad(props.amount())             == Approx(Ntot_n[i])  );
            CHECK( grad(props.mass())               == Approx(Mtot_n[i])  );
            CHECK( grad(props.gibbsEnergy())        == Approx(Gtot_n[i])  );
            CHECK( grad(props.enthalpy())           == Approx(Htot_n[i])  );
            CHECK( grad(props.volume())             == Approx(Vtot_n[i])  );
            CHECK( grad(props.volumeT())            == Approx(VTtot_n[i]) );
            CHECK( grad(props.volumeP())            == Approx(VPtot_n[i]) );
            CHECK( grad(props.entropy())            == Approx(Stot_n[i])  );
            CHECK( grad(props.internalEnergy())     == Approx(Utot_n[i])  );
            CHECK( grad(props.helmholtzEnergy())    == Approx(Atot_n[i])  );
            CHECK( grad(props.heatCapacityConstP()) == Approx(Cptot_n[i]) );
            CHECK( grad(props.heatCapacityConstV()) == Approx(Cvtot_n[i]) );
        }
    }

    SECTION("Testing when species have zero mole fractions")
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

    SECTION("Testing convenience methods")
    {
        const real T = 11.0;
        const real P = 13.0;
        const ArrayXr n = ArrayXr{{ 1.0, 2.0, 3.0 }};

        props.update(T, P, n);

        const auto A  = system.formulaMatrix();
        const auto Ae = system.formulaMatrixElements();
        const auto Az = system.formulaMatrixCharge();

        const auto idxElement = [&](auto symbol) { return system.elements().index(symbol); };
        const auto idxSpecies = [&](auto name) { return system.species().index(name); };

        const VectorXr b = A * n.matrix();                       // element and charge amounts in the system
        const VectorXr b0 = A.leftCols(2) * n.head(2).matrix();  // element and charge amounts in phase 0: H2O(g) CO2(g)
        const VectorXr b1 = A.rightCols(1) * n.tail(1).matrix(); // element and charge amounts in phase 1: CaCO3(s)

        const VectorXr be = Ae * n.matrix();                       // element amounts in the system
        const VectorXr be0 = Ae.leftCols(2) * n.head(2).matrix();  // element amounts in phase 0: H2O(g) CO2(g)
        const VectorXr be1 = Ae.rightCols(1) * n.tail(1).matrix(); // element amounts in phase 1: CaCO3(s)

        const real z  = dot(Az.row(0), n.matrix());                       // charge in the system
        const real z0 = dot(Az.row(0).leftCols(2), n.head(2).matrix());   // charge in phase 0: H2O(g) CO2(g)
        const real z1 = dot(Az.row(0).rightCols(1), n.tail(1).matrix());  // charge in phase 1: CaCO3(s)

        CHECK( props.charge() == Approx(z) );
        CHECK( props.chargeInPhase(0) == Approx(z0) );
        CHECK( props.chargeInPhase(1) == Approx(z1) );
        CHECK( props.chargeInPhase("SomeGas") == Approx(z0) );
        CHECK( props.chargeInPhase("SomeSolid") == Approx(z1) );
        CHECK( props.chargeAmongSpecies({0, 1}) == Approx(z0) );
        CHECK( props.chargeAmongSpecies({2}) == Approx(z1) );

        CHECK( props.elementAmounts().matrix().isApprox(be) );
        CHECK( props.elementAmountsInPhase(0).matrix().isApprox(be0) );
        CHECK( props.elementAmountsInPhase(1).matrix().isApprox(be1) );
        CHECK( props.elementAmountsInPhase("SomeGas").matrix().isApprox(be0) );
        CHECK( props.elementAmountsInPhase("SomeSolid").matrix().isApprox(be1) );
        CHECK( props.elementAmountsAmongSpecies({0, 1}).matrix().isApprox(be0) );
        CHECK( props.elementAmountsAmongSpecies({2}).matrix().isApprox(be1) );

        CHECK( props.componentAmounts().matrix().isApprox(b) );
        CHECK( props.componentAmountsInPhase(0).matrix().isApprox(b0) );
        CHECK( props.componentAmountsInPhase(1).matrix().isApprox(b1) );
        CHECK( props.componentAmountsInPhase("SomeGas").matrix().isApprox(b0) );
        CHECK( props.componentAmountsInPhase("SomeSolid").matrix().isApprox(b1) );
        CHECK( props.componentAmountsAmongSpecies({0, 1}).matrix().isApprox(b0) );
        CHECK( props.componentAmountsAmongSpecies({2}).matrix().isApprox(b1) );

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
            CHECK( props.speciesAmount(name)                     == Approx(props.speciesAmounts()[idx])                      );
            CHECK( props.speciesMass(name)                       == Approx(props.speciesMasses()[idx])                       );
            CHECK( props.speciesMoleFraction(name)               == Approx(props.speciesMoleFractions()[idx])                );
            CHECK( props.speciesActivityCoefficient(name)        == Approx(exp(props.speciesActivityCoefficientsLn()[idx]))  );
            CHECK( props.speciesActivityCoefficientLg(name)      == Approx(props.speciesActivityCoefficientsLn()[idx]/ln10)  );
            CHECK( props.speciesActivityCoefficientLn(name)      == Approx(props.speciesActivityCoefficientsLn()[idx])       );
            CHECK( props.speciesActivity(name)                   == Approx(exp(props.speciesActivitiesLn()[idx]))            );
            CHECK( props.speciesActivityLg(name)                 == Approx(props.speciesActivitiesLn()[idx]/ln10)            );
            CHECK( props.speciesActivityLn(name)                 == Approx(props.speciesActivitiesLn()[idx])                 );
            CHECK( props.speciesChemicalPotential(name)          == Approx(props.speciesChemicalPotentials()[idx])           );
            CHECK( props.speciesStandardVolume(name)             == Approx(props.speciesStandardVolumes()[idx])              );
            CHECK( props.speciesStandardVolumeT(name)            == Approx(props.speciesStandardVolumesT()[idx])             );
            CHECK( props.speciesStandardVolumeP(name)            == Approx(props.speciesStandardVolumesP()[idx])             );
            CHECK( props.speciesStandardGibbsEnergy(name)        == Approx(props.speciesStandardGibbsEnergies()[idx])        );
            CHECK( props.speciesStandardEnthalpy(name)           == Approx(props.speciesStandardEnthalpies()[idx])           );
            CHECK( props.speciesStandardEntropy(name)            == Approx(props.speciesStandardEntropies()[idx])            );
            CHECK( props.speciesStandardInternalEnergy(name)     == Approx(props.speciesStandardInternalEnergies()[idx])     );
            CHECK( props.speciesStandardHelmholtzEnergy(name)    == Approx(props.speciesStandardHelmholtzEnergies()[idx])    );
            CHECK( props.speciesStandardHeatCapacityConstP(name) == Approx(props.speciesStandardHeatCapacitiesConstP()[idx]) );
            CHECK( props.speciesStandardHeatCapacityConstV(name) == Approx(props.speciesStandardHeatCapacitiesConstV()[idx]) );

            CHECK( props.speciesAmount(i)                     == Approx(props.speciesAmounts()[idx])                         );
            CHECK( props.speciesMass(i)                       == Approx(props.speciesMasses()[idx])                          );
            CHECK( props.speciesMoleFraction(i)               == Approx(props.speciesMoleFractions()[idx])                   );
            CHECK( props.speciesActivityCoefficient(i)        == Approx(exp(props.speciesActivityCoefficientsLn()[idx]))     );
            CHECK( props.speciesActivityCoefficientLg(i)      == Approx(props.speciesActivityCoefficientsLn()[idx]/ln10)     );
            CHECK( props.speciesActivityCoefficientLn(i)      == Approx(props.speciesActivityCoefficientsLn()[idx])          );
            CHECK( props.speciesActivity(i)                   == Approx(exp(props.speciesActivitiesLn()[idx]))               );
            CHECK( props.speciesActivityLg(i)                 == Approx(props.speciesActivitiesLn()[idx]/ln10)               );
            CHECK( props.speciesActivityLn(i)                 == Approx(props.speciesActivitiesLn()[idx])                    );
            CHECK( props.speciesChemicalPotential(i)          == Approx(props.speciesChemicalPotentials()[idx])              );
            CHECK( props.speciesStandardVolume(i)             == Approx(props.speciesStandardVolumes()[idx])                 );
            CHECK( props.speciesStandardVolumeT(i)            == Approx(props.speciesStandardVolumesT()[idx])                );
            CHECK( props.speciesStandardVolumeP(i)            == Approx(props.speciesStandardVolumesP()[idx])                );
            CHECK( props.speciesStandardGibbsEnergy(i)        == Approx(props.speciesStandardGibbsEnergies()[idx])           );
            CHECK( props.speciesStandardEnthalpy(i)           == Approx(props.speciesStandardEnthalpies()[idx])              );
            CHECK( props.speciesStandardEntropy(i)            == Approx(props.speciesStandardEntropies()[idx])               );
            CHECK( props.speciesStandardInternalEnergy(i)     == Approx(props.speciesStandardInternalEnergies()[idx])        );
            CHECK( props.speciesStandardHelmholtzEnergy(i)    == Approx(props.speciesStandardHelmholtzEnergies()[idx])       );
            CHECK( props.speciesStandardHeatCapacityConstP(i) == Approx(props.speciesStandardHeatCapacitiesConstP()[idx])    );
            CHECK( props.speciesStandardHeatCapacityConstV(i) == Approx(props.speciesStandardHeatCapacitiesConstV()[idx])    );

            i += 1;
        }
    }

    SECTION("Testing indices methods")
    {
        real T = 3.0;
        real P = 5.0;
        ArrayXr n = ArrayXr{{ 4.0, 6.0, 5.0 }};

        props.update(T, P, n);

        CHECK( props.indicesPhasesWithState(StateOfMatter::Gas) == Indices{0} );
        CHECK( props.indicesPhasesWithState(StateOfMatter::Solid) == Indices{1} );
        CHECK( props.indicesPhasesWithStates({StateOfMatter::Solid, StateOfMatter::Gas}) == Indices{0, 1} );
        CHECK( props.indicesPhasesWithFluidState() == Indices{0} );
        CHECK( props.indicesPhasesWithSolidState() == Indices{1} );
    }

    SECTION("Testing correct increment of state id as a ChemicalProps object is updated")
    {
        ChemicalState state(system);
        state.temperature(345.6, "K");
        state.pressure(1.234, "bar");
        state.setSpeciesAmounts(0.1234);

        auto T = state.temperature();
        auto P = state.pressure();
        auto n = state.speciesAmounts();

        // Checking stateid with default constructor
        ChemicalProps props;
        CHECK(props.stateid() == 0);

        // Checking stateid with ChemicalProps(system) constructor
        props = ChemicalProps(system);
        CHECK(props.stateid() == 0);

        // Checking stateid with ChemicalProps(state) constructor
        props = ChemicalProps(state);
        CHECK(props.stateid() == 1);

        // Checking stateid with ChemicalProps::update(state) method
        props.update(state);
        CHECK(props.stateid() == 2);

        // Checking stateid with ChemicalProps::update(T, P, n) method
        props.update(T, P, n);
        CHECK(props.stateid() == 3);

        // Checking stateid with ChemicalProps::update(VectorXr) method
        props.update(VectorXr(props));
        CHECK(props.stateid() == 4);

        // Checking stateid with ChemicalProps::update(VectorXd) method
        props.update(VectorXd(props));
        CHECK(props.stateid() == 5);

        // Checking stateid with ChemicalProps::updateIdeal(state) method
        props.updateIdeal(state);
        CHECK(props.stateid() == 6);

        // Checking stateid with ChemicalProps::updateIdeal(T, P, n) method
        props.updateIdeal(T, P, n);
        CHECK(props.stateid() == 7);

        // Checking stateid with ChemicalProps::deserialize(ArrayStream<real> const&) method
        ArrayStream<real> rstream;
        props.serialize(rstream);
        props.deserialize(rstream);
        CHECK(props.stateid() == 8);

        // Checking stateid with ChemicalProps::deserialize(ArrayStream<double> const&) method
        ArrayStream<double> dstream;
        props.serialize(dstream);
        props.deserialize(dstream);
        CHECK(props.stateid() == 9);
    }
}
