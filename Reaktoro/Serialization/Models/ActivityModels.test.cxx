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
using namespace Catch;

// Reaktoro includes
#include <Reaktoro/Models/ActivityModels/ActivityModelExtendedUNIQUAC.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelPitzer.hpp>
#include <Reaktoro/Serialization/Core.hpp>
#include <Reaktoro/Serialization/Models/ActivityModels.hpp>
using namespace Reaktoro;

TEST_CASE("Testing serialization of ActivityModelParamsPitzer", "[Serialization][Models][ActivityModels][ActivityModelParamsPitzer]")
{
    ActivityModelParamsPitzer params;

    String yml = R"#(
        Beta0:
          - { Formulas: [ Ba+2, OH- ],        CorrectionModel: Phreeqc,  Parameters: [ 0.17175 ] }
          - { Formulas: [ Ba+2, Cl- ],        CorrectionModel: Phreeqc,  Parameters: [ 0.5268, 0.0, 0.0, 0.0, 0.0, 4.75e+4 ] }
          - { Formulas: [ Mg+2, SO4-2 ],      CorrectionModel: Phreeqc,  Parameters: [ 0.2135, -951, 0.0, -2.34e-2, 2.28e-5 ] }
          - { Formulas: [ Cl-, K+ ],          CorrectionModel: Phreeqc,  Parameters: [ 0.04808, -758.48, -4.7062, 0.010072, -3.7599e-6 ] }
          - { Formulas: [ Cl-, Na+ ],         CorrectionModel: Phreeqc,  Parameters: [ 7.534e-2, 9598.4, 35.48, -5.8731e-2, 1.798e-5, -5.0e+5 ] }
        Beta1:
          - { Formulas: [ Ba+2, OH- ],        CorrectionModel: Constant, Parameters: [ 1.2 ] }
          - { Formulas: [ Br-, H+ ],          CorrectionModel: Phreeqc,  Parameters: [ 0.3564, 0.0, 0.0, 4.467e-4 ] }
          - { Formulas: [ Cl-, Na+ ],         CorrectionModel: Phreeqc,  Parameters: [ 0.2769, 1.377e+4, 46.8, -6.9512e-2, 2.0e-5, -7.4823e+5 ] }
        Beta2:
          - { Formulas: [ Ca+2, OH- ],        CorrectionModel: Phreeqc,  Parameters: [ -5.72 ] }
          - { Formulas: [ Ca+2, Cl- ],        CorrectionModel: Phreeqc,  Parameters: [ -1.13, 0.0, 0.0, -0.0476 ] }
          - { Formulas: [ Mg+2, SO4-2 ],      CorrectionModel: Phreeqc,  Parameters: [ -32.45, 0.0, -3.236e+3, 21.812, -1.8859e-2 ] }
        Cphi:
          - { Formulas: [ B(OH)4-, Na+ ],     CorrectionModel: Constant, Parameters: [ 0.0114 ] }
          - { Formulas: [ Ba+2, Cl- ],        CorrectionModel: Phreeqc,  Parameters: [ -0.143, -114.5 ] }
          - { Formulas: [ Cl-, Na+ ],         CorrectionModel: Phreeqc,  Parameters: [ 1.48e-3, -120.5, -0.2081, 0.0, 1.166e-7, 11121 ] }
        Theta:
          - { Formulas: [ K+, Na+ ],          CorrectionModel: Constant, Parameters: [ -0.012 ] }
          - { Formulas: [ Ca+2, K+ ],         CorrectionModel: Phreeqc,  Parameters: [ -5.35e-3, 0.0, 0.0, 3.08e-4 ] }
          - { Formulas: [ Ca+2, Na+ ],        CorrectionModel: Phreeqc,  Parameters: [ 9.22e-2, 0.0, 0.0, -4.29e-4, 1.21e-6 ] }
        Lambda:
          - { Formulas: [ Cl-, CO2 ],         CorrectionModel: Constant, Parameters: [ -0.005 ] }
          - { Formulas: [ CO2, CO2 ],         CorrectionModel: Phreeqc,  Parameters: [ -1.34e-2, 348, 0.803 ] }
          - { Formulas: [ H4SiO4, Mg+2 ],     CorrectionModel: Phreeqc,  Parameters: [ 0.238, -1788, -9.023, 0.0103 ] }
        Zeta:
          - { Formulas: [ Cl-, H4SiO4, K+ ],  CorrectionModel: Constant, Parameters: [ -0.0153 ] }
          - { Formulas: [ CO2, Na+, SO4-2 ],  CorrectionModel: Constant, Parameters: [ -0.015 ] }
        Psi:
          - { Formulas: [ Br-, K+, Na+ ],     CorrectionModel: Constant, Parameters: [ -0.0022 ] }
          - { Formulas: [ Cl-, Mg+2, SO4-2 ], CorrectionModel: Phreeqc,  Parameters: [ -0.008, 32.63 ] }
          - { Formulas: [ Ca+2, Cl-, Na+ ],   CorrectionModel: Phreeqc,  Parameters: [ -1.48e-2, 0.0, 0.0, -5.2e-6 ] }
        Mu:
          - { Formulas: [ NH3, NH3, CO3-2 ],  CorrectionModel: Constant, Parameters: [ 0.1234 ] }
        Eta:
          - { Formulas: [ CO2, Na+, K+ ],     CorrectionModel: Constant, Parameters: [ -0.0022 ] }
        Alpha1:
          - { Formulas: [ Na+, Cl- ], Parameters: [ 2.0 ] }
        Alpha2:
          - { Formulas: [ Na+, Cl- ], Parameters: [ 12.0 ] }
    )#";

    params = Data::parse(yml).as<ActivityModelParamsPitzer>();

    CHECK( params.beta0.size()  == 5 );
    CHECK( params.beta1.size()  == 3 );
    CHECK( params.beta2.size()  == 3 );
    CHECK( params.Cphi.size()   == 3 );
    CHECK( params.theta.size()  == 3 );
    CHECK( params.lambda.size() == 3 );
    CHECK( params.zeta.size()   == 2 );
    CHECK( params.psi.size()    == 3 );
    CHECK( params.mu.size()     == 1 );
    CHECK( params.eta.size()    == 1 );
    CHECK( params.alpha1.size() == 1 );
    CHECK( params.alpha2.size() == 1 );

    CHECK( params.beta0[0].formulas == Vec<ChemicalFormula>{"Ba+2", "OH-"} );
    CHECK( params.beta0[0].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.beta0[0].parameters == Vec<real>{ 0.17175 } );

    CHECK( params.beta0[1].formulas == Vec<ChemicalFormula>{"Ba+2", "Cl-"} );
    CHECK( params.beta0[1].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.beta0[1].parameters == Vec<real>{ 0.5268, 0.0, 0.0, 0.0, 0.0, 4.75e+4 } );

    CHECK( params.beta0[2].formulas == Vec<ChemicalFormula>{"Mg+2", "SO4-2"} );
    CHECK( params.beta0[2].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.beta0[2].parameters == Vec<real>{ 0.2135, -951, 0.0, -2.34e-2, 2.28e-5 } );

    CHECK( params.beta0[3].formulas == Vec<ChemicalFormula>{"Cl-", "K+"} );
    CHECK( params.beta0[3].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.beta0[3].parameters == Vec<real>{ 0.04808, -758.48, -4.7062, 0.010072, -3.7599e-6 } );

    CHECK( params.beta0[4].formulas == Vec<ChemicalFormula>{"Cl-", "Na+"} );
    CHECK( params.beta0[4].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.beta0[4].parameters == Vec<real>{ 7.534e-2, 9598.4, 35.48, -5.8731e-2, 1.798e-5, -5.0e+5 } );

    CHECK( params.beta1[0].formulas == Vec<ChemicalFormula>{"Ba+2", "OH-"} );
    CHECK( params.beta1[0].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.beta1[0].parameters == Vec<real>{ 1.2 } );

    CHECK( params.beta1[1].formulas == Vec<ChemicalFormula>{"Br-", "H+"} );
    CHECK( params.beta1[1].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.beta1[1].parameters == Vec<real>{ 0.3564, 0.0, 0.0, 4.467e-4 } );

    CHECK( params.beta1[2].formulas == Vec<ChemicalFormula>{"Cl-", "Na+"} );
    CHECK( params.beta1[2].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.beta1[2].parameters == Vec<real>{ 0.2769, 1.377e+4, 46.8, -6.9512e-2, 2.0e-5, -7.4823e+5 } );

    CHECK( params.beta2[0].formulas == Vec<ChemicalFormula>{"Ca+2", "OH-"} );
    CHECK( params.beta2[0].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.beta2[0].parameters == Vec<real>{ -5.72 } );

    CHECK( params.beta2[1].formulas == Vec<ChemicalFormula>{"Ca+2", "Cl-"} );
    CHECK( params.beta2[1].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.beta2[1].parameters == Vec<real>{ -1.13, 0.0, 0.0, -0.0476 } );

    CHECK( params.beta2[2].formulas == Vec<ChemicalFormula>{"Mg+2", "SO4-2"} );
    CHECK( params.beta2[2].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.beta2[2].parameters == Vec<real>{ -32.45, 0.0, -3.236e+3, 21.812, -1.8859e-2 } );

    CHECK( params.Cphi[0].formulas == Vec<ChemicalFormula>{"B(OH)4-", "Na+"} );
    CHECK( params.Cphi[0].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.Cphi[0].parameters == Vec<real>{ 0.0114 } );

    CHECK( params.Cphi[1].formulas == Vec<ChemicalFormula>{"Ba+2", "Cl-"} );
    CHECK( params.Cphi[1].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.Cphi[1].parameters == Vec<real>{ -0.143, -114.5 } );

    CHECK( params.Cphi[2].formulas == Vec<ChemicalFormula>{"Cl-", "Na+"} );
    CHECK( params.Cphi[2].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.Cphi[2].parameters == Vec<real>{ 1.48e-3, -120.5, -0.2081, 0.0, 1.166e-7, 11121 } );

    CHECK( params.theta[0].formulas == Vec<ChemicalFormula>{"K+", "Na+"} );
    CHECK( params.theta[0].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.theta[0].parameters == Vec<real>{ -0.012 } );

    CHECK( params.theta[1].formulas == Vec<ChemicalFormula>{"Ca+2", "K+"} );
    CHECK( params.theta[1].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.theta[1].parameters == Vec<real>{ -5.35e-3, 0.0, 0.0, 3.08e-4 } );

    CHECK( params.theta[2].formulas == Vec<ChemicalFormula>{"Ca+2", "Na+"} );
    CHECK( params.theta[2].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.theta[2].parameters == Vec<real>{ 9.22e-2, 0.0, 0.0, -4.29e-4, 1.21e-6 } );

    CHECK( params.lambda[0].formulas == Vec<ChemicalFormula>{"Cl-", "CO2"} );
    CHECK( params.lambda[0].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.lambda[0].parameters == Vec<real>{ -0.005 } );

    CHECK( params.lambda[1].formulas == Vec<ChemicalFormula>{"CO2", "CO2"} );
    CHECK( params.lambda[1].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.lambda[1].parameters == Vec<real>{ -1.34e-2, 348, 0.803 } );

    CHECK( params.lambda[2].formulas == Vec<ChemicalFormula>{"H4SiO4", "Mg+2"} );
    CHECK( params.lambda[2].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.lambda[2].parameters == Vec<real>{ 0.238, -1788, -9.023, 0.0103 } );

    CHECK( params.zeta[0].formulas == Vec<ChemicalFormula>{"Cl-", "H4SiO4", "K+"} );
    CHECK( params.zeta[0].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.zeta[0].parameters == Vec<real>{ -0.0153 } );

    CHECK( params.zeta[1].formulas == Vec<ChemicalFormula>{"CO2", "Na+", "SO4-2"} );
    CHECK( params.zeta[1].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.zeta[1].parameters == Vec<real>{ -0.015 } );

    CHECK( params.psi[0].formulas == Vec<ChemicalFormula>{"Br-", "K+", "Na+"} );
    CHECK( params.psi[0].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.psi[0].parameters == Vec<real>{ -0.0022 } );

    CHECK( params.psi[1].formulas == Vec<ChemicalFormula>{"Cl-", "Mg+2", "SO4-2"} );
    CHECK( params.psi[1].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.psi[1].parameters == Vec<real>{ -0.008, 32.63 } );

    CHECK( params.psi[2].formulas == Vec<ChemicalFormula>{"Ca+2", "Cl-", "Na+"} );
    CHECK( params.psi[2].model == ActivityModelParamsPitzer::CorrectionModel::Phreeqc );
    CHECK( params.psi[2].parameters == Vec<real>{ -1.48e-2, 0.0, 0.0, -5.2e-6 } );

    CHECK( params.mu[0].formulas == Vec<ChemicalFormula>{"NH3", "NH3", "CO3-2"} );
    CHECK( params.mu[0].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.mu[0].parameters == Vec<real>{ 0.1234 } );

    CHECK( params.eta[0].formulas == Vec<ChemicalFormula>{"CO2", "Na+", "K+"} );
    CHECK( params.eta[0].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.eta[0].parameters == Vec<real>{ -0.0022 } );

    CHECK( params.alpha1[0].formulas == Vec<ChemicalFormula>{"Na+", "Cl-"} );
    CHECK( params.alpha1[0].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.alpha1[0].parameters == Vec<real>{ 2.0 } );

    CHECK( params.alpha2[0].formulas == Vec<ChemicalFormula>{"Na+", "Cl-"} );
    CHECK( params.alpha2[0].model == ActivityModelParamsPitzer::CorrectionModel::Constant );
    CHECK( params.alpha2[0].parameters == Vec<real>{ 12.0 } );
}

const String yml_euniquac_old_format = R"#(
q:
  - [ H2O, 1.4 ]
  - [ CO2, 2.4496 ]
r:
  - [ H2O, 0.92 ]
  - [ CO2, 0.74721 ]
u:
  - [ H2O, CO2, [ 8.838254, 0.86293 ] ]
  - [ H2O, HCO3-, [ 577.0502, -0.38795 ] ]
  - [ H2O, CO3-2, [ 361.3877, 3.3516 ] ]
  - [ CO2, CO2, [ 302.2475, 0.35871 ] ]
  - [ CO2, HCO3-, [ 526.305, -3.7342 ] ]
  - [ CO2, CO3-2, [ 2500, 0 ] ]
)#";

const String yml_euniquac_new_format = R"#(
q:
  H2O: [ 1.4, 0.0, 0.0 ]
  CO2: [ 2.4496, 0.0, 0.0 ]
r:
  H2O: [ 0.92, 0.0, 0.0 ]
  CO2: [ 0.74721, 0.0, 0.0 ]
u:
  H2O:
    CO2: [ 8.838254, 0.86293, 0.0 ]
    HCO3-: [ 577.0502, -0.38795, 0.0 ]
    CO3-2: [ 361.3877, 3.3516, 0.0 ]
  CO2:
    CO2: [ 302.2475, 0.35871, 0.0 ]
    HCO3-: [ 526.305, -3.7342, 0.0 ]
    CO3-2: [ 2500, 0.0, 0.0 ]
)#";

TEST_CASE("Testing serialization of ActivityModelParamsExtendedUNIQUAC", "[Serialization][Models][ActivityModels][ActivityModelParamsExtendedUNIQUAC]")
{
    ActivityModelParamsExtendedUNIQUAC params;

    auto yml = GENERATE(yml_euniquac_old_format, yml_euniquac_new_format);

    params = Data::parse(yml).as<ActivityModelParamsExtendedUNIQUAC>();

    CHECK( params.q.size() == 2 );
    CHECK( params.r.size() == 2 );
    CHECK( params.u.size() == 6 );

    using Q = Tuple<String, Vec<real>>;
    using R = Tuple<String, Vec<real>>;
    using U = Tuple<String, String, Vec<real>>;

    CHECK( params.q[0] == Q{"H2O", { 1.4, 0.0, 0.0 }} );
    CHECK( params.q[1] == Q{"CO2", { 2.4496, 0.0, 0.0 }} );

    CHECK( params.r[0] == R{"H2O", { 0.92, 0.0, 0.0 }} );
    CHECK( params.r[1] == R{"CO2", { 0.74721, 0.0, 0.0 }} );

    CHECK( params.u[0] == U{"H2O", "CO2", { 8.838254, 0.86293, 0.0 }} );
    CHECK( params.u[1] == U{"H2O", "HCO3-", { 577.0502, -0.38795, 0.0 }} );
    CHECK( params.u[2] == U{"H2O", "CO3-2", { 361.3877, 3.3516, 0.0 }} );
    CHECK( params.u[3] == U{"CO2", "CO2", { 302.2475, 0.35871, 0.0 }} );
    CHECK( params.u[4] == U{"CO2", "HCO3-", { 526.305, -3.7342, 0.0 }} );
    CHECK( params.u[5] == U{"CO2", "CO3-2", { 2500, 0.0, 0.0 }} );
}
