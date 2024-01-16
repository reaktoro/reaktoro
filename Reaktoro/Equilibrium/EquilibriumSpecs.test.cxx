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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumSpecs", "[EquilibriumSpecs]")
{
    ChemicalSystem system = test::createChemicalSystem();

    const auto Ne = system.elements().size();
    const auto Nb = 1 + Ne; // elements + charge components

    ChemicalState state(system);
    state.setTemperature(50.0, "celsius");
    state.setPressure(100.0, "bar");
    state.setSpeciesAmounts(1.0);

    state.props().update(state);

    const auto props = state.props();

    EquilibriumSpecs specs(system);

    WHEN("the EquilibriumSpecs object holds default state")
    {
        // Check both temperature and pressure are unknowns by default
        CHECK( specs.isTemperatureUnknown() );
        CHECK( specs.isPressureUnknown() );

        // Check both temperature and pressure are in the list of *p* control variables by default
        CHECK( contains(specs.namesControlVariablesP(), "T") );
        CHECK( contains(specs.namesControlVariablesP(), "P") );

        // Check both temperature and pressure are not in the list of *w* input variables by default
        CHECK( !contains(specs.namesInputs(), "T") );
        CHECK( !contains(specs.namesInputs(), "P") );
    }

    WHEN("temperature and pressure are input variables - the Gibbs energy minimization formulation")
    {
        specs.temperature();
        specs.pressure();

        CHECK( specs.numInputs()                              == 2 ); // T, P
        CHECK( specs.numControlVariables()                    == 0 );
        CHECK( specs.numControlVariablesP()                   == 0 );
        CHECK( specs.numControlVariablesQ()                   == 0 );
        CHECK( specs.numTitrants()                            == 0 );
        CHECK( specs.numTitrantsExplicit()                    == 0 );
        CHECK( specs.numTitrantsImplicit()                    == 0 );
        CHECK( specs.numConstraints()                         == 0 );
        CHECK( specs.numConservativeComponents()              == Nb );
        CHECK( specs.namesInputs()                            == Strings{"T", "P"} );
        CHECK( specs.isTemperatureUnknown()                   == false );
        CHECK( specs.isPressureUnknown()                      == false );
        CHECK( specs.indexTemperatureAmongControlVariablesP() == Index(-1) );
        CHECK( specs.indexPressureAmongControlVariablesP()    == Index(-1) );
    }

    WHEN("temperature and volume are input variables - the Helmholtz energy minimization formulation")
    {
        specs.temperature();
        specs.volume();

        CHECK( specs.numInputs()                              == 2 ); // T, V
        CHECK( specs.numControlVariables()                    == 1 ); // P
        CHECK( specs.numControlVariablesP()                   == 1 ); // P
        CHECK( specs.numControlVariablesQ()                   == 0 );
        CHECK( specs.numTitrants()                            == 0 );
        CHECK( specs.numTitrantsExplicit()                    == 0 );
        CHECK( specs.numTitrantsImplicit()                    == 0 );
        CHECK( specs.numConstraints()                         == 1 ); // V = V(given)
        CHECK( specs.numConservativeComponents()              == Nb );
        CHECK( specs.namesInputs()                            == Strings{"T", "V"} );
        CHECK( specs.namesControlVariables()                  == Strings{"P"} );
        CHECK( specs.namesConstraints()                       == Strings{"volume"} );
        CHECK( specs.isTemperatureUnknown()                   == false );
        CHECK( specs.isPressureUnknown()                      == true );
        CHECK( specs.indexTemperatureAmongControlVariablesP() == Index(-1) );
        CHECK( specs.indexPressureAmongControlVariablesP()    == 0 );
    }

    WHEN("volume and internal energy are input variables - the entropy maximization formulation")
    {
        specs.volume();
        specs.internalEnergy();

        CHECK( specs.numInputs()                              == 2 ); // V, U
        CHECK( specs.numControlVariables()                    == 2 ); // T, P
        CHECK( specs.numControlVariablesP()                   == 2 ); // T, P
        CHECK( specs.numControlVariablesQ()                   == 0 );
        CHECK( specs.numTitrants()                            == 0 );
        CHECK( specs.numTitrantsExplicit()                    == 0 );
        CHECK( specs.numTitrantsImplicit()                    == 0 );
        CHECK( specs.numConstraints()                         == 2 ); // V = V(given), U = U(given)
        CHECK( specs.numConservativeComponents()              == Nb );
        CHECK( specs.namesInputs()                            == Strings{"V", "U"} );
        CHECK( specs.namesControlVariables()                  == Strings{"T", "P"} );
        CHECK( specs.namesConstraints()                       == Strings{"volume", "internalEnergy"} );
        CHECK( specs.isTemperatureUnknown()                   == true );
        CHECK( specs.isPressureUnknown()                      == true );
        CHECK( specs.indexTemperatureAmongControlVariablesP() == 0 );
        CHECK( specs.indexPressureAmongControlVariablesP()    == 1 );
    }

    WHEN("temperature, pressure, and pH are input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.pH();

        CHECK( specs.numInputs()                              == 3 ); // T, P, pH
        CHECK( specs.numControlVariables()                    == 1 ); // n[H+]
        CHECK( specs.numControlVariablesP()                   == 0 );
        CHECK( specs.numControlVariablesQ()                   == 1 ); // n[H+]
        CHECK( specs.numTitrants()                            == 1 ); // [H+]
        CHECK( specs.numTitrantsExplicit()                    == 0 );
        CHECK( specs.numTitrantsImplicit()                    == 1 ); // [H+]
        CHECK( specs.numConstraints()                         == 1 ); // pH = pH(given)
        CHECK( specs.numConservativeComponents()              == Nb );
        CHECK( specs.namesInputs()                            == Strings{"T", "P", "pH"} );
        CHECK( specs.namesControlVariables()                  == Strings{"[H+]"} );
        CHECK( specs.namesTitrants()                          == Strings{"[H+]"} );
        CHECK( specs.namesTitrantsImplicit()                  == Strings{"[H+]"} );
        CHECK( specs.namesConstraints()                       == Strings{"pH"} );
        CHECK( specs.isTemperatureUnknown()                   == false );
        CHECK( specs.isPressureUnknown()                      == false );
        CHECK( specs.indexTemperatureAmongControlVariablesP() == Index(-1) );
        CHECK( specs.indexPressureAmongControlVariablesP()    == Index(-1) );
    }

    WHEN("volume, entropy, and activity[CO2(g)] are input variables")
    {
        specs.volume();
        specs.entropy();
        specs.activity("CO2(g)");

        CHECK( specs.numInputs()                              == 3 ); // V, S, ln(a[CO2(g)])
        CHECK( specs.numControlVariables()                    == 3 ); // T, P, n[CO2]
        CHECK( specs.numControlVariablesP()                   == 2 ); // T, P
        CHECK( specs.numControlVariablesQ()                   == 1 ); // n[CO2]
        CHECK( specs.numTitrants()                            == 1 ); // [CO2]
        CHECK( specs.numTitrantsExplicit()                    == 0 );
        CHECK( specs.numTitrantsImplicit()                    == 1 ); // [CO2]
        CHECK( specs.numConstraints()                         == 3 ); // V = V(given), S = S(given), a[CO2(g)] = a[CO2(g)](given)
        CHECK( specs.numConservativeComponents()              == Nb );
        CHECK( specs.namesInputs()                            == Strings{"V", "S", "ln(a[CO2(g)])"} );
        CHECK( specs.namesControlVariables()                  == Strings{"T", "P", "[CO2(g)]"} );
        CHECK( specs.namesTitrants()                          == Strings{"[CO2]"} );
        CHECK( specs.namesTitrantsImplicit()                  == Strings{"[CO2]"} );
        CHECK( specs.namesConstraints()                       == Strings{"volume", "entropy", "ln(a[CO2(g)])"} );
        CHECK( specs.isTemperatureUnknown()                   == true );
        CHECK( specs.isPressureUnknown()                      == true );
        CHECK( specs.indexTemperatureAmongControlVariablesP() == 0 );
        CHECK( specs.indexPressureAmongControlVariablesP()    == 1 );
    }

    WHEN("temperature, pressure, volume, internal energy, pH, and pE are input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.volume();
        specs.internalEnergy();
        specs.pH();
        specs.pE();
        specs.openTo("CO2");
        specs.openTo("CH4");

        CHECK( specs.numInputs()                              == 6 ); // T, P, V, U, pH, pE
        CHECK( specs.numControlVariables()                    == 4 ); // n[CO2], n[CH4], n[H+], n[e-]
        CHECK( specs.numControlVariablesP()                   == 2 ); // n[CO2], n[CH4]
        CHECK( specs.numControlVariablesQ()                   == 2 ); // n[H+], n[e-]
        CHECK( specs.numTitrants()                            == 4 ); // [CO2], [CH4], [H+], [e-]
        CHECK( specs.numTitrantsExplicit()                    == 2 ); // [CO2], [CH4]
        CHECK( specs.numTitrantsImplicit()                    == 2 ); // [H+], [e-]
        CHECK( specs.numConstraints()                         == 4 ); // V = V(given), U = U(given), pH = pH(given), pE = pE(given)
        CHECK( specs.numConservativeComponents()              == Nb );
        CHECK( specs.namesInputs()                            == Strings{"T", "P", "V", "U", "pH", "pE"} );
        CHECK( specs.namesControlVariables()                  == Strings{"[CO2]", "[CH4]", "[H+]", "[e-]"} );
        CHECK( specs.namesTitrants()                          == Strings{"[CO2]", "[CH4]", "[H+]", "[e-]"} );
        CHECK( specs.namesTitrantsExplicit()                  == Strings{"[CO2]", "[CH4]"} );
        CHECK( specs.namesTitrantsImplicit()                  == Strings{"[H+]", "[e-]"} );
        CHECK( specs.namesConstraints()                       == Strings{"volume", "internalEnergy", "pH", "pE"} );
        CHECK( specs.isTemperatureUnknown()                   == false );
        CHECK( specs.isPressureUnknown()                      == false );
        CHECK( specs.indexTemperatureAmongControlVariablesP() == Index(-1) );
        CHECK( specs.indexPressureAmongControlVariablesP()    == Index(-1) );
    }

    WHEN("temperature, pressure, volume, internal energy, pH, and pE are input variables and reactivity constraints, aka restricted reactions, are introduced")
    {
        const auto Nn = system.species().size();

        specs.temperature();
        specs.pressure();
        specs.volume();
        specs.internalEnergy();
        specs.pH();
        specs.pE();
        specs.openTo("CO2");
        specs.openTo("CH4");
        specs.addReactivityConstraint({ "xi1", VectorXd::Random(Nn), {} });
        specs.addReactivityConstraint({ "xi2", VectorXd::Random(Nn), {} });

        CHECK( specs.numInputs()                              == 6 ); // T, P, V, U, pH, pE
        CHECK( specs.numControlVariables()                    == 4 ); // n[CO2], n[CH4], n[H+], n[e-]
        CHECK( specs.numControlVariablesP()                   == 2 ); // n[CO2], n[CH4]
        CHECK( specs.numControlVariablesQ()                   == 2 ); // n[H+], n[e-]
        CHECK( specs.numTitrants()                            == 4 ); // [CO2], [CH4], [H+], [e-]
        CHECK( specs.numTitrantsExplicit()                    == 2 ); // [CO2], [CH4]
        CHECK( specs.numTitrantsImplicit()                    == 2 ); // [H+], [e-]
        CHECK( specs.numConstraints()                         == 6 ); // V = V(given), U = U(given), pH = pH(given), pE = pE(given), 2x reactivity constraints
        CHECK( specs.numConservativeComponents()              == Nb + 2 ); // elements, charge, and two reactivity constraints
        CHECK( specs.namesInputs()                            == Strings{"T", "P", "V", "U", "pH", "pE"} );
        CHECK( specs.namesControlVariables()                  == Strings{"[CO2]", "[CH4]", "[H+]", "[e-]"} );
        CHECK( specs.namesTitrants()                          == Strings{"[CO2]", "[CH4]", "[H+]", "[e-]"} );
        CHECK( specs.namesTitrantsExplicit()                  == Strings{"[CO2]", "[CH4]"} );
        CHECK( specs.namesTitrantsImplicit()                  == Strings{"[H+]", "[e-]"} );
        CHECK( specs.namesConstraints()                       == Strings{"volume", "internalEnergy", "pH", "pE", "xi1", "xi2"} );
        CHECK( specs.isTemperatureUnknown()                   == false );
        CHECK( specs.isPressureUnknown()                      == false );
        CHECK( specs.indexTemperatureAmongControlVariablesP() == Index(-1) );
        CHECK( specs.indexPressureAmongControlVariablesP()    == Index(-1) );
    }

    WHEN("model parameters are among the input variables")
    {
        specs.temperature();
        specs.pressure();
        specs.addInput("V");
        specs.addInput("G0[H2O]");

        CHECK( specs.numInputs()                              == 4 ); // T, P, V, G0[H2O]
        CHECK( specs.numControlVariables()                    == 0 );
        CHECK( specs.numControlVariablesP()                   == 0 );
        CHECK( specs.numControlVariablesQ()                   == 0 );
        CHECK( specs.numTitrants()                            == 0 );
        CHECK( specs.numTitrantsExplicit()                    == 0 );
        CHECK( specs.numTitrantsImplicit()                    == 0 );
        CHECK( specs.numConstraints()                         == 0 );
        CHECK( specs.numConservativeComponents()              == Nb );
        CHECK( specs.namesInputs()                            == Strings{"T", "P", "V", "G0[H2O]"} );
        CHECK( specs.namesControlVariables()                  == Strings{} );
        CHECK( specs.namesTitrants()                          == Strings{} );
        CHECK( specs.namesTitrantsExplicit()                  == Strings{} );
        CHECK( specs.namesTitrantsImplicit()                  == Strings{} );
        CHECK( specs.namesConstraints()                       == Strings{} );
        CHECK( specs.isTemperatureUnknown()                   == false );
        CHECK( specs.isPressureUnknown()                      == false );
        CHECK( specs.indexTemperatureAmongControlVariablesP() == Index(-1) );
        CHECK( specs.indexPressureAmongControlVariablesP()    == Index(-1) );
    }

    WHEN("temperature and pressure become unknowns and then inputs")
    {
        // Specify that temperature and pressure are unknowns
        specs.unknownTemperature();
        specs.unknownPressure();

        // Ensure temperature and pressure are unknowns
        CHECK( specs.isTemperatureUnknown() );
        CHECK( specs.isPressureUnknown() );

        // Ensure temperature and pressure are in the list of *p* control variables
        CHECK( contains(specs.namesControlVariablesP(), "T") );
        CHECK( contains(specs.namesControlVariablesP(), "P") );

        // Ensure temperature and pressure are not in the list of *w* input variables
        CHECK( !contains(specs.namesInputs(), "T") );
        CHECK( !contains(specs.namesInputs(), "P") );

        //=========================================================================================

        // Specify that temperature and pressure are now known and given
        specs.temperature();
        specs.pressure();

        // Ensure temperature and pressure are not unknowns
        CHECK( !specs.isTemperatureUnknown() );
        CHECK( !specs.isPressureUnknown() );

        // Ensure temperature and pressure are not in the list of *p* control variables
        CHECK( !contains(specs.namesControlVariablesP(), "T") );
        CHECK( !contains(specs.namesControlVariablesP(), "P") );

        // Ensure temperature and pressure are in the list of *w* input variables
        CHECK( contains(specs.namesInputs(), "T") );
        CHECK( contains(specs.namesInputs(), "P") );
    }

    SECTION("Checking lambda functions in equation constraints")
    {
        specs.volume();
        specs.internalEnergy();
        specs.enthalpy();
        specs.gibbsEnergy();
        specs.helmholtzEnergy();
        specs.entropy();

        const VectorXr p = {};
        const VectorXr w = random(specs.numInputs());

        const VectorXr v = specs.assembleEquationConstraints().fn(props, p, w);

        CHECK( v[0] == props.volume() - w[0] );
        CHECK( v[1] == props.internalEnergy() - w[1] );
        CHECK( v[2] == props.enthalpy() - w[2] );
        CHECK( v[3] == props.gibbsEnergy() - w[3] );
        CHECK( v[4] == props.helmholtzEnergy() - w[4] );
        CHECK( v[5] == props.entropy() - w[5] );
    }

    SECTION("Checking lambda functions in chemical potential constraints")
    {
        const auto constraintEh = GENERATE(true, false); // constrain Eh if true, pE if false (cannot be both constrained because of same titrant e-)

        specs.chemicalPotential("H2O(aq)");
        specs.lnActivity("CH4(g)");
        specs.lgActivity("CO2(g)");
        specs.activity("Ca++(aq)");
        specs.fugacity("O2");
        specs.pH();
        specs.pMg();

        if(constraintEh) specs.Eh(); else specs.pE(); // Not possible to constraint pE and Eh simultaneously

        const VectorXr p = {}; // the p control variables (empty in this example)
        const VectorXr w = random(specs.numInputs()).cwiseAbs(); // the input variables for the constrained properties above

        const auto T  = props.temperature();
        const auto P  = props.pressure();
        const auto RT = universalGasConstant * T;
        const auto F  = faradayConstant;

        const auto u0CH4  = system.species().get("CH4(g)").standardThermoProps(T, P).G0;
        const auto u0CO2  = system.species().get("CO2(g)").standardThermoProps(T, P).G0;
        const auto u0Capp = system.species().get("Ca++(aq)").standardThermoProps(T, P).G0;
        const auto u0O2   = system.species().get("O2(g)").standardThermoProps(T, P).G0;
        const auto u0Hp   = system.species().get("H+(aq)").standardThermoProps(T, P).G0;
        const auto u0Mgpp = system.species().get("Mg++(aq)").standardThermoProps(T, P).G0;

        auto const& qvars = specs.controlVariablesQ();

        CHECK( qvars[0].fn(props, p, w) == Approx(w[0]) );
        CHECK( qvars[1].fn(props, p, w) == Approx(u0CH4 + RT*w[1]) );
        CHECK( qvars[2].fn(props, p, w) == Approx(u0CO2 + RT*w[2]) );
        CHECK( qvars[3].fn(props, p, w) == Approx(u0Capp + RT*w[3]) );
        CHECK( qvars[4].fn(props, p, w) == Approx(u0O2 + RT*log(w[4])) );
        CHECK( qvars[5].fn(props, p, w) == Approx(u0Hp + RT*w[5] * (-ln10)) );
        CHECK( qvars[6].fn(props, p, w) == Approx(u0Mgpp + RT*w[6] * (-ln10)) );
        CHECK( qvars[7].fn(props, p, w) == Approx(constraintEh ? -F * w[7] : RT*w[7] * (-ln10)) );
    }

    SECTION("Checking when titrant unknowns are introduced")
    {
        WHEN("temperature and pressure are given inputs")
        {
            specs.temperature();
            specs.pressure();

            specs.addUnknownTitrantAmount("Ca++");
            specs.addUnknownTitrantAmount("H+");
            specs.addUnknownTitrantAmount("CaCO3");
            specs.addUnknownTitrantAmount("SiO2");

            CHECK( specs.numControlVariablesP() == 4 );

            auto pvars = specs.controlVariablesP();

            CHECK( pvars[0].substance == "Ca++" );
            CHECK( pvars[1].substance == "H+" );
            CHECK( pvars[2].substance == "CaCO3" );
            CHECK( pvars[3].substance == "SiO2" );

            CHECK( pvars[0].name == "[Ca++]" );
            CHECK( pvars[1].name == "[H+]" );
            CHECK( pvars[2].name == "[CaCO3]" );
            CHECK( pvars[3].name == "[SiO2]" );
        }
    }
}
