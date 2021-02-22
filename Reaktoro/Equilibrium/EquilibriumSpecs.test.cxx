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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumSpecs", "[EquilibriumSpecs]")
{
    ChemicalSystem system = test::createChemicalSystem();

    EquilibriumSpecs specs(system);

    ChemicalState state(system);
    state.setTemperature(50.0, "celsius");
    state.setPressure(100.0, "bar");
    state.setSpeciesAmounts(1.0);

    ChemicalProps props(state);

    WHEN("temperature and pressure are input parameters - the Gibbs energy minimization formulation")
    {
        specs.temperature();
        specs.pressure();

        CHECK( specs.numParameters()                       == 2 ); // T, P
        CHECK( specs.numControlVariables()                 == 0 );
        CHECK( specs.numTitrants()                         == 0 );
        CHECK( specs.numTitrantsExplicit()                 == 0 );
        CHECK( specs.numTitrantsImplicit()                 == 0 );
        CHECK( specs.numConstraints()                      == 0 );
        CHECK( specs.numConstraintsEquationType()          == 0 );
        CHECK( specs.numConstraintsChemicalPotentialType() == 0 );

        CHECK( specs.namesParameters() == Strings{"T", "P"} );

        CHECK( specs.isTemperatureUnknown() == false );
        CHECK( specs.isPressureUnknown()    == false );
    }

    WHEN("temperature and volume are input parameters - the Helmholtz energy minimization formulation")
    {
        specs.temperature();
        specs.volume();

        CHECK( specs.numParameters()                       == 2 ); // T, V
        CHECK( specs.numControlVariables()                 == 1 ); // P
        CHECK( specs.numTitrants()                         == 0 );
        CHECK( specs.numTitrantsExplicit()                 == 0 );
        CHECK( specs.numTitrantsImplicit()                 == 0 );
        CHECK( specs.numConstraints()                      == 1 ); // V = V(given)
        CHECK( specs.numConstraintsEquationType()          == 1 ); // V = V(given)
        CHECK( specs.numConstraintsChemicalPotentialType() == 0 );

        CHECK( specs.namesParameters()              == Strings{"T", "V"} );
        CHECK( specs.namesControlVariables()        == Strings{"P"} );
        CHECK( specs.namesConstraints()             == Strings{"volume"} );
        CHECK( specs.namesConstraintsEquationType() == Strings{"volume"} );

        CHECK( specs.isTemperatureUnknown() == false );
        CHECK( specs.isPressureUnknown()    == true );
    }

    WHEN("volume and internal energy are input parameters - the entropy maximization formulation")
    {
        specs.volume();
        specs.internalEnergy();

        CHECK( specs.numParameters()                       == 2 ); // V, U
        CHECK( specs.numControlVariables()                 == 2 ); // T, P
        CHECK( specs.numTitrants()                         == 0 );
        CHECK( specs.numTitrantsExplicit()                 == 0 );
        CHECK( specs.numTitrantsImplicit()                 == 0 );
        CHECK( specs.numConstraints()                      == 2 ); // V = V(given), U = U(given)
        CHECK( specs.numConstraintsEquationType()          == 2 ); // V = V(given), U = U(given)
        CHECK( specs.numConstraintsChemicalPotentialType() == 0 );

        CHECK( specs.namesParameters()              == Strings{"V", "U"} );
        CHECK( specs.namesControlVariables()        == Strings{"T", "P"} );
        CHECK( specs.namesConstraints()             == Strings{"volume", "internalEnergy"} );
        CHECK( specs.namesConstraintsEquationType() == Strings{"volume", "internalEnergy"} );

        CHECK( specs.isTemperatureUnknown() == true );
        CHECK( specs.isPressureUnknown()    == true );
    }

    WHEN("temperature, pressure, and pH are input parameters")
    {
        specs.temperature();
        specs.pressure();
        specs.pH();

        CHECK( specs.numParameters()                       == 3 ); // T, P, pH
        CHECK( specs.numControlVariables()                 == 1 ); // n[H+]
        CHECK( specs.numTitrants()                         == 1 ); // [H+]
        CHECK( specs.numTitrantsExplicit()                 == 0 );
        CHECK( specs.numTitrantsImplicit()                 == 1 ); // [H+]
        CHECK( specs.numConstraints()                      == 1 ); // pH = pH(given)
        CHECK( specs.numConstraintsEquationType()          == 0 );
        CHECK( specs.numConstraintsChemicalPotentialType() == 1 ); // pH = pH(given)

        CHECK( specs.namesParameters()                       == Strings{"T", "P", "pH"} );
        CHECK( specs.namesControlVariables()                 == Strings{"[H+]"} );
        CHECK( specs.namesTitrants()                         == Strings{"[H+]"} );
        CHECK( specs.namesTitrantsImplicit()                 == Strings{"[H+]"} );
        CHECK( specs.namesConstraints()                      == Strings{"pH"} );
        CHECK( specs.namesConstraintsChemicalPotentialType() == Strings{"pH"} );

        CHECK( specs.isTemperatureUnknown() == false );
        CHECK( specs.isPressureUnknown()    == false );
    }

    WHEN("volume, entropy, and activity[CO2(g)] are input parameters")
    {
        specs.volume();
        specs.entropy();
        specs.activity("CO2(g)");

        CHECK( specs.numParameters()                       == 3 ); // V, S, a[CO2(g)]
        CHECK( specs.numControlVariables()                 == 3 ); // T, P, n[CO2]
        CHECK( specs.numTitrants()                         == 1 ); // [CO2]
        CHECK( specs.numTitrantsExplicit()                 == 0 );
        CHECK( specs.numTitrantsImplicit()                 == 1 ); // [CO2]
        CHECK( specs.numConstraints()                      == 3 ); // V = V(given), S = S(given), a[CO2(g)] = a[CO2(g)](given)
        CHECK( specs.numConstraintsEquationType()          == 2 ); // V = V(given), S = S(given)
        CHECK( specs.numConstraintsChemicalPotentialType() == 1 ); // a[CO2(g)] = a[CO2(g)](given)

        CHECK( specs.namesParameters()                       == Strings{"V", "S", "lnActivity[CO2(g)]"} );
        CHECK( specs.namesControlVariables()                 == Strings{"T", "P", "[CO2]"} );
        CHECK( specs.namesTitrants()                         == Strings{"[CO2]"} );
        CHECK( specs.namesTitrantsImplicit()                 == Strings{"[CO2]"} );
        CHECK( specs.namesConstraints()                      == Strings{"volume", "entropy", "lnActivity[CO2(g)]"} );
        CHECK( specs.namesConstraintsEquationType()          == Strings{"volume", "entropy"} );
        CHECK( specs.namesConstraintsChemicalPotentialType() == Strings{"lnActivity[CO2(g)]"} );

        CHECK( specs.isTemperatureUnknown() == true );
        CHECK( specs.isPressureUnknown()    == true );
    }

    WHEN("temperature, pressure, volume, internal energy, pH, and pE are input parameters")
    {
        specs.temperature();
        specs.pressure();
        specs.volume();
        specs.internalEnergy();
        specs.pH();
        specs.pE();
        specs.openTo("CO2");
        specs.openTo("CH4");

        CHECK( specs.numParameters()                       == 6 ); // T, P, V, U, pH, pE
        CHECK( specs.numControlVariables()                 == 4 ); // n[CO2], n[CH4], n[H+], n[e-]
        CHECK( specs.numTitrants()                         == 4 ); // [CO2], [CH4], [H+], [e-]
        CHECK( specs.numTitrantsExplicit()                 == 2 ); // [CO2], [CH4]
        CHECK( specs.numTitrantsImplicit()                 == 2 ); // [H+], [e-]
        CHECK( specs.numConstraints()                      == 4 ); // V = V(given), U = U(given), pH = pH(given), pE = pE(given)
        CHECK( specs.numConstraintsEquationType()          == 2 ); // V = V(given), U = U(given)
        CHECK( specs.numConstraintsChemicalPotentialType() == 2 ); // pH = pH(given), pE = pE(given)

        CHECK( specs.namesParameters()                       == Strings{"T", "P", "V", "U", "pH", "pE"} );
        CHECK( specs.namesControlVariables()                 == Strings{"[CO2]", "[CH4]", "[H+]", "[e-]"} );
        CHECK( specs.namesTitrants()                         == Strings{"[CO2]", "[CH4]", "[H+]", "[e-]"} );
        CHECK( specs.namesTitrantsExplicit()                 == Strings{"[CO2]", "[CH4]"} );
        CHECK( specs.namesTitrantsImplicit()                 == Strings{"[H+]", "[e-]"} );
        CHECK( specs.namesConstraints()                      == Strings{"volume", "internalEnergy", "pH", "pE"} );
        CHECK( specs.namesConstraintsEquationType()          == Strings{"volume", "internalEnergy"} );
        CHECK( specs.namesConstraintsChemicalPotentialType() == Strings{"pH", "pE"} );

        CHECK( specs.isTemperatureUnknown() == false );
        CHECK( specs.isPressureUnknown()    == false );
    }

    SECTION("Checking lambda functions in equation constraints")
    {
        const auto& econstraints = specs.constraintsEquationType();

        Params params;
        const auto V = 1.0; params.append("V", V);
        const auto U = 2.0; params.append("U", U);
        const auto H = 3.0; params.append("H", H);
        const auto G = 4.0; params.append("G", G);
        const auto A = 5.0; params.append("A", A);
        const auto S = 6.0; params.append("S", S);

        specs.volume();
        specs.internalEnergy();
        specs.enthalpy();
        specs.gibbsEnergy();
        specs.helmholtzEnergy();
        specs.entropy();

        CHECK( econstraints[0].fn(props, params) == props.volume() - V );
        CHECK( econstraints[1].fn(props, params) == props.internalEnergy() - U );
        CHECK( econstraints[2].fn(props, params) == props.enthalpy() - H );
        CHECK( econstraints[3].fn(props, params) == props.gibbsEnergy() - G );
        CHECK( econstraints[4].fn(props, params) == props.helmholtzEnergy() - A );
        CHECK( econstraints[5].fn(props, params) == props.entropy() - S );
    }

    SECTION("Checking lambda functions in chemical potential constraints")
    {
        const auto& uconstraints = specs.constraintsChemicalPotentialType();

        Params params;
        const auto p0 = 1.0; params.append("u[H2O(aq)]", p0);
        const auto p1 = 2.0; params.append("lnActivity[CH4(g)]", p1);
        const auto p2 = 3.0; params.append("lnActivity[CO2(g)]", p2);
        const auto p3 = 4.0; params.append("lnActivity[Ca++(aq)]", p3);
        const auto p4 = 5.0; params.append("f[O2]", p4);
        const auto p5 = 6.0; params.append("pH", p5);
        const auto p6 = 7.0; params.append("pMg", p6);
        const auto p7 = 8.0; params.append("pE", p7);
        const auto p8 = 9.0; params.append("Eh", p8);

        specs.chemicalPotential("H2O(aq)");
        specs.lnActivity("CH4(g)");
        specs.lgActivity("CO2(g)");
        specs.activity("Ca++(aq)");
        specs.fugacity("O2");
        specs.pH();
        specs.pMg();

        const auto T = props.temperature();
        const auto P = props.pressure();
        const auto RT = universalGasConstant * T;

        const auto u0CH4  = system.species().get("CH4(g)").props(T, P).G0;
        const auto u0CO2  = system.species().get("CO2(g)").props(T, P).G0;
        const auto u0Capp = system.species().get("Ca++(aq)").props(T, P).G0;
        const auto u0O2   = system.species().get("O2(g)").props(T, P).G0;
        const auto u0Hp   = system.species().get("H+(aq)").props(T, P).G0;
        const auto u0Mgpp = system.species().get("Mg++(aq)").props(T, P).G0;

        CHECK( uconstraints[0].fn(props, params) == p0 );
        CHECK( uconstraints[1].fn(props, params) == u0CH4 + RT*p1 );
        CHECK( uconstraints[2].fn(props, params) == u0CO2 + RT*p2 );
        CHECK( uconstraints[3].fn(props, params) == u0Capp + RT*p3 );
        CHECK( uconstraints[4].fn(props, params) == u0O2 + RT*log(p4) );
        CHECK( uconstraints[5].fn(props, params) == u0Hp + RT*p5 * (-ln10) );
        CHECK( uconstraints[6].fn(props, params) == u0Mgpp + RT*p6 * (-ln10) );

        WHEN("pE is imposed instead of Eh")
        {
            specs.pE();

            CHECK( uconstraints.back().fn(props, params) == RT*p7 * (-ln10) );
        }

        WHEN("Eh is imposed instead of pE")
        {
            specs.Eh();

            const auto F = faradayConstant;

            CHECK( uconstraints.back().fn(props, params) == -F * p8 );
        }
    }
}
