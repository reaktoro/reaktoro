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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumSpecs", "[EquilibriumSpecs]")
{
    // The mock chemical system for the tests below.
    ChemicalSystem system = test::createChemicalSystem();

    EquilibriumSpecs specs(system);

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
}
