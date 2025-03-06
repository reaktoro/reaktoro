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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Extensions/Nasa/NasaDatabase.hpp>
using namespace Reaktoro;

TEST_CASE("Testing NasaDatabase module", "[NasaDatabase]")
{
    NasaDatabase db("nasa-cea");

    CHECK( db.species().size() == 1658 );

    CHECK( db.species().findWithName("AlF+") );
    CHECK( db.species().findWithName("B3H9") );
    CHECK( db.species().findWithName("C3S2") );
    CHECK( db.species().findWithName("Zr+") );
    CHECK( db.species().findWithName("Ti+") );
    CHECK( db.species().findWithName("CO2") );
    CHECK( db.species().findWithName("PbBr2(cr)") );
    CHECK( db.species().findWithName("Jet-A(g)") );

    const auto Tr = 298.15; // in K
    const auto Pr = 1.0e5;  // in Pa

    CHECK(db.species().getWithName("CO2").props(Tr, Pr).H0 == Approx( -393510.000) );
    CHECK(db.species().getWithName("Zr+").props(Tr, Pr).H0 == Approx( 1246246.292) );
    CHECK(db.species().getWithName("Ti+").props(Tr, Pr).H0 == Approx( 1137624.029) );
}

TEST_CASE("Testing NasaDatabase with EquilibriumSolver", "[NasaDatabase][EquilibriumSolver]")
{
    NasaDatabase db("nasa-cea");

    EquilibriumOptions options;
    options.optima.maxiters = 50;

    auto checkEquilibriumSolver = [&](auto reactantName, auto xiters)
    {
        auto oxidizerName = "Air";
        auto reactantAmount = 1.0; // in mol
        auto oxidizerAmount = 2.0; // in mol

        auto reactant = db.species().getWithName(reactantName);
        auto oxidizer = db.species().getWithName(oxidizerName);

        auto symbols = unique(concatenate(reactant.elements().symbols(), oxidizer.elements().symbols()));

        GaseousPhase gases(speciate(symbols));
        CondensedPhases condensed(speciate(symbols));

        ChemicalSystem system(db, gases, condensed);

        ChemicalState state0(system);
        state0.temperature(25.0, "celsius");
        state0.pressure(1.0, "atm");
        state0.setSpeciesAmounts(1e-16);
        state0.set(reactantName, reactantAmount, "mol");
        state0.set(oxidizerName, oxidizerAmount, "mol");

        ChemicalProps props0(state0);

        EquilibriumSpecs specs(system);
        specs.pressure();
        specs.enthalpy();

        EquilibriumConditions conditions(specs);
        conditions.pressure(props0.pressure());
        conditions.enthalpy(props0.enthalpy());
        conditions.setLowerBoundTemperature(100.0, "celsius");
        conditions.setUpperBoundTemperature(4000.0, "celsius");

        ChemicalState state(state0);

        EquilibriumSolver solver(specs);
        solver.setOptions(options);

        auto result = solver.solve(state, conditions);

        CHECK( result.succeeded() );
        CHECK( result.iterations() == xiters );
    };

    checkEquilibriumSolver("Mg(cd)", 43);
    checkEquilibriumSolver("CH4", 34);
}
