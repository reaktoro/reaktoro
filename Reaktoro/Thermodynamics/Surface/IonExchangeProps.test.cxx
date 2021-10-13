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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Thermodynamics/Surface/IonExchangeProps.hpp>
#include <Reaktoro/Singletons/Elements.hpp>

using namespace Reaktoro;

namespace test { extern auto createDatabase() -> Database; }

TEST_CASE("Testing IonExchangeProps class", "[IonExchangeProps]")
{
    // Extend the list of elements by the exchanger element
    const auto X = Element().withSymbol("X").withMolarMass(10.0);
    Elements::append(X);

    // Create custom database
    Database db = test::createDatabase();

    db.addSpecies( Species("X-"   ).withName("X-"   ).withAggregateState(AggregateState::IonExchange ).withStandardGibbsEnergy(0.0));
    db.addSpecies( Species("AlX3" ).withName("AlX3" ).withAggregateState(AggregateState::IonExchange ).withStandardGibbsEnergy(0.0));
    db.addSpecies( Species("CaX2" ).withName("CaX2" ).withAggregateState(AggregateState::IonExchange ).withStandardGibbsEnergy(0.0));
    db.addSpecies( Species("KX"   ).withName("KX"   ).withAggregateState(AggregateState::IonExchange ).withStandardGibbsEnergy(0.0));
    db.addSpecies( Species("MgX2" ).withName("MgX2" ).withAggregateState(AggregateState::IonExchange ).withStandardGibbsEnergy(0.0));
    db.addSpecies( Species("NaX"  ).withName("NaX"  ).withAggregateState(AggregateState::IonExchange ).withStandardGibbsEnergy(0.0));
    db.addSpecies( Species("NH4X" ).withName("NH4X" ).withAggregateState(AggregateState::IonExchange ).withStandardGibbsEnergy(0.0));

    Phases phases(db);

    phases.add( AqueousPhase(speciate("H O C Na Cl Ca Mg")) );
    phases.add( IonExchangePhase("NaX CaX2 KX AlX3 MgX2") );

    ChemicalSystem system(phases);

    EquilibriumSolver solver(system);

    IonExchangeProps exprops(system);

    auto phase = exprops.phase();

    const auto& species = system.species();
    const auto& elements = system.elements();
    const auto& exspecies = phase.species();
    const auto& exelements = phase.elements();

    auto num_species = species.size();

    SECTION("Testing correct initialization of the `IonExchangeProps` instance")
    {
        CHECK( phase.species().size()   == 5      );
        CHECK( phase.elements().size()  == 6      );
        CHECK( phase.species(0).name()  == "NaX"  );
        CHECK( phase.species(1).name()  == "CaX2" );
        CHECK( phase.species(2).name()  == "KX"   );
        CHECK( phase.species(3).name()  == "AlX3" );
        CHECK( phase.species(4).name()  == "MgX2" );
    }

    const real T = 25.0; // celsius
    const real P = 1.0;  // bar

    SECTION("Testing when species have zero amounts")
    {
        const ArrayXr n = ArrayXr::Zero(num_species);

        ChemicalState state(system);
        state.setTemperature(T);
        state.setPressure(P);
        state.setSpeciesAmounts(n);

        CHECK_THROWS(exprops.update(state));
    }

    SECTION("Testing when species have nonzero amounts")
    {
        ArrayXd n = ArrayXd::Ones(num_species);

        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesAmounts(n);

        exprops.update(state);

        // Check amounts of all species
        for(auto i = 0; i < exspecies.size(); i++)
            CHECK( exprops.speciesAmounts()[i] == Approx(1.00) );

        // Check amounts of all elements
        CHECK( exprops.elementAmount("Al") == Approx(1.0) );
        CHECK( exprops.elementAmount("Ca") == Approx(1.0) );
        CHECK( exprops.elementAmount("K")  == Approx(1.0) );
        CHECK( exprops.elementAmount("Na") == Approx(1.0) );
        CHECK( exprops.elementAmount("Mg") == Approx(1.0) );
        CHECK( exprops.elementAmount("X")  == Approx(9.0) );
    }

    SECTION("Testing when state is a brine")
    {
        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesMass("H2O(aq)"    , 1.00, "kg");
        state.setSpeciesAmount("Na+(aq)"  , 1.00, "mmol");
        state.setSpeciesAmount("Ca++(aq)" , 1.00, "mmol");
        state.setSpeciesAmount("Mg++(aq)" , 1.00, "mmol");
        state.setSpeciesAmount("NaX"      , 1.00, "umol");

        EquilibriumResult result = solver.solve(state);

        exprops.update(state);

        // Check molalities of all species
        CHECK( exprops.speciesAmount("NaX"  ) == Approx(1.36704e-08) );
        CHECK( exprops.speciesAmount("CaX2" ) == Approx(2.46581e-07) );
        CHECK( exprops.speciesAmount("KX"   ) == Approx(1e-16) );
        CHECK( exprops.speciesAmount("AlX3" ) == Approx(1e-16) );
        CHECK( exprops.speciesAmount("MgX2" ) == Approx(2.46581e-07) );

        // Check amounts of all elements
        CHECK( exprops.elementAmount("Na") == Approx(1.36704e-08 ) );
        CHECK( exprops.elementAmount("Mg") == Approx(2.46581e-07 ) );
        CHECK( exprops.elementAmount("Al") == Approx(1e-16       ) );
        CHECK( exprops.elementAmount("K" ) == Approx(1e-16       ) );
        CHECK( exprops.elementAmount("Ca") == Approx(2.46581e-07 ) );
        CHECK( exprops.elementAmount("X" ) == Approx(1e-06       ) );
    
        // Test convenience methods species and element amounts
        for(const auto& s : exspecies)
        {
            const auto name = s.name();
            const auto idx = exspecies.index(name);
            CHECK( exprops.speciesAmount(name) == Approx(exprops.speciesAmounts()[idx]) );
        }

        for(const auto& e : exelements)
        {
            const auto symbol = e.symbol();
            const auto idx = exelements.index(symbol);
            CHECK( exprops.elementAmount(symbol) == Approx(exprops.elementAmounts()[idx]) );
        }
    }
}
