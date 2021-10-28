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
#include <Reaktoro/Thermodynamics/Surface/ActivityModelIonExchange.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/ActivityModelHKF.hpp>
#include <Reaktoro/Singletons/Elements.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>

using namespace Reaktoro;

namespace test {

    extern auto createDatabase() -> Database;

    auto getPhreeqcDatabase(const String& name) -> PhreeqcDatabase;

} // namespace test

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
    phases.add( IonExchangePhase("NaX CaX2 KX AlX3 MgX2").setActivityModel(ActivityModelIonExchange()) );

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

        // Check equivalences of all species
        CHECK( exprops.speciesEquivalence("NaX"  ) == Approx(1.0) );
        CHECK( exprops.speciesEquivalence("CaX2" ) == Approx(2.0) );
        CHECK( exprops.speciesEquivalence("KX"   ) == Approx(1.0) );
        CHECK( exprops.speciesEquivalence("AlX3" ) == Approx(3.0) );
        CHECK( exprops.speciesEquivalence("MgX2" ) == Approx(2.0) );

        // Check equivalent fractions of all species
        CHECK( exprops.speciesEquivalentFraction("NaX"  ) == Approx(0.111111) );
        CHECK( exprops.speciesEquivalentFraction("CaX2" ) == Approx(0.222222) );
        CHECK( exprops.speciesEquivalentFraction("KX"   ) == Approx(0.111111) );
        CHECK( exprops.speciesEquivalentFraction("AlX3" ) == Approx(0.333333) );
        CHECK( exprops.speciesEquivalentFraction("MgX2" ) == Approx(0.222222) );

        // Check log10 gammas of all species
        for(auto i = 0; i < exspecies.size(); i++)
            CHECK( exprops.speciesActivityCoefficientsLg()[i] == Approx(0.00) );
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

        // Check equivalences of all species
        CHECK( exprops.speciesEquivalence("NaX"  ) == Approx(1.36704e-08 ) );
        CHECK( exprops.speciesEquivalence("CaX2" ) == Approx(4.93161e-07 ) );
        CHECK( exprops.speciesEquivalence("KX"   ) == Approx(1e-16       ) );
        CHECK( exprops.speciesEquivalence("AlX3" ) == Approx(3e-16       ) );
        CHECK( exprops.speciesEquivalence("MgX2" ) == Approx(4.93161e-07 ) );

        // Check equivalent fractions of all species
        CHECK( exprops.speciesEquivalentFraction("NaX"  ) == Approx(0.0136704 ) );
        CHECK( exprops.speciesEquivalentFraction("CaX2" ) == Approx(0.493161  ) );
        CHECK( exprops.speciesEquivalentFraction("KX"   ) == Approx(1e-10     ) );
        CHECK( exprops.speciesEquivalentFraction("AlX3" ) == Approx(3e-10     ) );
        CHECK( exprops.speciesEquivalentFraction("MgX2" ) == Approx(0.493161  ) );

        // Check log10 gammas of all species
        for(auto i = 0; i < exspecies.size(); i++)
            CHECK( exprops.speciesActivityCoefficientsLg()[i] == Approx(0.00) ); // because we are working on the custom library
    
        // Test convenience methods species amounts, equivalences, equivalent fractions, log10gamma and element amounts
        for(const auto& s : exspecies)
        {
            const auto name = s.name();
            const auto idx = exspecies.index(name);

            CHECK( exprops.speciesAmount(name)             == Approx(exprops.speciesAmounts()[idx]) );
            CHECK( exprops.speciesEquivalence(name)        == Approx(exprops.speciesEquivalences()[idx]) );
            CHECK( exprops.speciesEquivalentFraction(name) == Approx(exprops.speciesEquivalentFractions()[idx]) );
            CHECK(exprops.speciesActivityCoefficientLg(name) == Approx(exprops.speciesActivityCoefficientsLg()[idx]) );

            CHECK( exprops.speciesAmount(idx)             == Approx(exprops.speciesAmounts()[idx]) );
            CHECK( exprops.speciesEquivalence(idx)        == Approx(exprops.speciesEquivalences()[idx]) );
            CHECK( exprops.speciesEquivalentFraction(idx) == Approx(exprops.speciesEquivalentFractions()[idx]) );
            CHECK(exprops.speciesActivityCoefficientLg(idx) == Approx(exprops.speciesActivityCoefficientsLg()[idx]) );
        }

        for(const auto& e : exelements)
        {
            const auto symbol = e.symbol();
            const auto idx = exelements.index(symbol);

            CHECK( exprops.elementAmount(symbol) == Approx(exprops.elementAmounts()[idx]) );
            CHECK( exprops.elementAmount(idx) == Approx(exprops.elementAmounts()[idx]) );
        }
    }

    SECTION("Testing when state is a brine and database is initialized by phreeqc.dat")
    {
        // Initialize the database corresponding to the string `phreeqc.dat` has been already initialized
        auto dbphreeqc = test::getPhreeqcDatabase("phreeqc.dat");

        Phases phases(dbphreeqc);

        phases.add( AqueousPhase(speciate("H O C Na Cl Ca Mg")).setActivityModel(ActivityModelHKF()) );
        phases.add( IonExchangePhase("NaX CaX2 KX AlX3 MgX2").setActivityModel(ActivityModelIonExchange()) );

        ChemicalSystem system(phases);

        EquilibriumSolver solver(system);

        IonExchangeProps exprops(system);

        auto phase = exprops.phase();

        const auto& species = system.species();
        const auto& elements = system.elements();
        const auto& exspecies = phase.species();
        const auto& exelements = phase.elements();

        auto num_species = species.size();

        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesMass("H2O"   , 1.00, "kg");
        state.setSpeciesAmount("Na+" , 1.00, "mmol");
        state.setSpeciesAmount("Ca+2", 1.00, "mmol");
        state.setSpeciesAmount("Mg+2", 1.00, "mmol");
        state.setSpeciesAmount("NaX" , 1.00, "umol");

        EquilibriumResult result = solver.solve(state);

        exprops.update(state);

        // Check molalities of all species
        CHECK( exprops.speciesAmount("NaX" ) == Approx(9.84068e-09 ) );
        CHECK( exprops.speciesAmount("CaX2") == Approx(3.03997e-07 ) );
        CHECK( exprops.speciesAmount("KX"  ) == Approx(1e-16       ) );
        CHECK( exprops.speciesAmount("AlX3") == Approx(1e-16       ) );
        CHECK( exprops.speciesAmount("MgX2") == Approx(1.91083e-07 ) );

        // Check amounts of all elements
        CHECK( exprops.elementAmount("X" ) == Approx(1e-06       ) );
        CHECK( exprops.elementAmount("Na") == Approx(9.84068e-09 ) );
        CHECK( exprops.elementAmount("Mg") == Approx(1.91083e-07 ) );
        CHECK( exprops.elementAmount("Al") == Approx(1e-16       ) );
        CHECK( exprops.elementAmount("K" ) == Approx(1e-16       ) );
        CHECK( exprops.elementAmount("Ca") == Approx(3.03997e-07 ) );

        // Check equivalences of all species
        CHECK( exprops.speciesEquivalence("NaX" ) == Approx(9.84068e-09 ) );
        CHECK( exprops.speciesEquivalence("CaX2") == Approx(6.07993e-07 ) );
        CHECK( exprops.speciesEquivalence("KX"  ) == Approx(1e-16       ) );
        CHECK( exprops.speciesEquivalence("AlX3") == Approx(3e-16       ) );
        CHECK( exprops.speciesEquivalence("MgX2") == Approx(3.82166e-07 ) );

        // Check equivalent fractions of all species
        CHECK( exprops.speciesEquivalentFraction("NaX" ) == Approx(0.00984068 ) );
        CHECK( exprops.speciesEquivalentFraction("CaX2") == Approx(0.607993   ) );
        CHECK( exprops.speciesEquivalentFraction("KX"  ) == Approx(1e-10      ) );
        CHECK( exprops.speciesEquivalentFraction("AlX3") == Approx(3e-10      ) );
        CHECK( exprops.speciesEquivalentFraction("MgX2") == Approx(0.382166   ) );

        // Check log10 gammas of all species
        CHECK(exprops.speciesActivityCoefficientLg("NaX") == Approx(-0.0311024 ) );
        CHECK(exprops.speciesActivityCoefficientLg("CaX2") == Approx(-0.122843  ) );
        CHECK(exprops.speciesActivityCoefficientLg("KX") == Approx(-0.0317777 ) );
        CHECK(exprops.speciesActivityCoefficientLg("AlX3") == Approx(-0.257599  ) );
        CHECK(exprops.speciesActivityCoefficientLg("MgX2") == Approx(-0.12147   ) );

    }
}
