// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

/// Check implementation below for how to create a custom database.
auto createCustomDatabase() -> Database;

int main()
{
    // Create custom thermodynamic database
    Database db = createCustomDatabase();

    // Create phases with only aqueous phase
    Phases phases(db);
    phases.add( AqueousPhase("H2O H+ OH-") );

    // Construct the chemical system
    ChemicalSystem system(phases);

    // Define initial equilibrium state
    ChemicalState state(system);
    state.setTemperature(25.0, "celsius");
    state.setPressure(1.0, "bar");
    state.setSpeciesMass("H2O", 1.0, "kg");

    // Define equilibrium solver and equilibrate given initial state
    EquilibriumSolver solver(system);
    solver.solve(state);

    // Output chemical state to the txt-file
    state.output("state.txt");

    return 0;
}

auto createCustomDatabase() -> Database
{
    //================================================================================
    // Create the elements (the constituents of species)
    //================================================================================
    ElementList elements;

    elements.append(Element()
        .withName("Hydrogen")
        .withSymbol("H")
        .withMolarMass(1.0e-3));

    elements.append(Element()
        .withName("Oxygen")
        .withSymbol("O")
        .withMolarMass(16.0e-3));

    elements.append(Element()
        .withName("Electron")
        .withSymbol("E")
        .withMolarMass(5.48e-7));

    //================================================================================
    // Create the primary species (used to define the secondary species)
    //================================================================================
    SpeciesList primaryspecies;

    primaryspecies.append(Species()
        .withName("H+")
        .withFormula("H+")
        .withElements({
            {elements.get("H"),  1.0},
            {elements.get("E"), -1.0} })
        .withAggregateState(AggregateState::Aqueous)
        .withStandardGibbsEnergy(0.0));

    primaryspecies.append(Species()
        .withName("OH-")
        .withFormula("OH-")
        .withElements({
            {elements.get("O"), 1.0},
            {elements.get("H"), 1.0},
            {elements.get("E"), 1.0} })
        .withAggregateState(AggregateState::Aqueous)
        .withStandardGibbsEnergy(0.0));

    //================================================================================
    // Create the secondary species (using formation reactions from primary species)
    //================================================================================
    SpeciesList secondaryspecies;

    secondaryspecies.append(Species()
        .withName("H2O")
        .withFormula("H2O")
        .withElements({
            {elements.get("H"), 2.0},
            {elements.get("O"), 1.0} })
        .withAggregateState(AggregateState::Aqueous)
        .withFormationReaction(
            FormationReaction()
                .withReactants({
                    {primaryspecies.get("H+"),  1.0},
                    {primaryspecies.get("OH-"), 1.0} })
                .withEquilibriumConstant(14.0)));

    //================================================================================
    // Create the database containing all primary and secondary species
    //================================================================================
    SpeciesList species = concatenate(primaryspecies, secondaryspecies);

    return Database(species);
}
