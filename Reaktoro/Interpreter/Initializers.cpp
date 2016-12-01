// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "Initializers.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPath.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Kinetics/KineticPath.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralReaction.hpp>

// Reaktoro includes
#include "Keywords.hpp"
#include "Yaml.hpp"

namespace Reaktoro {

auto initializeChemicalPlot(ChemicalPlot& plot, const kwd::Plot& keyword) -> void
{
    // Set the other attributes of the plot
    plot.name(keyword.name);
    plot.x(keyword.x);
//    plot.y(StringList(keyword.y));
//    plot.legend(StringList(keyword.ytitles));
    plot.xlabel(keyword.xlabel);
    plot.ylabel(keyword.ylabel);
    plot.legend(keyword.key);
}

auto initializeMineralReaction(MineralReaction& reaction, const kwd::MineralReaction& node) -> void
{
    reaction.setMineral(node.mineral);
    for(auto mechanism : node.mechanisms)
        reaction.addMechanism(mechanism);
    reaction.setEquation(node.equation);
    reaction.setSpecificSurfaceArea(node.ssa.value, node.ssa.units);
}

auto initializeEquilibriumProblem(EquilibriumProblem& problem, const kwd::EquilibriumProblem& node) -> void
{
    // Auxiliary references
    const ChemicalSystem& system = problem.system();

    // Initialize temperature and pressure
    problem.setTemperature(node.temperature.value, node.temperature.units);
    problem.setPressure(node.pressure.value, node.pressure.units);

    // Collect the names of all species to be inert, including those in inert phases
    std::vector<std::string> inert_species;
    for(auto x : node.inert_species)
        inert_species.push_back(x.entity);
    for(auto x : node.inert_phases)
        for(auto s : system.phase(x).species())
            inert_species.push_back(s.name());

    // Initialize the partition with the inert species
    Partition partition(problem.system());
    partition.setInertSpecies(inert_species);
    problem.setPartition(partition);

    // Initialize the recipe conditions
    for(auto entry : node.recipe)
        problem.add(entry.entity, entry.value, entry.units);

//    // Initialize the pH constraint if existent
//    if(node.ph.size())
//    {
//        auto x = node.ph.back();
//        if(x.titrant1.size() && x.titrant2.size())
//            problem.pH(x.value, x.titrant1, x.titrant2);
//        else if(x.titrant1.size())
//            problem.pH(x.value, x.titrant1);
//        else problem.pH(x.value);
//    }
//
//    // Initialize the species amount constraints if existent
//    for(auto x : node.species_amounts)
//        if(x.titrant1.size())
//            problem.fixSpeciesAmount(x.entity, x.value, x.units, x.titrant1);
//        else problem.fixSpeciesAmount(x.entity, x.value, x.units);
//
//    // Initialize the species activity constraints if existent
//    for(auto x : node.species_activities)
//        if(x.titrant1.size() && x.titrant2.size())
//            problem.fixSpeciesActivity(x.entity, x.value, x.titrant1, x.titrant2);
//        else if(x.titrant1.size())
//            problem.fixSpeciesActivity(x.entity, x.value, x.titrant1);
//        else problem.fixSpeciesActivity(x.entity, x.value);
//
//    // Initialize the species fugacity constraints if existent
//    for(auto x : node.species_fugacities)
//        if(x.titrant1.size())
//            problem.setSpeciesFugacity(x.entity, x.value, x.units, x.titrant1);
//        else problem.setSpeciesFugacity(x.entity, x.value, x.units);
//
//    // Initialize the phase amount constraints if existent
//    for(auto x : node.phase_amounts)
//        if(x.titrant1.size())
//            problem.fixPhaseAmount(x.entity, x.value, x.units, x.titrant1);
//        else if(system.phase(x.entity).numSpecies() == 1.0)
//            problem.fixPhaseAmount(x.entity, x.value, x.units,
//                system.phase(x.entity).species(0).name());
//        else RuntimeError("Could not construct the equilibrium problem with "
//            "given `PhaseAmount` constraint for multi-component phase `" + x.entity + "`.",
//            "Expecting a titrant to control this phase amount, since a default titrant cannot "
//            "be determined like it happens for a single-component phase.");
//
//    // Initialize the phase volume constraints if existent
//    for(auto x : node.phase_volumes)
//        if(x.titrant1.size())
//            problem.fixPhaseVolume(x.entity, x.value, x.units, x.titrant1);
//        else if(system.phase(x.entity).numSpecies() == 1.0)
//            problem.fixPhaseVolume(x.entity, x.value, x.units,
//                system.phase(x.entity).species(0).name());
//        else RuntimeError("Could not construct the equilibrium problem with "
//            "given `PhaseVolume` constraint for multi-component phase `" + x.entity + "`.",
//            "Expecting a titrant to control this phase volume, since a default titrant cannot "
//            "be determined like it happens for a single-component phase.");
}

auto initializeEquilibriumPath(EquilibriumPath& path, const kwd::EquilibriumPath& keyword) -> void
{
    // Auxiliary references
    const ChemicalSystem& system = path.system();

    // Set the inert species in the equilibrium path calculation
    Partition partition(system);
    partition.setInertSpecies(keyword.inert_species);
    path.setPartition(partition);

    // Set the plots in the kinetic path problem
    for(auto p : keyword.plots)
    {
        auto plot = path.plot();
        initializeChemicalPlot(plot, p);
    }
}

auto initializeKineticPath(KineticPath& path, const kwd::KineticPath& keyword) -> void
{
    // Auxiliary references
    const ChemicalSystem& system = path.system();

    // Set the inert and kinetic species in the kinetic path calculation
    Partition partition(system);
    partition.setInertSpecies(keyword.inert_species);
    partition.setKineticSpecies(keyword.kinetic_species);
    path.setPartition(partition);

    // Set the plots in the kinetic path problem
    for(auto p : keyword.plots)
    {
        auto plot = path.plot();
        initializeChemicalPlot(plot, p);
    }
}

} // namespace Reaktoro
