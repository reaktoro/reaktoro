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

#include "ConvertUtils.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Kinetics/KineticPath.hpp>
#include <Reaktoro/Thermodynamics/Reactions/MineralReaction.hpp>

// Interpreter includes
#include "ParserOperators.hpp"
#include "YamlUtils.hpp"

namespace Reaktoro {

auto convertMineralReaction(const kwd::MineralReaction& node) -> MineralReaction
{
    MineralReaction reaction;
    reaction.setMineral(node.mineral);
    for(auto mechanism : node.mechanisms)
        reaction.addMechanism(mechanism);
    reaction.setEquation(node.equation);
    reaction.setSpecificSurfaceArea(node.ssa.value, node.ssa.units);
    return reaction;
}

auto convertEquilibriumProblem(const kwd::EquilibriumProblem& node, const ChemicalSystem& system) -> EquilibriumProblem
{
    // The equilibrium problem definition
    EquilibriumProblem problem(system);

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

    // Initialize the mixture conditions
    for(auto entry : node.mixture)
        problem.add(entry.entity, entry.value, entry.units);

    // Initialize the pH constraint if existent
    if(node.ph.size())
    {
        auto x = node.ph.back();
        if(x.titrant1.size() && x.titrant2.size())
            problem.pH(x.value, x.titrant1, x.titrant2);
        else if(x.titrant1.size())
            problem.pH(x.value, x.titrant1);
        else problem.pH(x.value);
    }

    // Initialize the species amount constraints if existent
    for(auto x : node.species_amounts)
        if(x.titrant1.size())
            problem.setSpeciesAmount(x.entity, x.value, x.units, x.titrant1);
        else problem.setSpeciesAmount(x.entity, x.value, x.units);

    // Initialize the species activity constraints if existent
    for(auto x : node.species_activities)
        if(x.titrant1.size() && x.titrant2.size())
            problem.setSpeciesActivity(x.entity, x.value, x.titrant1, x.titrant2);
        else if(x.titrant1.size())
            problem.setSpeciesActivity(x.entity, x.value, x.titrant1);
        else problem.setSpeciesActivity(x.entity, x.value);

    // Initialize the phase amount constraints if existent
    for(auto x : node.phase_amounts)
        if(x.titrant1.size())
            problem.setPhaseAmount(x.entity, x.value, x.units, x.titrant1);
        else if(system.phase(x.entity).numSpecies() == 1.0)
            problem.setPhaseAmount(x.entity, x.value, x.units,
                system.phase(x.entity).species(0).name());
        else RuntimeError("Could not construct the equilibrium problem with "
            "given `PhaseAmount` constraint for multi-component phase `" + x.entity + "`.",
            "Expecting a titrant to control this phase amount, since a default titrant cannot "
            "be determined like it happens for a single-component phase.");

    // Initialize the phase volume constraints if existent
    for(auto x : node.phase_volumes)
        if(x.titrant1.size())
            problem.setPhaseVolume(x.entity, x.value, x.units, x.titrant1);
        else if(system.phase(x.entity).numSpecies() == 1.0)
            problem.setPhaseVolume(x.entity, x.value, x.units,
                system.phase(x.entity).species(0).name());
        else RuntimeError("Could not construct the equilibrium problem with "
            "given `PhaseVolume` constraint for multi-component phase `" + x.entity + "`.",
            "Expecting a titrant to control this phase volume, since a default titrant cannot "
            "be determined like it happens for a single-component phase.");

    return problem;
}

auto convertKineticPath(const kwd::KineticPath& node, const ReactionSystem& reactions) -> KineticPath
{
    // Initialize kinetic path object
    KineticPath path(reactions);

    // Auxiliary references
    ChemicalSystem system = reactions.system();

    // Set the partition of the chemical system in the definition of the kinetic problem
    Partition partition(system);
    partition.setKineticSpecies(node.kinetic_species);
    path.setPartition(partition);

    // Set the plots in the kinetic path problem
    for(auto p : node.plots)
    {
        auto plot = path.plot();
        plot.xdata(p.x);
        plot.ydata(p.y);
        plot.xlabel(p.xlabel);
        plot.ylabel(p.ylabel);
        plot.legend(p.ytitles);
        plot.key(p.key);
    }

    return path;
}

} // namespace Reaktoro
