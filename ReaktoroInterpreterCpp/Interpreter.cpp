// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "Interpreter.hpp"

// ReaktoroInterpreter includes
#include <ReaktoroInterpreterCpp/AutomationUtils.hpp>
#include <ReaktoroInterpreterCpp/ParserUtils.hpp>

namespace Reaktoro {
namespace {

/// Initialize an EquilibriumProblem instance with an Equilibrium definition.
auto initEquilibriumProblem(EquilibriumProblem& problem, const Equilibrium& e) -> void
{
    // Auxiliary references
    ChemicalSystem system = problem.system();

    // Initialize temperature and pressure
    problem.setTemperature(e.temperature.value, e.temperature.units);
    problem.setPressure(e.pressure.value, e.pressure.units);

    // Collect the names of all species to be inert, including those in inert phases
    std::vector<std::string> inert_species;
    for(auto x : e.inert_species)
        inert_species.push_back(x.entity);
    for(auto x : e.inert_phases)
        for(auto s : system.phase(x).species())
            inert_species.push_back(s.name());

    // Initialize the partition with the inert species
    Partition partition(problem.system());
    partition.setInertSpecies(inert_species);
    problem.setPartition(partition);

    // Initialize the mixture conditions
    for(auto x : e.mixture)
        problem.add(x.entity, x.value, x.units);

    // Initialize the pH constraint if existent
    if(e.pH.size())
    {
        auto x = e.pH.back();
        if(x.titrant1.size() && x.titrant2.size())
            problem.pH(x.value, x.titrant1, x.titrant2);
        else if(x.titrant1.size())
            problem.pH(x.value, x.titrant1);
        else problem.pH(x.value);
    }

    // Initialize the species amount constraints if existent
    for(auto x : e.species_amounts)
        if(x.titrant1.size())
            problem.setSpeciesAmount(x.entity, x.value, x.units, x.titrant1);
        else problem.setSpeciesAmount(x.entity, x.value, x.units);

    // Initialize the species activity constraints if existent
    for(auto x : e.species_activities)
        if(x.titrant1.size() && x.titrant2.size())
            problem.setSpeciesActivity(x.entity, x.value, x.titrant1, x.titrant2);
        else if(x.titrant1.size())
            problem.setSpeciesActivity(x.entity, x.value, x.titrant1);
        else problem.setSpeciesActivity(x.entity, x.value);

    // Initialize the phase amount constraints if existent
    for(auto x : e.phase_amounts)
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
    for(auto x : e.phase_volumes)
        if(x.titrant1.size())
            problem.setPhaseVolume(x.entity, x.value, x.units, x.titrant1);
        else if(system.phase(x.entity).numSpecies() == 1.0)
            problem.setPhaseVolume(x.entity, x.value, x.units,
                system.phase(x.entity).species(0).name());
        else RuntimeError("Could not construct the equilibrium problem with "
            "given `PhaseVolume` constraint for multi-component phase `" + x.entity + "`.",
            "Expecting a titrant to control this phase volume, since a default titrant cannot "
            "be determined like it happens for a single-component phase.");
}

} // namespace

Interpreter::Interpreter()
{}

Interpreter::Interpreter(std::string str)
{
    execute(str);
}

Interpreter::Interpreter(std::istream& stream)
{
    execute(stream);
}

auto Interpreter::execute(std::string str) -> void
{
    str = preprocess(str);

    Node node = YAML::Load(str);
    Equilibrium equilibrium;
    node[0] >> equilibrium;

    auto compounds = collectCompounds(equilibrium);

    editor.addAqueousPhaseWithCompounds(compounds);

    system = editor;

    EquilibriumProblem problem(system);
    initEquilibriumProblem(problem, equilibrium);

    states[equilibrium.stateid] = equilibrate(problem);

    states[equilibrium.stateid].output(equilibrium.stateid + ".dat");
}

auto Interpreter::execute(std::istream& stream) -> void
{
    std::stringstream ss; ss << stream.rdbuf();
    execute(ss.str());
}

} // namespace Reaktoro
