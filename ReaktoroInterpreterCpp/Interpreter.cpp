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

/// Initialize a ChemicalState instance with an Equilibrium definition.
auto initChemicalState(ChemicalState& state, const Equilibrium& e) -> void
{
    // Initialize the amounts of the inert species
    for(auto x : e.inert_species)
        state.setSpeciesAmount(x.entity, x.value, x.units);
}


/// Initialize a KineticPath instance with a Kinetics definition.
auto initKineticPath(KineticPath& path, const Kinetics& k) -> void
{
    // Auxiliary references
    ChemicalSystem system = path.system();

    // Set the partition of the chemical system in the definition of the kinetic problem
    Partition partition(system);
    partition.setKineticSpecies(k.kinetic_species);
    path.setPartition(partition);

    // Set the plots in the kinetic path problem
    for(auto p : k.plots)
    {
        auto plot = path.plot();
        plot.xdata(p.x);
        plot.ydata(p.y);
        plot.xlabel(p.xlabel);
        plot.ylabel(p.ylabel);
        plot.legend(p.ytitles);
        plot.key(p.key);
    }
}

/// Return a list of compounds names that are present as mineral species in a database.
auto filterMineralSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>
{
    std::vector<std::string> minerals;
    for(auto compound : compounds)
        if(database.containsMineralSpecies(compound))
            minerals.push_back(compound);
    return minerals;
}

/// Return a list of compounds names that are present as gaseous species in a database.
auto filterGaseousSpecies(const std::vector<std::string>& compounds, const Database& database) -> std::vector<std::string>
{
    std::vector<std::string> gases;
    for(auto compound : compounds)
        if(database.containsGaseousSpecies(compound))
            gases.push_back(compound);
    return gases;
}

} // namespace

Interpreter::Interpreter()
: database("supcrt98")
{}

Interpreter::Interpreter(std::string str)
: Interpreter()
{
    execute(str);
}

Interpreter::Interpreter(std::istream& stream)
: Interpreter()
{
    execute(stream);
}

auto Interpreter::execute(std::string str) -> void
{
    str = preprocess(str);

    Node root = YAML::Load(str);

    ProcessFunction process_equilibrium = [&](const Node& node)
    {
        Equilibrium equilibrium;
        node >> equilibrium;

        // Collect all compounds that appears in the definition of the equilibrium problem
        auto compounds = collectCompounds(equilibrium);

        // Determine if there are gaseous and mineral species among those compounds
        auto mineral_species = filterMineralSpecies(compounds, database);
        auto gaseous_species = filterGaseousSpecies(compounds, database);

        // Add the aqueous phase using the compounds as the initializer
        editor.addAqueousPhaseWithCompounds(compounds);

        // Add a gaseous phase if there were gaseous species among the compound names
        if(gaseous_species.size())
            editor.addGaseousPhaseWithSpecies(gaseous_species);

        // Add a mineral phase for each mineral species among the compound names
        for(auto x : mineral_species)
            editor.addMineralPhaseWithSpecies({x});

        // Initialize the chemical system
        system = editor;

        EquilibriumProblem problem(system);
        initEquilibriumProblem(problem, equilibrium);

        ChemicalState state(system);
        initChemicalState(state, equilibrium);

        equilibrate(state, problem);

        state.output(equilibrium.stateid + ".dat");

        states[equilibrium.stateid] = state;
    };

    ProcessFunction process_kinetics = [&](const Node& node)
    {
        Kinetics kinetics;
        node >> kinetics;

        KineticPath path(reactions);
        initKineticPath(path, kinetics);

        auto& state = states[kinetics.initial_condition];

        const auto t = kinetics.duration.value;
        const auto units = kinetics.duration.units;

        path.solve(state, 0, t, units);

        state.output(kinetics.stateid + ".dat");

        states[kinetics.stateid] = state;
    };

    std::map<std::string, ProcessFunction> fmap = {
        {"equilibrium"     , process_equilibrium},
        {"kinetics"        , process_kinetics},
    };

    for(auto child : root)
    {
        std::cout << child << std::endl;
        std::string key = lowercase(keyword(child));
        auto it = fmap.find(key);
        Assert(it != fmap.end(), "Could not parse `" + child + "`.",
            "Expecting a valid keyword. Did you misspelled it?");
        it->second(child);
    }
}

auto Interpreter::execute(std::istream& stream) -> void
{
    std::stringstream ss; ss << stream.rdbuf();
    execute(ss.str());
}

} // namespace Reaktoro
