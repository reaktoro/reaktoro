// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "Interpreter.hpp"

// Reaktoro includes
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumUtils.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>
#include <Reaktoro/Thermodynamics/Core/Database.hpp>
#include <fstream>

namespace Reaktoro {

struct Interpreter::Impl
{
    ChemicalSystem system;

    std::map<std::string, ChemicalState> states;

    auto execute(json input) -> void
    {
        initializeChemicalSystem(input["system"]);
        executeCalculations(input["calculations"]);
    }

    auto initializeChemicalSystem(json node) -> void
    {
        auto dbname = node["database"].get<std::string>();
        auto elements = node["elements"].get<std::vector<std::string>>();

        Database database(dbname);
        ChemicalEditor editor(database);
        editor.initializePhasesWithElements(elements);

        system = ChemicalSystem(editor);
    }

    auto executeCalculations(json node) -> void
    {
        for(auto item : node) {
            if(item.count("equilibrium"))
                calculateEquilibrium(item["equilibrium"]);
        }
    }

    auto calculateEquilibrium(json node) -> void
    {
        EquilibriumProblem problem(system);
        problem.setTemperature(node["temperature"]["value"], node["temperature"]["units"]);
        problem.setPressure(node["pressure"]["value"], node["pressure"]["units"]);
        for(auto item : node["substances"])
            problem.add(item["substance"], item["quantity"], item["units"]);

        ChemicalState state = equilibrate(problem);

        std::string stateref = "default";
        if(node.count("stateReference"))
            stateref = node["stateReference"];

        states.insert({stateref, equilibrate(problem)});
    }
};

Interpreter::Interpreter()
    : pimpl(new Interpreter::Impl())
{}

Interpreter::Interpreter(const Interpreter& other)
    : pimpl(new Interpreter::Impl(*other.pimpl))
{}

Interpreter::~Interpreter()
{}

auto Interpreter::operator=(Interpreter other) -> Interpreter&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Interpreter::executeJsonObject(json input) -> void
{
    pimpl->execute(input);
}

auto Interpreter::executeJsonString(std::string input) -> void
{
    json jsoninput = json::parse(input);
    executeJsonObject(jsoninput);
}

auto Interpreter::executeJsonFile(std::string input) -> void
{
    json jsoninput;
    std::ifstream file(input, std::ios_base::in);
    file >> jsoninput;
    executeJsonObject(jsoninput);
}

auto Interpreter::system() -> const ChemicalSystem&
{
    return pimpl->system;
}

auto Interpreter::states() -> const std::map<std::string, ChemicalState>&
{
    return pimpl->states;
}

auto Interpreter::state(std::string reference) -> const ChemicalState&
{
    return pimpl->states.at(reference);
}

} // namespace Reaktoro
