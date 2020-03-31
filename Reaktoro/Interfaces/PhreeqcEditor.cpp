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

#include "PhreeqcEditor.hpp"

// C++ includes
#include <set>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Interfaces/Phreeqc.hpp>

namespace Reaktoro {

struct PhreeqcEditor::Impl
{
    // The name of the database file
    std::string database;

    // The names of the elements used to speciate the aqueous phase
    std::vector<std::string> elements;

    // The names of the gases used to set the gaseous phase
    std::vector<std::string> gases;

    // The names of the minerals used to set the pure mineral phases
    std::vector<std::string> minerals;
};

PhreeqcEditor::PhreeqcEditor()
: pimpl(new Impl())
{
}

PhreeqcEditor::PhreeqcEditor(std::string database)
: PhreeqcEditor()
{
    pimpl->database = database;
}

PhreeqcEditor::PhreeqcEditor(const PhreeqcEditor& other)
: pimpl(new Impl(*other.pimpl))
{}

PhreeqcEditor::~PhreeqcEditor()
{}

auto PhreeqcEditor::operator=(PhreeqcEditor other) -> PhreeqcEditor&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto PhreeqcEditor::setDatabase(std::string database) -> void
{
    pimpl->database = database;
}

auto PhreeqcEditor::setAqueousPhase(StringList elements) -> void
{
    // Remove H and O from the list, since PHREEQC always assume
    // these two elements in the aqueous phase, and keeping them
    // here will cause failures
    for(auto element : elements)
        if(element != "H" && element != "O")
            pimpl->elements.push_back(element);
}

auto PhreeqcEditor::setGaseousPhase(StringList gases) -> void
{
    pimpl->gases = gases;
}

auto PhreeqcEditor::setMineralPhases(StringList minerals) -> void
{
    pimpl->minerals = minerals;
}

PhreeqcEditor::operator ChemicalSystem() const
{
    // Assert the database has been given
    Assert(!pimpl->database.empty(),
        "Could not convert the PhreeqcEditor instance to a ChemicalSystem instance.",
        "No database was provided to PhreeqcEditor via its constructor or its "
        "method PhreeqcEditor::setDatabase.");

    Phreeqc phreeqc = *this;
    return phreeqc.system();
}

PhreeqcEditor::operator Phreeqc() const
{
    // Assert the database has been given
    Assert(!pimpl->database.empty(),
        "Could not convert the PhreeqcEditor instance to a Phreeqc instance.",
        "No database was provided to PhreeqcEditor via its constructor or its "
        "method PhreeqcEditor::setDatabase.");

    // The indentation string
    const std::string indent = "    ";

    // The set of ignored elements
    const std::set<std::string> ignored_elements = {"H"};

    // Create the auxiliary input string
    std::string input;

    // Define the SOLUTION block
    input = "SOLUTION\n"
            "    units   ppm\n";
    for(auto element : pimpl->elements)
        if(!ignored_elements.count(element))
            input += indent + element + " 1.0\n";

    // Define the EQUILIBRIUM_PHASES block containing the minerals
    if(pimpl->minerals.size())
    {
        input += "EQUILIBRIUM_PHASES\n";
        for(auto mineral : pimpl->minerals)
            input += indent + mineral + " 0.0\n";
    }

    // Define the GAS_PHASE block containing the gases
    if(pimpl->gases.size())
    {
        input += "GAS_PHASE\n";
        for(auto gas : pimpl->gases)
            input += indent + gas + " 0.0\n";
    }

    // End the input script
    input += "END\n";

    // Create the Phreeqc instance with available database and create input script
    Phreeqc phreeqc(pimpl->database);
    phreeqc.execute(input);

    return phreeqc;
}

} // namespace Reaktoro

