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

#include "ChemicalOutput.hpp"

// C++ includes
#include <fstream>
#include <iomanip>
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Units.hpp>
#include <Reaktoro/Core/ChemicalQuantity.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

struct ChemicalOutput::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The reaction system instance
    ReactionSystem reactions;

    /// The chemical quantity instance
    ChemicalQuantity quantity;

    /// The flag that indicates if output should be done at the terminal.
    bool terminal = false;

    /// The name of the file to which the output should be written.
    std::string filename;

    /// The names of the quantities to be output.
    std::vector<std::string> data;

    /// The names of the quantities to appear as column header in the output.
    std::vector<std::string> headings;

    /// The output stream of the data file.
    std::ofstream datafile;

    /// The counter for every update call
    Index counter = 0;

    Impl()
    {}

    Impl(const ChemicalSystem& system)
    : system(system), quantity(system)
    {}

    Impl(const ReactionSystem& reactions)
    : system(reactions.system()), reactions(reactions), quantity(reactions)
    {}

    ~Impl()
    {
        close();
    }

    auto open() -> void
    {
        // Ensure the output file is closed
        close();

        // Ensure output is done either to a file and/or terminal
        Assert(!filename.empty() || terminal,
            "Cannot open the ChemicalOutput instance for output.",
            "The instance has not been configured to output to the terminal or file.");

        // Make sure header is not empty
        if(headings.empty())
            headings = data;

        // Open the data file
        if(!filename.empty())
            datafile.open(filename, std::ofstream::out | std::ofstream::trunc);

        // Output the header of the data file
        for(auto word : headings)
        {
            if(datafile.is_open()) datafile << std::left << std::setw(20) << word;
            if(terminal) std::cout << std::left << std::setw(20) << word;
        }
        if(datafile.is_open()) datafile << std::endl;
        if(terminal) std::cout << std::endl;
    }

    auto close() -> void
    {
        datafile.close();
    }

    auto update(const ChemicalState& state, double t) -> void
    {
        // Update the counter
        ++counter;

        // Output the current chemical state to the data file.
        quantity.update(state, t);
        for(auto word : data)
        {
            auto val = word == "i" ? counter : quantity.value(word);
            if(datafile.is_open()) datafile << std::left << std::setw(20) << val;
            if(terminal) std::cout << std::left << std::setw(20) << val;
        }
        if(datafile.is_open()) datafile << std::endl;
        if(terminal) std::cout << std::endl;
    }
};

ChemicalOutput::ChemicalOutput()
: pimpl(new Impl())
{}

ChemicalOutput::ChemicalOutput(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ChemicalOutput::ChemicalOutput(const ReactionSystem& reactions)
: pimpl(new Impl(reactions))
{}

ChemicalOutput::~ChemicalOutput()
{}

auto ChemicalOutput::file(std::string filename) -> void
{
    pimpl->filename = filename;
}

auto ChemicalOutput::data(const StringList& quantities) -> void
{
    pimpl->data = quantities.strings();
}

auto ChemicalOutput::headings(const StringList& titles) -> void
{
    pimpl->headings = titles.strings();
}

auto ChemicalOutput::terminal(bool enabled) -> void
{
    pimpl->terminal = enabled;
}

auto ChemicalOutput::open() -> void
{
    pimpl->open();
}

auto ChemicalOutput::update(const ChemicalState& state, double t) -> void
{
    pimpl->update(state, t);
}

auto ChemicalOutput::close() -> void
{
    pimpl->close();
}

ChemicalOutput::operator bool() const
{
    return pimpl->terminal || pimpl->filename.size();
}

} // namespace Reaktoro
