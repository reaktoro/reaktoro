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
#include <cmath>

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

    /// The name of the output file.
    std::string filename;

    /// The base name of the output file.
    std::string basename;

    /// The extension name of the output file.
    std::string extension;

    /// The suffix word of the output file name.
    std::string suffix;

    /// The names of the quantities to be output.
    std::vector<std::string> data;

    /// The names of the quantities to appear as column header in the output.
    std::vector<std::string> headings;

    /// The floating-point precision in the output.
    int precision = 6;

    /// The flag that indicates if scientific format should be used.
    bool scientific = false;

    /// The output stream of the data file.
    std::ofstream datafile;

    /// The iteration number for every update call
    Index iteration = 0;

    /// The extra attachments to the output file.
    std::vector<std::string> attachments;

    /// The index of the active column
    Index icolumn = 0;

    /// The spacings between the columns
    std::vector<int> spacings;

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

    auto spacing(std::string word) const -> std::size_t
    {
        return word.size() + std::max(5, 20 - static_cast<int>(word.size()));
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

        // Check if scientific format should be used
        if(scientific)
            datafile << std::scientific;

        // Set the floating-point precision in the output.
        datafile << std::setprecision(precision);

        // Determine the spacings between the columns
        spacings.clear();
        for(auto word : headings)
            spacings.push_back(spacing(word));

        // Output the header of the data file
        icolumn = 0;
        for(auto word : headings)
        {
            auto space = spacings[icolumn];
            if(datafile.is_open()) datafile << std::left << std::setw(space) << word;
            if(terminal)
            {
                std::ios::fmtflags flags(std::cout.flags());
                if(scientific) std::cout << std::scientific;
                std::cout << std::setprecision(precision);
                std::cout << std::left << std::setw(space) << word;
                std::cout.flags(flags);
            }
            ++icolumn;
        }
    }

    auto close() -> void
    {
        datafile.close();
    }

    auto update(const ChemicalState& state, double t) -> void
    {
        // Output values on a new line
        if(datafile.is_open()) datafile << std::endl;
        if(terminal) std::cout << std::endl;

        // Output the current chemical state to the data file.
        quantity.update(state, t);

        // For each quantity, ouput its value on each column
        icolumn = 0;
        for(auto word : data)
        {
            auto space = spacings[icolumn];
            auto val = (word == "i") ? iteration : quantity.value(word);
            if(datafile.is_open()) datafile << std::left << std::setw(space) << val;
            if(terminal) std::cout << std::left << std::setw(space) << val;
            ++icolumn;
        }

        // Update the iteration number
        ++iteration;
    }

    template<typename ValueType>
    auto attach(ValueType value) -> void
    {
        auto space = spacings[icolumn];
        if(datafile.is_open()) datafile << std::left << std::setw(space) << value;
        if(terminal) std::cout << std::left << std::setw(space) << value;
        ++icolumn;
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

auto ChemicalOutput::filename(std::string filename) -> void
{
    pimpl->filename = filename;
    pimpl->basename = filename.substr(0, filename.rfind('.'));
    pimpl->extension = filename.substr(filename.rfind('.'));
}

auto ChemicalOutput::filename() const -> std::string
{
    return pimpl->filename;
}

auto ChemicalOutput::suffix(std::string word) -> void
{
    pimpl->suffix = word;
    pimpl->filename = pimpl->basename + pimpl->suffix + pimpl->extension;
}

auto ChemicalOutput::suffix() const -> std::string
{
    return pimpl->suffix;
}

auto ChemicalOutput::basename() const -> std::string
{
    return pimpl->basename;
}

auto ChemicalOutput::extension() const -> std::string
{
    return pimpl->extension;
}

auto ChemicalOutput::add(std::string quantity) -> void
{
    add(quantity, quantity);
}

auto ChemicalOutput::add(std::string quantity, std::string label) -> void
{
    pimpl->headings.insert(pimpl->headings.end() - pimpl->attachments.size(), label);
    pimpl->data.push_back(quantity);
}

auto ChemicalOutput::attachments(std::vector<std::string> titles) -> void
{
    pimpl->headings.insert(pimpl->headings.end(), titles.begin(), titles.end());
    pimpl->attachments.insert(pimpl->attachments.end(), titles.begin(), titles.end());
}

auto ChemicalOutput::attach(int value) -> void
{
    pimpl->attach(value);
}

auto ChemicalOutput::attach(double value) -> void
{
    pimpl->attach(value);
}

auto ChemicalOutput::attach(std::string value) -> void
{
    pimpl->attach(value);
}

auto ChemicalOutput::precision(int val) -> void
{
    pimpl->precision = std::abs(val);
}

auto ChemicalOutput::scientific(bool enable) -> void
{
    pimpl->scientific = enable;
}

auto ChemicalOutput::terminal(bool enabled) -> void
{
    pimpl->terminal = enabled;
}

auto ChemicalOutput::quantities() const -> std::vector<std::string>
{
    return pimpl->data;
}

auto ChemicalOutput::headings() const -> std::vector<std::string>
{
    return pimpl->headings;
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
