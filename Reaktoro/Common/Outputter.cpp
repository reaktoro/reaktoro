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

#include "Outputter.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {
namespace {

/// Return the width for std::setw to be used in Outputter
auto colwidth(unsigned width, const std::string& str) -> unsigned
{
    return width > str.size() + 4 ? width : str.size() + 4;
}

/// Return the bar string to be used in Outputter
auto barstr(unsigned width, const std::string& str) -> std::string
{
    return std::string(colwidth(width, str), '=');
}

} // namespace

auto OutputterOptions::operator=(bool act) -> OutputterOptions&
{
    active = act;
    return *this;
}

Outputter::Outputter()
{}

void Outputter::setOptions(const OutputterOptions& options_)
{
    options = options_;
}

void Outputter::addEntry(const std::string& name)
{
    entries.push_back(name);
}

void Outputter::addEntries(const std::string& prefix, unsigned size)
{
    for(unsigned i = 0; i < size; ++i)
    {
        std::stringstream ss;
        ss << prefix << "[" << i << "]";
        addEntry(ss.str());
    }
}

void Outputter::addEntries(const std::string& prefix, unsigned size, const std::vector<std::string>& names)
{
    if(names.size() != size)
        addEntries(prefix, size);
    else
    {
        for(std::string name : names)
        {
            std::stringstream ss;
            ss << prefix << "[" << name << "]";
            addEntry(ss.str());
        }
    }
}

void Outputter::addEntrySeparator()
{
    entries.push_back(options.separator);
}

void Outputter::addValueSeparator()
{
    values.push_back(options.separator);
}

void Outputter::outputHeader()
{
    if(options.active)
    {
        for(const std::string& entry : entries)
            std::cout << (entry == options.separator ? options.separator : barstr(options.width, entry));
        std::cout << std::endl;

        for(const std::string& entry : entries)
            if(entry == options.separator) std::cout << options.separator;
            else std::cout << std::setw(colwidth(options.width, entry)) << std::left << entry;
        std::cout << std::endl;

        for(const std::string& entry : entries)
            std::cout << (entry == options.separator ? options.separator : barstr(options.width, entry));
        std::cout << std::endl;
    }
}

void Outputter::outputState()
{
    if(options.active)
    {
        Assert(entries.size() == values.size(), "Could not output the state of the calculation.",
            "There are more entry names than values.");
        auto entry = entries.begin();
        for(const std::string& val : values)
            if(val == options.separator) std::cout << options.separator;
            else std::cout << std::setw(colwidth(options.width, *entry++)) << std::left << val;
        std::cout << std::endl;

        values.clear();
    }
}

void Outputter::outputMessage(const std::string& message)
{
    if(options.active) std::cout << message << std::endl;
}

} // namespace
