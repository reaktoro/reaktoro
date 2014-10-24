// Reaktor is a C++ library for computational reaction modelling.
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

#include "Outputter.hpp"

namespace Reaktor {

Outputter::Outputter()
{}

void Outputter::setOptions(const OutputOptions& options_)
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
        const std::string bar(options.width, '=');

        for(const std::string& entry : entries)
            std::cout << (entry == options.separator ? options.separator : bar);
        std::cout << std::endl;

        for(const std::string& entry : entries)
            if(entry == options.separator) std::cout << options.separator;
            else std::cout << std::setw(options.width) << std::left << entry;
        std::cout << std::endl;

        for(const std::string& entry : entries)
            std::cout << (entry == options.separator ? options.separator : bar);
        std::cout << std::endl;
    }
}

void Outputter::outputState()
{
    if(options.active)
    {
        for(const std::string& val : values)
            if(val == options.separator) std::cout << options.separator;
            else std::cout << std::setw(options.width) << std::left << val;
        std::cout << std::endl;

        values.clear();
    }
}

void Outputter::outputMessage(const std::string& message)
{
    if(options.active) std::cout << message << std::endl;
}

} // namespace
