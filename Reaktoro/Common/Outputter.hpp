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

#pragma once

// C++ includes
#include <iomanip>
#include <iostream>
#include <list>
#include <sstream>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

struct OutputterOptions
{
    bool active = false;

    bool fixed = false;

    bool scientific = false;

    unsigned precision = 6;

    unsigned width = 15;

    std::string separator = "|";

    /// Assign a boolean value to `active` member.
    auto operator=(bool active) -> OutputterOptions&;
};

/// A utility class for printing the progress of iterative calculations
class Outputter
{
  public:
    struct Options;

    Outputter();

    void setOptions(const OutputterOptions& options);

    void addEntry(const std::string& name);

    void addEntries(const std::string& prefix, unsigned size);

    void addEntries(const std::string& prefix, unsigned size, const std::vector<std::string>& names);

    template<typename Iter>
    void addEntries(const Iter& begin, const Iter& end);

    void addEntrySeparator();

    template<typename T>
    void addValue(const T& val);

    template<typename Iter>
    void addValues(const Iter& begin, const Iter& end);

    template<typename Vec>
    void addValues(const Vec& vec);

    void addValueSeparator();

    void outputHeader();

    void outputState();

    void outputMessage(const std::string& message);

    template<typename T>
    void outputMessage(const T& arg)
    {
        if(options.active)
            std::cout << arg;
    }

    template<typename T, typename... Args>
    void outputMessage(const T& arg, const Args&... args)
    {
        if(options.active) {
            std::cout << arg;
            outputMessage(args...);
        }
    }

  private:
    std::list<std::string> entries;

    std::list<std::string> values;

    OutputterOptions options;
};

template<typename Iter>
void Outputter::addEntries(const Iter& begin, const Iter& end)
{
    entries.insert(entries.end(), begin, end);
}

template<typename T>
void Outputter::addValue(const T& val)
{
    std::stringstream ss;
    ss << std::setprecision(options.precision);
    if(options.fixed)
        ss << std::fixed;
    if(options.scientific)
        ss << std::scientific;
    ss << val;
    values.push_back(ss.str());
}

template<typename Iter>
void Outputter::addValues(const Iter& begin, const Iter& end)
{
    for(Iter iter = begin; iter != end; ++iter)
        AddValue(*iter);
}

template<typename Vec>
void Outputter::addValues(const Vec& vec)
{
    for(unsigned i = 0; i < vec.size(); ++i)
        addValue(vec[i]);
}

} // namespace Reaktoro
