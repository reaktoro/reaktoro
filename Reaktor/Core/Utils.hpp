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

#pragma once

// Reaktor includes
#include <Reaktor/Common/Index.hpp>

namespace Reaktor {

/// Return the names of the entries in a container.
template<typename T>
auto names(const std::vector<T>& values) -> std::vector<std::string>
{
    std::vector<std::string> names;
    names.reserve(values.size());
    for(const auto& entry : values)
        names.push_back(entry.name());
    return names;
}

/// Return the index of an entry in a container.
template<typename T>
auto index(const std::string& name, const std::vector<T>& values) -> Index
{
    Index idx = 0;
    for(const auto& value : values)
        if(value.name() == name) return idx; else ++idx;
    return idx;
}

/// Return the index of an entry in a container.
template<typename T>
auto index(const T& value, const std::vector<T>& values) -> Index
{
    return index(value.name(), values);
}

/// Return the indices of some entries in a container.
template<typename T>
auto indices(const std::vector<std::string>& names, const std::vector<T>& values) -> Indices
{
    Indices idxs;
    idxs.reserve(names.size());
    for(const auto& name : names)
        idxs.push_back(index(name, values));
    return idxs;
}

/// Return the indices of some entries in a container.
template<typename T>
auto indices(const std::vector<T>& subvalues, const std::vector<T>& values) -> Indices
{
    Indices idxs;
    idxs.reserve(subvalues.size());
    for(const auto& value : subvalues)
        idxs.push_back(index(value, values));
    return idxs;
}

} // namespace Reaktor
