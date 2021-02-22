// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Param.hpp>

namespace Reaktoro {

/// The class used to store and retrieve model parameters.
/// @ingroup Core
class Params
{
public:
    /// Construct a default Params instance.
    Params();

    /// Append a new parameter to the list of parameters.
    /// @warning A runtime error is thrown if parameter id is not unique among other parameters already in the list.
    auto append(const Param& param) -> Param&;

    /// Append a new parameter to the list of parameters with given @p id and @p value.
    /// @warning A runtime error is thrown if parameter id is not unique among other parameters already in the list.
    auto append(const String& id, const real& value) -> Param&;

    /// Return the number of child parameters.
    auto size() const -> Index;

    /// Return the parameter with given index.
    auto operator[](Index i) -> Param&;

    /// Return the parameter with given index.
    auto operator[](Index i) const -> const Param&;

    /// Return the index of the first parameter with given identifier or the number of parameters if not found.
    auto find(const String& id) const -> Index;

    /// Return the index of the first parameter with given identifier or throw a runtime error if not found.
    auto index(const String& id) const -> Index;

    /// Return the first parameter with given identifier or throw a runtime error if not found.
    auto get(const String& id) -> Param&;

    /// Return the first parameter with given identifier or throw a runtime error if not found.
    auto get(const String& id) const -> const Param&;

    /// Return true if a parameter exists with given identifier.
    auto exists(const String& id) const -> bool;

private:
    Vec<Param> m_data;

public:
    /// Construct an Params object with given begin and end iterators.
    template<typename InputIterator>
    Params(InputIterator begin, InputIterator end) : m_data(begin, end) {}

    /// Return begin const iterator of this Params instance (for STL compatibility reasons).
    auto begin() const { return m_data.begin(); }

    /// Return begin iterator of this Params instance (for STL compatibility reasons).
    auto begin() { return m_data.begin(); }

    /// Return end const iterator of this Params instance (for STL compatibility reasons).
    auto end() const { return m_data.end(); }

    /// Return end iterator of this Params instance (for STL compatibility reasons).
    auto end() { return m_data.end(); }

    /// Append a new Param at the back of the container (for STL compatibility reasons).
    auto push_back(const Param& param) -> void { append(param); }

    /// Insert a container of Param objects into this Params instance (for STL compatibility reasons).
    template<typename Iterator, typename InputIterator>
    auto insert(Iterator pos, InputIterator begin, InputIterator end) -> void { m_data.insert(pos, begin, end); }

    /// The type of the value stored in a Params (for STL compatibility reasons).
    using value_type = Param;
};

} // namespace Reaktoro
