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

    /// Construct a default Params instance.
    Params(const std::initializer_list<Param>& params);

    /// Return a deep copy of this Params object.
    auto clone() const -> Params;

    /// Append a new parameter to the list of parameters.
    auto append(const Param& param) -> Param&;

    /// Append a new parameter to the list of parameters with given @p id and @p value.
    /// @warning A runtime error is thrown if @p id is not unique among other parameters already in the list.
    auto append(const String& id, const real& value) -> Param&;

    /// Resize this container of parameters.
    auto resize(Index size) -> void;

    /// Return the number of parameters.
    auto size() const -> Index;

    /// Return the parameter with given index.
    auto operator[](Index i) -> Param&;

    /// Return the parameter with given index.
    auto operator[](Index i) const -> const Param&;

    /// Return the parameter with given identifier.
    /// @warning A runtime error is thrown if there is no parameter with identifier @p id.
    auto operator[](const String& id) -> Param&;

    /// Return the parameter with given identifier.
    /// @warning A runtime error is thrown if there is no parameter with identifier @p id.
    auto operator[](const String& id) const -> const Param&;

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

//======================================================================
// CODE BELOW NEEDED FOR MEMOIZATION TECHNIQUE INVOLVING PARAMS
//======================================================================

namespace Reaktoro {
namespace detail {

template<typename T> struct SameValue;
template<typename T> struct AssignValue;
template<typename T> struct CloneValue;

template<>
struct SameValue<Params>
{
    static auto check(const Params& a, const Params& b)
    {
        errorif(a.size() != b.size(), "Expecting same size for both Params objects");
        for(auto i = 0; i < a.size(); ++i)
            if(a[i].value() != b[i].value()) // for memoization sake, a and b are equal if they Param objects with same real values (including seed numbers!)
                return false;
        return true;
    }
};

template<>
struct AssignValue<Params>
{
    static auto apply(Params& a, const Params& b)
    {
        a.resize(b.size());
        for(auto i = 0; i < b.size(); ++i)
            a[i].value() = b[i].value();
    }
};

template<>
struct CloneValue<Params>
{
    static auto apply(const Params& params) -> Params
    {
        return params.clone();
    }
};

} // namespace detail
} // namespace Reaktoro

//======================================================================
// CODE BELOW NEEDED FOR AUTOMATIC DIFFERENTIATION INVOLVING PARAMS
//======================================================================

namespace autodiff {
namespace detail {

/// Implementation of VectorTraits for Reaktoro::Params.
template<>
struct VectorTraits<Reaktoro::Params>
{
    using ValueType = Reaktoro::Param;

    template<typename NewValueType>
    using ReplaceValueType = std::vector<NewValueType>;
};

} // namespace autodiff
} // namespace detail


