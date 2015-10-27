// Reaktoro is a C++ library for computational reaction modelling.
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

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

// Boost includes
#include <boost/optional.hpp>

namespace Reaktoro {

/// A type for defining an optional instance
template<typename T>
class Optional
{
public:
    /// Construct a default, uninitialized, Optional instance
    Optional() {}

    /// Construct a default, initialized, Optional instance
    Optional(const T& value) : data(value) {}

    /// Retrieve the value of the Optional instance
    auto operator()() const -> const T& { return get(); }

    /// Initialize the value of the Optional instance
    auto set(const T& value) -> void { data.reset(value); }

    /// Retrieve the value of the Optional instance
    auto get() const -> const T&
    {
        if(!empty()) return data.get();
        else RuntimeError("Cannot get the value of the Optional instance.",
            "Its value has not been initialized.");
    }

    /// Retrieve the value of the Optional instance
    auto get() -> T&
    {
        if(!empty()) return data.get();
        else RuntimeError("Cannot get the value of the Optional instance.",
            "Its value has not been initialized.");
    }

    /// Check if the Optional instance is initialized
    auto empty() const -> bool { return !data; }

private:
    boost::optional<T> data;
};

} // namespace Reaktoro
