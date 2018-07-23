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
#include <memory>
#include <string>

namespace Reaktoro {

/// A type used to define a chemical element and its attributes
class Element
{
public:
    /// Construct a default Element instance
    Element();

    /// Set the name of the element
    auto setName(std::string name) -> void;

    /// Set the molar mass of the element (in units of kg/mol)
    auto setMolarMass(double value) -> void;

    /// Return the name of the element
    auto name() const -> std::string;

    /// Return the molar mass of the element (in units of kg/mol)
    auto molarMass() const -> double;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Compare two Element instances for less than
auto operator<(const Element& lhs, const Element& rhs) -> bool;

/// Compare two Element instances for equality
auto operator==(const Element& lhs, const Element& rhs) -> bool;

} // namespace Reaktoro
