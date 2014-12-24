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

// C++ includes
#include <memory>
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>

namespace Reaktor {

/// A type used to define the attributes of an Element instance
/// @see Element
/// @ingroup Core
struct ElementData
{
    /// The name of the chemical element
    std::string name;

    /// The molar mass of the chemical element (in units of kg/mol)
    double molar_mass = 0.0;
};

/// A type used to define a chemical element and its attributes
class Element
{
public:
	/// Construct a default Element instance
	Element();

	/// Construct an Element instance with all its attributes
    Element(const ElementData& data);

	/// Get the name of the element
    auto name() const -> std::string;

    /// Get the molar mass of the element (in units of kg/mol)
    auto molarMass() const -> double;

private:
	struct Impl;

	std::shared_ptr<Impl> pimpl;
};

/// Compare two Element instances for less than
auto operator<(const Element& lhs, const Element& rhs) -> bool;

/// Compare two Element instances for equality
auto operator==(const Element& lhs, const Element& rhs) -> bool;

} // namespace Reaktor
