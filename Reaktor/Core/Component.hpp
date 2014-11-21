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
#include <vector>
#include <string>

namespace Reaktor {

/// A type used to define a chemical component and its properties
class Component
{
public:
	/// Construct a default Component instance
	Component();

	/// Get the name of the component
	auto name() const -> std::string;
};

/// A type used to define a list of Component instances
typedef std::vector<Component> ComponentList;

} /* namespace Reaktor */
