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
#include <string>

namespace Reaktoro {

struct MineralCatalyst
{
    /// Construct a default MineralCatalyst instance
    MineralCatalyst();

    /// Construct a MineralCatalyst instance
    /// The options for the chemical quantities are:
    ///   - `"a"` or `"activity"` for activity of a species;
    ///   - `"p"` or `"pressure"` for partial pressure of a gaseous species.
    /// @param species The name of the species that participates as a catalyst in the mineral reaction
    /// @param quantity The name of the chemical quantity that acts as a catalyser in the mineral reaction
    /// @param power The power of the quantity that affects the rate of mineral reaction
    MineralCatalyst(std::string species, std::string quantity, double power);

    /// Construct a MineralCatalyst instance
    /// This constructor offers a convenient way to create a MineralCatalyst instance from a string.
    /// For example, in the code bellow, @c catalyst1 and @c catalyst2 are equivalent to @c catalyst3 and @c catalyst4
    /// respectively:
    /// ~~~~~~~~~~~~~~~
    /// MineralCatalyst catalyst1("a[H+]=1.0");     // "activity[H+]=1.0" is also allowed
    /// MineralCatalyst catalyst2("p[CO2(g)]=1.5"); // "pressure[CO2(g)]=1.0" is also allowed
    /// MineralCatalyst catalyst3("H+", "a", 1.0);
    /// MineralCatalyst catalyst4("CO2(g)", "p", 1.0);
    /// ~~~~~~~~~~~~~~~
    /// @param catalyst The catalyst definition as a string
    MineralCatalyst(std::string catalyst);

    /// The name of the species that participates as a catalyst in the mineral reaction
    std::string species;

    /// The name of the chemical quantity that acts as a catalyser in the mineral reaction
    std::string quantity = "activity";

    /// The power of the quantity that affects the rate of mineral reaction
    double power = 0.0;
};

/// Compare two MineralCatalyst instances for less than
auto operator<(const MineralCatalyst& lhs, const MineralCatalyst& rhs) -> bool;

/// Compare two MineralCatalyst instances for equality
auto operator==(const MineralCatalyst& lhs, const MineralCatalyst& rhs) -> bool;

} // namespace Reaktoro
