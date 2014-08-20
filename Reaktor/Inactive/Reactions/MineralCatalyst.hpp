/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// C++ includes
#include <string>

namespace Reaktor {

struct MineralCatalyst
{
    /**
     * Constructs a default MineralCatalyst instance
     */
    MineralCatalyst();

    /**
     * Constructs a MineralCatalyst instance
     *
     * * The options for the chemical quantities are:
     *   - `"a"` or `"activity"` for activity of a species;
     *   - `"p"` or `"pressure"` for partial pressure of a gaseous species.
     *
     * @param species The name of the species that participates as a catalyst in the mineral reaction
     * @param quantity The name of the chemical quantity that acts as a catalyser in the mineral reaction
     * @param power The power of the quantity that affects the rate of mineral reaction
     */
    MineralCatalyst(const std::string& species, const std::string& quantity, double power);

    /**
     * Constructs a MineralCatalyst instance
     *
     * This constructor offers a convenient way to create a MineralCatalyst instance from a string.
     * For example, in the code bellow, @c catalyst1 and @c catalyst2 are equivalent to @c catalyst3 and @c catalyst4
     * respectively:
     *
     * @code
     * MineralCatalyst catalyst1("a[H+]=1.0");     // "activity[H+]=1.0" is also allowed
     * MineralCatalyst catalyst2("p[CO2(g)]=1.5"); // "pressure[CO2(g)]=1.0" is also allowed
     * MineralCatalyst catalyst3("H+", "a", 1.0);
     * MineralCatalyst catalyst4("CO2(g)", "p", 1.0);
     * @endcode
     *
     * @param catalyst The catalyst definition as a string
     */
    MineralCatalyst(const std::string& catalyst);

    /// The name of the species that participates as a catalyst in the mineral reaction
    std::string species;

    /// The name of the chemical quantity that acts as a catalyser in the mineral reaction
    std::string quantity = "activity";

    /// The power of the quantity that affects the rate of mineral reaction
    double power = 0.0;
};

} /* namespace Reaktor */
