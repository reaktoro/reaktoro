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
#include <vector>

// Reaktor includes
#include <Reaktor/Reactions/MineralCatalyst.hpp>

namespace Reaktor {

struct MineralMechanism
{
    /**
     * Constructs a default instance of @ref MineralMechanism
     */
    MineralMechanism();

    /**
     * Constructs an instance of @ref MineralMechanism
     *
     * This constructor offers a convenient way to create a @ref MineralMechanism instance from a string.
     * For example, in the code bellow, @c mechanism1 and @c mechanism2 are respectively the neutral and
     * acidic mechanisms of the calcite reaction <tt>Calcite + H+ = HCO3- + Ca++</tt>:
     *
     * @code
     * MineralMechanism mechanism1("logk = -5.81 mol/(m2*s), Ea = 23.5 kJ/mol");
     * MineralMechanism mechanism2("logk = -0.30 mol/(m2*s), Ea = 14.4 kJ/mol, a[H+] = 1.0");
     * @endcode
     *
     * Note that the kinetic rate constant of the reaction @c logk is given in log scale, and its units
     * must be provided. The units of the Arrhenius activation energy @c Ea must also be provided.
     * @param mechanism The string representing the mineral mechanism
     */
    MineralMechanism(const std::string& mechanism);

    /**
     * Sets the kinetic rate constant of the mineral reaction at 298.15 K
     * @param value The value of the kinetic rate constant
     * @param unit The unit of the kinetic rate constant (must be convertible to m2/g)
     * @return A reference to this mineral mechanism instance
     */
    auto setRateConstant(double value, const std::string& unit) -> MineralMechanism&;

    /**
     * Sets the Arrhenius activation energy of the mineral reaction
     * @param value The value of the Arrhenius activation energy
     * @param unit The unit of the Arrhenius activation energy (must be convertible to kJ/mol)
     * @return A reference to this mineral mechanism instance
     */
    auto setActivationEnergy(double value, const std::string& unit) -> MineralMechanism&;

    /**
     * Sets the power parameter @a p of the mineral mechanism
     * @param value The value of the power parameter @a p
     * @return A reference to this mineral mechanism instance
     */
    auto setPowerP(double value) -> MineralMechanism&;

    /**
     * Sets the power parameter @a q of the mineral mechanism
     * @param value The value of the power parameter @a q
     * @return A reference to this mineral mechanism instance
     */
    auto setPowerQ(double value) -> MineralMechanism&;

    /**
     * Sets the mineral catalysts of the mineral mechanism
     * @param catalysts The string representing the catalysts of the mineral mechanism
     * @return A reference to this mineral mechanism instance
     * @see MineralCatalyst
     */
    auto setCatalysts(const std::string& catalysts) -> MineralMechanism&;

    /**
     * Sets the mineral catalysts of the mineral mechanism
     * @param catalyst The single catalyst instance of the mineral mechanism
     * @return A reference to this mineral mechanism instance
     * @see MineralCatalyst
     */
    auto setCatalysts(const MineralCatalyst& catalyst) -> MineralMechanism&;

    /**
     * Sets the mineral catalysts of the mineral mechanism
     * @param catalysts The vector of catalyst instances of the mineral mechanism
     * @return A reference to this mineral mechanism instance
     * @see MineralCatalyst
     */
    auto setCatalysts(const std::vector<MineralCatalyst>& catalysts) -> MineralMechanism&;

    /**
     * The kinetic rate constant of the mineral reaction at 298.15 K (in units of mol/(m2*s))
     */
    double kappa = 0.0;

    /**
     * The Arrhenius activation energy of the mineral reaction (in units of kJ/mol)
     */
    double Ea = 0.0;

    /**
     * The empirical and dimensionless power parameter @e p
     */
    double p = 1.0;

    /**
     * The empirical and dimensionless power parameter @e q
     */
    double q = 1.0;

    /**
     * The catalysts of the mineral reaction
     */
    std::vector<MineralCatalyst> catalysts;
};

} // namespace Reaktor
