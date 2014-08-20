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

// Reaktor includes
#include <Reaktor/Species/MineralSpecies.hpp>
#include <Reaktor/Mixtures/GeneralMixture.hpp>

namespace Reaktor {

/**
 * Provides a computational representation of a mineral mixture
 *
 * The MineralMixture class is defined as a collection of MineralSpecies objects,
 * representing, therefore, a mixture of mineral species. Its main purpose is to
 * provide the necessary operations in the calculation of activities of mineral
 * species.
 *
 * @see MineralSpecies
 *
 * @ingroup Mixtures
 */
class MineralMixture : public GeneralMixture<MineralSpecies>
{
public:
    /**
     * Constructs a default MineralMixture instance
     */
    MineralMixture();

    /**
     * Constructs a MineralMixture instance with given species
     * @param species The species that compose the mineral mixture
     */
    MineralMixture(const std::vector<MineralSpecies>& species);

    /**
     * Constructs a MineralMixture instance with a single species
     * @param species The species that compose the mineral mixture
     */
    MineralMixture(const MineralSpecies& species);

    /**
     * Destroys the instance
     */
    virtual ~MineralMixture();
};

} /* namespace Reaktor */
