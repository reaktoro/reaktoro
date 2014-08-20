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

#include "EquilibriumUtils.hpp"

// C++ includes
#include <set>
#include <vector>

// Reaktor includes
#include <Reaktor/Core/ChemicalState.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partitioning.hpp>

namespace Reaktor {

auto idxExistentPhases(const ChemicalState& state, const Partitioning& partitioning, double epsilon) -> Indices
{
    Indices idx_existent_phases;

    std::set<Index> idx_equilibrium_phases;
    for(const Index ispecies : partitioning.idxEquilibriumSpecies())
        idx_equilibrium_phases.insert(state.system().idxPhaseWithSpecies(ispecies));

    const double ntotal = state.equilibriumComposition(partitioning).sum();

    for(const Index& iphase : idx_equilibrium_phases)
    {
        const double nphase = state.phaseComposition(iphase).sum();

        if(nphase > ntotal * epsilon)
            idx_existent_phases.push_back(iphase);
    }

    return idx_existent_phases;
}

} /* namespace Reaktor */
