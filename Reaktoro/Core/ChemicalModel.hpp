// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <functional>

// Reaktoro includes
#include <Reaktoro/Common/ThermoVector.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>

namespace Reaktoro {

// Forward declarations
class PhaseChemicalModelResult;

/// The result of the chemical model function of a chemical system.
/// This class holds the chemical properties of the species
/// in a chemical system calculated by a ChemicalModel function.
/// @see ChemicalModel, ChemicalSystem
struct ChemicalModelResult
{
    /// Construct a default ChemicalModelResult instance
    ChemicalModelResult();

    /// Construct a ChemicalModelResult instance with given number of species.
    ChemicalModelResult(unsigned nspecies);

    /// Assign a vector of PhaseChemicalModelResult instances to this.
    auto operator=(const std::vector<PhaseChemicalModelResult>& results) -> ChemicalModelResult&;

    /// Return the ln activity constants of the species.
    auto lnActivityConstants() const -> ThermoVector;

    /// Return the ln activity coefficients of the species.
    auto lnActivityCoefficients() const -> ChemicalVector;

    /// Return the ln activities of the species.
    auto lnActivities() const -> ChemicalVector;

private:
    /// The number of species in the system.
    unsigned num_species = 0;

    /// The results of the phase chemical model function of each phase.
    std::vector<PhaseChemicalModelResult> results;
};

} // namespace Reaktoro
