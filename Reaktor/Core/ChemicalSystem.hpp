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

// Reaktor includes
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Common/Matrix.hpp>
#include <Reaktor/Common/ThermoVector.hpp>
#include <Reaktor/Core/ChemicalModels.hpp>
#include <Reaktor/Core/Element.hpp>
#include <Reaktor/Core/Multiphase.hpp>
#include <Reaktor/Core/Phase.hpp>
#include <Reaktor/Core/Species.hpp>

namespace Reaktor {

/// A type to describe the connectivity of elements, species, and phases in a chemical system
/// @see ChemicalSystem, Element, Species, Phase
/// @ingroup Core
struct Connectivity
{
    /// The mapping from the index of an element to the indices of the species that contains it.
    std::vector<Indices> element_to_species;

    /// The mapping from the index of a species to the indices of the elements that it contains.
    std::vector<Indices> species_to_elements;

    /// The mapping from the index of a species to the index of the phase that contains it.
    Indices species_to_phase;

    /// The mapping from the index of a phase to the indices of the species that it contains.
    std::vector<Indices> phase_to_species;

    /// The mapping from the index of an element to the indices of the phases that contains it.
    std::vector<Indices> element_to_phases;

    /// The mapping from the index of a phase to the indices of the elements that it contains.
    std::vector<Indices> phase_to_elements;
};

/// The type used to define a chemical system and its attributes
/// @see Multiphase, Phase, Species
/// @ingroup Core
class ChemicalSystem : public Multiphase
{
public:
    /// Construct a default ChemicalSystem instance
    ChemicalSystem();

    /// Construct a ChemicalSystem instance with given phases.
    explicit ChemicalSystem(const Multiphase& multiphase);

    /// Return the connectivity of the elements, species, and phases in the chemical system
    auto connectivity() const -> const Connectivity&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Output a ChemicalSystem instance
auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&;

} // namespace Reaktor
