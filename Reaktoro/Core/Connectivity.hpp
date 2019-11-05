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
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;

/// A type to describe the connectivity of elements, species, and phases in a chemical system
/// @see ChemicalSystem, Element, Species, Phase
/// @ingroup Core
class Connectivity
{
  public:
    /// Construct a default Connectivity instance
    Connectivity();

    /// Construct a Connectivity instance with given chemical system
    Connectivity(const ChemicalSystem& system);

    /// Construct a copy of a Connectivity instance
    Connectivity(const Connectivity& other);

    /// Destroy this instance
    virtual ~Connectivity();

    /// Assign a Connectivity instance to this instance
    auto operator=(Connectivity other) -> Connectivity&;

    /// Return the indices of the elements in a species.
    /// @param ispecies The index of the species.
    auto indicesElementsInSpecies(Index ispecies) const -> const Indices&;

    /// Return the indices of the elements in a phase.
    /// @param iphase The index of the phase.
    auto indicesElementsInPhase(Index iphase) const -> const Indices&;

    /// Return the indices of the species in a phase.
    /// @param iphase The index of the phase.
    auto indicesSpeciesInPhase(Index iphase) const -> const Indices&;

    /// Return the indices of the species with an element.
    /// @param ielement The index of the element.
    auto indicesSpeciesWithElement(Index ielement) const -> const Indices&;

    /// Return the indices of the phases with an element.
    /// @param ielement The index of the element.
    auto indicesPhasesWithElement(Index ielement) const -> const Indices&;

    /// Return the index of the phase with a species.
    /// @param ispecies The index of the species.
    auto indexPhaseWithSpecies(Index ispecies) const -> Index;

  private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
