// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <string>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;

/// The properties of a chemical system with respect to a chemical equilibrium state.
class EquilibriumProps
{
public:
    /// Construct a EquilibriumProps instance.
    EquilibriumProps(const ChemicalSystem& system);

    /// Construct a copy of a EquilibriumProps instance
    EquilibriumProps(const EquilibriumProps& other);

    /// Destroy this EquilibriumProps instance
    virtual ~EquilibriumProps();

    /// Assign a EquilibriumProps instance to this instance
    auto operator=(EquilibriumProps other) -> EquilibriumProps&;

    /// Set the indices of the equilibrium species partitioned as (primary, secondary).
    /// @param ips The indices of the equilibrium species ordered (primary, secondary)
    /// @param kp The number of primary species
    auto setIndicesEquilibriumSpecies(ArrayXlConstRef ips, Index kp) -> void;

    /// Set the indices of equilibrium elements whose amounts should be positive, but given amount was less or equal to zero.
    auto setIndicesStrictlyUnstableElements(ArrayXlConstRef isue) -> void;

    /// Set the indices of equilibrium species that contain one or more strictly unstable elements.
    /// @see setIndicesElementsStrictlyUnstable
    auto setIndicesStrictlyUnstableSpecies(ArrayXlConstRef isus) -> void;

    /// Set the chemical potentials of the species in the equilibrium state (in units of J/mol)
    auto setSpeciesChemicalPotentials(ArrayXrConstRef u) -> void;

    /// Set the chemical potentials of the elements in the equilibrium state (in units of J/mol)
    auto setElementChemicalPotentials(ArrayXrConstRef y) -> void;

    /// Set the stabilities of the species in the equilibrium state (in units of J/mol)
    auto setSpeciesStabilities(ArrayXrConstRef z) -> void;

    /// Return the number of equilibrium species.
    auto numEquilibriumSpecies() const -> Index;

    /// Return the number of primary equilibrium species.
    auto numPrimarySpecies() const -> Index;

    /// Return the number of secondary equilibrium species.
    auto numSecondarySpecies() const -> Index;

    /// Return the indices of the equilibrium species.
    auto indicesEquilibriumSpecies() const -> ArrayXlConstRef;

    /// Return the indices of the primary equilibrium species.
    auto indicesPrimarySpecies() const -> ArrayXlConstRef;

    /// Return the indices of the secondary equilibrium species.
    auto indicesSecondarySpecies() const -> ArrayXlConstRef;

    /// Return the indices of equilibrium elements whose amounts should be positive, but given amount was less or equal to zero.
    auto indicesStrictlyUnstableElements() const -> ArrayXlConstRef;

    /// Return the indices of equilibrium species that contain one or more strictly unstable elements.
    auto indicesStrictlyUnstableSpecies() const -> ArrayXlConstRef;

    /// Return the chemical potentials of the species (in units of J/mol).
    /// The chemical potentials of the species are the gradient of the Gibbs energy
    /// function with respect to species amounts.
    auto speciesChemicalPotentials() const -> ArrayXrConstRef;

    /// Return the chemical potentials of the elements (in units of J/mol).
    /// The chemical potentials of the elements are the Lagrange multipliers with
    /// respect to the mass conservation constraints on the amounts of the elements.
    auto elementChemicalPotentials() const -> ArrayXrConstRef;

    /// Return the stabilities of the species (in units of J/mol)
    /// The stabilities of the species are the slack variables with
    /// respect to the non-negative bound constraints on the amounts of the
    /// species in a chemical equilibrium calculation. They can be seen as
    /// measures of stability of a species at equilibrium, with values closer
    /// to zero meaning more stable.
    auto speciesStabilities() const -> ArrayXrConstRef;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
