// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalFormula;
class ChemicalSystem;
class ElementList;
class Phase;
class PhaseList;
class Reaction;
class ReactionList;
class SpeciesList;
class Surface;
class SurfaceList;

namespace detail {

/// Return the molar masses of the species.
auto molarMasses(SpeciesList const& species) -> ArrayXd;

/// Compute the amount of a species given a value in mass or amount unit.
/// @param system The chemical system in which the species is.
/// @param ispecies The index of the species in the chemical system.
/// @param value The quantity value of the species.
/// @param unit The quantity unit of the species.
auto computeSpeciesAmount(ChemicalSystem const& system, Index ispecies, real value, Chars unit) -> real;

/// Resolve the index of an element in a list of elements with given element symbol or its index.
auto resolveElementIndex(ElementList const& elementlist, StringOrIndex element) -> Index;

/// Resolve the index of an element in a chemical system with given element symbol or its index.
auto resolveElementIndex(ChemicalSystem const& system, StringOrIndex element) -> Index;

/// Resolve the index of an element in a phase with given element symbol or its index.
auto resolveElementIndex(Phase const& phase, StringOrIndex element) -> Index;

/// Resolve the index of a species in a list of species with given species name or its index.
auto resolveSpeciesIndex(SpeciesList const& specieslist, StringOrIndex species) -> Index;

/// Resolve the index of a species in a chemical system with given species name or its index.
auto resolveSpeciesIndex(ChemicalSystem const& system, StringOrIndex species) -> Index;

/// Resolve the index of a species in a phase with given species name or its index.
auto resolveSpeciesIndex(Phase const& phase, StringOrIndex species) -> Index;

/// Resolve the index of a phase in a list of phases with given phase name or its index.
auto resolvePhaseIndex(PhaseList const& phaselist, StringOrIndex phase) -> Index;

/// Resolve the index of a phase in a chemical system with given phase name or its index.
auto resolvePhaseIndex(ChemicalSystem const& system, StringOrIndex phase) -> Index;

/// Resolve the index of a surface in a list of surfaces with given surface name or its index.
auto resolveSurfaceIndex(SurfaceList const& surfacelist, StringOrIndex surface) -> Index;

/// Resolve the index of a surface in a chemical system with given surface name or its index.
auto resolveSurfaceIndex(ChemicalSystem const& system, StringOrIndex surface) -> Index;

/// Assemble the formula vector of a `substance` with respect to given list of `elements`.
auto assembleFormulaVector(ChemicalFormula const& susbtance, ElementList const& elements) -> VectorXd;

/// Assemble the formula matrix of the list of `species` with respect to given list of `elements`.
auto assembleFormulaMatrix(SpeciesList const& species, ElementList const& elements) -> MatrixXd;

/// Assemble the stoichiometric matrix of the list of `reactions` with respect to given list of `species`.
auto assembleStoichiometricMatrix(ReactionList const& reactions, SpeciesList const& species) -> MatrixXd;

/// Extract names of the species from the species' list
auto extractNames(SpeciesList const& list) -> Strings;

/// Extract names of the elements from the elements' list
auto extractNames(ElementList const& list) -> Strings;

/// Extract names of the phases from the phases' list
auto extractNames(PhaseList const& list) -> Strings;

/// Return pairs of indices of phases that participate in a reaction.
/// In case the reaction contains only species from a same phase, an empty
/// vector `{}` is returned. In case the reaction contains only one species, a
/// single-entry vector `{ {i, i} }` is returned where `i` is the index of the
/// phase in which the species exists.
auto determineReactingPhaseInterfacesInReaction(Reaction const& reaction, PhaseList const& phases) -> Pairs<Index, Index>;

/// Return the phase interfaces, as phase index pairs, across which reactions take place.
auto determineReactingPhaseInterfacesInReactions(Vec<Reaction> const& reactions, PhaseList const& phases) -> Pairs<Index, Index>;

/// Return the phase interfaces, as phase index pairs, across which reactions take place.
auto createSurfacesForReactingPhaseInterfacesInReactions(Vec<Reaction> const& reactions, PhaseList const& phases) -> Vec<Surface>;

} // namespace detail
} // namespace Reaktoro
