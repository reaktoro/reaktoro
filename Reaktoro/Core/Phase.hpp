// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Core/ActivityProps.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Core/StateOfMatter.hpp>

namespace Reaktoro {

/// A type used to define a phase and its attributes.
/// @see ChemicalSystem, Element, Species
/// @ingroup Core
class Phase
{
public:
    /// Construct a default Phase object.
    Phase();

    /// Return a deep copy of this Phase object.
    auto clone() const -> Phase;

    /// Return a copy of this Phase object with a new name.
    auto withName(String name) -> Phase;

    /// Return a copy of this Phase object with new list of species.
    auto withSpecies(SpeciesList species) -> Phase;

    /// Return a copy of this Phase object with a new state of matter.
    auto withStateOfMatter(StateOfMatter state) -> Phase;

    /// Return a copy of this Phase object with a new activity model function.
    auto withActivityPropsFn(ActivityPropsFn fn) -> Phase;

    /// Return a copy of this Phase object with a new ideal activity model function.
    auto withIdealActivityPropsFn(ActivityPropsFn fn) -> Phase;

    /// Return the name of the phase.
    auto name() const -> String;

    /// Return the state of matter of the phase.
    auto stateOfMatter() const -> StateOfMatter;

    /// Return the common aggregate state of the species in the phase.
    auto aggregateState() const -> AggregateState;

    /// Return the elements of the phase.
    auto elements() const -> const ElementList&;

    /// Return the element in the phase with given index.
    auto element(Index idx) const -> const Element&;

    /// Return the species of the phase.
    auto species() const -> const SpeciesList&;

    /// Return the species in the phase with given index.
    auto species(Index idx) const -> const Species&;

    /// Return the function that computes activity properties of the phase.
    auto activityPropsFn() const -> const ActivityPropsFn&;

    /// Return the function that computes ideal activity properties of the phase.
    auto idealActivityPropsFn() const -> const ActivityPropsFn&;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

/// Compare two Phase instances for less than
auto operator<(const Phase& lhs, const Phase& rhs) -> bool;

/// Compare two Phase instances for equality
auto operator==(const Phase& lhs, const Phase& rhs) -> bool;

} // namespace Reaktoro
