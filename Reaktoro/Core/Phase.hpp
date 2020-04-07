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
#include <string>

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Core/StateOfMatter.hpp>

namespace Reaktoro {

/// A type used to define a phase and its attributes.
/// @see ChemicalSystem, Element, Species
/// @ingroup Core
class Phase
{
public:
    /// Construct a Phase object.
    Phase();

    /// Construct a copy of a Phase object.
    Phase(const Phase& other);

    /// Return a copy of this Phase object with a new name.
    auto withName(std::string name) -> Phase;

    /// Return a copy of this Phase object with new list of species.
    auto withSpecies(SpeciesList species) -> Phase;

    /// Return a copy of this Phase object with a new physical state.
    auto withPhysicalState(StateOfMatter state) -> Phase;

    /// Return a copy of this Phase object with a new activity model.
    auto withActivityModel(ActivityModelFn model) -> Phase;

    /// Return the name of the phase.
    auto name() const -> std::string;

    /// Return the physical state of the phase.
    auto physicalState() const -> StateOfMatter;

    /// Return the species of the phase.
    auto species() const -> const SpeciesList&;

    /// Return the species in the phase with given index.
    auto species(Index idx) const -> const Species&;

    /// Return the activity model of the phase.
    /// @see ActivityModel
    auto activityModel() const -> const ActivityModelFn&;

    /// Calculate the activity and excess thermodynamic properties of the phase.
    /// @param T The temperature of the system (in units of K)
    /// @param P The pressure of the system (in units of Pa)
    /// @param n The amounts of the species in the phase (in units of mol)
    auto activityProps(real T, real P, VectorXrConstRef n) const -> ActivityProps;

    /// Return a deep copy of this Phase object.
    auto clone() const -> Phase;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Compare two Phase instances for less than
auto operator<(const Phase& lhs, const Phase& rhs) -> bool;

/// Compare two Phase instances for equality
auto operator==(const Phase& lhs, const Phase& rhs) -> bool;

} // namespace Reaktoro
