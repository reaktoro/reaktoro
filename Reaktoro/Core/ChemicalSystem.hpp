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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Core/ReactionList.hpp>
#include <Reaktoro/Core/Reactions.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// The class used to represent a chemical system and its attributes and properties.
/// @see Species, Phase
/// @ingroup Core
class ChemicalSystem
{
public:
    /// Construct a default uninitialized ChemicalSystem instance.
    ChemicalSystem();

    /// Construct a ChemicalSystem instance with given database and phases.
    ChemicalSystem(const Database& database, const Vec<Phase>& phases);

    /// Construct a ChemicalSystem instance with given database, phases, and reactions.
    ChemicalSystem(const Database& database, const Vec<Phase>& phases, const Vec<Reaction>& reactions);

    /// Construct a ChemicalSystem instance with given phases.
    explicit ChemicalSystem(const Phases& phases);

    /// Construct a ChemicalSystem instance with given phases and reactions.
    explicit ChemicalSystem(const Phases& phases, const Reactions& reactions);

    /// Construct a ChemicalSystem instance with given database and one or more generic phases.
    template<typename... GenericPhases, EnableIf<areGenericPhases<GenericPhases...>>...>
    ChemicalSystem(const Database& database, const GenericPhases&... genericPhases)
    : ChemicalSystem(Phases(database, genericPhases...)) {}

    /// Return the database used to construct the chemical system.
    auto database() const -> const Database&;

    /// Return the element in the system with given index.
    auto element(Index index) const -> const Element&;

    /// Return the list of elements in the system.
    auto elements() const -> const ElementList&;

    /// Return the species in the system with given index.
    auto species(Index index) const -> const Species&;

    /// Return the list of species in the system.
    auto species() const -> const SpeciesList&;

    /// Return the phase in the system with given index.
    auto phase(Index index) const -> const Phase&;

    /// Return the list of phases in the system.
    auto phases() const -> const PhaseList&;

    /// Return the reaction in the system with given index.
    auto reaction(Index index) const -> const Reaction&;

    /// Return the list of reactions in the system.
    auto reactions() const -> const ReactionList&;

    /// Return the formula matrix of the system.
    /// The formula matrix is defined as the matrix whose entry *(j, i)* is
    /// given by the coefficient of the *j*th element in the *i*th species. The
    /// last row of this matrix contains the electric charge of each species.
    auto formulaMatrix() const -> MatrixXdConstRef;

    /// Return the top rows of the formula matrix corresponding to elements.
    auto formulaMatrixElements() const -> MatrixXdConstRef;

    /// Return the bottom row of the formula matrix corresponding to electric charge.
    auto formulaMatrixCharge() const -> MatrixXdConstRef;

    /// Return the stoichiometric matrix of the reactions corresponding to the species in the system.
    /// The stoichiometric matrix is defined as the matrix whose entry *(i, j)*
    /// is given by the coefficient of the *i*th species in the *j*th reaction.
    auto stoichiometricMatrix() const -> MatrixXdConstRef;

    /// Return the phase interfaces across which kinetically controlled
    /// reactions take place. These interfaces are automatically detected from
    /// the given reactions composing the chemical system. They are implemented
    /// as phase index pairs `{(i, j)}`. If a kinetically controlled reaction
    /// is defined with a single species, both indices `i` and `j` are
    /// identical, and this is the index of the phase in which the species
    /// exists. This is needed when, for example, mineral-aqueous reactions are
    /// defined and all aqueous species assumed to be in equilibrium.
    auto reactingPhaseInterfaces() const -> const Pairs<Index, Index>&;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Output a ChemicalSystem instance
auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&;

} // namespace Reaktoro
