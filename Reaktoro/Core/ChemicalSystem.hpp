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
#include <Reaktoro/Common/Meta.hpp>
#include <Reaktoro/Common/TraitsUtils.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/ElementList.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/PhaseList.hpp>
#include <Reaktoro/Core/Phases.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/ReactionList.hpp>
#include <Reaktoro/Core/Reactions.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>
#include <Reaktoro/Core/Surface.hpp>
#include <Reaktoro/Core/SurfaceList.hpp>
#include <Reaktoro/Core/Surfaces.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalSystem;

template<typename T, typename... Ts>
constexpr auto _arePhaseReactionOrSurfaceConvertible()
{
    constexpr auto isPhaseConvertible = isBaseOf<GeneralPhase, T> || isBaseOf<GeneralPhasesGenerator, T>;
    constexpr auto isReactionConvertible = isConvertible<T, Reaction> || isConvertible<T, GeneralReaction> || isConvertible<T, ReactionGenerator>;
    constexpr auto isSurfaceConvertible = isConvertible<T, Surface> || isConvertible<T, GeneralSurface> || isConvertible<T, SurfaceGenerator>;
    constexpr auto aux = isPhaseConvertible || isReactionConvertible || isSurfaceConvertible;

    if constexpr (sizeof...(Ts))
        return aux && _arePhaseReactionOrSurfaceConvertible<Ts...>();
    else return aux;
}

/// Used to determine if `T` and all types in `Ts` are either GeneralPhase or GeneralPhaseGenerator.
template<typename T, typename... Ts>
constexpr auto arePhaseReactionOrSurfaceConvertible = _arePhaseReactionOrSurfaceConvertible<T, Ts...>();

/// Create a ChemicalSystem object with given database and a list of phase and reaction convertible objects.
template<typename... Args>
auto createChemicalSystem(Database const& db, Args const&... args) -> ChemicalSystem;

/// The class used to represent a chemical system and its attributes and properties.
/// @see Species, Phase
/// @ingroup Core
class ChemicalSystem
{
public:
    /// Construct a default uninitialized ChemicalSystem object.
    ChemicalSystem();

    /// Construct a ChemicalSystem object with given database and phases.
    explicit ChemicalSystem(Database const& database, PhaseList const& phases);

    /// Construct a ChemicalSystem object with given database, phases, and surfaces.
    explicit ChemicalSystem(Database const& database, PhaseList const& phases, SurfaceList const& surfaces);

    /// Construct a ChemicalSystem object with given database, phases, and reactions.
    explicit ChemicalSystem(Database const& database, PhaseList const& phases, ReactionList const& reactions);

    /// Construct a ChemicalSystem object with given database, phases, reactions, and surfaces.
    explicit ChemicalSystem(Database const& database, PhaseList const& phases, ReactionList const& reactions, SurfaceList const& surfaces);

    /// Construct a ChemicalSystem object with given phases.
    explicit ChemicalSystem(Phases const& phases);

    /// Construct a ChemicalSystem object with given phases and surfaces.
    explicit ChemicalSystem(Phases const& phases, Surfaces const& surfaces);

    /// Construct a ChemicalSystem object with given phases and reactions.
    explicit ChemicalSystem(Phases const& phases, Reactions const& reactions);

    /// Construct a ChemicalSystem object with given phases, reactions, and surfaces.
    explicit ChemicalSystem(Phases const& phases, Reactions const& reactions, Surfaces const& surfaces);

    /// Construct a ChemicalSystem object with given database and a list of phase and reaction convertible objects.
    template<typename... Args, EnableIf<arePhaseReactionOrSurfaceConvertible<Args...>>...>
    explicit ChemicalSystem(Database const& db, Args const&... args)
    : ChemicalSystem(createChemicalSystem(db, args...)) {}

    /// Return the unique identification number of this ChemicalSystem object.
    /// ChemicalSystem objects are guaranteed to be the same if they have the same id.
    auto id() const -> Index;

    /// Return the database used to construct the chemical system.
    auto database() const -> Database const&;

    /// Return the element in the system with given index.
    auto element(Index index) const -> Element const&;

    /// Return the list of elements in the system.
    auto elements() const -> ElementList const&;

    /// Return the species in the system with given index.
    auto species(Index index) const -> Species const&;

    /// Return the list of species in the system.
    auto species() const -> SpeciesList const&;

    /// Return the phase in the system with given index.
    auto phase(Index index) const -> Phase const&;

    /// Return the list of phases in the system.
    auto phases() const -> PhaseList const&;

    /// Return the reaction in the system with given index.
    auto reaction(Index index) const -> Reaction const&;

    /// Return the list of reactions in the system.
    auto reactions() const -> ReactionList const&;

    /// Return the surface in the system with given index.
    auto surface(Index index) const -> Surface const&;

    /// Return the list of surfaces in the system.
    auto surfaces() const -> SurfaceList const&;

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

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

/// Output a ChemicalSystem object
auto operator<<(std::ostream& out, ChemicalSystem const& system) -> std::ostream&;

template<typename... Args>
auto createChemicalSystem(Database const& db, Args const&... args) -> ChemicalSystem
{
    Phases phases(db);
    Reactions reactions;
    Surfaces surfaces;

    ForEach([&](auto arg) constexpr {
        using T = Decay<decltype(arg)>;
        constexpr auto isPhaseConvertible = isBaseOf<GeneralPhase, T> || isBaseOf<GeneralPhasesGenerator, T>;
        constexpr auto isReactionConvertible = isConvertible<T, Reaction> || isConvertible<T, GeneralReaction> || isConvertible<T, ReactionGenerator>;
        constexpr auto isSurfaceConvertible = isConvertible<T, Surface> || isConvertible<T, GeneralSurface> || isConvertible<T, SurfaceGenerator>;

        static_assert(isPhaseConvertible || isReactionConvertible || isSurfaceConvertible, "One of the arguments in your list of arguments for the construction of a ChemicalSystem object has a non-convertible type to either phase, reaction, or surface.");

        if constexpr(isPhaseConvertible) phases.add(arg);
        else if constexpr(isReactionConvertible) reactions.add(arg);
        else surfaces.add(arg);
    }, args...);

    return ChemicalSystem(phases, reactions, surfaces);
}

} // namespace Reaktoro
