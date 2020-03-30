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
#include <deque>
#include <optional>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalFormula.hpp>

namespace Reaktoro {

/// A type used to represent a dissociation reaction of a neutral substance into ions.
struct DissociationReaction
{
    /// The chemical formula of the complex substance that dissociates into ions.
    ChemicalFormula complex;

    /// The dissociated ions and their stoichiometric coefficients.
    std::vector<std::pair<double, ChemicalFormula>> ions;
};

/// A type used store a collection of dissociation reactions.
/// @see Element
class DissociationReactions
{
public:
    /// Construct a copy of a DissociationReactions object [deleted].
    DissociationReactions(const DissociationReactions&) = delete;

    /// Assign a DissociationReactions object to this [deleted].
    auto operator=(const DissociationReactions&) -> DissociationReactions& = delete;

    /// Return the single DissociationReactions object.
    static auto instance() -> DissociationReactions&;

    /// Return the dissociation reactions in the database.
    static auto reactions() -> const std::deque<DissociationReaction>&;

    /// Append a dissociation reaction in to the database.
    static auto append(DissociationReaction reaction) -> void;

    /// Return the number of dissociation reactions in the database.
    static auto size() -> std::size_t;

    /// Return the dissociation reaction of the substance with given chemical formula.
    static auto get(const ChemicalFormula& complex) -> std::optional<DissociationReaction>;

    /// Return begin const iterator of this DissociationReactions instance
    auto begin() const;

    /// Return begin iterator of this DissociationReactions instance
    auto begin();

    /// Return end const iterator of this DissociationReactions instance
    auto end() const;

    /// Return end iterator of this DissociationReactions instance
    auto end();

private:
    /// The dissociation reactions stored in the database.
    std::deque<DissociationReaction> m_reactions;

private:
    /// Construct a default DissociationReactions object [private].
    DissociationReactions();

    /// Destroy this DissociationReactions object [private].
    ~DissociationReactions();
};

} // namespace Reaktoro
