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
#include <Reaktoro/Core/Reaction.hpp>

namespace Reaktoro {

/// The class used to configure mineral dissolution/precipitation reactions.
class ReactionGenerator
{
public:
    /// Destroy this ReactionGenerator object.
    virtual ~ReactionGenerator() = default;

    /// Convert a ReactionGenerator object into a vector of Reaction objects.
    /// In case the derived class always produce a single reaction, return a vector with a single reaction.
    virtual auto convert(ChemicalSystem const& system) const -> Vec<Reaction> = 0;
};

/// The function type for the generation of reactions with a given chemical system.
using ReactionGeneratorFn = Fn<Vec<Reaction>(ChemicalSystem const&)>;

} // namespace Reaktoro
