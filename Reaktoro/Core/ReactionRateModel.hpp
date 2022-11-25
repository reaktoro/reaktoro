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
#include <Reaktoro/Core/Model.hpp>
#include <Reaktoro/Core/ReactionRate.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalProps;
class Database;
class PhaseList;
class ReactionEquation;
class SpeciesList;
class SurfaceList;

/// The type of functions for calculation of reaction rates (in mol/s).
/// @param props The chemical properties of the chemical system
/// @return The rate of the reaction (in mol/s)
/// @see Reaction
/// @ingroup Core
using ReactionRateModel = Model<ReactionRate(ChemicalProps const& props)>;

/// The data provided to a ReactionRateModelGenerator to construct the ReactionRateModel of a reaction.
/// @see ReactionRateModelGenerator, ReactionRateModel
/// @ingroup Core
struct ReactionRateModelGeneratorArgs
{
    /// The name of the reaction for which the rate model is generated.
    String const& name;

    /// The equation of the reaction for which the rate model is generated.
    ReactionEquation const& equation;

    /// The thermodynamic database used to construct the chemical system where the reaction belongs to.
    Database const& database;

    /// The species in the chemical system where the reaction belongs to.
    SpeciesList const& species;

    /// The phases in the chemical system where the reaction belongs to.
    PhaseList const& phases;

    /// The surfaces in the chemical system where the reaction belongs to.
    SurfaceList const& surfaces;
};

/// The function signature for functions that generates a ReactionRateModel for a reaction.
using ReactionRateModelGenerator = Fn<ReactionRateModel(ReactionRateModelGeneratorArgs)>;

} // namespace Reaktoro
