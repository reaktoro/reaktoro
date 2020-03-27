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
#include <string>

namespace Reaktoro {

/// The states of aggregation of chemical species according to IUPAC.
/// Cox, J. D. (1982). Notation for states and processes, significance of the
/// word standard in chemical thermodynamics, and remarks on commonly tabulated
/// forms of thermodynamic functions, Pure and Applied Chemistry, 54(6),
/// 1239-1250. doi: https://doi.org/10.1351/pac198254061239
enum class AggregateState
{
    Gas,              ///< for a gas or a vapour (symbol g)
    Liquid,           ///< for a liquid (symbol l)
    Solid,            ///< for a solid (symbol s)
    CondensedPhase,   ///< for either the solid or the liquid state (symbol cd)
    Fluid,            ///< for either the gaseous or the liquid state) (symbol fl)
    LiquidCrystal,    ///< for a liquid crystal (crystalline liquid) (symbol lc)
    CrystallineSolid, ///< for a crystalline solid (symbol cr)
    AmorphousSolid,   ///< for an amorphous solid (symbol am)
    Vitreous,         ///< for a vitreous substance (a glass) (symbol vit)
    Adsorbed,         ///< for a species adsorbed on a substrate (an adsorbate) (symbol ads)
    Monomeric,        ///< for a monomeric form (symbol mon)
    Polymeric,        ///< for a polymeric form (symbol pol)
    SolidSolution,    ///< for a solid solution (symbol ss)
    Aqueous,          ///< for a species in a solution in which water is the solvent (symbol aq)
    Undefined         ///< when aggregate state is not explicitly provided (default)
};

/// Return the AggregateState value from given aggregate state symbol.
/// The following table maps a symbol to an aggregate state value.
///
/// | Symbol | Aggregate State                  |
/// |:------ |:-------------------------------- |
/// | `g`    | AggregateState::Gas              |
/// | `l`    | AggregateState::Liquid           |
/// | `s`    | AggregateState::Solid            |
/// | `cd`   | AggregateState::CondensedPhase   |
/// | `fl`   | AggregateState::Fluid            |
/// | `lc`   | AggregateState::LiquidCrystal    |
/// | `cr`   | AggregateState::CrystallineSolid |
/// | `am`   | AggregateState::AmorphousSolid   |
/// | `vit`  | AggregateState::Vitreous         |
/// | `ads`  | AggregateState::Adsorbed         |
/// | `mon`  | AggregateState::Monomeric        |
/// | `pol`  | AggregateState::Polymeric        |
/// | `ss`   | AggregateState::SolidSolution    |
/// | `aq`   | AggregateState::Aqueous          |
///
/// AggregateState::Undefined is returned if symbol is none of above.
auto parseAggregateState(const std::string& symbol) -> AggregateState;

} // namespace Reaktoro
