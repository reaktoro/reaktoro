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

namespace Reaktoro {

// Forward declarations
class ChemicalState;
class EquilibriumRestrictions;
struct EquilibriumOptions;
struct EquilibriumResult;

/// Perform a chemical equilibrium calculation on a given chemical state.
/// The calculation is performed with fixed temperature and pressure obtained
/// from the chemical state, and the chemical system is closed, so chemical
/// elements and electric charge are conserved.
///@{
auto equilibrate(ChemicalState& state) -> EquilibriumResult;
auto equilibrate(ChemicalState& state, const EquilibriumOptions& options) -> EquilibriumResult;
auto equilibrate(ChemicalState& state, const EquilibriumRestrictions& restrictions) -> EquilibriumResult;
auto equilibrate(ChemicalState& state, const EquilibriumRestrictions& restrictions, const EquilibriumOptions& options) -> EquilibriumResult;
///@}

} // namespace Reaktoro
