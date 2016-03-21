// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// A type that describes a solver for many chemical calculations.
class ChemicalSolver
{
public:
    /// Construct a default ChemicalSolver instance.
    ChemicalSolver();

    /// Destroy a ChemicalSolver instance.
    virtual ~ChemicalSolver();

    /// Set the partitioning of the chemical system.
    auto setPartition(const Partition& partition) -> void;

    /// Calculate the porosity field.
    /// @param values The input array containing the output calculated data.
    auto porosity(const double* values) const -> void;

    /// Calculate the equilibrium-sensitivity of the porosity field with respect to temperature.
    /// @param values The input array containing the output calculated data.
    auto porositySensitivityT(const double* values) const -> void;

    /// Calculate the equilibrium-sensitivity of the porosity field with respect to pressure.
    /// @param values The input array containing the output calculated data.
    auto porositySensitivityP(const double* values) const -> void;

    /// Calculate the equilibrium-sensitivity of the porosity field with respect to the molar amount of an equilibrium element.
    /// @param values The input array containing the output calculated data.
    /// @param iequilibriumelement The index of the equilibrium element.
    auto porositySensitivityB(const double* values, unsigned iequilibriumelement) const -> void;

    /// Calculate the saturation field of a fluid phase.
    /// @param values The input array containing the output calculated data.
    /// @param ifluidphase The index of the fluid phase.
    auto saturation(const double* values, unsigned ifluidphase) const -> void;

    /// Calculate the equilibrium-sensitivity of the saturation field of a fluid phase with respect to temperature.
    /// @param values The input array containing the output calculated data.
    /// @param ifluidphase The index of the fluid phase.
    auto saturationSensitivityT(const double* values, unsigned ifluidphase) const -> void;

    /// Calculate the equilibrium-sensitivity of the saturation field of a fluid phase with respect to pressure.
    /// @param values The input array containing the output calculated data.
    /// @param ifluidphase The index of the fluid phase.
    auto saturationSensitivityP(const double* values, unsigned ifluidphase) const -> void;

    /// Calculate the equilibrium-sensitivity of the saturation field of a fluid phase with respect to the molar amount of an element.
    /// @param values The input array containing the output calculated data.
    /// @param ifluidphase The index of the fluid phase.
    /// @param iequilibriumelement The index of the equilibrium element.
    auto saturationSensitivityB(const double* values, unsigned ifluidphase, unsigned iequilibriumelement) const -> void;

    /// Calculate the density field of a fluid phase.
    /// @param values The input array containing the output calculated data.
    /// @param ifluidphase The index of the fluid phase.
    auto density(const double* vals, const double* ddt, const double* ddp, const double* ddbe) const -> void;

    /// Calculate the equilibrium-sensitivity of the density field of a fluid phase with respect to temperature.
    /// @param values The input array containing the output calculated data.
    /// @param ifluidphase The index of the fluid phase.
    auto densitySensitivityT(const double* values, unsigned ifluidphase) const -> void;

    /// Calculate the equilibrium-sensitivity of the density field of a fluid phase with respect to pressure.
    /// @param values The input array containing the output calculated data.
    /// @param ifluidphase The index of the fluid phase.
    auto densitySensitivityP(const double* values, unsigned ifluidphase) const -> void;

    /// Calculate the equilibrium-sensitivity of the density field of a fluid phase with respect to the molar amount of an element.
    /// @param values The input array containing the output calculated data.
    /// @param ifluidphase The index of the fluid phase.
    /// @param iequilibriumelement The index of the equilibrium element.
    auto densitySensitivityBe(const double* values, unsigned ifluidphase, unsigned iequilibriumelement) const -> void;

    auto densitySensitivityNk(const double* values, unsigned ifluidphase, unsigned iequilibriumelement) const -> void;
};

} // namespace Reaktoro

