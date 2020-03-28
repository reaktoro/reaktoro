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

// Reaktoro includes

namespace Reaktoro {

/// Describe the thermodynamic state of a species
struct SpeciesThermoState
{
    /// The apparent standard molar Gibbs free energy @f$\Delta G_{f}^{\circ}@f$ of the species (in units of J/mol)
    real gibbs_energy;

    /// The apparent standard molar Helmholtz free energy @f$\Delta A_{f}^{\circ}@f$ of the species (in units of J/mol)
    real helmholtz_energy;

    /// The apparent standard molar internal energy @f$\Delta U_{f}^{\circ}@f$ of the species (in units of J/mol)
    real internal_energy;

    /// The apparent standard molar enthalpy @f$\Delta H_{f}^{\circ}@f$ of the species (in units of J/mol)
    real enthalpy;

    /// The standard molar entropy @f$ S^{\circ}@f$ of the species (in units of J/K)
    real entropy;

    /// The standard molar volume @f$ V^{\circ}@f$ of the species (in units of m3/mol)
    real volume;

    /// The standard molar isobaric heat capacity @f$ C_{P}^{\circ}@f$ of the species (in units of J/(mol K))
    real heat_capacity_cp;

    /// The standard molar isochoric heat capacity @f$ C_{V}^{\circ}@f$ of the species (in units of J/(mol K))
    real heat_capacity_cv;
};

} // namespace Reaktoro
