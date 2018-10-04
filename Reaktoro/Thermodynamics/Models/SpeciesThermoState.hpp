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
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

/// Describe the thermodynamic state of a species
struct SpeciesThermoState
{
    /// The apparent standard molar Gibbs free energy @f$\Delta G_{f}^{\circ}@f$ of the species (in units of J/mol)
    ThermoScalar gibbs_energy;

    /// The apparent standard molar Helmholtz free energy @f$\Delta A_{f}^{\circ}@f$ of the species (in units of J/mol)
    ThermoScalar helmholtz_energy;

    /// The apparent standard molar internal energy @f$\Delta U_{f}^{\circ}@f$ of the species (in units of J/mol)
    ThermoScalar internal_energy;

    /// The apparent standard molar enthalpy @f$\Delta H_{f}^{\circ}@f$ of the species (in units of J/mol)
    ThermoScalar enthalpy;

    /// The standard molar entropy @f$ S^{\circ}@f$ of the species (in units of J/K)
    ThermoScalar entropy;

    /// The standard molar volume @f$ V^{\circ}@f$ of the species (in units of m3/mol)
    ThermoScalar volume;

    /// The standard molar isobaric heat capacity @f$ C_{P}^{\circ}@f$ of the species (in units of J/(mol K))
    ThermoScalar heat_capacity_cp;

    /// The standard molar isochoric heat capacity @f$ C_{V}^{\circ}@f$ of the species (in units of J/(mol K))
    ThermoScalar heat_capacity_cv;
};

} // namespace Reaktoro
