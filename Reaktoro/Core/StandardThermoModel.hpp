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
#include <functional>

// Reaktoro includes
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

/// Forward declarations
class Species;

/// The standard thermodynamic properties of a chemical species.
struct StandardThermoProps
{
    /// The standard molar Gibbs energy @f$G^{\circ}@f$ of the species (in unit of J/mol)
    ThermoScalar G0;

    /// The standard molar Helmholtz energy @f$A^{\circ}@f$ of the species (in unit of J/mol)
    ThermoScalar A0;

    /// The standard molar internal energy @f$U^{\circ}@f$ of the species (in unit of J/mol)
    ThermoScalar U0;

    /// The standard molar enthalpy @f$H^{\circ}@f$ of the species (in unit of J/mol)
    ThermoScalar H0;

    /// The standard molar entropy @f$S^{\circ}@f$ of the species (in unit of J/K)
    ThermoScalar S0;

    /// The standard molar volume @f$V^{\circ}@f$ of the species (in unit of m3/mol)
    ThermoScalar V0;

    /// The standard molar isobaric heat capacity @f$C_{P}^{\circ}@f$ of the species (in unit of J/(mol·K))
    ThermoScalar Cp0;

    /// The standard molar isochoric heat capacity @f$C_{V}^{\circ}@f$ of the species (in unit of J/(mol·K))
    ThermoScalar Cv0;
};

/// The function type for the standard thermodynamic model of a chemical species.
using StandardThermoModelFn = std::function<StandardThermoProps(Temperature, Pressure, const Species&)>;

} // namespace Reaktoro
