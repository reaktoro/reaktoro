// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro {

/// The primary standard thermodynamic properties of a chemical species.
struct StandardThermoProps
{
    /// The standard molar Gibbs energy @f$G^{\circ}@f$ of the species (in J/mol)
    real G0 = {};

    /// The standard molar enthalpy @f$H^{\circ}@f$ of the species (in J/mol)
    real H0 = {};

    /// The standard molar volume @f$V^{\circ}@f$ of the species (in m3/mol)
    real V0 = {};

    /// The standard molar isobaric heat capacity @f$C_{P}^{\circ}@f$ of the species (in J/(mol·K))
    real Cp0 = {};

    /// The standard molar isochoric heat capacity @f$C_{V}^{\circ}@f$ of the species (in J/(mol·K))
    real Cv0 = {};
};

/// The function type for calculation of standard thermodynamic properties of species.
using StandardThermoPropsFn = std::function<StandardThermoProps(real, real, const Species&)>;

} // namespace Reaktoro
