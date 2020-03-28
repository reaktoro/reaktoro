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

namespace Reaktoro {

/// The activity and excess thermodynamic properties of a phase.
/// @see ActivityModelFn, StandardThermoModelFn, StandardThermoProps
struct ActivityProps
{
    /// The activity coefficients (natural log) of the species in the phase.
    VectorXd ln_activity_coefficients;

    /// The activities (natural log) of the species in the phase.
    VectorXd ln_activities;

    /// The partial molar volumes of the species in the phase (in units of m3/mol).
    VectorXd partial_molar_volumes;

    /// The molar volume of the phase (in units of m3/mol).
    real molar_volume;

    /// The excess molar Gibbs energy of the phase (in units of J/mol).
    real excess_molar_gibbs_energy;

    /// The excess molar enthalpy of the phase (in units of J/mol).
    real excess_molar_enthalpy;

    /// The excess molar isobaric heat capacity of the phase (in units of J/(mol*K)).
    real excess_molar_heat_capacity_cp;

    /// The excess molar isochoric heat capacity of the phase (in units of J/(mol*K)).
    real excess_molar_heat_capacity_cv;
};

/// The function type for the activity model of a phase.
using ActivityModelFn = std::function<ActivityProps(Temperature, Pressure, VectorXrConstRef)>;

} // namespace Reaktoro
