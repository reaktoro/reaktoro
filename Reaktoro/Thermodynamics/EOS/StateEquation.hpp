// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

// C++ includes
#include <functional>

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>

namespace Reaktoro {

/// The result of the equation of state function that calculates the thermodynamic properties of a phase.
/// @see StateEquation
struct StateEquationResult
{
    // Construct a default StateEquationResult instance
    StateEquationResult();

    // Construct a StateEquationResult instance
    explicit StateEquationResult(unsigned nspecies);

    /// The natural log of the activity coefficients of the species.
    ChemicalVector ln_activity_coefficients;

    /// The natural log of the activities of the species.
    ChemicalVector ln_activities;

    /// The residual molar Gibbs energy of the phase w.r.t. to its ideal state (in units of J/mol).
    ChemicalScalar residual_molar_gibbs_energy;

    /// The residual molar enthalpy of the phase w.r.t. to its ideal state (in units of J/mol).
    ChemicalScalar residual_molar_enthalpy;

    /// The residual molar volume of the phase w.r.t. to its ideal state (in units of m3/mol).
    ChemicalScalar residual_molar_volume;

    /// The residual molar isobaric heat capacity of the phase w.r.t. to its ideal state (in units of J/(mol*K)).
    ChemicalScalar residual_molar_heat_capacity_cp;

    /// The residual molar isochoric heat capacity of the phase w.r.t. to its ideal state (in units of J/(mol*K)).
    ChemicalScalar residual_molar_heat_capacity_cv;
};

/// The signature for an equation of state function that calculates the thermodynamic properties of a phase.
/// @param T The temperature of the phase (in unints of K)
/// @param P The pressure of the phase (in unints of Pa)
/// @param n The composition of the phase (in unints of mol)
/// @return The thermodynamic properties of the phase.
/// @see StateEquationResult
using StateEquation = std::function<StateEquationResult(double, double, const Vector&)>;

} // namespace Reaktoro
