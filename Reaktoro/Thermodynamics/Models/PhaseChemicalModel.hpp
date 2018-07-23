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
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>

namespace Reaktoro {

/// The result of the chemical model function that calculates the chemical properties of a phase.
struct PhaseChemicalModelResult
{
    /// Construct a default PhaseChemicalModelResult instance
    PhaseChemicalModelResult();

    /// Construct a PhaseChemicalModelResult instance with allocated memory
    explicit PhaseChemicalModelResult(unsigned nspecies);

    /// Resize this PhaseChemicalModelResult with a given number of species
    auto resize(unsigned nspecies) -> void;

    /// The number of species in the phase.
    unsigned num_species = 0;

    /// The natural log of the activity coefficients of the species.
    ChemicalVector ln_activity_coefficients;

    /// The natural log of the activity constants of the species.
    ThermoVector ln_activity_constants;

    /// The natural log of the activities of the species.
    ChemicalVector ln_activities;

    /// The molar volume of the phase (in units of m3/mol).
    ChemicalScalar molar_volume;

    /// The residual molar Gibbs energy of the phase w.r.t. to its ideal state (in units of J/mol).
    ChemicalScalar residual_molar_gibbs_energy;

    /// The residual molar enthalpy of the phase w.r.t. to its ideal state (in units of J/mol).
    ChemicalScalar residual_molar_enthalpy;

    /// The residual molar isobaric heat capacity of the phase w.r.t. to its ideal state (in units of J/(mol*K)).
    ChemicalScalar residual_molar_heat_capacity_cp;

    /// The residual molar isochoric heat capacity of the phase w.r.t. to its ideal state (in units of J/(mol*K)).
    ChemicalScalar residual_molar_heat_capacity_cv;
};

/// The signature of the chemical model function that calculates the chemical properties of a phase.
using PhaseChemicalModel = std::function<PhaseChemicalModelResult(double, double, const Vector&)>;

} // namespace Reaktoro
