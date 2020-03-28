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
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

/// The result of a chemical model function that calculates the chemical properties of species.
template<typename ScalarType, typename VectorType>
struct PhaseChemicalModelResultBase
{
    /// The natural log of the activity coefficients of the species in the phase.
    VectorType ln_activity_coefficients;

    /// The natural log of the activities of the species in the phase.
    VectorType ln_activities;

    /// The partial molar volumes of the species in the phase (in units of m3/mol).
    VectorType partial_molar_volumes;

    /// The molar volume of the phase (in units of m3/mol).
    ScalarType molar_volume;

    /// The residual molar Gibbs energy of the phase w.r.t. to its ideal state (in units of J/mol).
    ScalarType residual_molar_gibbs_energy;

    /// The residual molar enthalpy of the phase w.r.t. to its ideal state (in units of J/mol).
    ScalarType residual_molar_enthalpy;

    /// The residual molar isobaric heat capacity of the phase w.r.t. to its ideal state (in units of J/(mol*K)).
    ScalarType residual_molar_heat_capacity_cp;

    /// The residual molar isochoric heat capacity of the phase w.r.t. to its ideal state (in units of J/(mol*K)).
    ScalarType residual_molar_heat_capacity_cv;
};

/// The chemical properties of the species in a phase.
using PhaseChemicalModelResult = PhaseChemicalModelResultBase<ChemicalScalarRef, VectorXdRef>;

/// The chemical properties of the species in a phase (constant).
using PhaseChemicalModelResultConst = PhaseChemicalModelResultBase<ChemicalScalarConstRef, VectorXdConstRef>;

/// The signature of the chemical model function that calculates the chemical properties of the species in a phase.
using PhaseChemicalModel = std::function<void(PhaseChemicalModelResult&, Temperature, Pressure, VectorConstRef)>;

} // namespace Reaktoro
