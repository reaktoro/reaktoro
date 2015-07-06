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

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

class PhaseThermoProperties
{
private:

    /// The temperature of the phase (in units of K)
    ThermoScalar temperature;

    /// The pressure of the phase (in units of Pa)
    ThermoScalar pressure;

    /// The composition of the species of the phase (in units of mol)
    ChemicalVector composition;

    /// The concentrations of the species of the phase ()
    ChemicalVector concentrations;

    /// The activity_factors
    ThermoVector activity_factors;

    /// The activity_coefficients
    ChemicalVector activity_coefficients;

    /// The activities
    ChemicalVector activities;

    /// The partial_molar_gibbs_energies
    ChemicalVector partial_molar_gibbs_energies;

    /// The partial_molar_enthalpies
    ChemicalVector partial_molar_enthalpies;

    /// The partial_molar_volumes
    ChemicalVector partial_molar_volumes;

public:
    auto temperature() const -> ThermoScalar;

    auto pressure() const -> ThermoScalar;

    auto composition() const -> ChemicalVector;

    auto concentrations() const -> ChemicalVector;

    auto activityFactors() const -> ThermoVector;

    auto activityCoefficients() const -> ChemicalVector;

    auto activities() const -> ChemicalVector;

    auto partialMolarGibbsEnergies() const -> ChemicalVector;

    auto partialMolarEnthalpies() const -> ChemicalVector;

    auto partialMolarVolumes() const -> ChemicalVector;

    auto molarGibbsEnergy() const -> ChemicalScalar;

    auto molarEnthalpy() const -> ChemicalScalar;

    auto molarVolume() const -> ChemicalScalar;

    auto molarHeatCapacityP() const -> ChemicalScalar;

    auto molarHeatCapacityV() const -> ChemicalScalar;

    auto molarDensity() const -> ChemicalScalar;

    auto density() const -> ChemicalScalar;

    auto moles() const -> ChemicalScalar;

    auto mass() const -> ChemicalScalar;

    auto gibbsEnergy() const -> ChemicalScalar;

    auto volume() const -> ChemicalScalar;

    auto enthalpy() const -> ChemicalScalar;

    auto heatCapacityP() const -> ChemicalScalar;

    auto heatCapacityV() const -> ChemicalScalar;
};

} // namespace Reaktoro
