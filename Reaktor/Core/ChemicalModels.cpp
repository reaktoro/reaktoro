// Reaktor is a C++ library for computational reaction modelling.
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

#include "ChemicalModels.hpp"

namespace Reaktor {

struct ChemicalModels::Impl
{
    /// The function for the apparent standard molar Gibbs free energies of the species (in units of J/mol).
    ThermoVectorFunction standard_gibbs_energy_fn;

    /// The function for the apparent standard molar Helmholtz free energies of the species (in units of J/mol).
    ThermoVectorFunction standard_helmholtz_energy_fn;

    /// The function for the apparent standard molar internal energies of the species (in units of J/mol).
    ThermoVectorFunction standard_internal_energy_fn;

    /// The function for the apparent standard molar enthalpies of the species (in units of J/mol).
    ThermoVectorFunction standard_enthalpy_fn;

    /// The function for the standard molar entropies of the species (in units of J/K).
    ThermoVectorFunction standard_entropy_fn;

    /// The function for the standard molar volumes of the species (in units of m3/mol).
    ThermoVectorFunction standard_volume_fn;

    /// The function for the standard molar isobaric heat capacity of the species (in units of J/(mol*K)).
    ThermoVectorFunction standard_heat_capacity_fn;

    /// The function for the concentrations of the species (no uniform units).
    ChemicalVectorFunction concentration_fn;

    /// The function for the natural log of the activity coefficients of the species.
    ChemicalVectorFunction activity_coefficient_fn;

    /// The function for the natural log of the activities of the species.
    ChemicalVectorFunction activity_fn;

    /// The function for the chemical potentials of the species (in units of J/mol).
    ChemicalVectorFunction chemical_potential_fn;

    /// The function for the molar volumes of the phases (in units of m3/mol).
    ChemicalVectorFunction phase_molar_volume_fn;
};

ChemicalModels::ChemicalModels()
: pimpl(new Impl())
{}

ChemicalModels::ChemicalModels(const ChemicalModels& other)
: pimpl(new Impl(*other.pimpl))
{}

ChemicalModels::~ChemicalModels()
{}

auto ChemicalModels::operator=(ChemicalModels other) -> ChemicalModels&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ChemicalModels::setStandardGibbsEnergyFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_gibbs_energy_fn = function;
}

auto ChemicalModels::setStandardEnthalpyFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_helmholtz_energy_fn = function;
}

auto ChemicalModels::setStandardHelmholtzEnergyFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_internal_energy_fn = function;
}

auto ChemicalModels::setStandardInternalEnergyFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_enthalpy_fn = function;
}

auto ChemicalModels::setStandardEntropyFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_entropy_fn = function;
}

auto ChemicalModels::setStandardVolumeFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_volume_fn = function;
}

auto ChemicalModels::setStandardHeatCapacityFunction(const ThermoVectorFunction& function) -> void
{
    pimpl->standard_heat_capacity_fn = function;
}

auto ChemicalModels::setConcentrationFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->concentration_fn = function;
}

auto ChemicalModels::setActivityCoefficientFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->activity_coefficient_fn = function;
}

auto ChemicalModels::setActivityFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->activity_fn = function;
}

auto ChemicalModels::setChemicalPotentialFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->chemical_potential_fn = function;
}

auto ChemicalModels::setPhaseMolarVolumeFunction(const ChemicalVectorFunction& function) -> void
{
    pimpl->phase_molar_volume_fn = function;
}

auto ChemicalModels::standardGibbsEnergyFunction() const -> const ThermoVectorFunction&
{
    return pimpl->standard_gibbs_energy_fn;
}

auto ChemicalModels::standardEnthalpyFunction() const -> const ThermoVectorFunction&
{
    return pimpl->standard_helmholtz_energy_fn;
}

auto ChemicalModels::standardHelmholtzEnergyFunction() const -> const ThermoVectorFunction&
{
    return pimpl->standard_internal_energy_fn;
}

auto ChemicalModels::standardInternalEnergyFunction() const -> const ThermoVectorFunction&
{
    return pimpl->standard_enthalpy_fn;
}

auto ChemicalModels::standardEntropyFunction() const -> const ThermoVectorFunction&
{
    return pimpl->standard_entropy_fn;
}

auto ChemicalModels::standardVolumeFunction() const -> const ThermoVectorFunction&
{
    return pimpl->standard_volume_fn;
}

auto ChemicalModels::standardHeatCapacityFunction() const -> const ThermoVectorFunction&
{
    return pimpl->standard_heat_capacity_fn;
}

auto ChemicalModels::concentrationFunction() const -> const ChemicalVectorFunction&
{
    return pimpl->concentration_fn;
}

auto ChemicalModels::activityCoefficientFunction() const -> const ChemicalVectorFunction&
{
    return pimpl->activity_coefficient_fn;
}

auto ChemicalModels::activityFunction() const -> const ChemicalVectorFunction&
{
    return pimpl->activity_fn;
}

auto ChemicalModels::chemicalPotentialFunction() const -> const ChemicalVectorFunction&
{
    return pimpl->chemical_potential_fn;
}

auto ChemicalModels::phaseMolarVolumeFunction() const -> const ChemicalVectorFunction&
{
    return pimpl->phase_molar_volume_fn;
}

} // namespace Reaktor
