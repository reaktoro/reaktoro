// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "ThermoProperties.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Thermodynamics/Models/ThermoModel.hpp>

namespace Reaktoro {

struct ThermoProperties::Impl
{
    /// The chemical system
    ChemicalSystem system;

    /// The number of species in the system
    Index num_species = 0;

    /// The number of phases in the system
    Index num_phases = 0;

    /// The results of the evaluation of the PhaseThermoModel functions of each phase.
    ThermoModelResult tres;

    /// The temperature of the system (in units of K)
    Temperature T;

    /// The pressure of the system (in units of Pa)
    Pressure P;

    /// Construct a default Impl instance
    Impl()
    {}

    /// Construct a Impl instance with given ChemicalSystem
    Impl(const ChemicalSystem& system)
    : system(system), num_species(system.numSpecies()),
      num_phases(system.numPhases()), tres(num_species), T(298.15), P(1e-5)
    {}

    /// Update the thermodynamic properties of the chemical system.
    auto update(double T_, double P_) -> void
    {
        // Update both temperature and pressure
        T = T_;
        P = P_;

        // Update the thermodynamic properties of each phase
        Index ispecies = 0;
        for(Index iphase = 0; iphase < num_phases; ++iphase)
        {
            const auto nspecies = system.numSpeciesInPhase(iphase);
            auto tp = tres.phaseProperties(ispecies, nspecies);
            system.phase(iphase).thermoModel()(tp, T, P);
            ispecies += nspecies;
        }
    }

    /// Return the standard partial molar Gibbs energies of the species (in units of J/mol).
    auto standardPartialMolarGibbsEnergies() const -> ThermoVector
    {
        return tres.standardPartialMolarGibbsEnergies();
    }

    /// Return the standard partial molar enthalpies of the species (in units of J/mol).
    auto standardPartialMolarEnthalpies() const -> ThermoVector
    {
        return tres.standardPartialMolarEnthalpies();
    }

    /// Return the standard partial molar volumes of the species (in units of m3/mol).
    auto standardPartialMolarVolumes() const -> ThermoVector
    {
        return tres.standardPartialMolarVolumes();
    }

    /// Return the standard partial molar entropies of the species (in units of J/(mol*K)).
    auto standardPartialMolarEntropies() const -> ThermoVector
    {
        const auto& G = standardPartialMolarGibbsEnergies();
        const auto& H = standardPartialMolarEnthalpies();
        return (H - G)/T;
    }

    /// Return the standard partial molar internal energies of the species (in units of J/mol).
    auto standardPartialMolarInternalEnergies() const -> ThermoVector
    {
        const auto& H = standardPartialMolarEnthalpies();
        const auto& V = standardPartialMolarVolumes();
        return H - P*V;
    }

    /// Return the standard partial molar Helmholtz energies of the species (in units of J/mol).
    auto standardPartialMolarHelmholtzEnergies() const -> ThermoVector
    {
        const auto& G = standardPartialMolarGibbsEnergies();
        const auto& V = standardPartialMolarVolumes();
        return G - P*V;
    }

    /// Return the standard partial molar isobaric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
    {
        return tres.standardPartialMolarHeatCapacitiesConstP();
    }

    /// Return the standard partial molar isochoric heat capacities of the species (in units of J/(mol*K)).
    auto standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
    {
        return tres.standardPartialMolarHeatCapacitiesConstV();
    }
};

ThermoProperties::ThermoProperties()
: pimpl(new Impl())
{}

ThermoProperties::ThermoProperties(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

ThermoProperties::ThermoProperties(const ThermoProperties& other)
: pimpl(new Impl(*other.pimpl))
{}

ThermoProperties::~ThermoProperties()
{}

auto ThermoProperties::operator=(ThermoProperties other) -> ThermoProperties&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ThermoProperties::update(double T, double P) -> void
{
    pimpl->update(T, P);
}

auto ThermoProperties::temperature() const -> double
{
    return pimpl->T.val;
}

auto ThermoProperties::pressure() const -> double
{
    return pimpl->P.val;
}

auto ThermoProperties::standardPartialMolarGibbsEnergies() const -> ThermoVector
{
    return pimpl->standardPartialMolarGibbsEnergies();
}

auto ThermoProperties::standardPartialMolarEnthalpies() const -> ThermoVector
{
    return pimpl->standardPartialMolarEnthalpies();
}

auto ThermoProperties::standardPartialMolarVolumes() const -> ThermoVector
{
    return pimpl->standardPartialMolarVolumes();
}

auto ThermoProperties::standardPartialMolarEntropies() const -> ThermoVector
{
    return pimpl->standardPartialMolarEntropies();
}

auto ThermoProperties::standardPartialMolarInternalEnergies() const -> ThermoVector
{
    return pimpl->standardPartialMolarInternalEnergies();
}

auto ThermoProperties::standardPartialMolarHelmholtzEnergies() const -> ThermoVector
{
    return pimpl->standardPartialMolarHelmholtzEnergies();
}

auto ThermoProperties::standardPartialMolarHeatCapacitiesConstP() const -> ThermoVector
{
    return pimpl->standardPartialMolarHeatCapacitiesConstP();
}

auto ThermoProperties::standardPartialMolarHeatCapacitiesConstV() const -> ThermoVector
{
    return pimpl->standardPartialMolarHeatCapacitiesConstV();
}

} // namespace Reaktoro
