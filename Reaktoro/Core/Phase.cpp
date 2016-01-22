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

#include "Phase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

struct Phase::Impl
{
    /// The name of the phase
    std::string name;

    /// The list of Species instances defining the phase
    std::vector<Species> species;

    /// The list of Element instances in the phase
    std::vector<Element> elements;

    /// The function that calculates the standard thermodynamic properties of the phase and its species
    PhaseThermoModel thermo_model;

    /// The function that calculates the chemical properties of the phase and its species
    PhaseChemicalModel chemical_model;

    // The molar masses of the species
    Vector molar_masses;

    auto properties(double T, double P) const -> ThermoProperties
    {
        // The thermodynamic properties of the species
        ThermoProperties prop;

        // Set temperature, pressure and composition
        prop.T = T;
        prop.P = P;

        // Calculate the standard thermodynamic properties of the phase
        auto res = thermo_model(T, P);

        // Set the standard thermodynamic properties of the species in the phase
        prop.standard_partial_molar_gibbs_energies = res.standard_partial_molar_gibbs_energies;
        prop.standard_partial_molar_enthalpies = res.standard_partial_molar_enthalpies;
        prop.standard_partial_molar_volumes = res.standard_partial_molar_volumes;
        prop.standard_partial_molar_heat_capacities_cp = res.standard_partial_molar_heat_capacities_cp;
        prop.standard_partial_molar_heat_capacities_cv = res.standard_partial_molar_heat_capacities_cv;

        return prop;
    }

    auto properties(double T, double P, const Vector& n) const -> PhaseChemicalProperties
    {
        // The chemical properties of the phase and its species
        PhaseChemicalProperties prop;

        // Set temperature, pressure and composition
        prop.T = T;
        prop.P = P;
        prop.n = n;

        // Calculate the standard thermodynamic properties of the species
        ThermoProperties tp = properties(T, P);

        // Set the standard thermodynamic properties of the species in the phase
        prop.standard_partial_molar_gibbs_energies = tp.standard_partial_molar_gibbs_energies;
        prop.standard_partial_molar_enthalpies = tp.standard_partial_molar_enthalpies;
        prop.standard_partial_molar_volumes = tp.standard_partial_molar_volumes;
        prop.standard_partial_molar_heat_capacities_cp = tp.standard_partial_molar_heat_capacities_cp;
        prop.standard_partial_molar_heat_capacities_cv = tp.standard_partial_molar_heat_capacities_cv;

        // Calculate the molar fractions of the species
        prop.molar_fractions = molarFractions(prop.n);

        // Calculate the ideal contribution for the thermodynamic properties of the phase
        const ChemicalVector& x = prop.molar_fractions;
        prop.phase_molar_gibbs_energy     = sum(x % tp.standardPartialMolarGibbsEnergies());
        prop.phase_molar_enthalpy         = sum(x % tp.standardPartialMolarEnthalpies());
        prop.phase_molar_volume           = sum(x % tp.standardPartialMolarVolumes());
        prop.phase_molar_heat_capacity_cp = sum(x % tp.standardPartialMolarHeatCapacitiesConstP());
        prop.phase_molar_heat_capacity_cv = sum(x % tp.standardPartialMolarHeatCapacitiesConstV());

        // Calculate the chemical properties of the phase, otherwise leave it with the ideal gas/solution volume
        auto res = chemical_model(T, P, n);

        // Check if the molar volume of the phase was calculated
        if(res.molar_volume.val > 0.0)
            prop.phase_molar_volume = res.molar_volume;

        // Add the non-ideal residual contribution to the thermodynamic properties of the phase
        prop.phase_molar_gibbs_energy     += res.residual_molar_gibbs_energy;
        prop.phase_molar_enthalpy         += res.residual_molar_enthalpy;
        prop.phase_molar_heat_capacity_cp += res.residual_molar_heat_capacity_cp;
        prop.phase_molar_heat_capacity_cv += res.residual_molar_heat_capacity_cv;

        // Set the thermodynamic properties of the species
        prop.ln_activity_coefficients = res.ln_activity_coefficients;
        prop.ln_activity_constants    = res.ln_activity_constants;
        prop.ln_activities            = res.ln_activities;

        // Set the number of moles and mass of the phase
        prop.phase_moles = sum(prop.n);
        prop.phase_mass = sum(molar_masses % prop.n);

        return prop;
    }
};

Phase::Phase()
: pimpl(new Impl())
{}

Phase::Phase(const Phase& other)
: pimpl(new Impl(*other.pimpl))
{}

Phase::~Phase()
{}

auto Phase::operator=(Phase other) -> Phase&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Phase::setName(std::string name) -> void
{
    pimpl->name = name;
}

auto Phase::setSpecies(const std::vector<Species>& species) -> void
{
    pimpl->species = species;
    pimpl->molar_masses = molarMasses(species);
}

auto Phase::setThermoModel(const PhaseThermoModel& model) -> void
{
    pimpl->thermo_model = model;
}

auto Phase::setChemicalModel(const PhaseChemicalModel& model) -> void
{
    pimpl->chemical_model = model;
}

auto Phase::numElements() const -> unsigned
{
    return elements().size();
}

auto Phase::numSpecies() const -> unsigned
{
    return species().size();
}

auto Phase::name() const -> std::string
{
    return pimpl->name;
}

auto Phase::elements() const -> const std::vector<Element>&
{
    return pimpl->elements;
}

auto Phase::elements() -> std::vector<Element>&
{
    return pimpl->elements;
}

auto Phase::species() const -> const std::vector<Species>&
{
    return pimpl->species;
}

auto Phase::species() -> std::vector<Species>&
{
    return pimpl->species;
}

auto Phase::species(Index index) const -> const Species&
{
    return pimpl->species[index];
}

auto Phase::thermoModel() const -> const PhaseThermoModel&
{
    return pimpl->thermo_model;
}

auto Phase::chemicalModel() const -> const PhaseChemicalModel&
{
    return pimpl->chemical_model;
}

auto Phase::indexSpecies(std::string name) const -> Index
{
    return index(name, species());
}

auto Phase::indexSpeciesWithError(std::string name) const -> Index
{
    const Index index = indexSpecies(name);
    Assert(index < numSpecies(),
        "Could not get the index of species `" + name + "`.",
        "There is no species called `" + name + "` in the phase.");
    return index;
}

auto Phase::indexSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    return indexAny(names, species());
}

auto Phase::indexSpeciesAnyWithError(const std::vector<std::string>& names) const -> Index
{
    const Index index = indexSpeciesAny(names);
    Assert(index < numSpecies(),
        "Could not get the index of the species with "
        "any of the following names `" + join(names, ", ") + "`.",
        "There is no species in phase `" + name() + "` with any of these names.");
    return index;
}

auto Phase::properties(double T, double P) const -> ThermoProperties
{
    return pimpl->properties(T, P);
}

auto Phase::properties(double T, double P, const Vector& n) const -> PhaseChemicalProperties
{
    return pimpl->properties(T, P, n);
}

auto operator<(const Phase& lhs, const Phase& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Phase& lhs, const Phase& rhs) -> bool
{
    return lhs.name() == rhs.name();
}

} // namespace Reaktoro
