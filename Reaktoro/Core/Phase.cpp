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

#include "Phase.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
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

    /// The standard reference state type of the phase
    PhaseReferenceState reftype = PhaseReferenceState::IdealGas;

    // The molar masses of the species
    Vector molar_masses;

    auto molarFractions(const Vector& n) const -> ChemicalVector
    {
        const unsigned nspecies = species.size();
        if(nspecies == 1)
        {
            ChemicalVector x(1, 1);
            x.val[0] = 1.0;
            return x;
        }
        ChemicalVector x(nspecies, nspecies);
        const double nt = n.sum();
        if(nt == 0.0) return x;
        x.val = n/nt;
        for(unsigned i = 0; i < nspecies; ++i)
        {
            x.ddn.row(i).fill(-x.val[i]/nt);
            x.ddn(i, i) += 1.0/nt;
        }
        return x;
    }

    auto properties(double T, double P) const -> ThermoProperties
    {
        // The thermodynamic properties of the species
        ThermoProperties prop;

        // Get a reference to the internal members of ThermoProperties
        auto& inter = prop.internal;

        // Set temperature, pressure and composition
        inter.T = ThermoScalar::Temperature(T);
        inter.P = ThermoScalar::Pressure(P);

        // Calculate the standard thermodynamic properties of the phase
        auto res = thermo_model(T, P);

        // Set the standard thermodynamic properties of the species in the phase
        inter.standard_partial_molar_gibbs_energies = res.standard_partial_molar_gibbs_energies;
        inter.standard_partial_molar_enthalpies = res.standard_partial_molar_enthalpies;
        inter.standard_partial_molar_volumes = res.standard_partial_molar_volumes;
        inter.standard_partial_molar_heat_capacities_cp = res.standard_partial_molar_heat_capacities_cp;
        inter.standard_partial_molar_heat_capacities_cv = res.standard_partial_molar_heat_capacities_cv;

        return prop;
    }

    auto properties(double T, double P, const Vector& n) const -> PhaseChemicalProperties
    {
        // The chemical properties of the phase and its species
        PhaseChemicalProperties prop;

        // Get a reference to the internal members of PhaseChemicalProperties
        auto& inter = prop.internal;

        // Set temperature, pressure and composition
        inter.T = ThermoScalar::Temperature(T);
        inter.P = ThermoScalar::Pressure(P);
        inter.n = ChemicalVector::Composition(n);

        // Calculate the molar fractions of the species
        inter.molar_fractions = molarFractions(n);

        // Calculate the standard thermodynamic properties of the species
        ThermoProperties tp = properties(T, P);

        // Calculate the ideal contribution for the thermodynamic properties of the phase
        const ChemicalVector& x = inter.molar_fractions;
        inter.phase_molar_gibbs_energy     = sum(x % tp.standardPartialMolarGibbsEnergies());
        inter.phase_molar_enthalpy         = sum(x % tp.standardPartialMolarEnthalpies());
        inter.phase_molar_volume           = sum(x % tp.standardPartialMolarVolumes());
        inter.phase_molar_heat_capacity_cp = sum(x % tp.standardPartialMolarHeatCapacitiesConstP());
        inter.phase_molar_heat_capacity_cv = sum(x % tp.standardPartialMolarHeatCapacitiesConstV());

        // Calculate the chemical properties of the phase
        auto res = chemical_model(T, P, n);

        // Add the non-ideal residual contribution to the thermodynamic properties of the phase
        inter.phase_molar_gibbs_energy     += res.residual_molar_gibbs_energy;
        inter.phase_molar_enthalpy         += res.residual_molar_enthalpy;
        inter.phase_molar_volume           += res.residual_molar_volume;
        inter.phase_molar_heat_capacity_cp += res.residual_molar_heat_capacity_cp;
        inter.phase_molar_heat_capacity_cv += res.residual_molar_heat_capacity_cv;

        // Set the thermodynamic properties of the species
        inter.ln_activity_coefficients = res.ln_activity_coefficients;
        inter.ln_activities            = res.ln_activities;

        // Set the mass of the phase
        inter.phase_mass = sum(molar_masses % inter.n);

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

auto Phase::setReferenceState(PhaseReferenceState reftype) -> void
{
    pimpl->reftype = reftype;
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

auto Phase::species() const -> const std::vector<Species>&
{
    return pimpl->species;
}

auto Phase::species(Index index) const -> const Species&
{
    return pimpl->species[index];
}

auto Phase::indexSpecies(std::string name) const -> Index
{
    return index(name, species());
}

auto Phase::referenceState() const -> PhaseReferenceState
{
    return pimpl->reftype;
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
