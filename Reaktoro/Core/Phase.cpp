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
#include <Reaktoro/Core/PhaseProperties.hpp>

namespace Reaktoro {

struct Phase::Impl
{
    /// The name of the phase
    std::string name;

    /// The list of Species instances defining the phase
    std::vector<Species> species;

    /// The list of Element instances in the phase
    std::vector<Element> elements;

    /// The function that calculates the thermodynamic properties of the phase and its species
    PhaseThermoModel model;

    auto standardPartialMolarGibbsEnergies(double T, double P) const -> ThermoVector
    {
        const unsigned num_species = species.size();
        ThermoVector res(num_species);
        for(unsigned i = 0; i < num_species; ++i)
            res.row(i) = species[i].standardPartialMolarGibbsEnergy(T, P);
        return res;
    }

    auto standardPartialMolarEnthalpies(double T, double P) const -> ThermoVector
    {
        const unsigned num_species = species.size();
        ThermoVector res(num_species);
        for(unsigned i = 0; i < num_species; ++i)
            res.row(i) = species[i].standardPartialMolarEnthalpy(T, P);
        return res;
    }

    auto standardPartialMolarVolumes(double T, double P) const -> ThermoVector
    {
        const unsigned num_species = species.size();
        ThermoVector res(num_species);
        for(unsigned i = 0; i < num_species; ++i)
            res.row(i) = species[i].standardPartialMolarVolume(T, P);
        return res;
    }

    auto standardPartialMolarHeatCapacitiesConstP(double T, double P) const -> ThermoVector
    {
        const unsigned num_species = species.size();
        ThermoVector res(num_species);
        for(unsigned i = 0; i < num_species; ++i)
            res.row(i) = species[i].standardPartialMolarHeatCapacityConstP(T, P);
        return res;
    }

    auto standardPartialMolarHeatCapacitiesConstV(double T, double P) const -> ThermoVector
    {
        const unsigned num_species = species.size();
        ThermoVector res(num_species);
        for(unsigned i = 0; i < num_species; ++i)
            res.row(i) = species[i].standardPartialMolarHeatCapacityConstV(T, P);
        return res;
    }

    auto molarFractions(const Vector& n) const -> ChemicalVector
    {
        const unsigned nspecies = numSpecies();
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

    auto properties(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& n) const -> PhaseProperties
    {
        // The thermodynamic properties of the phase
        PhaseProperties prop;

        // Set temperature, pressure and composition
        prop.T = T;
        prop.P = P;
        prop.n = n;

        // Calculate the molar fractions of the species
        prop.x = molarFractions(n.val);

        // Calculate the standard thermodynamic properties of the species
        prop.standard_partial_molar_gibbs_energies     = standardPartialMolarGibbsEnergies(T.val, P.val);
        prop.standard_partial_molar_enthalpies         = standardPartialMolarEnthalpies(T.val, P.val);
        prop.standard_partial_molar_volumes            = standardPartialMolarVolumes(T.val, P.val);
        prop.standard_partial_molar_heat_capacities_cp = standardPartialMolarHeatCapacitiesConstP(T.val, P.val);
        prop.standard_partial_molar_heat_capacities_cv = standardPartialMolarHeatCapacitiesConstV(T.val, P.val);

        // Calculate the thermodynamic properties of mixing
        auto res = model(T, P, n);

        // Set the thermodynamic properties of the phase and its species
        prop.molar_gibbs_energy       = res.molar_gibbs_energy;
        prop.molar_enthalpy           = res.molar_enthalpy;
        prop.molar_volume             = res.molar_volume;
        prop.molar_heat_capacity_cp   = res.molar_heat_capacity_cp;
        prop.molar_heat_capacity_cv   = res.molar_heat_capacity_cv;
        prop.ln_activity_constants    = res.ln_activity_constants;
        prop.ln_activity_coefficients = res.ln_activity_coefficients;
        prop.ln_activities            = res.ln_activities;

        return res;
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
}

auto Phase::setThermoModel(const PhaseThermoModel& model) -> void
{
    pimpl->model = model;
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

auto Phase::properties(const ThermoScalar& T, const ThermoScalar& P, const ChemicalVector& n) const -> PhaseProperties
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
