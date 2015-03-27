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

#include "ChemicalSystem.hpp"

// C++ includes
#include <iostream>
#include <iomanip>
#include <set>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/Utils.hpp>

namespace Reaktor {
namespace {

auto connectivity(const Multiphase& multiphase) -> Connectivity
{
    const auto& elements = multiphase.elements();
    const auto& species  = multiphase.species();
    const auto& phases   = multiphase.phases();

    const unsigned num_elements = elements.size();
    const unsigned num_species = species.size();
    const unsigned num_phases = phases.size();

    Connectivity c;

    c.element_to_species.resize(num_elements);
    c.species_to_elements.resize(num_species);
    for(unsigned j = 0; j < num_elements; ++j)
        for(unsigned i = 0; i < num_species; ++i)
            if(species[i].elements().count(elements[j])) {
                c.element_to_species[j].push_back(i);
                c.species_to_elements[i].push_back(j); }

    c.species_to_phase.resize(num_species);
    c.phase_to_species.resize(num_phases);
    for(unsigned k = 0; k < num_phases; ++k)
        for(unsigned i = 0; i < num_species; ++i)
            if(contains(species[i], phases[k].species())) {
                c.species_to_phase[i] = k;
                c.phase_to_species[k].push_back(i);
            }

    c.element_to_phases.resize(num_elements);
    c.phase_to_elements.resize(num_phases);
    for(unsigned k = 0; k < num_phases; ++k)
        for(unsigned j = 0; j < num_elements; ++j)
            if(contains(elements[j], phases[k].elements())) {
                c.element_to_phases[j].push_back(k);
                c.phase_to_elements[k].push_back(j);
            }

    return c;
}

} // namespace

struct ChemicalSystem::Impl
{
    /// The connectivity of the chemical system
    Connectivity connectivity;
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const Multiphase& multiphase)
: Multiphase(multiphase), pimpl(new Impl())
{
    pimpl->connectivity = Reaktor::connectivity(multiphase);
}

auto ChemicalSystem::connectivity() const -> const Connectivity&
{
    return pimpl->connectivity;
}

auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&
{
    const auto& phases = system.phases();
    for(unsigned i = 0; i < phases.size(); ++i)
    {
        out << "Phase(" << i << "): " << phases[i].name() << std::endl;
        const auto& species = phases[i].species();
        for(unsigned i = 0; i < species.size(); ++i)
        {
            const auto name = species[i].name();
            const auto idx  = system.indexSpecies(name);
            out << std::setw(5) << std::left << idx;
            out << std::setw(30) << std::left << name;
            out << std::endl;
        }
    }

    out << std::endl;

    out << "Elements:" << std::endl;
    const auto& elements = system.elements();
    for(unsigned i = 0; i < elements.size(); ++i)
        out << (i > 0 ? ", " : "") << i << ":" << elements[i].name();

    out << std::endl;

    return out;
}

} // namespace Reaktor
