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

#include "Connectivity.hpp"

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

Connectivity::Connectivity()
{}

Connectivity::Connectivity(const ChemicalSystem& system)
{
    const auto& elements = system.elements();
    const auto& species  = system.species();
    const auto& phases   = system.phases();

    const unsigned num_elements = elements.size();
    const unsigned num_species = species.size();
    const unsigned num_phases = phases.size();

    element_to_species.resize(num_elements);
    species_to_elements.resize(num_species);
    for(unsigned j = 0; j < num_elements; ++j)
        for(unsigned i = 0; i < num_species; ++i)
            if(species[i].elements().count(elements[j])) {
                element_to_species[j].push_back(i);
                species_to_elements[i].push_back(j); }

    species_to_phase.resize(num_species);
    phase_to_species.resize(num_phases);
    for(unsigned k = 0; k < num_phases; ++k)
        for(unsigned i = 0; i < num_species; ++i)
            if(contains(species[i], phases[k].species())) {
                species_to_phase[i] = k;
                phase_to_species[k].push_back(i);
            }

    element_to_phases.resize(num_elements);
    phase_to_elements.resize(num_phases);
    for(unsigned k = 0; k < num_phases; ++k)
        for(unsigned j = 0; j < num_elements; ++j)
            if(contains(elements[j], phases[k].elements())) {
                element_to_phases[j].push_back(k);
                phase_to_elements[k].push_back(j);
            }
}

} // namespace Reaktoro
