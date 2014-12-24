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

// Boost includes
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

using namespace boost::python;

// Reaktor includes
#include <Reaktor/Core/ChemicalSystem.hpp>

// PyReator includes
#include <python/Utils/Converters.hpp>

namespace Reaktor {

auto exportChemicalSystem() -> void
{
    class_<ChemicalSystemData>("ChemicalSystemData")
        .def_readwrite("phases", &ChemicalSystemData::phases)
        .def_readwrite("gibbs_energies", &ChemicalSystemData::gibbs_energies)
        .def_readwrite("enthalpies", &ChemicalSystemData::enthalpies)
        .def_readwrite("helmholtz_energies", &ChemicalSystemData::helmholtz_energies)
        .def_readwrite("entropies", &ChemicalSystemData::entropies)
        .def_readwrite("volumes", &ChemicalSystemData::volumes)
        .def_readwrite("internal_energies", &ChemicalSystemData::internal_energies)
        .def_readwrite("heat_capacities_cp", &ChemicalSystemData::heat_capacities_cp)
        .def_readwrite("concentrations", &ChemicalSystemData::concentrations)
        .def_readwrite("ln_activity_coefficients", &ChemicalSystemData::ln_activity_coefficients)
        .def_readwrite("ln_activities", &ChemicalSystemData::ln_activities)
        .def_readwrite("chemical_potentials", &ChemicalSystemData::chemical_potentials)
        .def_readwrite("densities", &ChemicalSystemData::densities)
        ;

    class_<ChemicalSystem>("ChemicalSystem")
        .def(init<>())
        .def(init<const ChemicalSystemData&>())
        .def("elements", &ChemicalSystem::elements, return_value_policy<copy_const_reference>())
        .def("species", &ChemicalSystem::species, return_value_policy<copy_const_reference>())
        .def("phases", &ChemicalSystem::phases, return_value_policy<copy_const_reference>())
        .def("gibbsEnergies", &ChemicalSystem::gibbsEnergies)
        .def("enthalpies", &ChemicalSystem::enthalpies)
        .def("helmholtzEnergies", &ChemicalSystem::helmholtzEnergies)
        .def("entropies", &ChemicalSystem::entropies)
        .def("volumes", &ChemicalSystem::volumes)
        .def("internalEnergies", &ChemicalSystem::internalEnergies)
        .def("heatCapacitiesCp", &ChemicalSystem::heatCapacitiesCp)
        .def("concentrations", &ChemicalSystem::concentrations)
        .def("lnActivityCoefficients", &ChemicalSystem::lnActivityCoefficients)
        .def("lnActivities", &ChemicalSystem::lnActivities)
        .def("chemicalPotentials", &ChemicalSystem::chemicalPotentials)
        .def("densities", &ChemicalSystem::densities)
        ;

    def("formulaMatrix", formulaMatrix);
    def("balanceMatrix", balanceMatrix);
}

} // namespace Reaktor
