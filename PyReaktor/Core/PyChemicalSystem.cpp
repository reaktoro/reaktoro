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

#include "PyChemicalSystem.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Core/ChemicalSystem.cpp>
#include <Reaktor/Interfaces/Gems.hpp>

namespace Reaktor {
namespace {

auto createChemicalSystem(const Gems& gems) -> boost::shared_ptr<ChemicalSystem>
{
	return boost::make_shared<ChemicalSystem>(gems);
}

} // namespace

auto export_ChemicalSystem() -> void
{
    py::class_<ChemicalSystemData>("ChemicalSystemData")
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

    using element_ftype1 = const Element&(ChemicalSystem::*)(Index) const;
    using element_ftype2 = const Element&(ChemicalSystem::*)(std::string) const;

    using species_ftype1 = const std::vector<Species>&(ChemicalSystem::*)() const;
    using species_ftype2 = const Species&(ChemicalSystem::*)(Index) const;
    using species_ftype3 = const Species&(ChemicalSystem::*)(std::string) const;

    auto species_all =
    py::class_<ChemicalSystem>("ChemicalSystem")
        .def(py::init<>())
        .def(py::init<const ChemicalSystemData&>())
        .def("__init__", py::make_constructor(createChemicalSystem))
        .def("elements", &ChemicalSystem::elements, py::return_value_policy<py::copy_const_reference>())
        .def("element", static_cast<element_ftype1>(&ChemicalSystem::element), py::return_value_policy<py::copy_const_reference>())
        .def("element", static_cast<element_ftype2>(&ChemicalSystem::element), py::return_value_policy<py::copy_const_reference>())
        .def("species", static_cast<species_ftype1>(&ChemicalSystem::species), py::return_value_policy<py::copy_const_reference>())
        .def("species", static_cast<species_ftype2>(&ChemicalSystem::species), py::return_value_policy<py::copy_const_reference>())
        .def("species", static_cast<species_ftype3>(&ChemicalSystem::species), py::return_value_policy<py::copy_const_reference>())
        .def("phases", &ChemicalSystem::phases, py::return_value_policy<py::copy_const_reference>())
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

    py::def("formulaMatrix", formulaMatrix);
    py::def("balanceMatrix", balanceMatrix);
}

} // namespace Reaktor
