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

#include "PySpecies.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Core/Element.hpp>
#include <Reaktor/Core/Species.hpp>

// PyReator includes
#include <PyReaktor/Utils/PyConverters.hpp>

namespace Reaktor {

auto export_Species() -> void
{
    py::class_<SpeciesData>("SpeciesData")
        .def_readwrite("name", &SpeciesData::name)
        .def_readwrite("formula", &SpeciesData::formula)
        .def_readwrite("elements", &SpeciesData::elements)
        .def_readwrite("atoms", &SpeciesData::atoms)
        .def_readwrite("charge", &SpeciesData::charge)
        .def_readwrite("molar_mass", &SpeciesData::molar_mass)
        .def_readwrite("standard_gibbs_energy", &SpeciesData::standard_gibbs_energy)
        .def_readwrite("standard_enthalpy", &SpeciesData::standard_enthalpy)
        .def_readwrite("standard_helmholtz_energy", &SpeciesData::standard_helmholtz_energy)
        .def_readwrite("standard_entropy", &SpeciesData::standard_entropy)
        .def_readwrite("standard_volume", &SpeciesData::standard_volume)
        .def_readwrite("standard_internal_energy", &SpeciesData::standard_internal_energy)
        .def_readwrite("standard_heat_capacity", &SpeciesData::standard_heat_capacity)
        ;

    py::class_<Species>("Species")
        .def(py::init<>())
        .def(py::init<const SpeciesData&>())
        .def("name", &Species::name, py::return_value_policy<py::copy_const_reference>())
        .def("formula", &Species::formula, py::return_value_policy<py::copy_const_reference>())
        .def("elements", &Species::elements, py::return_value_policy<py::copy_const_reference>())
        .def("atoms", &Species::atoms, py::return_value_policy<py::copy_const_reference>())
        .def("charge", &Species::charge)
        .def("molarMass", &Species::molarMass)
        .def("data", &Species::data, py::return_value_policy<py::copy_const_reference>())
        .def("standardGibbsEnergy", &Species::standardGibbsEnergy)
        .def("standardHelmholtzEnergy", &Species::standardHelmholtzEnergy)
        .def("standardInternalEnergy", &Species::standardInternalEnergy)
        .def("standardEnthalpy", &Species::standardEnthalpy)
        .def("standardEntropy", &Species::standardEntropy)
        .def("standardVolume", &Species::standardVolume)
        .def("standardHeatCapacity", &Species::standardHeatCapacity)
        ;

    py::def("atoms", atoms);
    py::def("collectElements", collectElements);
    py::def("charges", charges);
    py::def("molarMasses", molarMasses);

    export_std_vector<Species>("SpeciesVector");
}

} // namespace Reaktor
