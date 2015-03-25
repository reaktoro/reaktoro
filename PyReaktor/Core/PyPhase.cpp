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

#include "PyPhase.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/Phase.hpp>

// PyReator includes
#include <PyReaktor/Utils/PyConverters.hpp>

namespace Reaktor {

auto export_Phase() -> void
{
    py::class_<PhaseData>("PhaseData")
        .def_readwrite("name", &PhaseData::name)
        .def_readwrite("species", &PhaseData::species)
        .def_readwrite("concentrations", &PhaseData::concentrations)
        .def_readwrite("ln_activity_coefficients", &PhaseData::ln_activity_coefficients)
        .def_readwrite("ln_activities", &PhaseData::ln_activities)
        .def_readwrite("chemical_potentials", &PhaseData::chemical_potentials)
        .def_readwrite("molar_volume", &PhaseData::molar_volume)
        ;

    py::class_<Phase>("Phase")
        .def(py::init<>())
        .def(py::init<const PhaseData&>())
        .def("numElements", &Phase::numElements)
        .def("numSpecies", &Phase::numSpecies)
        .def("name", &Phase::name, py::return_value_policy<py::copy_const_reference>())
        .def("elements", &Phase::elements, py::return_value_policy<py::copy_const_reference>())
        .def("species", &Phase::species, py::return_value_policy<py::copy_const_reference>())
        .def("data", &Phase::data, py::return_value_policy<py::copy_const_reference>())
        .def("concentrations", &Phase::concentrations)
        .def("activityCoefficients", &Phase::activityCoefficients)
        .def("activities", &Phase::activities)
        .def("chemicalPotentials", &Phase::chemicalPotentials)
        .def("molarVolume", &Phase::molarVolume)
        ;

    export_std_vector<Phase>("PhaseVector");
}

} // namespace Reaktor
