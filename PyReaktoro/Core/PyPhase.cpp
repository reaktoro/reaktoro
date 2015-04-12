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

#include "PyPhase.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/Phase.hpp>

// PyReator includes
#include <PyReaktoro/Utils/PyConverters.hpp>

namespace Reaktoro {

auto export_Phase() -> void
{
    auto species1 = static_cast<const std::vector<Species>&(Phase::*)() const>(&Phase::species);
    auto species2 = static_cast<const Species&(Phase::*)(Index) const>(&Phase::species);

    py::class_<Phase>("Phase")
        .def(py::init<>())
        .def("numElements", &Phase::numElements)
        .def("numSpecies", &Phase::numSpecies)
        .def("name", &Phase::name, py::return_value_policy<py::copy_const_reference>())
        .def("elements", &Phase::elements, py::return_value_policy<py::copy_const_reference>())
        .def("species", species1, py::return_value_policy<py::copy_const_reference>())
        .def("species", species2, py::return_value_policy<py::copy_const_reference>())
        .def("concentrations", &Phase::concentrations)
        .def("activityCoefficients", &Phase::activityCoefficients)
        .def("activities", &Phase::activities)
        .def("chemicalPotentials", &Phase::chemicalPotentials)
        .def("molarVolume", &Phase::molarVolume)
        ;

    export_std_vector<Phase>("PhaseVector");
}

} // namespace Reaktoro
