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

#include "PySpecies.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/Element.hpp>
#include <Reaktoro/Core/Species.hpp>

// PyReator includes
#include <PyReaktoro/Utils/PyConverters.hpp>

namespace Reaktoro {

auto export_Species() -> void
{
    using return_const_ref = py::return_value_policy<py::copy_const_reference>;

    py::class_<Species>("Species")
        .def(py::init<>())
        .def("setName", &Species::setName)
        .def("setFormula", &Species::setFormula)
        .def("setElements", &Species::setElements)
        .def("numElements", &Species::numElements)
        .def("name", &Species::name)
        .def("formula", &Species::formula)
        .def("elements", &Species::elements, return_const_ref())
        .def("molarMass", &Species::molarMass)
        .def("elementCoefficient", &Species::elementCoefficient)
        ;

    export_std_vector<Species>("SpeciesVector");
    export_std_map<Element, double>("ElementDoubleMap");
}

} // namespace Reaktoro
