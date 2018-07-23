// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "PyElement.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/Element.hpp>

// PyReator includes
#include <PyReaktoro/Common/PyConverters.hpp>

namespace Reaktoro {

auto export_Element() -> void
{
    py::class_<Element>("Element")
        .def(py::init<>())

        .def("setName", &Element::setName,
            "Set the name of the element.",
            py::args("self", "name"))

        .def("setMolarMass", &Element::setMolarMass,
            "Set the molar mass of the element (in units of kg/mol).")

        .def("name", &Element::name,
            "Return the name of the element.")

        .def("molarMass", &Element::molarMass,
            "Return the molar mass of the element (in units of kg/mol).")
        ;

    export_std_vector<Element>("ElementVector");
}

} // namespace Reaktoro
