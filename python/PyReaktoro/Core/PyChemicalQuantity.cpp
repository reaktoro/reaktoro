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

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalQuantity.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

void exportChemicalQuantity(py::module& m)
{
    auto update1 = static_cast<ChemicalQuantity&(ChemicalQuantity::*)(const ChemicalState&)>(&ChemicalQuantity::update);
    auto update2 = static_cast<ChemicalQuantity&(ChemicalQuantity::*)(const ChemicalState&,double)>(&ChemicalQuantity::update);

    py::class_<ChemicalQuantity>(m, "ChemicalQuantity")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ReactionSystem&>())
        .def(py::init<const ChemicalState&>())
        .def("system", &ChemicalQuantity::system, py::return_value_policy::reference_internal)
        .def("reactions", &ChemicalQuantity::reactions, py::return_value_policy::reference_internal)
        .def("state", &ChemicalQuantity::state, py::return_value_policy::reference_internal)
        .def("properties", &ChemicalQuantity::properties, py::return_value_policy::reference_internal)
        .def("rates", &ChemicalQuantity::rates, py::return_value_policy::reference_internal)
        .def("tag", &ChemicalQuantity::tag)
        .def("update", update1, py::return_value_policy::reference_internal)
        .def("update", update2, py::return_value_policy::reference_internal)
        .def("value", &ChemicalQuantity::value)
        .def("__call__", &ChemicalQuantity::value)
        ;
}

} // namespace Reaktoro
