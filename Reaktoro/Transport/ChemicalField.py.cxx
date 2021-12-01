// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Transport/ChemicalField.hpp>
using namespace Reaktoro;

auto ChemicalField_setitem(ChemicalField& self, Index i, const ChemicalState& state) -> void
{
    self[i] = state;
}

auto ChemicalField_getitem(const ChemicalField& self, Index i) -> const ChemicalState&
{
    return self[i];
}

void exportChemicalField(py::module& m)
{
    py::class_<ChemicalField>(m, "ChemicalField")
        .def(py::init<Index, const ChemicalSystem&>())
        .def(py::init<Index, const ChemicalState&>())
        .def("size", &ChemicalField::size)
        .def("set", &ChemicalField::set)
        .def("temperature", &ChemicalField::temperature)
        .def("pressure", &ChemicalField::pressure)
        .def("elementAmounts", &ChemicalField::elementAmounts)
        .def("output", &ChemicalField::output)
        .def("__setitem__", ChemicalField_setitem)
        .def("__getitem__", ChemicalField_getitem, py::return_value_policy::reference_internal)
        ;
}
