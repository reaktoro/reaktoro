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

#include "PyChemicalVector.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalVector.hpp>

namespace Reaktoro {

auto export_ChemicalVector() -> void
{
    py::class_<ChemicalVector>("ChemicalVector")
        .def(py::init<>())
        .def(py::init<unsigned, unsigned>())
        .def(py::init<const Vector&, const Vector&, const Vector&, const Matrix&>())
        .def_readwrite("val", &ChemicalVector::val)
        .def_readwrite("ddt", &ChemicalVector::ddt)
        .def_readwrite("ddp", &ChemicalVector::ddp)
        .def_readwrite("ddn", &ChemicalVector::ddn)
        ;

//    py::class_<ChemicalVectorRow>("ChemicalVectorRow", py::no_init)
//        .def(py::init<ChemicalVector&, unsigned>())
//        .def("assign", assignChemicalVectorRow)
//        .add_property("val", &ChemicalVectorRow::val)
//        .add_property("ddt", &ChemicalVectorRow::ddt)
//        .add_property("ddp", &ChemicalVectorRow::ddp)
//        .add_property("ddn", &ChemicalVectorRow::ddn)
//        ;
//
//    py::class_<ChemicalVectorConstRow>("ChemicalVectorConstRow", py::no_init)
//        .def(py::init<const ChemicalVector&, unsigned>())
//        .def_readonly("val", &ChemicalVectorConstRow::val)
//        .def_readonly("ddt", &ChemicalVectorConstRow::ddt)
//        .def_readonly("ddp", &ChemicalVectorConstRow::ddp)
//        .def_readonly("ddn", &ChemicalVectorConstRow::ddn)
//        ;
}

} // namespace Reaktoro

