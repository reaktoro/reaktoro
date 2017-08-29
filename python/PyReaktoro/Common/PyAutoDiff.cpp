// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "PyAutoDiff.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

auto export_ThermoScalar() -> void
{
    py::class_<ThermoScalar>("ThermoScalar")
        .def(py::init<>())
        .def(py::init<double>())
        .def(py::init<double, double, double>())
        .def_readwrite("val", &ThermoScalar::val)
        .def_readwrite("ddT", &ThermoScalar::ddT)
        .def_readwrite("ddP", &ThermoScalar::ddP)
        ;
}

auto export_ThermoVector() -> void
{
    py::class_<ThermoVector>("ThermoVector")
        .def(py::init<>())
        .def(py::init<Index>())
        .def(py::init<Index, double>())
        .def(py::init<const Vector&, const Vector&, const Vector&>())
        .def_readwrite("val", &ThermoVector::val)
        .def_readwrite("ddT", &ThermoVector::ddT)
        .def_readwrite("ddP", &ThermoVector::ddP)
        ;
}

auto export_ChemicalScalar() -> void
{
    py::class_<ChemicalScalar>("ChemicalScalar")
        .def(py::init<>())
        .def(py::init<Index>())
        .def(py::init<Index, double>())
        .def(py::init<double, double, double, const Vector&>())
        .def_readwrite("val", &ChemicalScalar::val)
        .def_readwrite("ddT", &ChemicalScalar::ddT)
        .def_readwrite("ddP", &ChemicalScalar::ddP)
        .def_readwrite("ddn", &ChemicalScalar::ddn)
        ;
}

auto export_ChemicalVector() -> void
{
    py::class_<ChemicalVector>("ChemicalVector")
        .def(py::init<>())
        .def(py::init<Index>())
        .def(py::init<Index, Index>())
        .def(py::init<const Vector&, const Vector&, const Vector&, const Matrix&>())
        .def_readwrite("val", &ChemicalVector::val)
        .def_readwrite("ddT", &ChemicalVector::ddT)
        .def_readwrite("ddP", &ChemicalVector::ddP)
        .def_readwrite("ddn", &ChemicalVector::ddn)
        ;
}

auto export_Temperature() -> void
{
    py::class_<Temperature, py::bases<ThermoScalar>>("Temperature")
        .def(py::init<>())
        .def(py::init<double>())
        ;

    py::implicitly_convertible<Temperature, double>();
    py::implicitly_convertible<double, Temperature>();
}

auto export_Pressure() -> void
{
    py::class_<Pressure, py::bases<ThermoScalar>>("Pressure")
        .def(py::init<>())
        .def(py::init<double>())
        ;

    py::implicitly_convertible<Pressure, double>();
    py::implicitly_convertible<double, Pressure>();
}

auto export_AutoDiff() -> void
{
    export_ThermoScalar();
    export_ThermoVector();
    export_ChemicalScalar();
    export_ChemicalVector();
    export_Temperature();
    export_Pressure();
}

} // namespace Reaktoro

