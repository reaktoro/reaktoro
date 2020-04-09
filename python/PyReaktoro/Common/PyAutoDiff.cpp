// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Common/ThermoVector.hpp>

namespace Reaktoro {

void exportThermoScalar(py::module& m)
{
    py::class_<ThermoScalar>(m, "ThermoScalar")
        .def(py::init<>())
        .def(py::init<double>())
        .def(py::init<double, double, double>())
        .def_readwrite("val", &ThermoScalar::val)
        .def_readwrite("ddT", &ThermoScalar::ddT)
        .def_readwrite("ddP", &ThermoScalar::ddP)
        ;
}

void exportThermoVector(py::module& m)
{
    py::class_<ThermoVector>(m, "ThermoVector")
        .def(py::init<>())
        .def(py::init<Index>())
        .def(py::init<Index, double>())
        .def(py::init<VectorConstRef, VectorConstRef, VectorConstRef>())
        .def_readwrite("val", &ThermoVector::val)
        .def_readwrite("ddT", &ThermoVector::ddT)
        .def_readwrite("ddP", &ThermoVector::ddP)
        ;

    py::class_<ThermoVectorRef>(m, "ThermoVectorRef")
        .def_readwrite("val", &ThermoVectorRef::val)
        .def_readwrite("ddT", &ThermoVectorRef::ddT)
        .def_readwrite("ddP", &ThermoVectorRef::ddP)
        ;

    py::class_<ThermoVectorConstRef>(m, "ThermoVectorConstRef")
        .def_readonly("val", &ThermoVectorConstRef::val)
        .def_readonly("ddT", &ThermoVectorConstRef::ddT)
        .def_readonly("ddP", &ThermoVectorConstRef::ddP)
        ;
}

void exportChemicalScalar(py::module& m)
{
    py::class_<ChemicalScalar>(m, "ChemicalScalar")
        .def(py::init<>())
        .def(py::init<Index>())
        .def(py::init<Index, double>())
        .def(py::init<double, double, double, VectorConstRef>())
        .def_readwrite("val", &ChemicalScalar::val)
        .def_readwrite("ddT", &ChemicalScalar::ddT)
        .def_readwrite("ddP", &ChemicalScalar::ddP)
        .def_readwrite("ddn", &ChemicalScalar::ddn)
        ;
}

void exportChemicalVector(py::module& m)
{
    py::class_<ChemicalVector>(m, "ChemicalVector")
        .def(py::init<>())
        .def(py::init<Index>())
        .def(py::init<Index, Index>())
        .def(py::init<VectorConstRef, VectorConstRef, VectorConstRef, MatrixConstRef>())
        .def_readwrite("val", &ChemicalVector::val)
        .def_readwrite("ddT", &ChemicalVector::ddT)
        .def_readwrite("ddP", &ChemicalVector::ddP)
        .def_readwrite("ddn", &ChemicalVector::ddn)
        ;

    py::class_<ChemicalVectorRef>(m, "ChemicalVectorRef")
        .def_readwrite("val", &ChemicalVectorRef::val)
        .def_readwrite("ddT", &ChemicalVectorRef::ddT)
        .def_readwrite("ddP", &ChemicalVectorRef::ddP)
        .def_readwrite("ddn", &ChemicalVectorRef::ddn)
        ;

    py::class_<ChemicalVectorConstRef>(m, "ChemicalVectorConstRef")
        .def_readonly("val", &ChemicalVectorConstRef::val)
        .def_readonly("ddT", &ChemicalVectorConstRef::ddT)
        .def_readonly("ddP", &ChemicalVectorConstRef::ddP)
        .def_readonly("ddn", &ChemicalVectorConstRef::ddn)
        ;
}

void exportTemperature(py::module& m)
{
    py::class_<Temperature, ThermoScalar>(m, "Temperature")
        .def(py::init<>())
        .def(py::init<double>())
        ;

    py::implicitly_convertible<double, Temperature>();
}

void exportPressure(py::module& m)
{
    py::class_<Pressure, ThermoScalar>(m, "Pressure")
        .def(py::init<>())
        .def(py::init<double>())
        ;

    py::implicitly_convertible<double, Pressure>();
}

void exportAutoDiff(py::module& m)
{
    exportThermoScalar(m);
    exportThermoVector(m);
    exportChemicalScalar(m);
    exportChemicalVector(m);
    exportTemperature(m);
    exportPressure(m);
}

} // namespace Reaktoro

