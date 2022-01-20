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

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/SpeciesThermoProps.hpp>
#include <Reaktoro/Core/StandardThermoProps.hpp>
using namespace Reaktoro;

void exportSpeciesThermoProps(py::module& m)
{
    py::class_<SpeciesThermoProps>(m, "SpeciesThermoProps")
        .def(py::init<>())
        .def(py::init<const real&, const real&, const StandardThermoProps&>())
        .def_readwrite("T"  , &SpeciesThermoProps::T)
        .def_readwrite("P"  , &SpeciesThermoProps::P)
        .def_readwrite("G0" , &SpeciesThermoProps::G0)
        .def_readwrite("H0" , &SpeciesThermoProps::H0)
        .def_readwrite("V0" , &SpeciesThermoProps::V0)
        .def_readwrite("VT0", &SpeciesThermoProps::VT0)
        .def_readwrite("VP0", &SpeciesThermoProps::VP0)
        .def_readwrite("Cp0", &SpeciesThermoProps::Cp0)
        .def_readwrite("Cv0", &SpeciesThermoProps::Cv0)
        .def_readwrite("U0" , &SpeciesThermoProps::U0)
        .def_readwrite("S0" , &SpeciesThermoProps::S0)
        .def_readwrite("A0" , &SpeciesThermoProps::A0)
        .def("__repr__", [](const SpeciesThermoProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(double() * py::self)
        .def(py::self * double())
        .def(py::self / double())
        ;
}
