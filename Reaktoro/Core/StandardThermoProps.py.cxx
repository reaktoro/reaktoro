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
#include <Reaktoro/Core/Model.py.hxx>
#include <Reaktoro/Core/StandardThermoProps.hpp>
using namespace Reaktoro;

void exportStandardThermoProps(py::module& m)
{
    py::class_<StandardThermoProps>(m, "StandardThermoProps")
        .def(py::init<>())
        .def_readwrite("G0", &StandardThermoProps::G0)
        .def_readwrite("H0", &StandardThermoProps::H0)
        .def_readwrite("V0", &StandardThermoProps::V0)
        .def_readwrite("Cp0", &StandardThermoProps::Cp0)
        .def_readwrite("Cv0", &StandardThermoProps::Cv0)
        ;

    exportModel<StandardThermoProps, real, real>(m, "StandardThermoModel");
}
