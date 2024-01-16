// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Models/StandardThermoModels/StandardThermoModelExtendedUNIQUAC.hpp>
using namespace Reaktoro;

void exportStandardThermoModelExtendedUNIQUAC(py::module& m)
{
    py::class_<StandardThermoModelParamsExtendedUNIQUAC>(m, "StandardThermoModelParamsExtendedUNIQUAC")
        .def(py::init<>())
        .def_readwrite("Gr"   , &StandardThermoModelParamsExtendedUNIQUAC::Gr   , "The standard molar Gibbs energy of formation of the substance at reference temperature 298.15 K (in kJ/mol).")
        .def_readwrite("Hr"   , &StandardThermoModelParamsExtendedUNIQUAC::Hr   , "The standard molar enthalpy of formation of the substance at reference temperature 298.15 K (in kJ/mol).")
        .def_readwrite("Sr"   , &StandardThermoModelParamsExtendedUNIQUAC::Sr   , "The standard molar entropy of the substance at reference temperature 298.15 K (in J/(mol·K)).")
        .def_readwrite("Vr"   , &StandardThermoModelParamsExtendedUNIQUAC::Vr   , "The standard molar volume of the substance at reference temperature 298.15 K (in m³/mol).")
        .def_readwrite("Cp"   , &StandardThermoModelParamsExtendedUNIQUAC::Cp   , "The standard molar heat capacity (constant pressure) of the substance at reference temperature 298.15 K (in J/(mol·K)).")
        .def_readwrite("a"    , &StandardThermoModelParamsExtendedUNIQUAC::a    , "The parameter `a` in the standard heat capacity model (in J/(mol·K)).")
        .def_readwrite("b"    , &StandardThermoModelParamsExtendedUNIQUAC::b    , "The parameter `b` in the standard heat capacity model (in J/(mol·K²)).")
        .def_readwrite("c"    , &StandardThermoModelParamsExtendedUNIQUAC::c    , "The parameter `c` in the standard heat capacity model (in J/mol).")
        .def_readwrite("alpha", &StandardThermoModelParamsExtendedUNIQUAC::alpha, "The parameter `α` in the standard heat capacity model (in 1/bar).")
        .def_readwrite("beta" , &StandardThermoModelParamsExtendedUNIQUAC::beta , "The parameter `β` in the standard heat capacity model (in 1/bar²).")
        .def_readwrite("Theta", &StandardThermoModelParamsExtendedUNIQUAC::Theta, "The parameter `Θ` in the standard heat capacity model (in K).")
        ;

    m.def("StandardThermoModelExtendedUNIQUAC", StandardThermoModelExtendedUNIQUAC);
}
