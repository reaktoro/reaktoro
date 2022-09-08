// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/ReactionStandardThermoProps.hpp>
using namespace Reaktoro;

void exportReactionStandardThermoProps(py::module& m)
{
    py::class_<ReactionStandardThermoProps>(m, "ReactionStandardThermoProps")
        .def(py::init<>())
        .def_readwrite("dG0" , &ReactionStandardThermoProps::dG0 , "The standard molar Gibbs energy change of the reaction (in J/mol)")
        .def_readwrite("dH0" , &ReactionStandardThermoProps::dH0 , "The standard molar enthalpy change of the reaction (in J/mol)")
        .def_readwrite("dCp0", &ReactionStandardThermoProps::dCp0, "The standard molar isobaric heat capacity change of the reaction (in J/(mol·K))")
        ;

    py::class_<ReactionStandardThermoModelArgs>(m, "ReactionStandardThermoModelArgs")
        .def_property_readonly("T", [](const ReactionStandardThermoModelArgs& self) { return self.T; })
        .def_property_readonly("P", [](const ReactionStandardThermoModelArgs& self) { return self.P; })
        .def_property_readonly("dV0", [](const ReactionStandardThermoModelArgs& self) { return self.dV0; })
        ;

    exportModel<ReactionStandardThermoProps, ReactionStandardThermoModelArgs>(m, "ReactionStandardThermoModel");
}
