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
#include <Reaktoro/Equilibrium/SmartEquilibriumOptions.hpp>
using namespace Reaktoro;

void exportSmartEquilibriumOptions(py::module& m)
{
    py::class_<SmartEquilibriumOptions>(m, "SmartEquilibriumOptions")
        .def(py::init<>())
        .def_readwrite("learning", &SmartEquilibriumOptions::learning, "The options for the chemical equilibrium calculations during learning operations.")
        .def_readwrite("reltol_negative_amounts", &SmartEquilibriumOptions::reltol_negative_amounts, "The relative tolerance for negative species amounts when predicting with first-order Taylor approximation.")
        .def_readwrite("reltol", &SmartEquilibriumOptions::reltol, "The relative tolerance used in the acceptance test for the predicted chemical equilibrium state.")
        .def_readwrite("abstol", &SmartEquilibriumOptions::abstol, "The absolute tolerance used in the acceptance test for the predicted chemical equilibrium state.")
        ;
}

