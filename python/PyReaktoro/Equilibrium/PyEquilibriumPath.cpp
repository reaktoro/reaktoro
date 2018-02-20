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

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPath.hpp>

namespace Reaktoro {

void exportEquilibriumPath(py::module& m)
{
    py::class_<EquilibriumPathOptions>(m, "EquilibriumPathOptions")
        .def_readwrite("equilibrium", &EquilibriumPathOptions::equilibrium)
        .def_readwrite("ode", &EquilibriumPathOptions::ode)
        ;

    py::class_<EquilibriumPathResult>(m, "EquilibriumPathResult")
        .def_readwrite("equilibrium", &EquilibriumPathResult::equilibrium)
        ;

    py::class_<EquilibriumPath>(m, "EquilibriumPath")
        .def(py::init<const ChemicalSystem&>())
        .def("setOptions", &EquilibriumPath::setOptions)
        .def("setPartition", &EquilibriumPath::setPartition)
        .def("solve", &EquilibriumPath::solve)
        .def("output", &EquilibriumPath::output)
        .def("plot", &EquilibriumPath::plot)
        .def("plots", &EquilibriumPath::plots)
        .def("system", &EquilibriumPath::system, py::return_value_policy::reference_internal)
        .def("partition", &EquilibriumPath::partition, py::return_value_policy::reference_internal)
        ;
}

} // namespace Reaktoro
