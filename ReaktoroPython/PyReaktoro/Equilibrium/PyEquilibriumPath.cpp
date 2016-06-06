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

#include "PyEquilibriumPath.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalOutput.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPath.hpp>
#include <Reaktoro/Equilibrium/EquilibriumState.hpp>

namespace Reaktoro {

auto export_EquilibriumPath() -> void
{
    py::class_<EquilibriumPathOptions>("EquilibriumPathOptions")
        .def_readwrite("equilibrium", &EquilibriumPathOptions::equilibrium)
        .def_readwrite("ode", &EquilibriumPathOptions::ode)
        ;

    py::class_<EquilibriumPathResult>("EquilibriumPathResult")
        .def_readwrite("equilibrium", &EquilibriumPathResult::equilibrium)
        ;

    py::class_<EquilibriumPath>("EquilibriumPath", py::no_init)
        .def(py::init<const ChemicalSystem&>())
        .def("setOptions", &EquilibriumPath::setOptions)
        .def("setPartition", &EquilibriumPath::setPartition)
        .def("solve", &EquilibriumPath::solve)
        .def("output", &EquilibriumPath::output)
        .def("plot", &EquilibriumPath::plot)
        .def("plots", &EquilibriumPath::plots)
        .def("system", &EquilibriumPath::system, py::return_internal_reference<>())
        .def("partition", &EquilibriumPath::partition, py::return_internal_reference<>())
        ;
}

} // namespace Reaktoro
