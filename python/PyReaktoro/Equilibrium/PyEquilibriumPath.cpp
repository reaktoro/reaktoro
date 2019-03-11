// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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
