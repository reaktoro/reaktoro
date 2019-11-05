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
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Kinetics/KineticOptions.hpp>
#include <Reaktoro/Kinetics/KineticPath.hpp>

namespace Reaktoro {

void exportKineticPath(py::module& m)
{
    py::class_<KineticPath>(m, "KineticPath")
        .def(py::init<const ReactionSystem&>())
        .def("setOptions", &KineticPath::setOptions)
        .def("setPartition", &KineticPath::setPartition)
        .def("addSource", &KineticPath::addSource)
        .def("addPhaseSink", &KineticPath::addPhaseSink)
        .def("addFluidSink", &KineticPath::addFluidSink)
        .def("addSolidSink", &KineticPath::addSolidSink)
        .def("solve", &KineticPath::solve)
        .def("output", &KineticPath::output)
        .def("plot", &KineticPath::plot)
        .def("plots", &KineticPath::plots)
        .def("system", &KineticPath::system, py::return_value_policy::reference_internal)
        .def("reactions", &KineticPath::reactions, py::return_value_policy::reference_internal)
        .def("partition", &KineticPath::partition, py::return_value_policy::reference_internal);
}

} // namespace Reaktoro
