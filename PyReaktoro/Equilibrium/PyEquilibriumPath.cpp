// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPath.hpp>

namespace Reaktoro {

auto export_EquilibriumPath() -> void
{
    py::class_<EquilibriumPathPlotOptions>("EquilibriumPathPlotOptions")
        .def_readwrite("active", &EquilibriumPathPlotOptions::active)
        .def_readwrite("period", &EquilibriumPathPlotOptions::period)
        .def_readwrite("execute", &EquilibriumPathPlotOptions::execute)
        ;

    py::class_<EquilibriumPathOptions>("EquilibriumPathOptions")
        .def_readwrite("equilibrium", &EquilibriumPathOptions::equilibrium)
        .def_readwrite("plot", &EquilibriumPathOptions::plot)
        ;

    auto setPartition1 = static_cast<void(EquilibriumPath::*)(const Partition&)>(&EquilibriumPath::setPartition);
    auto setPartition2 = static_cast<void(EquilibriumPath::*)(std::string)>(&EquilibriumPath::setPartition);

    py::class_<EquilibriumPath>("EquilibriumPath")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def("setOptions", &EquilibriumPath::setOptions)
        .def("setPartition", setPartition1)
        .def("setPartition", setPartition2)
        .def("solve", &EquilibriumPath::solve)
        .def("output", &EquilibriumPath::output)
        .def("plot", &EquilibriumPath::plot)
        ;
}

} // namespace Reaktoro
