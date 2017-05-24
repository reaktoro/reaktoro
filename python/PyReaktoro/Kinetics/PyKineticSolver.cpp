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

#include "PyKineticSolver.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Kinetics/KineticOptions.hpp>
#include <Reaktoro/Kinetics/KineticSolver.hpp>

namespace Reaktoro {

auto export_KineticSolver() -> void
{
    auto step1 = static_cast<void(KineticSolver::*)(ChemicalState&, double&)>(&KineticSolver::step);
    auto step2 = static_cast<void(KineticSolver::*)(ChemicalState&, double&, double)>(&KineticSolver::step);

    py::class_<KineticSolver>("KineticSolver", py::no_init)
        .def(py::init<const ReactionSystem&>())
        .def("setOptions", &KineticSolver::setOptions)
        .def("setPartition", &KineticSolver::setPartition)
        .def("addSource", &KineticSolver::addSource)
        .def("addPhaseSink", &KineticSolver::addPhaseSink)
        .def("addFluidSink", &KineticSolver::addFluidSink)
        .def("addSolidSink", &KineticSolver::addSolidSink)
        .def("initialize", &KineticSolver::initialize)
        .def("step", step1)
        .def("step", step2)
        .def("solve", &KineticSolver::solve)
        ;
}

} // namespace Reaktoro
