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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Kinetics/KineticOptions.hpp>
#include <Reaktoro/Kinetics/KineticSolver.hpp>
#include <Reaktoro/Kinetics/KineticResult.hpp>

namespace Reaktoro {

void exportKineticSolver(py::module& m)
{
    auto step1 = static_cast<double(KineticSolver::*)(ChemicalState&, double)>(&KineticSolver::step);
    auto step2 = static_cast<double(KineticSolver::*)(ChemicalState&, double, double)>(&KineticSolver::step);

    auto initialize1 = static_cast<void(KineticSolver::*)(ChemicalState&, double)>(&KineticSolver::initialize);
    auto initialize2 = static_cast<void(KineticSolver::*)(ChemicalState&, double, VectorConstRef)>(&KineticSolver::initialize);

    auto solve1 = static_cast<double(KineticSolver::*)(ChemicalState&, double, double)>(&KineticSolver::solve);
    auto solve2 = static_cast<void(KineticSolver::*)(ChemicalState&, double, double, VectorConstRef)>(&KineticSolver::solve);

    py::class_<KineticSolver>(m, "KineticSolver")
        .def(py::init<const ReactionSystem&, const Partition&>())
        .def("setOptions", &KineticSolver::setOptions)
        .def("addSource", &KineticSolver::addSource)
        .def("addPhaseSink", &KineticSolver::addPhaseSink)
        .def("addFluidSink", &KineticSolver::addFluidSink)
        .def("addSolidSink", &KineticSolver::addSolidSink)
        .def("initialize", initialize1)
        .def("initialize", initialize2)
        .def("step", step1)
        .def("step", step2)
        .def("solve", solve1)
        .def("solve", solve2)
        .def("result", &KineticSolver::result, py::return_value_policy::reference_internal)

        // DEPRECATED METHODS: TO BE REMOVED
        .def("setPartition", &KineticSolver::setPartition)
        ;
}

} // namespace Reaktoro
