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

#include "PyChemicalSolver.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>
#include <Reaktoro/Util/ChemicalField.hpp>
#include <Reaktoro/Util/ChemicalSolver.hpp>

namespace Reaktoro {

auto export_ChemicalSolver() -> void
{
    py::class_<ChemicalSolver>("ChemicalSolver")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&, Index>())
        .def(py::init<const ReactionSystem&, Index>())
        .def("setPartition", &ChemicalSolver::setPartition)
        .def("setState", &ChemicalSolver::setState)
        .def("setStates", &ChemicalSolver::setStates)
        .def("state", &ChemicalSolver::state, py::return_internal_reference<>())
        .def("states", &ChemicalSolver::states, py::return_internal_reference<>())
        .def("ne", &ChemicalSolver::ne, py::return_internal_reference<>())
        .def("porosity", &ChemicalSolver::porosity, py::return_internal_reference<>())
        .def("saturations", &ChemicalSolver::saturations, py::return_internal_reference<>())
        .def("densities", &ChemicalSolver::densities, py::return_internal_reference<>())
        .def("equilibrate", &ChemicalSolver::equilibrate)
        .def("react", &ChemicalSolver::react)
        ;
}

} // namespace Reaktoro
