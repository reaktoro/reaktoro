// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <Reaktoro/Transport/Mesh.hpp>
#include <Reaktoro/Transport/TransportOptions.hpp>
#include <Reaktoro/Transport/TransportResult.hpp>
#include <Reaktoro/Transport/TransportSolver.hpp>
using namespace Reaktoro;

void exportTransportSolver(py::module& m)
{
    py::class_<TransportSolver>(m, "TransportSolver")
        .def(py::init<>())
        .def("setMesh", &TransportSolver::setMesh)
        .def("setVelocity", &TransportSolver::setVelocity)
        .def("setDiffusionCoeff", &TransportSolver::setDiffusionCoeff)
        .def("setBoundaryValue", &TransportSolver::setBoundaryValue)
        .def("setTimeStep", &TransportSolver::setTimeStep)
        .def("mesh", &TransportSolver::mesh, py::return_value_policy::reference_internal)
        .def("timeStep", &TransportSolver::timeStep)
        .def("initialize", &TransportSolver::initialize)
        .def("step", py::overload_cast<ArrayXrRef, ArrayXrConstRef>(&TransportSolver::step))
        .def("step", py::overload_cast<ArrayXrRef>(&TransportSolver::step))
        ;
}
