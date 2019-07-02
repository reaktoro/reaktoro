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
#include <Reaktoro/Transport/Mesh.hpp>
#include <Reaktoro/Transport/TransportOptions.hpp>
#include <Reaktoro/Transport/TransportResult.hpp>
#include <Reaktoro/Transport/TransportSolver.hpp>

namespace Reaktoro {

void exportTransportSolver(py::module& m)
{
    auto step1 = static_cast<TransportResult(TransportSolver::*)(VectorRef, VectorConstRef)>(&TransportSolver::step);
    auto step2 = static_cast<TransportResult(TransportSolver::*)(VectorRef)>(&TransportSolver::step);

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
        .def("step", step1)
        .def("step", step2)
        ;
}

} // namespace Reaktoro
