// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Transport/TransportSolver.hpp>

namespace Reaktoro {

void exportMesh(py::module& m)
{
    py::class_<Mesh>(m, "Mesh")
        .def(py::init<>())
        .def(py::init<Index, double, double>(), py::arg("num_cells"), py::arg("xl") = 0.0, py::arg("xr") = 1.0)
        .def("setDiscretization", &Mesh::setDiscretization, py::arg("num_cells"), py::arg("xl") = 0.0, py::arg("xr") = 1.0)
        .def("numCells", &Mesh::numCells)
        .def("xl", &Mesh::xl)
        .def("xr", &Mesh::xr)
        .def("dx", &Mesh::dx)
        .def("xcells", &Mesh::xcells, py::return_value_policy::reference_internal)
        ;
}

auto ChemicalField_setitem(ChemicalField& self, Index i, const ChemicalState& state) -> void
{
    self[i] = state;
}

auto ChemicalField_getitem(const ChemicalField& self, Index i) -> const ChemicalState&
{
    return self[i];
}

void exportChemicalField(py::module& m)
{
    py::class_<ChemicalField>(m, "ChemicalField")
        .def(py::init<Index, const ChemicalSystem&>())
        .def(py::init<Index, const ChemicalState&>())
        .def("size", &ChemicalField::size)
        .def("set", &ChemicalField::set)
        .def("temperature", &ChemicalField::temperature)
        .def("pressure", &ChemicalField::pressure)
        .def("elementAmounts", &ChemicalField::elementAmounts)
        .def("output", &ChemicalField::output)
        .def("__setitem__", ChemicalField_setitem)
        .def("__getitem__", ChemicalField_getitem, py::return_value_policy::reference_internal)
        ;
}

void exportTransportSolver(py::module& m)
{
    auto step1 = static_cast<void(TransportSolver::*)(VectorRef, VectorConstRef)>(&TransportSolver::step);
    auto step2 = static_cast<void(TransportSolver::*)(VectorRef)>(&TransportSolver::step);

    py::class_<TransportSolver>(m, "TransportSolver")
        .def(py::init<>())
        .def("setMesh", &TransportSolver::setMesh)
        .def("setVelocity", &TransportSolver::setVelocity)
        .def("setDiffusionCoeff", &TransportSolver::setDiffusionCoeff)
        .def("setBoundaryValue", &TransportSolver::setBoundaryValue)
        .def("setTimeStep", &TransportSolver::setTimeStep)
        .def("mesh", &TransportSolver::mesh, py::return_value_policy::reference_internal)
        .def("initialize", &TransportSolver::initialize)
        .def("step", step1)
        .def("step", step2)
        ;
}

void exportReactiveTransportSolver(py::module& m)
{
    py::class_<ReactiveTransportSolver>(m, "ReactiveTransportSolver")
        .def(py::init<const ChemicalSystem&>())
        .def("setMesh", &ReactiveTransportSolver::setMesh)
        .def("setVelocity", &ReactiveTransportSolver::setVelocity)
        .def("setDiffusionCoeff", &ReactiveTransportSolver::setDiffusionCoeff)
        .def("setBoundaryState", &ReactiveTransportSolver::setBoundaryState)
        .def("setTimeStep", &ReactiveTransportSolver::setTimeStep)
        .def("system", &ReactiveTransportSolver::system, py::return_value_policy::reference_internal)
        .def("output", &ReactiveTransportSolver::output)
        .def("initialize", &ReactiveTransportSolver::initialize)
        .def("step", &ReactiveTransportSolver::step)
        ;
}

} // namespace Reaktoro
