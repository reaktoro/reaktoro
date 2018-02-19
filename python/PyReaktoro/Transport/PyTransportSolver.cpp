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

#include "PyTransportSolver.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Transport/TransportSolver.hpp>

namespace Reaktoro {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(setDiscretizationOverloads, setDiscretization, 1, 3);

auto export_Mesh() -> void
{
    py::class_<Mesh>("Mesh")
        .def(py::init<>())
        .def(py::init<Index, py::optional<double, double>>())
        .def("setDiscretization", &Mesh::setDiscretization, setDiscretizationOverloads())
        .def("numCells", &Mesh::numCells)
        .def("xl", &Mesh::xl)
        .def("xr", &Mesh::xr)
        .def("dx", &Mesh::dx)
//        .def("xcells", &Mesh::xcells, py::return_internal_reference<>())
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

auto export_ChemicalField() -> void
{
    py::class_<ChemicalField>("ChemicalField", py::no_init)
        .def(py::init<Index, const ChemicalSystem&>())
        .def(py::init<Index, const ChemicalState&>())
        .def("size", &ChemicalField::size)
        .def("set", &ChemicalField::set)
        .def("temperature", &ChemicalField::temperature)
        .def("pressure", &ChemicalField::pressure)
        .def("elementAmounts", &ChemicalField::elementAmounts)
        .def("output", &ChemicalField::output)
        .def("__setitem__", ChemicalField_setitem)
        .def("__getitem__", ChemicalField_getitem, py::return_internal_reference<>())
        ;
}

auto export_TransportSolver() -> void
{
    auto step1 = static_cast<void(TransportSolver::*)(VectorRef, VectorConstRef)>(&TransportSolver::step);
    auto step2 = static_cast<void(TransportSolver::*)(VectorRef)>(&TransportSolver::step);

    py::class_<TransportSolver>("TransportSolver")
        .def(py::init<>())
        .def("setMesh", &TransportSolver::setMesh)
        .def("setVelocity", &TransportSolver::setVelocity)
        .def("setDiffusionCoeff", &TransportSolver::setDiffusionCoeff)
        .def("setBoundaryValue", &TransportSolver::setBoundaryValue)
        .def("setTimeStep", &TransportSolver::setTimeStep)
        .def("mesh", &TransportSolver::mesh, py::return_internal_reference<>())
        .def("initialize", &TransportSolver::initialize)
        .def("step", step1)
        .def("step", step2)
        ;
}

auto export_ReactiveTransportSolver() -> void
{
    py::class_<ReactiveTransportSolver>("ReactiveTransportSolver", py::no_init)
        .def(py::init<const ChemicalSystem&>())
        .def("setMesh", &ReactiveTransportSolver::setMesh)
        .def("setVelocity", &ReactiveTransportSolver::setVelocity)
        .def("setDiffusionCoeff", &ReactiveTransportSolver::setDiffusionCoeff)
        .def("setBoundaryState", &ReactiveTransportSolver::setBoundaryState)
        .def("setTimeStep", &ReactiveTransportSolver::setTimeStep)
        .def("system", &ReactiveTransportSolver::system, py::return_internal_reference<>())
        .def("initialize", &ReactiveTransportSolver::initialize)
        .def("step", &ReactiveTransportSolver::step)
        ;
}

} // namespace Reaktoro
