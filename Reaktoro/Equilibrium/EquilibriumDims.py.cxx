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

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

void exportEquilibriumDims(py::module& m)
{
    py::class_<EquilibriumDims>(m, "EquilibriumDims")
        .def(py::init<>())
        .def(py::init<const EquilibriumSpecs&>())
        .def_readwrite("Ne", &EquilibriumDims::Ne)
        .def_readwrite("Nc", &EquilibriumDims::Nc)
        .def_readwrite("Nn", &EquilibriumDims::Nn)
        .def_readwrite("Np", &EquilibriumDims::Np)
        .def_readwrite("Nq", &EquilibriumDims::Nq)
        .def_readwrite("Nt", &EquilibriumDims::Nt)
        .def_readwrite("Nx", &EquilibriumDims::Nx)
        .def_readwrite("Nu", &EquilibriumDims::Nu)
        ;
}
