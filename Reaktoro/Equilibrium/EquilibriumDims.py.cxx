// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

void exportEquilibriumDims(py::module& m)
{
    py::class_<EquilibriumDims>(m, "EquilibriumDims")
        .def(py::init<>())
        .def(py::init<EquilibriumSpecs const&>())
        .def_readwrite("Ne", &EquilibriumDims::Ne, "The number of elements in the chemical system.")
        .def_readwrite("Nn", &EquilibriumDims::Nn, "The number of species in the chemical system.")
        .def_readwrite("Np", &EquilibriumDims::Np, "The number of *p* control variables (temperature, pressure, amounts of explicit titrants, and custom variables).")
        .def_readwrite("Nq", &EquilibriumDims::Nq, "The number of *q* control variables (amounts of implicit titrants).")
        .def_readwrite("Nv", &EquilibriumDims::Nv, "The number of equations constraints in the chemical equilibrium problem.")
        .def_readwrite("Nr", &EquilibriumDims::Nr, "The number of reactivity constraints (i.e., *restricted reactions*) in the chemical equilibrium problem.")
        .def_readwrite("Nb", &EquilibriumDims::Nb, "The number of elements and charge in the chemical system (equivalent to 1 + Ne).")
        .def_readwrite("Nc", &EquilibriumDims::Nc, "The number of components (electric charge, chemical elements, extent of restricted reactions) in the chemical equilibrium problem (equivalent to `1 + Ne + Nr`).")
        .def_readwrite("Nt", &EquilibriumDims::Nt, "The number of substances for which the chemical system is open to (the number of explicit and implicit titrants).")
        .def_readwrite("Nx", &EquilibriumDims::Nx, "The number of variables *x* in *x = (n, q)* (equivalent to `Nn + Nq`).")
        .def_readwrite("Nu", &EquilibriumDims::Nu, "The number of unknown variables in the chemical equilibrium problem (equivalent to `Nn + Np + Nq`).")
        .def_readwrite("Nw", &EquilibriumDims::Nw, "The number of input variables *w* in the chemical equilibrium problem.")
        ;
}
