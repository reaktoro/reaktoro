// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Kinetics/KineticsSensitivity.hpp>
using namespace Reaktoro;

void exportKineticsSensitivity(py::module& m)
{
    py::class_<KineticsSensitivity, EquilibriumSensitivity>(m, "KineticsSensitivity")
        .def(py::init<>())
        .def(py::init<EquilibriumSpecs const&>())
        .def(py::init<EquilibriumSensitivity const&>())
        .def("initialize", &KineticsSensitivity::initialize, "Initialize this KineticsSensitivity object with given equilibrium problem specifications.")
        .def("dnddt", &KineticsSensitivity::dnddt, return_internal_ref, "Return the derivatives of the species amounts n with respect to time step dt.")
        .def("dpddt", &KineticsSensitivity::dpddt, return_internal_ref, "Return the derivatives of the p control variables with respect to time step dt.")
        .def("dqddt", &KineticsSensitivity::dqddt, return_internal_ref, "Return the derivatives of the q control variables with respect to time step dt.")
        .def("duddt", &KineticsSensitivity::duddt, return_internal_ref, "Return the sensitivity derivatives of the chemical properties u with respect to time step dt.")
        ;
}
