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
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Param.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <pybind11/detail/common.h>
using namespace Reaktoro;

void exportEquilibriumSensitivity(py::module& m)
{
    py::class_<EquilibriumSensitivity>(m, "EquilibriumSensitivity")
        .def(py::init<>())
        .def(py::init<EquilibriumSpecs const&>())
        .def("initialize", &EquilibriumSensitivity::initialize, "Initialize this EquilibriumSensitivity object with given equilibrium problem specifications.")
        .def("dndw", py::overload_cast<String const&>(&EquilibriumSensitivity::dndw, py::const_), return_internal_ref, "Return the derivatives of the species amounts n with respect to an input variable in w.")
        .def("dndw", py::overload_cast<Param const&>(&EquilibriumSensitivity::dndw, py::const_), return_internal_ref, "Return the derivatives of the species amounts n with respect to an input variable in w.")
        .def("dndw", py::overload_cast<>(&EquilibriumSensitivity::dndw, py::const_), return_internal_ref, "Return the derivatives of the species amounts n with respect to the input variables w.")
        .def("dndw", py::overload_cast<MatrixXdConstRef>(&EquilibriumSensitivity::dndw), "Set the derivatives of the species amounts n with respect to the input variables w.")
        .def("dpdw", py::overload_cast<String const&>(&EquilibriumSensitivity::dpdw, py::const_), return_internal_ref, "Return the derivatives of the p control variables with respect to an input variable in w.")
        .def("dpdw", py::overload_cast<Param const&>(&EquilibriumSensitivity::dpdw, py::const_), return_internal_ref, "Return the derivatives of the p control variables with respect to an input variable in w.")
        .def("dpdw", py::overload_cast<>(&EquilibriumSensitivity::dpdw, py::const_), return_internal_ref, "Return the derivatives of the p control variables with respect to the input variables w.")
        .def("dpdw", py::overload_cast<MatrixXdConstRef>(&EquilibriumSensitivity::dpdw), "Set the derivatives of the p control variables with respect to the input variables w.")
        .def("dqdw", py::overload_cast<String const&>(&EquilibriumSensitivity::dqdw, py::const_), return_internal_ref, "Return the derivatives of the q control variables with respect to an input variable in w.")
        .def("dqdw", py::overload_cast<Param const&>(&EquilibriumSensitivity::dqdw, py::const_), return_internal_ref, "Return the derivatives of the q control variables with respect to an input variable in w.")
        .def("dqdw", py::overload_cast<>(&EquilibriumSensitivity::dqdw, py::const_), return_internal_ref, "Return the derivatives of the q control variables with respect to the input variables w.")
        .def("dqdw", py::overload_cast<MatrixXdConstRef>(&EquilibriumSensitivity::dqdw), "Set the derivatives of the q control variables with respect to the input variables w.")
        .def("dndb", py::overload_cast<>(&EquilibriumSensitivity::dndb, py::const_), return_internal_ref, "Return the derivatives of the species amounts n with respect to component amounts b.")
        .def("dndb", py::overload_cast<MatrixXdConstRef>(&EquilibriumSensitivity::dndb), "Set the derivatives of the species amounts n with respect to component amounts b.")
        .def("dpdb", py::overload_cast<>(&EquilibriumSensitivity::dpdb, py::const_), return_internal_ref, "Return the derivatives of the control variables p with respect to component amounts b.")
        .def("dpdb", py::overload_cast<MatrixXdConstRef>(&EquilibriumSensitivity::dpdb), "Set the derivatives of the control variables p with respect to component amounts b.")
        .def("dqdb", py::overload_cast<>(&EquilibriumSensitivity::dqdb, py::const_), return_internal_ref, "Return the derivatives of the control variables q with respect to component amounts b.")
        .def("dqdb", py::overload_cast<MatrixXdConstRef>(&EquilibriumSensitivity::dqdb), "Set the derivatives of the control variables q with respect to component amounts b.")
        .def("dudw", py::overload_cast<>(&EquilibriumSensitivity::dudw, py::const_), return_internal_ref, "Return the total derivatives of the chemical properties u with respect to input variables w.")
        .def("dudb", py::overload_cast<>(&EquilibriumSensitivity::dudb, py::const_), return_internal_ref, "Return the total derivatives of the chemical properties u with respect to component amounts b.")
        .def("dudw", py::overload_cast<MatrixXdConstRef>(&EquilibriumSensitivity::dudw), "Set the total derivatives of the chemical properties u with respect to input variables w.")
        .def("dudb", py::overload_cast<MatrixXdConstRef>(&EquilibriumSensitivity::dudb), "Set the total derivatives of the chemical properties u with respect to component amounts b.")
        ;
}
