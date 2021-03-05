// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <pybind11/functional.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/Params.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

void exportEquilibriumConditions(py::module& m)
{
    const auto return_internal_ref = py::return_value_policy::reference_internal;

    py::class_<EquilibriumConditions>(m, "EquilibriumConditions")
        .def(py::init<const EquilibriumSpecs>())
        .def("temperature", &EquilibriumConditions::temperature)
        .def("pressure", &EquilibriumConditions::pressure)
        .def("volume", &EquilibriumConditions::volume)
        .def("internalEnergy", &EquilibriumConditions::internalEnergy)
        .def("enthalpy", &EquilibriumConditions::enthalpy)
        .def("gibbsEnergy", &EquilibriumConditions::gibbsEnergy)
        .def("helmholtzEnergy", &EquilibriumConditions::helmholtzEnergy)
        .def("entropy", &EquilibriumConditions::entropy)
        .def("chemicalPotential", &EquilibriumConditions::chemicalPotential)
        .def("lnActivity", &EquilibriumConditions::lnActivity)
        .def("lgActivity", &EquilibriumConditions::lgActivity)
        .def("activity", &EquilibriumConditions::activity)
        .def("fugacity", &EquilibriumConditions::fugacity)
        .def("pH", &EquilibriumConditions::pH)
        .def("pMg", &EquilibriumConditions::pMg)
        .def("pE", &EquilibriumConditions::pE)
        .def("Eh", &EquilibriumConditions::Eh)
        .def("initialComponentAmounts", py::overload_cast<VectorXrConstRef>(&EquilibriumConditions::initialComponentAmounts))
        .def("initialComponentAmounts", py::overload_cast<>(&EquilibriumConditions::initialComponentAmounts, py::const_))
        .def("initialComponentAmountsCompute", &EquilibriumConditions::initialComponentAmountsCompute)
        .def("initialComponentAmountsComputeOrRetrieve", &EquilibriumConditions::initialComponentAmountsComputeOrRetrieve)
        .def("system", &EquilibriumConditions::system)
        .def("params", &EquilibriumConditions::params)
        ;
}
