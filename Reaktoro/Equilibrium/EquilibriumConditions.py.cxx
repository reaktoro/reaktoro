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
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
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
        .def("startWith", py::overload_cast<String, real, String>(&EquilibriumConditions::startWith))
        .def("startWith", py::overload_cast<Index, real, String>(&EquilibriumConditions::startWith))
        .def("startWith", py::overload_cast<const ChemicalState&>(&EquilibriumConditions::startWith))
        .def("startWithComponentAmounts", &EquilibriumConditions::startWithComponentAmounts)
        .def("chemicalPotential", &EquilibriumConditions::chemicalPotential)
        .def("lnActivity", &EquilibriumConditions::lnActivity)
        .def("lgActivity", &EquilibriumConditions::lgActivity)
        .def("activity", &EquilibriumConditions::activity)
        .def("fugacity", &EquilibriumConditions::fugacity)
        .def("pH", &EquilibriumConditions::pH)
        .def("pMg", &EquilibriumConditions::pMg)
        .def("pE", &EquilibriumConditions::pE)
        .def("Eh", &EquilibriumConditions::Eh)
        .def("set", &EquilibriumConditions::set)
        .def("initialSpeciesAmounts", &EquilibriumConditions::initialSpeciesAmounts, return_internal_ref)
        .def("initialComponentAmounts", &EquilibriumConditions::initialComponentAmounts, return_internal_ref)
        .def("system", &EquilibriumConditions::system, return_internal_ref)
        .def("inputNames", &EquilibriumConditions::inputNames, return_internal_ref)
        .def("inputValues", &EquilibriumConditions::inputValues, return_internal_ref)
        ;
}
