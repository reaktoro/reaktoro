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
#include <pybind11/functional.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
using namespace Reaktoro;

void exportEquilibriumConditions(py::module& m)
{
    const auto return_internal_ref = py::return_value_policy::reference_internal;

    py::class_<EquilibriumConstraintEquation>(m, "EquilibriumConstraintEquation")
        .def_readwrite("id", &EquilibriumConstraintEquation::id)
        .def_readwrite("fn", &EquilibriumConstraintEquation::fn)
        ;

    py::class_<EquilibriumConstraintConstantProperty>(m, "EquilibriumConstraintConstantProperty")
        .def_readwrite("id", &EquilibriumConstraintConstantProperty::id)
        .def_readwrite("fn", &EquilibriumConstraintConstantProperty::fn)
        ;

    py::class_<EquilibriumConstraintChemicalPotential>(m, "EquilibriumConstraintChemicalPotential")
        .def_readwrite("substance", &EquilibriumConstraintChemicalPotential::substance)
        .def_readwrite("fn", &EquilibriumConstraintChemicalPotential::fn)
        ;

    py::class_<EquilibriumConditionsDetails>(m, "EquilibriumConditionsDetails")
        .def_readwrite("unknownT", &EquilibriumConditionsDetails::unknownT)
        .def_readwrite("unknownP", &EquilibriumConditionsDetails::unknownP)
        .def_readwrite("constantT", &EquilibriumConditionsDetails::constantT)
        .def_readwrite("constantP", &EquilibriumConditionsDetails::constantP)
        .def_readwrite("T", &EquilibriumConditionsDetails::T)
        .def_readwrite("P", &EquilibriumConditionsDetails::P)
        .def_readwrite("b0", &EquilibriumConditionsDetails::b0)
        .def_readwrite("titrants", &EquilibriumConditionsDetails::titrants)
        .def_readwrite("econstraints", &EquilibriumConditionsDetails::econstraints)
        .def_readwrite("pconstraints", &EquilibriumConditionsDetails::pconstraints)
        .def_readwrite("uconstraints", &EquilibriumConditionsDetails::uconstraints)
        ;

    py::class_<EquilibriumConditions>(m, "EquilibriumConditions")
        .def("temperature", &EquilibriumConditions::temperature)
        .def("pressure", &EquilibriumConditions::pressure)
        .def("volume", &EquilibriumConditions::volume)
        .def("internalEnergy", &EquilibriumConditions::internalEnergy)
        .def("enthalpy", &EquilibriumConditions::enthalpy)
        .def("gibbsEnergy", &EquilibriumConditions::gibbsEnergy)
        .def("helmholtzEnergy", &EquilibriumConditions::helmholtzEnergy)
        .def("entropy", &EquilibriumConditions::entropy)
        .def("constantTemperature", &EquilibriumConditions::constantTemperature)
        .def("constantPressure", &EquilibriumConditions::constantPressure)
        .def("constantVolume", &EquilibriumConditions::constantVolume)
        .def("constantInternalEnergy", &EquilibriumConditions::constantInternalEnergy)
        .def("constantEnthalpy", &EquilibriumConditions::constantEnthalpy)
        .def("constantGibbsEnergy", &EquilibriumConditions::constantGibbsEnergy)
        .def("constantHelmholtzEnergy", &EquilibriumConditions::constantHelmholtzEnergy)
        .def("constantEntropy", &EquilibriumConditions::constantEntropy)
        .def("constantProperty", &EquilibriumConditions::constantProperty)
        .def("chemicalPotential", py::overload_cast<const ChemicalFormula&, const Fn<real(real,real)>&>(&EquilibriumConditions::chemicalPotential))
        .def("chemicalPotential", py::overload_cast<String, real, String>(&EquilibriumConditions::chemicalPotential))
        .def("lnActivity", py::overload_cast<const Species&, real>(&EquilibriumConditions::lnActivity))
        .def("lnActivity", py::overload_cast<String, real>(&EquilibriumConditions::lnActivity))
        .def("lgActivity", &EquilibriumConditions::lgActivity)
        .def("activity", &EquilibriumConditions::activity)
        .def("fugacity", &EquilibriumConditions::fugacity)
        .def("pH", &EquilibriumConditions::pH)
        .def("pMg", &EquilibriumConditions::pMg)
        .def("pE", &EquilibriumConditions::pE)
        .def("Eh", &EquilibriumConditions::Eh)
        .def("titrate", &EquilibriumConditions::titrate)
        .def("titrateEither", &EquilibriumConditions::titrateEither)
        .def("initialElementAmounts", &EquilibriumConditions::initialElementAmounts)
        .def("enforce", &EquilibriumConditions::enforce)
        .def("details", &EquilibriumConditions::details, return_internal_ref)
        ;
}
