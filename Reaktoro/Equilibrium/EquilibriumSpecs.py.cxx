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
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

void exportEquilibriumSpecs(py::module& m)
{
    const auto return_internal_ref = py::return_value_policy::reference_internal;

    py::class_<EquilibriumConstraintEquation>(m, "EquilibriumConstraintEquation")
        .def_readwrite("name", &EquilibriumConstraintEquation::name)
        .def_readwrite("fn", &EquilibriumConstraintEquation::fn)
        ;

    py::class_<EquilibriumConstraintChemicalPotential>(m, "EquilibriumConstraintChemicalPotential")
        .def_readwrite("name", &EquilibriumConstraintChemicalPotential::name)
        .def_readwrite("substance", &EquilibriumConstraintChemicalPotential::substance)
        .def_readwrite("fn", &EquilibriumConstraintChemicalPotential::fn)
        ;

    py::class_<EquilibriumSpecs>(m, "EquilibriumSpecs")
        .def(py::init<const ChemicalSystem&>())
        .def("temperature", &EquilibriumSpecs::temperature)
        .def("pressure", &EquilibriumSpecs::pressure)
        .def("volume", &EquilibriumSpecs::volume)
        .def("internalEnergy", &EquilibriumSpecs::internalEnergy)
        .def("enthalpy", &EquilibriumSpecs::enthalpy)
        .def("gibbsEnergy", &EquilibriumSpecs::gibbsEnergy)
        .def("helmholtzEnergy", &EquilibriumSpecs::helmholtzEnergy)
        .def("entropy", &EquilibriumSpecs::entropy)
        .def("chemicalPotential", &EquilibriumSpecs::chemicalPotential)
        .def("lnActivity", py::overload_cast<const Species&>(&EquilibriumSpecs::lnActivity))
        .def("lnActivity", py::overload_cast<String>(&EquilibriumSpecs::lnActivity))
        .def("lgActivity", &EquilibriumSpecs::lgActivity)
        .def("activity", &EquilibriumSpecs::activity)
        .def("fugacity", &EquilibriumSpecs::fugacity)
        .def("pH", &EquilibriumSpecs::pH)
        .def("pMg", &EquilibriumSpecs::pMg)
        .def("pE", &EquilibriumSpecs::pE)
        .def("Eh", &EquilibriumSpecs::Eh)
        .def("openTo", &EquilibriumSpecs::openTo)
        .def("numInputs", &EquilibriumSpecs::numInputs)
        .def("numParams", &EquilibriumSpecs::numParams)
        .def("numControlVariables", &EquilibriumSpecs::numControlVariables)
        .def("numControlVariablesP", &EquilibriumSpecs::numControlVariablesP)
        .def("numControlVariablesQ", &EquilibriumSpecs::numControlVariablesQ)
        .def("numTitrants", &EquilibriumSpecs::numTitrants)
        .def("numTitrantsExplicit", &EquilibriumSpecs::numTitrantsExplicit)
        .def("numTitrantsImplicit", &EquilibriumSpecs::numTitrantsImplicit)
        .def("numConstraints", &EquilibriumSpecs::numConstraints)
        .def("numConstraintsEquationType", &EquilibriumSpecs::numConstraintsEquationType)
        .def("numConstraintsChemicalPotentialType", &EquilibriumSpecs::numConstraintsChemicalPotentialType)
        .def("namesInputs", &EquilibriumSpecs::namesInputs)
        .def("namesParams", &EquilibriumSpecs::namesParams)
        .def("namesControlVariables", &EquilibriumSpecs::namesControlVariables)
        .def("namesControlVariablesP", &EquilibriumSpecs::namesControlVariablesP)
        .def("namesControlVariablesQ", &EquilibriumSpecs::namesControlVariablesQ)
        .def("namesTitrants", &EquilibriumSpecs::namesTitrants)
        .def("namesTitrantsExplicit", &EquilibriumSpecs::namesTitrantsExplicit)
        .def("namesTitrantsImplicit", &EquilibriumSpecs::namesTitrantsImplicit)
        .def("namesConstraints", &EquilibriumSpecs::namesConstraints)
        .def("namesConstraintsEquationType", &EquilibriumSpecs::namesConstraintsEquationType)
        .def("namesConstraintsChemicalPotentialType", &EquilibriumSpecs::namesConstraintsChemicalPotentialType)
        .def("addConstraint", py::overload_cast<const EquilibriumConstraintEquation&>(&EquilibriumSpecs::addConstraint))
        .def("addConstraint", py::overload_cast<const EquilibriumConstraintChemicalPotential&>(&EquilibriumSpecs::addConstraint))
        .def("addInput", py::overload_cast<const String&>(&EquilibriumSpecs::addInput), return_internal_ref)
        .def("addInput", py::overload_cast<const Param&>(&EquilibriumSpecs::addInput), return_internal_ref)
        .def("system", &EquilibriumSpecs::system, return_internal_ref)
        .def("inputs", &EquilibriumSpecs::inputs, return_internal_ref)
        .def("params", &EquilibriumSpecs::params, return_internal_ref)
        .def("indicesParams", &EquilibriumSpecs::indicesParams, return_internal_ref)
        .def("isTemperatureUnknown", &EquilibriumSpecs::isTemperatureUnknown)
        .def("isPressureUnknown", &EquilibriumSpecs::isPressureUnknown)
        .def("titrants", &EquilibriumSpecs::titrants)
        .def("titrantsExplicit", &EquilibriumSpecs::titrantsExplicit)
        .def("titrantsImplicit", &EquilibriumSpecs::titrantsImplicit)
        .def("constraintsEquationType", &EquilibriumSpecs::constraintsEquationType, return_internal_ref)
        .def("constraintsChemicalPotentialType", &EquilibriumSpecs::constraintsChemicalPotentialType, return_internal_ref)
        ;
}
