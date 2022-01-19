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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

void exportEquilibriumSpecs(py::module& m)
{
    py::class_<ControlVariableQ>(m, "ControlVariableQ")
        .def(py::init<>())
        .def_readwrite("name", &ControlVariableQ::name)
        .def_readwrite("substance", &ControlVariableQ::substance)
        .def_readwrite("id", &ControlVariableQ::id)
        .def_readwrite("fn", &ControlVariableQ::fn)
        ;

    py::class_<ControlVariableP>(m, "ControlVariableP")
        .def(py::init<>())
        .def_readwrite("name", &ControlVariableP::name)
        .def_readwrite("substance", &ControlVariableP::substance)
        .def_readwrite("ispecies", &ControlVariableP::ispecies)
        .def_readwrite("fn", &ControlVariableP::fn)
        ;

    py::class_<ConstraintEquation>(m, "ConstraintEquation")
        .def(py::init<>())
        .def_readwrite("id", &ConstraintEquation::id)
        .def_readwrite("fn", &ConstraintEquation::fn)
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
        .def("charge", &EquilibriumSpecs::charge)
        .def("elementAmount", &EquilibriumSpecs::elementAmount)
        .def("elementAmountInPhase", &EquilibriumSpecs::elementAmountInPhase)
        .def("elementMass", &EquilibriumSpecs::elementMass)
        .def("elementMassInPhase", &EquilibriumSpecs::elementMassInPhase)
        .def("phaseAmount", &EquilibriumSpecs::phaseAmount)
        .def("phaseMass", &EquilibriumSpecs::phaseMass)
        .def("phaseVolume", &EquilibriumSpecs::phaseVolume)
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
        .def("addUnknownTitrantAmount", &EquilibriumSpecs::addUnknownTitrantAmount)
        .def("addUnknownChemicalPotential", &EquilibriumSpecs::addUnknownChemicalPotential)
        .def("addUnknownStandardChemicalPotential", &EquilibriumSpecs::addUnknownStandardChemicalPotential)
        .def("addUnknownActivity", &EquilibriumSpecs::addUnknownActivity)
        .def("addUnknownActivityCoefficient", &EquilibriumSpecs::addUnknownActivityCoefficient)
        .def("numInputs", &EquilibriumSpecs::numInputs)
        .def("numParams", &EquilibriumSpecs::numParams)
        .def("numControlVariables", &EquilibriumSpecs::numControlVariables)
        .def("numControlVariablesP", &EquilibriumSpecs::numControlVariablesP)
        .def("numControlVariablesQ", &EquilibriumSpecs::numControlVariablesQ)
        .def("numTitrants", &EquilibriumSpecs::numTitrants)
        .def("numTitrantsExplicit", &EquilibriumSpecs::numTitrantsExplicit)
        .def("numTitrantsImplicit", &EquilibriumSpecs::numTitrantsImplicit)
        .def("numConstraints", &EquilibriumSpecs::numConstraints)
        .def("namesInputs", &EquilibriumSpecs::namesInputs)
        .def("namesParams", &EquilibriumSpecs::namesParams)
        .def("namesControlVariables", &EquilibriumSpecs::namesControlVariables)
        .def("namesControlVariablesP", &EquilibriumSpecs::namesControlVariablesP)
        .def("namesControlVariablesQ", &EquilibriumSpecs::namesControlVariablesQ)
        .def("namesTitrants", &EquilibriumSpecs::namesTitrants)
        .def("namesTitrantsExplicit", &EquilibriumSpecs::namesTitrantsExplicit)
        .def("namesTitrantsImplicit", &EquilibriumSpecs::namesTitrantsImplicit)
        .def("namesConstraints", &EquilibriumSpecs::namesConstraints)
        .def("addControlVariableQ", &EquilibriumSpecs::addControlVariableQ)
        .def("addControlVariableP", &EquilibriumSpecs::addControlVariableP)
        .def("addConstraint", &EquilibriumSpecs::addConstraint)
        .def("addInput", py::overload_cast<const String&>(&EquilibriumSpecs::addInput), return_internal_ref)
        .def("addInput", py::overload_cast<const Param&>(&EquilibriumSpecs::addInput), return_internal_ref)
        .def("system", &EquilibriumSpecs::system, return_internal_ref)
        .def("inputs", &EquilibriumSpecs::inputs, return_internal_ref)
        .def("params", &EquilibriumSpecs::params, return_internal_ref)
        .def("indicesParams", &EquilibriumSpecs::indicesParams, return_internal_ref)
        .def("isTemperatureUnknown", &EquilibriumSpecs::isTemperatureUnknown)
        .def("isPressureUnknown", &EquilibriumSpecs::isPressureUnknown)
        .def("indexControlVariableTemperature", &EquilibriumSpecs::indexControlVariableTemperature)
        .def("indexControlVariablePressure", &EquilibriumSpecs::indexControlVariablePressure)
        .def("controlVariablesQ", &EquilibriumSpecs::controlVariablesQ, return_internal_ref)
        .def("controlVariablesP", &EquilibriumSpecs::controlVariablesP, return_internal_ref)
        .def("titrants", &EquilibriumSpecs::titrants)
        .def("titrantsExplicit", &EquilibriumSpecs::titrantsExplicit)
        .def("titrantsImplicit", &EquilibriumSpecs::titrantsImplicit)
        .def("constraintsEquationType", &EquilibriumSpecs::constraintsEquationType, return_internal_ref)
        ;
}
