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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

void exportEquilibriumConditions(py::module& m)
{
    py::class_<EquilibriumConditions>(m, "EquilibriumConditions")
        .def(py::init<ChemicalSystem const&>())
        .def(py::init<EquilibriumSpecs const&>())

        .def("temperature", &EquilibriumConditions::temperature, "Specify the temperature of the system at chemical equilibrium", py::arg("value"), py::arg("unit")="K")
        .def("pressure", &EquilibriumConditions::pressure, "Specify the pressure of the system at chemical equilibrium", py::arg("value"), py::arg("unit")="Pa")
        .def("volume", &EquilibriumConditions::volume, "Specify the volume of the system at chemical equilibrium", py::arg("value"), py::arg("unit")="m3")
        .def("internalEnergy", &EquilibriumConditions::internalEnergy, "Specify the internal energy of the system at chemical equilibrium", py::arg("value"), py::arg("unit")="J")
        .def("enthalpy", &EquilibriumConditions::enthalpy, "Specify the enthalpy of the system at chemical equilibrium", py::arg("value"), py::arg("unit")="J")
        .def("gibbsEnergy", &EquilibriumConditions::gibbsEnergy, "Specify the Gibbs energy of the system at chemical equilibrium", py::arg("value"), py::arg("unit")="J")
        .def("helmholtzEnergy", &EquilibriumConditions::helmholtzEnergy, "Specify the Helmholtz energy of the system at chemical equilibrium", py::arg("value"), py::arg("unit")="J")
        .def("entropy", &EquilibriumConditions::entropy, "Specify the entropy of the system at chemical equilibrium", py::arg("value"), py::arg("unit")="J/K")
        .def("charge", &EquilibriumConditions::charge, "Specify the electric charge at chemical equilibrium", py::arg("value"), py::arg("unit")="mol")
        .def("elementAmount", &EquilibriumConditions::elementAmount, "Specify the amount of an element at chemical equilibrium", py::arg("element"), py::arg("value"), py::arg("unit")="mol")
        .def("elementAmountInPhase", &EquilibriumConditions::elementAmountInPhase, "Specify the amount of an element in a phase at chemical equilibrium", py::arg("element"), py::arg("phase"), py::arg("value"), py::arg("unit")="mol")
        .def("elementMass", &EquilibriumConditions::elementMass, "Specify the mass of an element at chemical equilibrium", py::arg("element"), py::arg("value"), py::arg("unit")="kg")
        .def("elementMassInPhase", &EquilibriumConditions::elementMassInPhase, "Specify the mass of an element in a phase at chemical equilibrium", py::arg("element"), py::arg("phase"), py::arg("value"), py::arg("unit")="kg")
        .def("phaseAmount", &EquilibriumConditions::phaseAmount, "Specify the amount of a phase at chemical equilibrium", py::arg("phase"), py::arg("value"), py::arg("unit")="mol")
        .def("phaseMass", &EquilibriumConditions::phaseMass, "Specify the mass of a phase at chemical equilibrium", py::arg("phase"), py::arg("value"), py::arg("unit")="kg")
        .def("phaseVolume", &EquilibriumConditions::phaseVolume, "Specify the volume of a phase at chemical equilibrium", py::arg("phase"), py::arg("value"), py::arg("unit")="m3")
        .def("chemicalPotential", &EquilibriumConditions::chemicalPotential, "Specify the chemical potential of a substance at chemical equilibrium", py::arg("species"), py::arg("value"), py::arg("unit")="J/mol")
        .def("lnActivity", &EquilibriumConditions::lnActivity, "Specify the ln activity of a species at chemical equilibrium", py::arg("species"), py::arg("value"))
        .def("lgActivity", &EquilibriumConditions::lgActivity, "Specify the lg activity of a species at chemical equilibrium", py::arg("species"), py::arg("value"))
        .def("activity", &EquilibriumConditions::activity, "Specify the activity of a species at chemical equilibrium", py::arg("species"), py::arg("value"))
        .def("fugacity", &EquilibriumConditions::fugacity, "Specify the fugacity of a gaseous species at chemical equilibrium", py::arg("species"), py::arg("value"), py::arg("unit")="Pa")
        .def("pH", &EquilibriumConditions::pH, "Specify the pH at chemical equilibrium", py::arg("value"))
        .def("pMg", &EquilibriumConditions::pMg, "Specify the pMg at chemical equilibrium", py::arg("value"))
        .def("pE", &EquilibriumConditions::pE, "Specify the pE at chemical equilibrium", py::arg("value"))
        .def("Eh", &EquilibriumConditions::Eh, "Specify the Eh at chemical equilibrium", py::arg("value"), py::arg("unit")="V")

        .def("setLowerBoundTemperature", &EquilibriumConditions::setLowerBoundTemperature, "Set the lower bound for temperature during the equilibrium calculation", py::arg("value"), py::arg("unit")="K")
        .def("setUpperBoundTemperature", &EquilibriumConditions::setUpperBoundTemperature, "Set the upper bound for temperature during the equilibrium calculation", py::arg("value"), py::arg("unit")="K")
        .def("setLowerBoundPressure", &EquilibriumConditions::setLowerBoundPressure, "Set the lower bound for pressure during the equilibrium calculation", py::arg("value"), py::arg("unit")="Pa")
        .def("setUpperBoundPressure", &EquilibriumConditions::setUpperBoundPressure, "Set the upper bound for pressure during the equilibrium calculation", py::arg("value"), py::arg("unit")="Pa")
        .def("setLowerBoundTitrant", &EquilibriumConditions::setLowerBoundTitrant, "Set the lower bound for the amount of a titrant during the equilibrium calculation", py::arg("substance"), py::arg("value"), py::arg("unit")="mol")
        .def("setUpperBoundTitrant", &EquilibriumConditions::setUpperBoundTitrant, "Set the upper bound for the amount of a titrant during the equilibrium calculation", py::arg("substance"), py::arg("value"), py::arg("unit")="mol")
        .def("setLowerBoundsControlVariablesP", &EquilibriumConditions::setLowerBoundsControlVariablesP, "Set the values of the specified lower bounds for the *p* control variables.")
        .def("setUpperBoundsControlVariablesP", &EquilibriumConditions::setUpperBoundsControlVariablesP, "Set the values of the specified upper bounds for the *p* control variables.")
        .def("lowerBoundsControlVariablesP", &EquilibriumConditions::lowerBoundsControlVariablesP, return_internal_ref, "Return the specified lower bounds for the *p* control variables.")
        .def("upperBoundsControlVariablesP", &EquilibriumConditions::upperBoundsControlVariablesP, return_internal_ref, "Return the specified upper bounds for the *p* control variables.")

        .def("set", &EquilibriumConditions::set, "Set the value of an input variable with given name.")
        .def("setInputVariable", py::overload_cast<String const&, real const&>(&EquilibriumConditions::setInputVariable), "Set the value of an input variable with given name.")
        .def("setInputVariable", py::overload_cast<Index, real const&>(&EquilibriumConditions::setInputVariable), "Set the value of an input variable with given index.")
        .def("setInputVariables", &EquilibriumConditions::setInputVariables, "Set the input variables with given vector of input values.")
        .def("inputNames", &EquilibriumConditions::inputNames, return_internal_ref, "Return the names of the input variables associated with the equilibrium conditions.")
        .def("inputValues", &EquilibriumConditions::inputValues, return_internal_ref, "Return the values of the input variables associated with the equilibrium conditions.")
        .def("inputValuesGetOrCompute", &EquilibriumConditions::inputValuesGetOrCompute, "Get the values of the input variables associated with the equilibrium conditions if specified, otherwise fetch them from given initial state.")
        .def("inputValue", &EquilibriumConditions::inputValue, return_internal_ref, "Return the values of the input variables associated with the equilibrium conditions.")

        .def("setInitialComponentAmounts", &EquilibriumConditions::setInitialComponentAmounts, "Set the initial amounts of the conservative components c0 before the chemical system reacts.")
        .def("setInitialComponentAmountsFromSpeciesAmounts", &EquilibriumConditions::setInitialComponentAmountsFromSpeciesAmounts, "Set the initial amounts of the conservative components c0 before the chemical system reacts.")
        .def("setInitialComponentAmountsFromState", &EquilibriumConditions::setInitialComponentAmountsFromState, "Set the initial amounts of the conservative components c0 before the chemical system reacts.")
        .def("initialComponentAmounts", &EquilibriumConditions::initialComponentAmounts, return_internal_ref, "Get the initial amounts of the conservative components c0 before the chemical system reacts.")
        .def("initialComponentAmountsGetOrCompute", py::overload_cast<VectorXdConstRef const&>(&EquilibriumConditions::initialComponentAmountsGetOrCompute, py::const_), "Get the initial amounts of the conservative components c0 before the chemical system reacts if available, otherwise compute it.")
        .def("initialComponentAmountsGetOrCompute", py::overload_cast<ChemicalState const&>(&EquilibriumConditions::initialComponentAmountsGetOrCompute, py::const_), "Get the initial amounts of the conservative components c0 before the chemical system reacts if available, otherwise compute it.")

        .def("system", &EquilibriumConditions::system, return_internal_ref, "Return the chemical system associated with the equilibrium conditions.")
        ;
}
