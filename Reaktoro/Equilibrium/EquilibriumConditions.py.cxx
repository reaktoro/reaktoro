// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
        .def(py::init<const EquilibriumSpecs>())
        .def("temperature", &EquilibriumConditions::temperature)
        .def("pressure", &EquilibriumConditions::pressure)
        .def("volume", &EquilibriumConditions::volume)
        .def("internalEnergy", &EquilibriumConditions::internalEnergy)
        .def("enthalpy", &EquilibriumConditions::enthalpy)
        .def("gibbsEnergy", &EquilibriumConditions::gibbsEnergy)
        .def("helmholtzEnergy", &EquilibriumConditions::helmholtzEnergy)
        .def("entropy", &EquilibriumConditions::entropy)
        .def("charge", &EquilibriumConditions::charge)
        .def("elementAmount", &EquilibriumConditions::elementAmount)
        .def("elementAmountInPhase", &EquilibriumConditions::elementAmountInPhase)
        .def("elementMass", &EquilibriumConditions::elementMass)
        .def("elementMassInPhase", &EquilibriumConditions::elementMassInPhase)
        .def("phaseAmount", &EquilibriumConditions::phaseAmount)
        .def("phaseMass", &EquilibriumConditions::phaseMass)
        .def("phaseVolume", &EquilibriumConditions::phaseVolume)
        .def("chemicalPotential", &EquilibriumConditions::chemicalPotential)
        .def("lnActivity", &EquilibriumConditions::lnActivity)
        .def("lgActivity", &EquilibriumConditions::lgActivity)
        .def("activity", &EquilibriumConditions::activity)
        .def("fugacity", &EquilibriumConditions::fugacity)
        .def("pH", &EquilibriumConditions::pH)
        .def("pMg", &EquilibriumConditions::pMg)
        .def("pE", &EquilibriumConditions::pE)
        .def("Eh", &EquilibriumConditions::Eh)
        .def("setLowerBoundTemperature", &EquilibriumConditions::setLowerBoundTemperature)
        .def("setUpperBoundTemperature", &EquilibriumConditions::setUpperBoundTemperature)
        .def("setLowerBoundPressure", &EquilibriumConditions::setLowerBoundPressure)
        .def("setUpperBoundPressure", &EquilibriumConditions::setUpperBoundPressure)
        .def("setLowerBoundTitrant", &EquilibriumConditions::setLowerBoundTitrant)
        .def("setUpperBoundTitrant", &EquilibriumConditions::setUpperBoundTitrant)
        .def("set", &EquilibriumConditions::set)
        .def("system", &EquilibriumConditions::system, return_internal_ref)
        .def("inputNames", &EquilibriumConditions::inputNames, return_internal_ref)
        .def("inputValues", &EquilibriumConditions::inputValues, return_internal_ref)
        .def("lowerBoundsControlVariablesP", &EquilibriumConditions::lowerBoundsControlVariablesP, return_internal_ref)
        .def("upperBoundsControlVariablesP", &EquilibriumConditions::upperBoundsControlVariablesP, return_internal_ref)
        ;
}
