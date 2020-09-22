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
#include <pybind11/eigen.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
using namespace Reaktoro;

void exportChemicalState(py::module& m)
{
    py::class_<ChemicalState>(m, "ChemicalState")
        .def(py::init<const ChemicalSystem&>())
        .def("setTemperature", py::overload_cast<real>(&ChemicalState::setTemperature))
        .def("setTemperature", py::overload_cast<real, String>(&ChemicalState::setTemperature))
        .def("setPressure", py::overload_cast<real>(&ChemicalState::setPressure))
        .def("setPressure", py::overload_cast<real, String>(&ChemicalState::setPressure))
        .def("setSpeciesAmounts", py::overload_cast<real>(&ChemicalState::setSpeciesAmounts))
        .def("setSpeciesAmounts", py::overload_cast<ArrayXrConstRef>(&ChemicalState::setSpeciesAmounts))
        .def("setSpeciesAmounts", py::overload_cast<ArrayXdConstRef>(&ChemicalState::setSpeciesAmounts))
        .def("setSpeciesAmount", py::overload_cast<Index, real>(&ChemicalState::setSpeciesAmount))
        .def("setSpeciesAmount", py::overload_cast<Index, real, String>(&ChemicalState::setSpeciesAmount))
        .def("setSpeciesAmount", py::overload_cast<String, real>(&ChemicalState::setSpeciesAmount))
        .def("setSpeciesAmount", py::overload_cast<String, real, String>(&ChemicalState::setSpeciesAmount))
        .def("setSpeciesMass", py::overload_cast<Index, real>(&ChemicalState::setSpeciesMass))
        .def("setSpeciesMass", py::overload_cast<Index, real, String>(&ChemicalState::setSpeciesMass))
        .def("setSpeciesMass", py::overload_cast<String, real>(&ChemicalState::setSpeciesMass))
        .def("setSpeciesMass", py::overload_cast<String, real, String>(&ChemicalState::setSpeciesMass))
        .def("system", &ChemicalState::system, py::return_value_policy::reference_internal)
        .def("temperature", &ChemicalState::temperature)
        .def("pressure", &ChemicalState::pressure)
        .def("speciesAmounts", py::overload_cast<>(&ChemicalState::speciesAmounts, py::const_), py::return_value_policy::reference_internal)
        .def("speciesAmounts", py::overload_cast<>(&ChemicalState::speciesAmounts), py::return_value_policy::reference_internal)
        .def("elementAmounts", &ChemicalState::elementAmounts)
        .def("speciesAmount", py::overload_cast<Index>(&ChemicalState::speciesAmount, py::const_))
        .def("speciesAmount", py::overload_cast<Index, String>(&ChemicalState::speciesAmount, py::const_))
        .def("speciesAmount", py::overload_cast<String>(&ChemicalState::speciesAmount, py::const_))
        .def("speciesAmount", py::overload_cast<String, String>(&ChemicalState::speciesAmount, py::const_))
        .def("speciesMass", py::overload_cast<Index>(&ChemicalState::speciesMass, py::const_))
        .def("speciesMass", py::overload_cast<Index, String>(&ChemicalState::speciesMass, py::const_))
        .def("speciesMass", py::overload_cast<String>(&ChemicalState::speciesMass, py::const_))
        .def("speciesMass", py::overload_cast<String, String>(&ChemicalState::speciesMass, py::const_))
        .def("phaseProps", &ChemicalState::phaseProps, py::return_value_policy::reference_internal)
        .def("props", py::overload_cast<>(&ChemicalState::props, py::const_), py::return_value_policy::reference_internal)
        .def("props", py::overload_cast<>(&ChemicalState::props), py::return_value_policy::reference_internal)
        .def("equilibrium", py::overload_cast<>(&ChemicalState::equilibrium, py::const_), py::return_value_policy::reference_internal)
        .def("equilibrium", py::overload_cast<>(&ChemicalState::equilibrium), py::return_value_policy::reference_internal)
        ;

    py::class_<ChemicalState::Equilibrium>(m, "_ChemicalStateEquilibrium")
        .def(py::init<const ChemicalSystem&>())
        .def("setIndicesPrimarySecondarySpecies", &ChemicalState::Equilibrium::setIndicesPrimarySecondarySpecies)
        .def("setIndicesStrictlyUnstableElements", &ChemicalState::Equilibrium::setIndicesStrictlyUnstableElements)
        .def("setIndicesStrictlyUnstableSpecies", &ChemicalState::Equilibrium::setIndicesStrictlyUnstableSpecies)
        .def("setLagrangeMultipliers", &ChemicalState::Equilibrium::setLagrangeMultipliers)
        .def("setComplementarityVariables", &ChemicalState::Equilibrium::setComplementarityVariables)
        .def("setControlVariables", &ChemicalState::Equilibrium::setControlVariables)
        .def("numPrimarySpecies", &ChemicalState::Equilibrium::numPrimarySpecies)
        .def("numSecondarySpecies", &ChemicalState::Equilibrium::numSecondarySpecies)
        .def("indicesPrimarySpecies", &ChemicalState::Equilibrium::indicesPrimarySpecies, py::return_value_policy::reference_internal)
        .def("indicesSecondarySpecies", &ChemicalState::Equilibrium::indicesSecondarySpecies, py::return_value_policy::reference_internal)
        .def("indicesStrictlyUnstableElements", &ChemicalState::Equilibrium::indicesStrictlyUnstableElements, py::return_value_policy::reference_internal)
        .def("indicesStrictlyUnstableSpecies", &ChemicalState::Equilibrium::indicesStrictlyUnstableSpecies, py::return_value_policy::reference_internal)
        .def("lagrangeMultipliers", py::overload_cast<>(&ChemicalState::Equilibrium::lagrangeMultipliers, py::const_), py::return_value_policy::reference_internal)
        .def("lagrangeMultipliers", py::overload_cast<>(&ChemicalState::Equilibrium::lagrangeMultipliers), py::return_value_policy::reference_internal)
        .def("complementarityVariables", py::overload_cast<>(&ChemicalState::Equilibrium::complementarityVariables, py::const_), py::return_value_policy::reference_internal)
        .def("complementarityVariables", py::overload_cast<>(&ChemicalState::Equilibrium::complementarityVariables), py::return_value_policy::reference_internal)
        .def("controlVariables", py::overload_cast<>(&ChemicalState::Equilibrium::controlVariables, py::const_), py::return_value_policy::reference_internal)
        .def("controlVariables", py::overload_cast<>(&ChemicalState::Equilibrium::controlVariables), py::return_value_policy::reference_internal)
        ;
}
