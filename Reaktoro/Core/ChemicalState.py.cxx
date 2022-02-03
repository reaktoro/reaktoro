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

// Optima includes
#include <Optima/State.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
using namespace Reaktoro;

void exportChemicalState(py::module& m)
{
    const auto return_internal_ref = py::return_value_policy::reference_internal;

    py::class_<ChemicalState>(m, "ChemicalState")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalState&>())
        .def("temperature", py::overload_cast<real>(&ChemicalState::temperature))
        .def("temperature", py::overload_cast<real, String>(&ChemicalState::temperature))
        .def("pressure", py::overload_cast<real>(&ChemicalState::pressure))
        .def("pressure", py::overload_cast<real, String>(&ChemicalState::pressure))
        .def("charge", &ChemicalState::charge)
        .def("add", py::overload_cast<String, real, String>(&ChemicalState::add))
        .def("add", py::overload_cast<Index, real, String>(&ChemicalState::add))
        .def("set", py::overload_cast<String, real, String>(&ChemicalState::set))
        .def("set", py::overload_cast<Index, real, String>(&ChemicalState::set))
        .def("setTemperature", py::overload_cast<real>(&ChemicalState::setTemperature))
        .def("setTemperature", py::overload_cast<real, String>(&ChemicalState::setTemperature))
        .def("setPressure", py::overload_cast<real>(&ChemicalState::setPressure))
        .def("setPressure", py::overload_cast<real, String>(&ChemicalState::setPressure))
        .def("setSpeciesAmounts", [](ChemicalState& s, double val) { s.setSpeciesAmounts(val); })
        .def("setSpeciesAmounts", [](ChemicalState& s, const real& val) { s.setSpeciesAmounts(val); })
        .def("setSpeciesAmounts", [](ChemicalState& s, ArrayXrConstRef vals) { s.setSpeciesAmounts(vals); })
        .def("setSpeciesAmounts", [](ChemicalState& s, const py::array_t<double>& vals) { s.setSpeciesAmounts(ArrayXd::Map(vals.data(), vals.size())); })
        // .def("setSpeciesAmounts", [](ChemicalState& s, VectorXdConstRef vals) { s.setSpeciesAmounts(vals.array()); })
        // .def("setSpeciesAmounts", py::overload_cast<real>(&ChemicalState::setSpeciesAmounts))
        // .def("setSpeciesAmounts", py::overload_cast<ArrayXrConstRef>(&ChemicalState::setSpeciesAmounts))
        // .def("setSpeciesAmounts", py::overload_cast<ArrayXdConstRef>(&ChemicalState::setSpeciesAmounts))
        .def("setSpeciesAmount", &ChemicalState::setSpeciesAmount)
        .def("setSpeciesMass", &ChemicalState::setSpeciesMass)
        .def("scaleSpeciesAmounts", &ChemicalState::scaleSpeciesAmounts)
        .def("scaleSpeciesAmountsInPhase", &ChemicalState::scaleSpeciesAmountsInPhase)
        .def("scaleVolume", &ChemicalState::scaleVolume)
        .def("scalePhaseVolume", &ChemicalState::scalePhaseVolume)
        .def("scaleFluidVolume", &ChemicalState::scaleFluidVolume)
        .def("scaleSolidVolume", &ChemicalState::scaleSolidVolume)
        .def("scaleMass", &ChemicalState::scaleMass)
        .def("scalePhaseMass", &ChemicalState::scalePhaseMass)
        .def("scaleFluidMass", &ChemicalState::scaleFluidMass)
        .def("scaleSolidMass", &ChemicalState::scaleSolidMass)

        .def("setSurfaceArea", py::overload_cast<const StringOrIndex&, const StringOrIndex&, real, Chars>(&ChemicalState::setSurfaceArea))
        .def("setSurfaceArea", py::overload_cast<Index, real, Chars>(&ChemicalState::setSurfaceArea))
        .def("surfaceArea", py::overload_cast<const StringOrIndex&, const StringOrIndex&, real, Chars>(&ChemicalState::surfaceArea))
        .def("surfaceArea", py::overload_cast<const StringOrIndex&, const StringOrIndex&>(&ChemicalState::surfaceArea, py::const_))
        .def("surfaceArea", py::overload_cast<Index>(&ChemicalState::surfaceArea, py::const_))
        .def("surfaceAreas", &ChemicalState::surfaceAreas, return_internal_ref)
        .def("surfaces", &ChemicalState::surfaces, return_internal_ref)
        .def("surfaceIndex", &ChemicalState::surfaceIndex)

        .def("system", &ChemicalState::system, return_internal_ref)
        .def("temperature", py::overload_cast<>(&ChemicalState::temperature, py::const_))
        .def("pressure", py::overload_cast<>(&ChemicalState::pressure, py::const_))
        .def("speciesAmounts", &ChemicalState::speciesAmounts, return_internal_ref)
        .def("speciesAmountsInPhase", &ChemicalState::speciesAmountsInPhase, return_internal_ref)
        .def("componentAmounts", &ChemicalState::componentAmounts)
        .def("elementAmounts", &ChemicalState::elementAmounts)
        .def("speciesAmount", &ChemicalState::speciesAmount)
        .def("speciesMass", &ChemicalState::speciesMass)
        .def("equilibrium", py::overload_cast<>(&ChemicalState::equilibrium, py::const_), return_internal_ref)
        .def("equilibrium", py::overload_cast<>(&ChemicalState::equilibrium), return_internal_ref)
        .def("props", py::overload_cast<>(&ChemicalState::props, py::const_), return_internal_ref)
        .def("props", py::overload_cast<>(&ChemicalState::props), return_internal_ref)
        .def("output", py::overload_cast<std::ostream&>(&ChemicalState::output, py::const_))
        .def("output", py::overload_cast<const String&>(&ChemicalState::output, py::const_))
        .def("__repr__", [](const ChemicalState& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;

    py::class_<ChemicalState::Equilibrium>(m, "_ChemicalStateEquilibrium")
        .def("setInputNames", &ChemicalState::Equilibrium::setInputNames)
        .def("setInputValues", &ChemicalState::Equilibrium::setInputValues)
        .def("setInitialComponentAmounts", &ChemicalState::Equilibrium::setInitialComponentAmounts)
        .def("setControlVariablesP", &ChemicalState::Equilibrium::setControlVariablesP)
        .def("setControlVariablesQ", &ChemicalState::Equilibrium::setControlVariablesQ)
        .def("setOptimaState", &ChemicalState::Equilibrium::setOptimaState)
        .def("setIndicesPrimarySecondarySpecies", &ChemicalState::Equilibrium::setIndicesPrimarySecondarySpecies)
        .def("setIndicesStrictlyUnstableElements", &ChemicalState::Equilibrium::setIndicesStrictlyUnstableElements)
        .def("setIndicesStrictlyUnstableSpecies", &ChemicalState::Equilibrium::setIndicesStrictlyUnstableSpecies)
        .def("numPrimarySpecies", &ChemicalState::Equilibrium::numPrimarySpecies)
        .def("numSecondarySpecies", &ChemicalState::Equilibrium::numSecondarySpecies)
        .def("indicesPrimarySpecies", &ChemicalState::Equilibrium::indicesPrimarySpecies, return_internal_ref)
        .def("indicesSecondarySpecies", &ChemicalState::Equilibrium::indicesSecondarySpecies, return_internal_ref)
        .def("indicesStrictlyUnstableElements", &ChemicalState::Equilibrium::indicesStrictlyUnstableElements, return_internal_ref)
        .def("indicesStrictlyUnstableSpecies", &ChemicalState::Equilibrium::indicesStrictlyUnstableSpecies, return_internal_ref)
        .def("elementChemicalPotentials", &ChemicalState::Equilibrium::elementChemicalPotentials, return_internal_ref)
        .def("speciesStabilities", &ChemicalState::Equilibrium::speciesStabilities, return_internal_ref)
        .def("explicitTitrantAmounts", &ChemicalState::Equilibrium::explicitTitrantAmounts, return_internal_ref)
        .def("implicitTitrantAmounts", &ChemicalState::Equilibrium::implicitTitrantAmounts, return_internal_ref)
        .def("inputNames", &ChemicalState::Equilibrium::inputNames, return_internal_ref)
        .def("inputValues", &ChemicalState::Equilibrium::inputValues, return_internal_ref)
        .def("initialComponentAmounts", &ChemicalState::Equilibrium::initialComponentAmounts, return_internal_ref)
        .def("controlVariablesP", &ChemicalState::Equilibrium::controlVariablesP, return_internal_ref)
        .def("controlVariablesQ", &ChemicalState::Equilibrium::controlVariablesQ, return_internal_ref)
        .def("p", &ChemicalState::Equilibrium::p, return_internal_ref)
        .def("q", &ChemicalState::Equilibrium::q, return_internal_ref)
        .def("w", &ChemicalState::Equilibrium::w, return_internal_ref)
        .def("b", &ChemicalState::Equilibrium::b, return_internal_ref)
        .def("optimaState", &ChemicalState::Equilibrium::optimaState, return_internal_ref)
        ;
}
