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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
using namespace Reaktoro;

void exportChemicalState(py::module& m)
{
    py::class_<ChemicalState>(m, "ChemicalState")
        .def(py::init<ChemicalSystem const&>())
        .def(py::init<ChemicalState const&>())

        .def("setTemperature", py::overload_cast<real const&>(&ChemicalState::setTemperature))
        .def("setTemperature", py::overload_cast<real, Chars>(&ChemicalState::setTemperature))
        .def("temperature", py::overload_cast<real const&>(&ChemicalState::temperature))
        .def("temperature", py::overload_cast<real, Chars>(&ChemicalState::temperature))
        .def("temperature", py::overload_cast<>(&ChemicalState::temperature, py::const_))

        .def("setPressure", py::overload_cast<real const&>(&ChemicalState::setPressure))
        .def("setPressure", py::overload_cast<real, Chars>(&ChemicalState::setPressure))
        .def("pressure", py::overload_cast<real const&>(&ChemicalState::pressure))
        .def("pressure", py::overload_cast<real, Chars>(&ChemicalState::pressure))
        .def("pressure", py::overload_cast<>(&ChemicalState::pressure, py::const_))

        .def("setSpeciesAmounts", [](ChemicalState& s, double val) { s.setSpeciesAmounts(val); })
        .def("setSpeciesAmounts", [](ChemicalState& s, real val) { s.setSpeciesAmounts(val); })
        .def("setSpeciesAmounts", [](ChemicalState& s, ArrayXrConstRef const& vals) { s.setSpeciesAmounts(vals); })
        .def("setSpeciesAmounts", [](ChemicalState& s, py::array_t<double> const& vals) { s.setSpeciesAmounts(ArrayXd::Map(vals.data(), vals.size())); })
        .def("setSpeciesAmount", py::overload_cast<Index, real const&>(&ChemicalState::setSpeciesAmount))
        .def("setSpeciesAmount", py::overload_cast<StringOrIndex const&, real, Chars>(&ChemicalState::setSpeciesAmount))
        .def("setSpeciesMass", &ChemicalState::setSpeciesMass)
        .def("set", &ChemicalState::set)
        .def("add", &ChemicalState::add)

        .def("speciesAmounts", &ChemicalState::speciesAmounts, return_internal_ref)
        .def("speciesAmountsInPhase", &ChemicalState::speciesAmountsInPhase, return_internal_ref)
        .def("speciesAmount", &ChemicalState::speciesAmount)
        .def("speciesMass", &ChemicalState::speciesMass)
        .def("componentAmounts", &ChemicalState::componentAmounts)
        .def("elementAmounts", &ChemicalState::elementAmounts)
        .def("charge", &ChemicalState::charge)

        .def("scaleSpeciesAmounts", py::overload_cast<real const&>(&ChemicalState::scaleSpeciesAmounts))
        .def("scaleSpeciesAmounts", py::overload_cast<real const&, Indices const&>(&ChemicalState::scaleSpeciesAmounts))
        .def("scaleSpeciesAmountsInPhase", &ChemicalState::scaleSpeciesAmountsInPhase)

        .def("scaleVolume", &ChemicalState::scaleVolume)
        .def("scalePhaseVolume", &ChemicalState::scalePhaseVolume)
        .def("scaleFluidVolume", &ChemicalState::scaleFluidVolume)
        .def("scaleSolidVolume", &ChemicalState::scaleSolidVolume)

        .def("scaleMass", &ChemicalState::scaleMass)
        .def("scalePhaseMass", &ChemicalState::scalePhaseMass)
        .def("scaleFluidMass", &ChemicalState::scaleFluidMass)
        .def("scaleSolidMass", &ChemicalState::scaleSolidMass)

        .def("setSurfaceAreas", [](ChemicalState& s, ArrayXrConstRef const& vals) { s.setSurfaceAreas(vals); })
        .def("setSurfaceAreas", [](ChemicalState& s, py::array_t<double> const& vals) { s.setSurfaceAreas(ArrayXd::Map(vals.data(), vals.size())); })
        .def("setSurfaceArea", py::overload_cast<StringOrIndex const&, StringOrIndex const&, real, Chars>(&ChemicalState::setSurfaceArea))
        .def("setSurfaceArea", py::overload_cast<StringOrIndex const&, real, Chars>(&ChemicalState::setSurfaceArea))
        .def("surfaceArea", py::overload_cast<StringOrIndex const&, StringOrIndex const&, real, Chars>(&ChemicalState::surfaceArea))
        .def("surfaceArea", py::overload_cast<StringOrIndex const&, real, Chars>(&ChemicalState::surfaceArea))
        .def("surfaceArea", py::overload_cast<StringOrIndex const&, StringOrIndex const&>(&ChemicalState::surfaceArea, py::const_))
        .def("surfaceArea", py::overload_cast<StringOrIndex const&>(&ChemicalState::surfaceArea, py::const_))
        .def("surfaceAreas", &ChemicalState::surfaceAreas, return_internal_ref)

        .def("update", py::overload_cast<real const&, real const&, ArrayXrConstRef const&>(&ChemicalState::update))
        .def("update", py::overload_cast<real const&, real const&, ArrayXrConstRef const&, ArrayXrConstRef const&>(&ChemicalState::update))

        .def("updateIdeal", py::overload_cast<real const&, real const&, ArrayXrConstRef const&>(&ChemicalState::updateIdeal))
        .def("updateIdeal", py::overload_cast<real const&, real const&, ArrayXrConstRef const&, ArrayXrConstRef const&>(&ChemicalState::updateIdeal))

        .def("system", &ChemicalState::system, return_internal_ref)
        .def("props", py::overload_cast<>(&ChemicalState::props, py::const_), return_internal_ref)
        .def("props", py::overload_cast<>(&ChemicalState::props), return_internal_ref)
        .def("equilibrium", py::overload_cast<>(&ChemicalState::equilibrium, py::const_), return_internal_ref)
        .def("equilibrium", py::overload_cast<>(&ChemicalState::equilibrium), return_internal_ref)
        .def("output", py::overload_cast<std::ostream&>(&ChemicalState::output, py::const_))
        .def("output", py::overload_cast<String const&>(&ChemicalState::output, py::const_))
        .def("__repr__", [](ChemicalState const& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;

    py::class_<ChemicalState::Equilibrium>(m, "_ChemicalStateEquilibrium")
        .def("reset", &ChemicalState::Equilibrium::reset)
        .def("setNamesInputVariables", &ChemicalState::Equilibrium::setNamesInputVariables)
        .def("setNamesControlVariablesP", &ChemicalState::Equilibrium::setNamesControlVariablesP)
        .def("setNamesControlVariablesQ", &ChemicalState::Equilibrium::setNamesControlVariablesQ)
        .def("setInputVariables", &ChemicalState::Equilibrium::setInputVariables)
        .def("setControlVariablesP", &ChemicalState::Equilibrium::setControlVariablesP)
        .def("setControlVariablesQ", &ChemicalState::Equilibrium::setControlVariablesQ)
        .def("setInitialComponentAmounts", &ChemicalState::Equilibrium::setInitialComponentAmounts)
        .def("setOptimaState", &ChemicalState::Equilibrium::setOptimaState)
        .def("numPrimarySpecies", &ChemicalState::Equilibrium::numPrimarySpecies)
        .def("numSecondarySpecies", &ChemicalState::Equilibrium::numSecondarySpecies)
        .def("indicesPrimarySpecies", &ChemicalState::Equilibrium::indicesPrimarySpecies, return_internal_ref)
        .def("indicesSecondarySpecies", &ChemicalState::Equilibrium::indicesSecondarySpecies, return_internal_ref)
        .def("elementChemicalPotentials", &ChemicalState::Equilibrium::elementChemicalPotentials, return_internal_ref)
        .def("speciesStabilities", &ChemicalState::Equilibrium::speciesStabilities, return_internal_ref)
        .def("explicitTitrantAmount", &ChemicalState::Equilibrium::explicitTitrantAmount)
        .def("implicitTitrantAmount", &ChemicalState::Equilibrium::implicitTitrantAmount)
        .def("titrantAmount", &ChemicalState::Equilibrium::titrantAmount)
        .def("namesInputVariables", &ChemicalState::Equilibrium::namesInputVariables, return_internal_ref)
        .def("namesControlVariablesP", &ChemicalState::Equilibrium::namesControlVariablesP, return_internal_ref)
        .def("namesControlVariablesQ", &ChemicalState::Equilibrium::namesControlVariablesQ, return_internal_ref)
        .def("inputVariables", &ChemicalState::Equilibrium::inputVariables, return_internal_ref)
        .def("controlVariablesP", &ChemicalState::Equilibrium::controlVariablesP, return_internal_ref)
        .def("controlVariablesQ", &ChemicalState::Equilibrium::controlVariablesQ, return_internal_ref)
        .def("initialComponentAmounts", &ChemicalState::Equilibrium::initialComponentAmounts, return_internal_ref)
        .def("w", &ChemicalState::Equilibrium::w, return_internal_ref)
        .def("p", &ChemicalState::Equilibrium::p, return_internal_ref)
        .def("q", &ChemicalState::Equilibrium::q, return_internal_ref)
        .def("c", &ChemicalState::Equilibrium::c, return_internal_ref)
        .def("optimaState", &ChemicalState::Equilibrium::optimaState, return_internal_ref)
        ;
}
