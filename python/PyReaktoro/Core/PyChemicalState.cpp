// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>
#include <Reaktoro/Interfaces/Gems.hpp>
#include <Reaktoro/Interfaces/Phreeqc.hpp>

namespace Reaktoro {
namespace {

auto assignChemicalState(ChemicalState& state, const ChemicalState& other) -> void
{
    state = other;
}

auto cloneChemicalState(ChemicalState& state) -> ChemicalState
{
    return state;
}

}  // namespace

void exportChemicalState(py::module& m)
{
    auto setTemperature1 = static_cast<void(ChemicalState::*)(double)>(&ChemicalState::setTemperature);
    auto setTemperature2 = static_cast<void(ChemicalState::*)(double, std::string)>(&ChemicalState::setTemperature);

    auto setPressure1 = static_cast<void(ChemicalState::*)(double)>(&ChemicalState::setPressure);
    auto setPressure2 = static_cast<void(ChemicalState::*)(double, std::string)>(&ChemicalState::setPressure);

    auto setSpeciesAmounts1 = static_cast<void(ChemicalState::*)(double)>(&ChemicalState::setSpeciesAmounts);
    auto setSpeciesAmounts2 = static_cast<void(ChemicalState::*)(VectorConstRef)>(&ChemicalState::setSpeciesAmounts);
    auto setSpeciesAmounts3 = static_cast<void(ChemicalState::*)(VectorConstRef, const Indices&)>(&ChemicalState::setSpeciesAmounts);

    auto setSpeciesAmount1 = static_cast<void(ChemicalState::*)(Index, double)>(&ChemicalState::setSpeciesAmount);
    auto setSpeciesAmount2 = static_cast<void(ChemicalState::*)(std::string, double)>(&ChemicalState::setSpeciesAmount);
    auto setSpeciesAmount3 = static_cast<void(ChemicalState::*)(Index, double, std::string)>(&ChemicalState::setSpeciesAmount);
    auto setSpeciesAmount4 = static_cast<void(ChemicalState::*)(std::string, double, std::string)>(&ChemicalState::setSpeciesAmount);

    auto setSpeciesMass1 = static_cast<void(ChemicalState::*)(Index, double)>(&ChemicalState::setSpeciesMass);
    auto setSpeciesMass2 = static_cast<void(ChemicalState::*)(std::string, double)>(&ChemicalState::setSpeciesMass);
    auto setSpeciesMass3 = static_cast<void(ChemicalState::*)(Index, double, std::string)>(&ChemicalState::setSpeciesMass);
    auto setSpeciesMass4 = static_cast<void(ChemicalState::*)(std::string, double, std::string)>(&ChemicalState::setSpeciesMass);

    auto scalePhaseVolume1 = static_cast<void(ChemicalState::*)(Index, double)>(&ChemicalState::scalePhaseVolume);
    auto scalePhaseVolume2 = static_cast<void(ChemicalState::*)(Index, double, std::string)>(&ChemicalState::scalePhaseVolume);
    auto scalePhaseVolume3 = static_cast<void(ChemicalState::*)(std::string, double)>(&ChemicalState::scalePhaseVolume);
    auto scalePhaseVolume4 = static_cast<void(ChemicalState::*)(std::string, double, std::string)>(&ChemicalState::scalePhaseVolume);

    auto scaleFluidVolume1 = static_cast<void(ChemicalState::*)(double)>(&ChemicalState::scaleFluidVolume);
    auto scaleFluidVolume2 = static_cast<void(ChemicalState::*)(double, std::string)>(&ChemicalState::scaleFluidVolume);

    auto scaleSolidVolume1 = static_cast<void(ChemicalState::*)(double)>(&ChemicalState::scaleSolidVolume);
    auto scaleSolidVolume2 = static_cast<void(ChemicalState::*)(double, std::string)>(&ChemicalState::scaleSolidVolume);

    auto scaleVolume1 = static_cast<void(ChemicalState::*)(double)>(&ChemicalState::scaleVolume);
    auto scaleVolume2 = static_cast<void(ChemicalState::*)(double, std::string)>(&ChemicalState::scaleVolume);

    auto speciesAmounts1 = static_cast<VectorConstRef(ChemicalState::*)() const>(&ChemicalState::speciesAmounts);
    auto speciesAmounts2 = static_cast<Vector(ChemicalState::*)(const Indices&) const>(&ChemicalState::speciesAmounts);

    auto speciesAmount1 = static_cast<double(ChemicalState::*)(Index) const>(&ChemicalState::speciesAmount);
    auto speciesAmount2 = static_cast<double(ChemicalState::*)(std::string) const>(&ChemicalState::speciesAmount);
    auto speciesAmount3 = static_cast<double(ChemicalState::*)(Index, std::string) const>(&ChemicalState::speciesAmount);
    auto speciesAmount4 = static_cast<double(ChemicalState::*)(std::string, std::string) const>(&ChemicalState::speciesAmount);

    auto elementAmount1 = static_cast<double(ChemicalState::*)(Index) const>(&ChemicalState::elementAmount);
    auto elementAmount2 = static_cast<double(ChemicalState::*)(std::string) const>(&ChemicalState::elementAmount);
    auto elementAmount3 = static_cast<double(ChemicalState::*)(Index, std::string) const>(&ChemicalState::elementAmount);
    auto elementAmount4 = static_cast<double(ChemicalState::*)(std::string, std::string) const>(&ChemicalState::elementAmount);

    auto elementAmountInPhase1 = static_cast<double(ChemicalState::*)(Index, Index) const>(&ChemicalState::elementAmountInPhase);
    auto elementAmountInPhase2 = static_cast<double(ChemicalState::*)(std::string, std::string) const>(&ChemicalState::elementAmountInPhase);
    auto elementAmountInPhase3 = static_cast<double(ChemicalState::*)(Index, Index, std::string) const>(&ChemicalState::elementAmountInPhase);
    auto elementAmountInPhase4 = static_cast<double(ChemicalState::*)(std::string, std::string, std::string) const>(&ChemicalState::elementAmountInPhase);

    auto elementAmountInSpecies1 = static_cast<double(ChemicalState::*)(Index, const Indices&) const>(&ChemicalState::elementAmountInSpecies);
    auto elementAmountInSpecies2 = static_cast<double(ChemicalState::*)(Index, const Indices&, std::string) const>(&ChemicalState::elementAmountInSpecies);

    auto phaseAmount1 = static_cast<double(ChemicalState::*)(Index) const>(&ChemicalState::phaseAmount);
    auto phaseAmount2 = static_cast<double(ChemicalState::*)(std::string) const>(&ChemicalState::phaseAmount);
    auto phaseAmount3 = static_cast<double(ChemicalState::*)(Index, std::string) const>(&ChemicalState::phaseAmount);
    auto phaseAmount4 = static_cast<double(ChemicalState::*)(std::string, std::string) const>(&ChemicalState::phaseAmount);

    auto output1 = static_cast<void(ChemicalState::*)(std::ostream&, int) const>(&ChemicalState::output);
    auto output2 = static_cast<void(ChemicalState::*)(std::string const&, int) const>(&ChemicalState::output);

    py::class_<ChemicalState>(m, "ChemicalState")
        .def(py::init<const ChemicalSystem&>())
        .def("assign", assignChemicalState)
        .def("clone", cloneChemicalState)
        .def("setTemperature", setTemperature1)
        .def("setTemperature", setTemperature2)
        .def("setPressure", setPressure1)
        .def("setPressure", setPressure2)
        .def("setSpeciesAmounts", setSpeciesAmounts1)
        .def("setSpeciesAmounts", setSpeciesAmounts2)
        .def("setSpeciesAmounts", setSpeciesAmounts3)
        .def("setSpeciesAmount", setSpeciesAmount1)
        .def("setSpeciesAmount", setSpeciesAmount2)
        .def("setSpeciesAmount", setSpeciesAmount3)
        .def("setSpeciesAmount", setSpeciesAmount4)
        .def("setSpeciesMass", setSpeciesMass1)
        .def("setSpeciesMass", setSpeciesMass2)
        .def("setSpeciesMass", setSpeciesMass3)
        .def("setSpeciesMass", setSpeciesMass4)
        .def("setSpeciesDualPotentials", &ChemicalState::setSpeciesDualPotentials)
        .def("setElementDualPotentials", &ChemicalState::setElementDualPotentials)
        .def("scaleSpeciesAmounts", &ChemicalState::scaleSpeciesAmounts)
        .def("scaleSpeciesAmountsInPhase", &ChemicalState::scaleSpeciesAmountsInPhase)
        .def("scalePhaseVolume", scalePhaseVolume1)
        .def("scalePhaseVolume", scalePhaseVolume2)
        .def("scalePhaseVolume", scalePhaseVolume3)
        .def("scalePhaseVolume", scalePhaseVolume4)
        .def("scaleFluidVolume", scaleFluidVolume1)
        .def("scaleFluidVolume", scaleFluidVolume2)
        .def("scaleSolidVolume", scaleSolidVolume1)
        .def("scaleSolidVolume", scaleSolidVolume2)
        .def("scaleVolume", scaleVolume1)
        .def("scaleVolume", scaleVolume2)
        .def("system", &ChemicalState::system, py::return_value_policy::reference_internal)
        .def("temperature", &ChemicalState::temperature)
        .def("pressure", &ChemicalState::pressure)
        .def("speciesAmounts", speciesAmounts1, py::return_value_policy::reference_internal)
        .def("speciesAmounts", speciesAmounts2)
        .def("speciesAmount", speciesAmount1)
        .def("speciesAmount", speciesAmount2)
        .def("speciesAmount", speciesAmount3)
        .def("speciesAmount", speciesAmount4)
        .def("speciesDualPotentials", &ChemicalState::speciesDualPotentials, py::return_value_policy::reference_internal)
        .def("elementAmounts", &ChemicalState::elementAmounts)
        .def("elementAmountsInPhase", &ChemicalState::elementAmountsInPhase)
        .def("elementAmountsInSpecies", &ChemicalState::elementAmountsInSpecies)
        .def("elementAmount", elementAmount1)
        .def("elementAmount", elementAmount2)
        .def("elementAmount", elementAmount3)
        .def("elementAmount", elementAmount4)
        .def("elementAmountInPhase", elementAmountInPhase1)
        .def("elementAmountInPhase", elementAmountInPhase2)
        .def("elementAmountInPhase", elementAmountInPhase3)
        .def("elementAmountInPhase", elementAmountInPhase4)
        .def("elementAmountInSpecies", elementAmountInSpecies1)
        .def("elementAmountInSpecies", elementAmountInSpecies2)
        .def("elementDualPotentials", &ChemicalState::elementDualPotentials, py::return_value_policy::reference_internal)
        .def("phaseAmount", phaseAmount1)
        .def("phaseAmount", phaseAmount2)
        .def("phaseAmount", phaseAmount3)
        .def("phaseAmount", phaseAmount4)
        .def("phaseStabilityIndices", &ChemicalState::phaseStabilityIndices)
        .def("properties", &ChemicalState::properties, py::keep_alive<1, 0>()) // keep returned ChemicalProperties object alive until ChemicalState object is garbage collected!
        .def("output", output1, py::arg("out"), py::arg("precision"))
        .def("output", output2, py::arg("out"), py::arg("precision"))
        .def("__repr__", [](const ChemicalState& self) { std::stringstream ss; ss << self; return ss.str(); })
        .def(py::self + py::self)
        .def(double() * py::self)
        .def(py::self * double())
        ;
}

} // namespace Reaktoro
