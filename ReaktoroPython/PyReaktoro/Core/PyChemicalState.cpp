// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "PyChemicalState.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalPropertiesAqueousPhase.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>
#include <Reaktoro/Interfaces/Gems.hpp>
#include <Reaktoro/Interfaces/Phreeqc.hpp>

namespace Reaktoro {
namespace {

auto createChemicalStateGems(Gems& gems) -> boost::shared_ptr<ChemicalState>
{
    return boost::make_shared<ChemicalState>(gems.operator ChemicalState());
}

auto createChemicalStatePhreeqc(Phreeqc& phreeqc) -> boost::shared_ptr<ChemicalState>
{
    return boost::make_shared<ChemicalState>(phreeqc.operator ChemicalState());
}

auto assignChemicalState(ChemicalState& state, const ChemicalState& other) -> void
{
    state = other;
}

auto cloneChemicalState(ChemicalState& state) -> ChemicalState
{
    return state;
}

}  // namespace

auto export_ChemicalState() -> void
{
    auto setTemperature1 = static_cast<void(ChemicalState::*)(double)>(&ChemicalState::setTemperature);
    auto setTemperature2 = static_cast<void(ChemicalState::*)(double, std::string)>(&ChemicalState::setTemperature);

    auto setPressure1 = static_cast<void(ChemicalState::*)(double)>(&ChemicalState::setPressure);
    auto setPressure2 = static_cast<void(ChemicalState::*)(double, std::string)>(&ChemicalState::setPressure);

    auto setSpeciesAmounts1 = static_cast<void(ChemicalState::*)(double)>(&ChemicalState::setSpeciesAmounts);
    auto setSpeciesAmounts2 = static_cast<void(ChemicalState::*)(const Vector&)>(&ChemicalState::setSpeciesAmounts);

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

    auto speciesAmounts1 = static_cast<const Vector&(ChemicalState::*)() const>(&ChemicalState::speciesAmounts);
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

    py::class_<ChemicalState>("ChemicalState")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def("__init__", py::make_constructor(createChemicalStateGems))
        .def("__init__", py::make_constructor(createChemicalStatePhreeqc))
        .def("assign", assignChemicalState)
        .def("clone", cloneChemicalState)
        .def("setTemperature", setTemperature1)
        .def("setTemperature", setTemperature2)
        .def("setPressure", setPressure1)
        .def("setPressure", setPressure2)
        .def("setSpeciesAmounts", setSpeciesAmounts1)
        .def("setSpeciesAmounts", setSpeciesAmounts2)
        .def("setSpeciesAmount", setSpeciesAmount1)
        .def("setSpeciesAmount", setSpeciesAmount2)
        .def("setSpeciesAmount", setSpeciesAmount3)
        .def("setSpeciesAmount", setSpeciesAmount4)
        .def("setSpeciesMass", setSpeciesMass1)
        .def("setSpeciesMass", setSpeciesMass2)
        .def("setSpeciesMass", setSpeciesMass3)
        .def("setSpeciesMass", setSpeciesMass4)
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
        .def("scaleVolume", &ChemicalState::scaleVolume)
        .def("system", &ChemicalState::system, py::return_internal_reference<>())
        .def("temperature", &ChemicalState::temperature)
        .def("pressure", &ChemicalState::pressure)
        .def("speciesAmounts", speciesAmounts1, py::return_internal_reference<>())
        .def("speciesAmounts", speciesAmounts2)
        .def("speciesAmount", speciesAmount1)
        .def("speciesAmount", speciesAmount2)
        .def("speciesAmount", speciesAmount3)
        .def("speciesAmount", speciesAmount4)
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
        .def("phaseAmount", phaseAmount1)
        .def("phaseAmount", phaseAmount2)
        .def("phaseAmount", phaseAmount3)
        .def("phaseAmount", phaseAmount4)
        .def("properties", &ChemicalState::properties)
        .def("output", &ChemicalState::output)
        .def(py::self + py::self)
        .def(double() * py::self)
        .def(py::self * double())
        .def(py::self_ns::str(py::self_ns::self))
        ;
}

} // namespace Reaktoro
