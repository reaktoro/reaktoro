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

#include "PyChemicalSystem.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ThermoProperties.hpp>
#include <Reaktoro/Interfaces/Gems.hpp>
#include <Reaktoro/Interfaces/Phreeqc.hpp>
#include <Reaktoro/Interfaces/PhreeqcEditor.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>

namespace Reaktoro {
namespace {

auto createChemicalSystemFromChemicalEditor(const ChemicalEditor& editor) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(editor);
}

auto createChemicalSystemFromGems(Gems& gems) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(gems);
}

auto createChemicalSystemFromPhreeqc(Phreeqc& phreeqc) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(phreeqc);
}

auto createChemicalSystemFromPhreeqcEditor(PhreeqcEditor& editor) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(editor);
}

} // namespace

auto export_ChemicalSystem() -> void
{
    auto element1 = static_cast<const Element&(ChemicalSystem::*)(Index) const>(&ChemicalSystem::element);
    auto element2 = static_cast<const Element&(ChemicalSystem::*)(std::string) const>(&ChemicalSystem::element);

    auto species1 = static_cast<const std::vector<Species>&(ChemicalSystem::*)() const>(&ChemicalSystem::species);
    auto species2 = static_cast<const Species&(ChemicalSystem::*)(Index) const>(&ChemicalSystem::species);
    auto species3 = static_cast<const Species&(ChemicalSystem::*)(std::string) const>(&ChemicalSystem::species);

    auto phase1 = static_cast<const Phase&(ChemicalSystem::*)(Index) const>(&ChemicalSystem::phase);
    auto phase2 = static_cast<const Phase&(ChemicalSystem::*)(std::string) const>(&ChemicalSystem::phase);

    auto indicesElementsInSpecies1 = static_cast<Indices(ChemicalSystem::*)(Index) const>(&ChemicalSystem::indicesElementsInSpecies);
    auto indicesElementsInSpecies2 = static_cast<Indices(ChemicalSystem::*)(const Indices&) const>(&ChemicalSystem::indicesElementsInSpecies);

    auto properties1 = static_cast<ThermoProperties(ChemicalSystem::*)(double,double) const>(&ChemicalSystem::properties);
    auto properties2 = static_cast<ChemicalProperties(ChemicalSystem::*)(double,double,const Vector&) const>(&ChemicalSystem::properties);

    py::class_<ChemicalSystem>("ChemicalSystem")
        .def(py::init<>())
        .def(py::init<const std::vector<Phase>&>())
        .def("__init__", py::make_constructor(createChemicalSystemFromChemicalEditor))
        .def("__init__", py::make_constructor(createChemicalSystemFromGems))
        .def("__init__", py::make_constructor(createChemicalSystemFromPhreeqc))
        .def("__init__", py::make_constructor(createChemicalSystemFromPhreeqcEditor))
        .def("numElements", &ChemicalSystem::numElements)
        .def("numSpecies", &ChemicalSystem::numSpecies)
        .def("numSpeciesInPhase", &ChemicalSystem::numSpeciesInPhase)
        .def("numPhases", &ChemicalSystem::numPhases)
        .def("elements", &ChemicalSystem::elements, py::return_internal_reference<>())
        .def("species", species1, py::return_internal_reference<>())
        .def("phases", &ChemicalSystem::phases, py::return_internal_reference<>())
        .def("formulaMatrix", &ChemicalSystem::formulaMatrix, py::return_internal_reference<>())
        .def("element", element1, py::return_internal_reference<>())
        .def("element", element2, py::return_internal_reference<>())
        .def("species", species2, py::return_internal_reference<>())
        .def("species", species3, py::return_internal_reference<>())
        .def("phase", phase1, py::return_internal_reference<>())
        .def("phase", phase2, py::return_internal_reference<>())
        .def("indexElement", &ChemicalSystem::indexElement)
        .def("indexElementWithError", &ChemicalSystem::indexElementWithError)
        .def("indexSpecies", &ChemicalSystem::indexSpecies)
        .def("indexSpeciesWithError", &ChemicalSystem::indexSpeciesWithError)
        .def("indexSpeciesAny", &ChemicalSystem::indexSpeciesAny)
        .def("indexSpeciesAnyWithError", &ChemicalSystem::indexSpeciesAnyWithError)
        .def("indexPhase", &ChemicalSystem::indexPhase)
        .def("indexPhaseWithError", &ChemicalSystem::indexPhaseWithError)
        .def("indexPhaseWithSpecies", &ChemicalSystem::indexPhaseWithSpecies)
        .def("indexFirstSpeciesInPhase", &ChemicalSystem::indexFirstSpeciesInPhase)
        .def("indicesElements", &ChemicalSystem::indicesElements)
        .def("indicesElementsInSpecies", indicesElementsInSpecies1)
        .def("indicesElementsInSpecies", indicesElementsInSpecies2)
        .def("indicesSpecies", &ChemicalSystem::indicesSpecies)
        .def("indicesSpeciesInPhases", &ChemicalSystem::indicesSpeciesInPhases)
        .def("indicesPhases", &ChemicalSystem::indicesPhases)
        .def("indicesPhasesWithSpecies", &ChemicalSystem::indicesPhasesWithSpecies)
        .def("indicesFluidPhases", &ChemicalSystem::indicesFluidPhases)
        .def("indicesFluidSpecies", &ChemicalSystem::indicesFluidSpecies)
        .def("indicesSolidPhases", &ChemicalSystem::indicesSolidPhases)
        .def("indicesSolidSpecies", &ChemicalSystem::indicesSolidSpecies)
        .def("elementAmounts", &ChemicalSystem::elementAmounts)
        .def("elementAmountsInPhase", &ChemicalSystem::elementAmountsInPhase)
        .def("elementAmountsInSpecies", &ChemicalSystem::elementAmountsInSpecies)
        .def("elementAmount", &ChemicalSystem::elementAmount)
        .def("elementAmountInPhase", &ChemicalSystem::elementAmountInPhase)
        .def("elementAmountInSpecies", &ChemicalSystem::elementAmountInSpecies)
        .def("properties", properties1)
        .def("properties", properties2)
        .def(py::self_ns::str(py::self_ns::self));
        ;
}

} // namespace Reaktoro
