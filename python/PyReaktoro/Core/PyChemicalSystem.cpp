// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Interfaces/Gems.hpp>
#include <Reaktoro/Interfaces/Phreeqc.hpp>
#include <Reaktoro/Thermodynamics/Core/ChemicalEditor.hpp>

namespace Reaktoro {
namespace {

auto createChemicalSystemFromChemicalEditor(const ChemicalEditor& editor) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(editor);
}

auto createChemicalSystemFromGems(const Gems& gems) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(gems);
}

auto createChemicalSystemFromPhreeqc(const Phreeqc& phreeqx) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(phreeqx);
}

} // namespace

auto export_ChemicalSystem() -> void
{
    py::class_<ChemicalSystemModel>("ChemicalSystemModel")
        .def_readwrite("standard_gibbs_energies", &ChemicalSystemModel::standard_gibbs_energies)
        .def_readwrite("standard_helmholtz_energies", &ChemicalSystemModel::standard_helmholtz_energies)
        .def_readwrite("standard_internal_energies", &ChemicalSystemModel::standard_internal_energies)
        .def_readwrite("standard_enthalpies", &ChemicalSystemModel::standard_enthalpies)
        .def_readwrite("standard_entropies", &ChemicalSystemModel::standard_entropies)
        .def_readwrite("standard_volumes", &ChemicalSystemModel::standard_volumes)
        .def_readwrite("standard_heat_capacities_cp", &ChemicalSystemModel::standard_heat_capacities_cp)
        .def_readwrite("concentrations", &ChemicalSystemModel::concentrations)
        .def_readwrite("activity_coefficients", &ChemicalSystemModel::activity_coefficients)
        .def_readwrite("activities", &ChemicalSystemModel::activities)
        .def_readwrite("chemical_potentials", &ChemicalSystemModel::chemical_potentials)
        .def_readwrite("phase_molar_volumes", &ChemicalSystemModel::phase_molar_volumes)
        ;

    using return_const_ref = py::return_value_policy<py::copy_const_reference>;

    auto element1 = static_cast<const Element&(ChemicalSystem::*)(Index) const>(&ChemicalSystem::element);
    auto element2 = static_cast<const Element&(ChemicalSystem::*)(std::string) const>(&ChemicalSystem::element);

    auto species1 = static_cast<const std::vector<Species>&(ChemicalSystem::*)() const>(&ChemicalSystem::species);
    auto species2 = static_cast<const Species&(ChemicalSystem::*)(Index) const>(&ChemicalSystem::species);
    auto species3 = static_cast<const Species&(ChemicalSystem::*)(std::string) const>(&ChemicalSystem::species);

    auto phase1 = static_cast<const Phase&(ChemicalSystem::*)(Index) const>(&ChemicalSystem::phase);
    auto phase2 = static_cast<const Phase&(ChemicalSystem::*)(std::string) const>(&ChemicalSystem::phase);

    auto indicesElementsInSpecies1 = static_cast<Indices(ChemicalSystem::*)(Index) const>(&ChemicalSystem::indicesElementsInSpecies);
    auto indicesElementsInSpecies2 = static_cast<Indices(ChemicalSystem::*)(const Indices&) const>(&ChemicalSystem::indicesElementsInSpecies);

    py::class_<ChemicalSystem>("ChemicalSystem")
        .def(py::init<>())
        .def(py::init<const std::vector<Phase>&>())
        .def(py::init<const std::vector<Phase>&, const ChemicalSystemModel&>())
        .def("__init__", py::make_constructor(createChemicalSystemFromChemicalEditor))
        .def("__init__", py::make_constructor(createChemicalSystemFromGems))
        .def("__init__", py::make_constructor(createChemicalSystemFromPhreeqc))
        .def("numElements", &ChemicalSystem::numElements)
        .def("numSpecies", &ChemicalSystem::numSpecies)
        .def("numSpeciesInPhase", &ChemicalSystem::numSpeciesInPhase)
        .def("numPhases", &ChemicalSystem::numPhases)
        .def("elements", &ChemicalSystem::elements, return_const_ref())
        .def("species", species1, return_const_ref())
        .def("phases", &ChemicalSystem::phases, return_const_ref())
        .def("formulaMatrix", &ChemicalSystem::formulaMatrix, return_const_ref())
        .def("element", element1, return_const_ref())
        .def("element", element2, return_const_ref())
        .def("species", species2, return_const_ref())
        .def("species", species3, return_const_ref())
        .def("phase", phase1, return_const_ref())
        .def("phase", phase2, return_const_ref())
        .def("indexElement", &ChemicalSystem::indexElement)
        .def("indexElementWithError", &ChemicalSystem::indexElementWithError)
        .def("indexSpecies", &ChemicalSystem::indexSpecies)
        .def("indexSpeciesWithError", &ChemicalSystem::indexSpeciesWithError)
        .def("indexPhase", &ChemicalSystem::indexPhase)
        .def("indexPhaseWithError", &ChemicalSystem::indexPhaseWithError)
        .def("indexPhaseWithSpecies", &ChemicalSystem::indexPhaseWithSpecies)
        .def("indicesElements", &ChemicalSystem::indicesElements)
        .def("indicesSpecies", &ChemicalSystem::indicesSpecies)
        .def("indicesPhases", &ChemicalSystem::indicesPhases)
        .def("indicesPhasesWithSpecies", &ChemicalSystem::indicesPhasesWithSpecies)
        .def("indicesElementsInSpecies", indicesElementsInSpecies1)
        .def("indicesElementsInSpecies", indicesElementsInSpecies2)
        .def("indexFirstSpeciesInPhase", &ChemicalSystem::indexFirstSpeciesInPhase)
        .def("standardGibbsEnergies", &ChemicalSystem::standardGibbsEnergies)
        .def("standardEnthalpies", &ChemicalSystem::standardEnthalpies)
        .def("standardHelmholtzEnergies", &ChemicalSystem::standardHelmholtzEnergies)
        .def("standardEntropies", &ChemicalSystem::standardEntropies)
        .def("standardVolumes", &ChemicalSystem::standardVolumes)
        .def("standardInternalEnergies", &ChemicalSystem::standardInternalEnergies)
        .def("standardHeatCapacities", &ChemicalSystem::standardHeatCapacities)
        .def("molarFractions", &ChemicalSystem::molarFractions)
        .def("concentrations", &ChemicalSystem::concentrations)
        .def("activityCoefficients", &ChemicalSystem::activityCoefficients)
        .def("activityConstants", &ChemicalSystem::activityConstants)
        .def("activities", &ChemicalSystem::activities)
        .def("chemicalPotentials", &ChemicalSystem::chemicalPotentials)
        .def("phaseMolarVolumes", &ChemicalSystem::phaseMolarVolumes)
        .def("phaseDensities", &ChemicalSystem::phaseDensities)
        .def("phaseMolarAmounts", &ChemicalSystem::phaseMolarAmounts)
        .def("phaseMassAmounts", &ChemicalSystem::phaseMassAmounts)
        .def("phaseVolumes", &ChemicalSystem::phaseVolumes)
        .def("elementAmounts", &ChemicalSystem::elementAmounts)
        .def("elementAmountsInPhase", &ChemicalSystem::elementAmountsInPhase)
        .def("elementAmountsInSpecies", &ChemicalSystem::elementAmountsInSpecies)
        .def("elementAmount", &ChemicalSystem::elementAmount)
        .def("elementAmountInPhase", &ChemicalSystem::elementAmountInPhase)
        .def("elementAmountInSpecies", &ChemicalSystem::elementAmountInSpecies)
        .def(py::self_ns::str(py::self_ns::self));
        ;
}

} // namespace Reaktoro
