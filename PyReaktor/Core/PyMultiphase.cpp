// Reaktor is a C++ library for computational reaction modelling.
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

#include "PyMultiphase.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Core/Multiphase.hpp>
#include <Reaktor/Interfaces/Gems.hpp>
#include <Reaktor/Interfaces/Phreeqx.hpp>

namespace Reaktor {

auto export_Multiphase() -> void
{
    py::class_<MultiphaseModel>("MultiphaseModel")
        .def_readwrite("standard_gibbs_energy_fn", &MultiphaseModel::standard_gibbs_energy_fn)
        .def_readwrite("standard_helmholtz_energy_fn", &MultiphaseModel::standard_helmholtz_energy_fn)
        .def_readwrite("standard_internal_energy_fn", &MultiphaseModel::standard_internal_energy_fn)
        .def_readwrite("standard_enthalpy_fn", &MultiphaseModel::standard_enthalpy_fn)
        .def_readwrite("standard_entropy_fn", &MultiphaseModel::standard_entropy_fn)
        .def_readwrite("standard_volume_fn", &MultiphaseModel::standard_volume_fn)
        .def_readwrite("standard_heat_capacity_fn", &MultiphaseModel::standard_heat_capacity_fn)
        .def_readwrite("concentration_fn", &MultiphaseModel::concentration_fn)
        .def_readwrite("activity_coefficient_fn", &MultiphaseModel::activity_coefficient_fn)
        .def_readwrite("activity_fn", &MultiphaseModel::activity_fn)
        .def_readwrite("chemical_potential_fn", &MultiphaseModel::chemical_potential_fn)
        .def_readwrite("phase_molar_volume_fn", &MultiphaseModel::phase_molar_volume_fn)
        ;

    auto element1 = static_cast<const Element&(Multiphase::*)(Index) const>(&Multiphase::element);
    auto element2 = static_cast<const Element&(Multiphase::*)(std::string) const>(&Multiphase::element);

    auto species1 = static_cast<const std::vector<Species>&(Multiphase::*)() const>(&Multiphase::species);
    auto species2 = static_cast<const Species&(Multiphase::*)(Index) const>(&Multiphase::species);
    auto species3 = static_cast<const Species&(Multiphase::*)(std::string) const>(&Multiphase::species);

    auto phase1 = static_cast<const Phase&(Multiphase::*)(Index) const>(&Multiphase::phase);
    auto phase2 = static_cast<const Phase&(Multiphase::*)(std::string) const>(&Multiphase::phase);

    auto indicesElementsInSpecies1 = static_cast<Indices(Multiphase::*)(Index) const>(&Multiphase::indicesElementsInSpecies);
    auto indicesElementsInSpecies2 = static_cast<Indices(Multiphase::*)(const Indices&) const>(&Multiphase::indicesElementsInSpecies);

    py::class_<Multiphase>("Multiphase")
        .def(py::init<>())
        .def(py::init<const std::vector<Phase>&>())
        .def(py::init<const std::vector<Phase>&, const MultiphaseModel&>())
        .def("numElements", &Multiphase::numElements)
        .def("numSpecies", &Multiphase::numSpecies)
        .def("numSpeciesInPhase", &Multiphase::numSpeciesInPhase)
        .def("numPhases", &Multiphase::numPhases)
        .def("elements", &Multiphase::elements, py::return_value_policy<py::copy_const_reference>())
        .def("species", species1, py::return_value_policy<py::copy_const_reference>())
        .def("phases", &Multiphase::phases, py::return_value_policy<py::copy_const_reference>())
        .def("formulaMatrix", &Multiphase::formulaMatrix, py::return_value_policy<py::copy_const_reference>())
        .def("element", element1, py::return_value_policy<py::copy_const_reference>())
        .def("element", element2, py::return_value_policy<py::copy_const_reference>())
        .def("species", species2, py::return_value_policy<py::copy_const_reference>())
        .def("species", species3, py::return_value_policy<py::copy_const_reference>())
        .def("phase", phase1, py::return_value_policy<py::copy_const_reference>())
        .def("phase", phase2, py::return_value_policy<py::copy_const_reference>())
        .def("indexElement", &Multiphase::indexElement)
        .def("indexElementWithError", &Multiphase::indexElementWithError)
        .def("indexSpecies", &Multiphase::indexSpecies)
        .def("indexSpeciesWithError", &Multiphase::indexSpeciesWithError)
        .def("indexPhase", &Multiphase::indexPhase)
        .def("indexPhaseWithError", &Multiphase::indexPhaseWithError)
        .def("indexPhaseWithSpecies", &Multiphase::indexPhaseWithSpecies)
        .def("indicesElements", &Multiphase::indicesElements)
        .def("indicesSpecies", &Multiphase::indicesSpecies)
        .def("indicesPhases", &Multiphase::indicesPhases)
        .def("indicesPhasesWithSpecies", &Multiphase::indicesPhasesWithSpecies)
        .def("indicesElementsInSpecies", indicesElementsInSpecies1)
        .def("indicesElementsInSpecies", indicesElementsInSpecies2)
        .def("indexFirstSpeciesInPhase", &Multiphase::indexFirstSpeciesInPhase)
        .def("standardGibbsEnergies", &Multiphase::standardGibbsEnergies)
        .def("standardEnthalpies", &Multiphase::standardEnthalpies)
        .def("standardHelmholtzEnergies", &Multiphase::standardHelmholtzEnergies)
        .def("standardEntropies", &Multiphase::standardEntropies)
        .def("standardVolumes", &Multiphase::standardVolumes)
        .def("standardInternalEnergies", &Multiphase::standardInternalEnergies)
        .def("standardHeatCapacities", &Multiphase::standardHeatCapacities)
        .def("concentrations", &Multiphase::concentrations)
        .def("activityCoefficients", &Multiphase::activityCoefficients)
        .def("activities", &Multiphase::activities)
        .def("chemicalPotentials", &Multiphase::chemicalPotentials)
        .def("phaseMolarVolumes", &Multiphase::phaseMolarVolumes)
        .def("phaseDensities", &Multiphase::phaseDensities)
        .def("phaseMolarAmounts", &Multiphase::phaseMolarAmounts)
        .def("phaseMassAmounts", &Multiphase::phaseMassAmounts)
        .def("phaseVolumes", &Multiphase::phaseVolumes)
        .def("elementAmounts", &Multiphase::elementAmounts)
        .def("elementAmountsInPhase", &Multiphase::elementAmountsInPhase)
        .def("elementAmountsInSpecies", &Multiphase::elementAmountsInSpecies)
        .def("elementAmount", &Multiphase::elementAmount)
        .def("elementAmountInPhase", &Multiphase::elementAmountInPhase)
        .def("elementAmountsInSpecies", &Multiphase::elementAmountsInSpecies)
        .def(py::self_ns::str(py::self_ns::self));
        ;
}

} // namespace Reaktor
