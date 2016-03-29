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

#include "PyChemicalProperties.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

auto export_ChemicalProperties() -> void
{
    py::class_<ChemicalProperties>("ChemicalProperties")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def("temperature", &ChemicalProperties::temperature)
        .def("pressure", &ChemicalProperties::pressure)
        .def("composition", &ChemicalProperties::composition, py::return_internal_reference<>())
        .def("molarFractions", &ChemicalProperties::molarFractions)
        .def("lnActivityCoefficients", &ChemicalProperties::lnActivityCoefficients)
        .def("lnActivities", &ChemicalProperties::lnActivities)
        .def("chemicalPotentials", &ChemicalProperties::chemicalPotentials)
        .def("standardPartialMolarGibbsEnergies", &ChemicalProperties::standardPartialMolarGibbsEnergies)
        .def("standardPartialMolarEnthalpies", &ChemicalProperties::standardPartialMolarEnthalpies)
        .def("standardPartialMolarVolumes", &ChemicalProperties::standardPartialMolarVolumes)
        .def("standardPartialMolarEntropies", &ChemicalProperties::standardPartialMolarEntropies)
        .def("standardPartialMolarInternalEnergies", &ChemicalProperties::standardPartialMolarInternalEnergies)
        .def("standardPartialMolarHelmholtzEnergies", &ChemicalProperties::standardPartialMolarHelmholtzEnergies)
        .def("standardPartialMolarHeatCapacitiesConstP", &ChemicalProperties::standardPartialMolarHeatCapacitiesConstP)
        .def("standardPartialMolarHeatCapacitiesConstV", &ChemicalProperties::standardPartialMolarHeatCapacitiesConstV)
        .def("phaseMolarGibbsEnergies", &ChemicalProperties::phaseMolarGibbsEnergies)
        .def("phaseMolarEnthalpies", &ChemicalProperties::phaseMolarEnthalpies)
        .def("phaseMolarVolumes", &ChemicalProperties::phaseMolarVolumes)
        .def("phaseMolarEntropies", &ChemicalProperties::phaseMolarEntropies)
        .def("phaseMolarInternalEnergies", &ChemicalProperties::phaseMolarInternalEnergies)
        .def("phaseMolarHelmholtzEnergies", &ChemicalProperties::phaseMolarHelmholtzEnergies)
        .def("phaseMolarHeatCapacitiesConstP", &ChemicalProperties::phaseMolarHeatCapacitiesConstP)
        .def("phaseMolarHeatCapacitiesConstV", &ChemicalProperties::phaseMolarHeatCapacitiesConstV)
        .def("phaseSpecificGibbsEnergies", &ChemicalProperties::phaseSpecificGibbsEnergies)
        .def("phaseSpecificEnthalpies", &ChemicalProperties::phaseSpecificEnthalpies)
        .def("phaseSpecificVolumes", &ChemicalProperties::phaseSpecificVolumes)
        .def("phaseSpecificEntropies", &ChemicalProperties::phaseSpecificEntropies)
        .def("phaseSpecificInternalEnergies", &ChemicalProperties::phaseSpecificInternalEnergies)
        .def("phaseSpecificHelmholtzEnergies", &ChemicalProperties::phaseSpecificHelmholtzEnergies)
        .def("phaseSpecificHeatCapacitiesConstP", &ChemicalProperties::phaseSpecificHeatCapacitiesConstP)
        .def("phaseSpecificHeatCapacitiesConstV", &ChemicalProperties::phaseSpecificHeatCapacitiesConstV)
        .def("phaseDensities", &ChemicalProperties::phaseDensities)
        .def("phaseMasses", &ChemicalProperties::phaseMasses)
        .def("phaseAmounts", &ChemicalProperties::phaseAmounts)
        .def("phaseVolumes", &ChemicalProperties::phaseVolumes)
        .def("volume", &ChemicalProperties::volume)
        ;

    py::class_<PhaseChemicalProperties>("PhaseChemicalProperties")
        .def(py::init<>())
        .def(py::init<unsigned>())
        .def("temperature", &PhaseChemicalProperties::temperature)
        .def("pressure", &PhaseChemicalProperties::pressure)
        .def("composition", &PhaseChemicalProperties::composition)
        .def("molarFractions", &PhaseChemicalProperties::molarFractions)
        .def("lnActivityCoefficients", &PhaseChemicalProperties::lnActivityCoefficients)
        .def("lnActivities", &PhaseChemicalProperties::lnActivities)
        .def("chemicalPotentials", &PhaseChemicalProperties::chemicalPotentials)
        .def("standardPartialMolarGibbsEnergies", &PhaseChemicalProperties::standardPartialMolarGibbsEnergies)
        .def("standardPartialMolarEnthalpies", &PhaseChemicalProperties::standardPartialMolarEnthalpies)
        .def("standardPartialMolarVolumes", &PhaseChemicalProperties::standardPartialMolarVolumes)
        .def("standardPartialMolarEntropies", &PhaseChemicalProperties::standardPartialMolarEntropies)
        .def("standardPartialMolarInternalEnergies", &PhaseChemicalProperties::standardPartialMolarInternalEnergies)
        .def("standardPartialMolarHelmholtzEnergies", &PhaseChemicalProperties::standardPartialMolarHelmholtzEnergies)
        .def("standardPartialMolarHeatCapacitiesConstP", &PhaseChemicalProperties::standardPartialMolarHeatCapacitiesConstP)
        .def("standardPartialMolarHeatCapacitiesConstV", &PhaseChemicalProperties::standardPartialMolarHeatCapacitiesConstV)
        .def("molarGibbsEnergy", &PhaseChemicalProperties::molarGibbsEnergy)
        .def("molarEnthalpy", &PhaseChemicalProperties::molarEnthalpy)
        .def("molarVolume", &PhaseChemicalProperties::molarVolume)
        .def("molarEntropy", &PhaseChemicalProperties::molarEntropy)
        .def("molarInternalEnergy", &PhaseChemicalProperties::molarInternalEnergy)
        .def("molarHelmholtzEnergy", &PhaseChemicalProperties::molarHelmholtzEnergy)
        .def("molarHeatCapacityConstP", &PhaseChemicalProperties::molarHeatCapacityConstP)
        .def("molarHeatCapacityConstV", &PhaseChemicalProperties::molarHeatCapacityConstV)
        .def("specificGibbsEnergy", &PhaseChemicalProperties::specificGibbsEnergy)
        .def("specificEnthalpy", &PhaseChemicalProperties::specificEnthalpy)
        .def("specificVolume", &PhaseChemicalProperties::specificVolume)
        .def("specificEntropy", &PhaseChemicalProperties::specificEntropy)
        .def("specificInternalEnergy", &PhaseChemicalProperties::specificInternalEnergy)
        .def("specificHelmholtzEnergy", &PhaseChemicalProperties::specificHelmholtzEnergy)
        .def("specificHeatCapacityConstP", &PhaseChemicalProperties::specificHeatCapacityConstP)
        .def("specificHeatCapacityConstV", &PhaseChemicalProperties::specificHeatCapacityConstV)
        .def("density", &PhaseChemicalProperties::density)
        .def("mass", &PhaseChemicalProperties::mass)
        .def("amount", &PhaseChemicalProperties::amount)
        .def("volume", &PhaseChemicalProperties::volume)
        ;
}

} // namespace Reaktoro
