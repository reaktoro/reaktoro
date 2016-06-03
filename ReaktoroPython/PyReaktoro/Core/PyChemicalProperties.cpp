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
#include <Reaktoro/Core/ChemicalPropertiesAqueousPhase.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

auto export_ChemicalProperties() -> void
{
    auto update1 = static_cast<void (ChemicalProperties::*)(double, double)>(&ChemicalProperties::update);
    auto update2 = static_cast<void (ChemicalProperties::*)(double, double, const Vector&)>(&ChemicalProperties::update);

    py::class_<ChemicalProperties>("ChemicalProperties")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def("update", update1)
        .def("update", update2)
        .def("temperature", &ChemicalProperties::temperature)
        .def("pressure", &ChemicalProperties::pressure)
        .def("composition", &ChemicalProperties::composition, py::return_internal_reference<>())
        .def("system", &ChemicalProperties::system, py::return_internal_reference<>())
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
        .def("subvolume", &ChemicalProperties::subvolume)
        .def("fluidVolume", &ChemicalProperties::fluidVolume)
        .def("solidVolume", &ChemicalProperties::solidVolume)
        .def("aqueous", &ChemicalProperties::aqueous)
        ;
}

} // namespace Reaktoro
