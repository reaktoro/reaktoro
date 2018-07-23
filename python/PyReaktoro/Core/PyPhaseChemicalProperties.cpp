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

#include "PyPhaseChemicalProperties.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/PhaseChemicalProperties.hpp>

namespace Reaktoro {

auto export_PhaseChemicalProperties() -> void
{
    auto update1 = static_cast<void (PhaseChemicalProperties::*)(double, double)>(&PhaseChemicalProperties::update);
    auto update2 = static_cast<void (PhaseChemicalProperties::*)(double, double, const Vector&)>(&PhaseChemicalProperties::update);

    py::class_<PhaseChemicalProperties>("PhaseChemicalProperties")
        .def(py::init<>())
        .def(py::init<const Phase&>())
        .def("update", update1)
        .def("update", update2)
        .def("temperature", &PhaseChemicalProperties::temperature)
        .def("pressure", &PhaseChemicalProperties::pressure)
        .def("composition", &PhaseChemicalProperties::composition, py::return_internal_reference<>())
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
