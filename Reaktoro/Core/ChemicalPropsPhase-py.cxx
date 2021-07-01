// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
using namespace Reaktoro;

void exportChemicalPropsPhase(py::module& m)
{
    auto update1 = [](ChemicalPropsPhase& self, const real& T, const real& P, ArrayXrConstRef n)
    {
        self.update(T, P, n);
    };

    auto updateIdeal1 = [](ChemicalPropsPhase& self, const real& T, const real& P, ArrayXrConstRef n)
    {
        self.updateIdeal(T, P, n);
    };

    py::class_<ChemicalPropsPhase>(m, "ChemicalPropsPhase")
        .def(py::init<const Phase&>())
        .def("update", update1)
        .def("updateIdeal", updateIdeal1)
        .def("phase", &ChemicalPropsPhase::phase, py::return_value_policy::reference_internal)
        .def("data", &ChemicalPropsPhase::data, py::return_value_policy::reference_internal)
        .def("temperature", &ChemicalPropsPhase::temperature)
        .def("pressure", &ChemicalPropsPhase::pressure)
        .def("speciesAmounts", &ChemicalPropsPhase::speciesAmounts, py::return_value_policy::reference_internal)
        .def("moleFractions", &ChemicalPropsPhase::moleFractions, py::return_value_policy::reference_internal)
        .def("lnActivityCoefficients", &ChemicalPropsPhase::lnActivityCoefficients, py::return_value_policy::reference_internal)
        .def("lnActivities", &ChemicalPropsPhase::lnActivities, py::return_value_policy::reference_internal)
        .def("chemicalPotentials", &ChemicalPropsPhase::chemicalPotentials, py::return_value_policy::reference_internal)
        .def("standardGibbsEnergies", &ChemicalPropsPhase::standardGibbsEnergies, py::return_value_policy::reference_internal)
        .def("standardEnthalpies", &ChemicalPropsPhase::standardEnthalpies, py::return_value_policy::reference_internal)
        .def("standardVolumes", &ChemicalPropsPhase::standardVolumes, py::return_value_policy::reference_internal)
        .def("standardEntropies", &ChemicalPropsPhase::standardEntropies)
        .def("standardInternalEnergies", &ChemicalPropsPhase::standardInternalEnergies)
        .def("standardHelmholtzEnergies", &ChemicalPropsPhase::standardHelmholtzEnergies)
        .def("standardHeatCapacitiesConstP", &ChemicalPropsPhase::standardHeatCapacitiesConstP, py::return_value_policy::reference_internal)
        .def("standardHeatCapacitiesConstV", &ChemicalPropsPhase::standardHeatCapacitiesConstV, py::return_value_policy::reference_internal)
        .def("molarGibbsEnergy", &ChemicalPropsPhase::molarGibbsEnergy)
        .def("molarEnthalpy", &ChemicalPropsPhase::molarEnthalpy)
        .def("molarVolume", &ChemicalPropsPhase::molarVolume)
        .def("molarEntropy", &ChemicalPropsPhase::molarEntropy)
        .def("molarInternalEnergy", &ChemicalPropsPhase::molarInternalEnergy)
        .def("molarHelmholtzEnergy", &ChemicalPropsPhase::molarHelmholtzEnergy)
        .def("molarHeatCapacityConstP", &ChemicalPropsPhase::molarHeatCapacityConstP)
        .def("molarHeatCapacityConstV", &ChemicalPropsPhase::molarHeatCapacityConstV)
        .def("molarDensity", &ChemicalPropsPhase::molarDensity)
        .def("amount", &ChemicalPropsPhase::amount)
        .def("mass", &ChemicalPropsPhase::mass)
        .def("gibbsEnergy", &ChemicalPropsPhase::gibbsEnergy)
        .def("enthalpy", &ChemicalPropsPhase::enthalpy)
        .def("volume", &ChemicalPropsPhase::volume)
        .def("entropy", &ChemicalPropsPhase::entropy)
        .def("internalEnergy", &ChemicalPropsPhase::internalEnergy)
        .def("helmholtzEnergy", &ChemicalPropsPhase::helmholtzEnergy)
        ;
}
