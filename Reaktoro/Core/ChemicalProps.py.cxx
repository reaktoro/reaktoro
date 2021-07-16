// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
using namespace Reaktoro;

void exportChemicalProps(py::module& m)
{
    py::class_<ChemicalProps>(m, "ChemicalProps")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalState&>())
        .def("update", py::overload_cast<const ChemicalState&>(&ChemicalProps::update))
        .def("update", py::overload_cast<const real&, const real&, ArrayXrConstRef>(&ChemicalProps::update))
        .def("update", py::overload_cast<ArrayXrConstRef>(&ChemicalProps::update))
        .def("update", py::overload_cast<ArrayXdConstRef>(&ChemicalProps::update))
        .def("updateIdeal", py::overload_cast<const ChemicalState&>(&ChemicalProps::updateIdeal))
        .def("updateIdeal", py::overload_cast<const real&, const real&, ArrayXrConstRef>(&ChemicalProps::updateIdeal))
        .def("system", &ChemicalProps::system, py::return_value_policy::reference_internal)
        .def("phaseProps", &ChemicalProps::phaseProps, py::return_value_policy::reference_internal)
        .def("temperature", &ChemicalProps::temperature, py::return_value_policy::reference_internal)
        .def("pressure", &ChemicalProps::pressure, py::return_value_policy::reference_internal)
        .def("speciesAmounts", &ChemicalProps::speciesAmounts, py::return_value_policy::reference_internal)
        .def("moleFractions", &ChemicalProps::moleFractions, py::return_value_policy::reference_internal)
        .def("lnActivityCoefficients", &ChemicalProps::lnActivityCoefficients, py::return_value_policy::reference_internal)
        .def("lnActivities", &ChemicalProps::lnActivities, py::return_value_policy::reference_internal)
        .def("chemicalPotentials", &ChemicalProps::chemicalPotentials, py::return_value_policy::reference_internal)
        .def("standardVolumes", &ChemicalProps::standardVolumes, py::return_value_policy::reference_internal)
        .def("standardGibbsEnergies", &ChemicalProps::standardGibbsEnergies, py::return_value_policy::reference_internal)
        .def("standardEnthalpies", &ChemicalProps::standardEnthalpies, py::return_value_policy::reference_internal)
        .def("standardEntropies", &ChemicalProps::standardEntropies)
        .def("standardInternalEnergies", &ChemicalProps::standardInternalEnergies)
        .def("standardHelmholtzEnergies", &ChemicalProps::standardHelmholtzEnergies)
        .def("standardHeatCapacitiesConstP", &ChemicalProps::standardHeatCapacitiesConstP, py::return_value_policy::reference_internal)
        .def("standardHeatCapacitiesConstV", &ChemicalProps::standardHeatCapacitiesConstV, py::return_value_policy::reference_internal)
        .def("amount", &ChemicalProps::amount)
        .def("mass", &ChemicalProps::mass)
        .def("volume", &ChemicalProps::volume)
        .def("gibbsEnergy", &ChemicalProps::gibbsEnergy)
        .def("enthalpy", &ChemicalProps::enthalpy)
        .def("entropy", &ChemicalProps::entropy)
        .def("internalEnergy", &ChemicalProps::internalEnergy)
        .def("helmholtzEnergy", &ChemicalProps::helmholtzEnergy)
        ;
}
