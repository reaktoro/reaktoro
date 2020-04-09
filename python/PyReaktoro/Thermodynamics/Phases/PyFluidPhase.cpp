// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Thermodynamics/Mixtures/FluidMixture.hpp>
#include <Reaktoro/Thermodynamics/Phases/FluidPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/LiquidPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/FluidSpecies.hpp>

namespace Reaktoro {

void exportFluidPhase(py::module& m)
{
    py::class_<FluidPhase, Phase>(m, "_FluidPhase")
        .def(py::init<const std::string&, PhaseType>())
        .def(py::init<const FluidMixture&, const std::string&, PhaseType>())
        .def("setChemicalModelIdeal", &FluidPhase::setChemicalModelIdeal, py::return_value_policy::reference_internal)
        .def("setChemicalModelVanDerWaals", &FluidPhase::setChemicalModelVanDerWaals, py::return_value_policy::reference_internal)
        .def("setChemicalModelRedlichKwong", &FluidPhase::setChemicalModelRedlichKwong, py::return_value_policy::reference_internal)
        .def("setChemicalModelSoaveRedlichKwong", &FluidPhase::setChemicalModelSoaveRedlichKwong, py::return_value_policy::reference_internal)
        .def("setChemicalModelPengRobinson", &FluidPhase::setChemicalModelPengRobinson,
            py::arg("params") = CubicEOS::Params{},
            py::return_value_policy::reference_internal)
        .def("setChemicalModelCubicEOS", &FluidPhase::setChemicalModelCubicEOS,
            py::arg("params") = CubicEOS::Params{},
            py::return_value_policy::reference_internal
        )
        .def("setChemicalModelSpycherPruessEnnis", &FluidPhase::setChemicalModelSpycherPruessEnnis, py::return_value_policy::reference_internal)
        .def("setChemicalModelSpycherReed", &FluidPhase::setChemicalModelSpycherReed, py::return_value_policy::reference_internal)
        .def("mixture", (FluidMixture& (FluidPhase::*)(void)) &FluidPhase::mixture, py::return_value_policy::reference_internal)
        ;

    py::class_<GaseousPhase, FluidPhase>(m, "GaseousPhase")
        .def(py::init<>())
        ;

    py::class_<LiquidPhase, FluidPhase>(m, "LiquidPhase")
        .def(py::init<>())
        ;
}

} // namespace Reaktoro