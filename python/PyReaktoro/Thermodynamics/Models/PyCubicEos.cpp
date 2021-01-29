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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/FluidMixture.hpp>
#include <Reaktoro/Thermodynamics/Phases/FluidPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/GaseousPhase.hpp>
#include <Reaktoro/Thermodynamics/Phases/LiquidPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/FluidSpecies.hpp>

namespace Reaktoro {

void exportCubicEos(py::module& m)
{
    py::enum_<CubicEOS::Model>(m, "CubicEOSModel")
        .value("None", CubicEOS::VanDerWaals)
        .value("RedlichKwong", CubicEOS::RedlichKwong)
        .value("SoaveRedlichKwong", CubicEOS::SoaveRedlichKwong)
        .value("PengRobinson", CubicEOS::PengRobinson)
        ;

    py::enum_<PhaseIdentificationMethod>(m, "PhaseIdentificationMethod")
        .value("None", PhaseIdentificationMethod::None)
        .value("VolumeMethod", PhaseIdentificationMethod::VolumeMethod)
        .value("IsothermalCompressibilityMethods", PhaseIdentificationMethod::IsothermalCompressibilityMethods)
        .value("GibbsEnergyAndEquationOfStateMethod", PhaseIdentificationMethod::GibbsEnergyAndEquationOfStateMethod)
        ;

    py::class_<CubicEOS::Params>(m, "CubicEOSParams")
        .def(
            py::init<CubicEOS::Model, PhaseIdentificationMethod>(),
            py::arg("model") = CubicEOS::PengRobinson,
            py::arg("phase_identification_method") = PhaseIdentificationMethod::None
        )
        .def_readwrite("model", &CubicEOS::Params::model)
        .def_readwrite("phase_identification_method", &CubicEOS::Params::phase_identification_method)
        ;
}

} // namespace Reaktoro