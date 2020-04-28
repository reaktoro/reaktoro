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
#include <Reaktoro/Common/ReactionEquation.hpp>
#include <Reaktoro/Extensions/Interfaces/Phreeqc.hpp>

namespace Reaktoro {

void exportPhreeqc(py::module& m)
{
	auto execute1 = static_cast<void(Phreeqc::*)(std::string,std::string)>(&Phreeqc::execute);
	auto execute2 = static_cast<void(Phreeqc::*)(std::string)>(&Phreeqc::execute);

    py::class_<Phreeqc, Interface>(m, "Phreeqc")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("load", &Phreeqc::load)
        .def("execute", execute1)
        .def("execute", execute2)
        .def("reset", &Phreeqc::reset)
        .def("reactions", &Phreeqc::reactions)
        .def("stoichiometricMatrix", &Phreeqc::stoichiometricMatrix)
        .def("standardMolarGibbsEnergies", &Phreeqc::standardMolarGibbsEnergies)
        .def("standardMolarEnthalpies", &Phreeqc::standardMolarEnthalpies)
        .def("standardMolarVolumes", &Phreeqc::standardMolarVolumes)
        .def("standardMolarHeatCapacitiesConstP", &Phreeqc::standardMolarHeatCapacitiesConstP)
        .def("standardMolarHeatCapacitiesConstV", &Phreeqc::standardMolarHeatCapacitiesConstV)
        .def("lnActivityCoefficients", &Phreeqc::lnActivityCoefficients)
        .def("lnActivityConstants", &Phreeqc::lnActivityConstants)
        .def("lnActivities", &Phreeqc::lnActivities)
        .def("lnEquilibriumConstants", &Phreeqc::lnEquilibriumConstants)
        .def("phaseMolarVolumes", &Phreeqc::phaseMolarVolumes)
        ;
}

} // namespace Reaktoro
