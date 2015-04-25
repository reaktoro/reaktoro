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

#include "PyMineralPhase.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Mixtures/MineralMixture.hpp>
#include <Reaktoro/Thermodynamics/Phases/MineralPhase.hpp>
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

namespace Reaktoro {

auto export_MineralPhase() -> void
{
    py::class_<MineralPhase>("MineralPhase")
        .def(py::init<>())
        .def(py::init<const std::vector<MineralSpecies>&>())
        .def("setActivityModel", &MineralPhase::setActivityModel)
        .def("setActivityModelIdeal", &MineralPhase::setActivityModelIdeal)
        .def("concentrations", &MineralPhase::concentrations)
        .def("activityCoefficients", &MineralPhase::activityCoefficients)
        .def("activities", &MineralPhase::activities)
        ;
}

} // namespace Reaktoro
