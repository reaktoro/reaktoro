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

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Species/MineralSpecies.hpp>

// PyReator includes
#include <PyReaktoro/Common/PyConverters.hpp>

namespace Reaktoro {

auto export_MineralSpecies() -> void
{
    py::class_<MineralSpecies, py::bases<Species>>("MineralSpecies")
        .def(py::init<>())
        .def("setThermoData", &MineralSpecies::setThermoData)
        .def("thermoData", &MineralSpecies::thermoData, py::return_internal_reference<>())
        ;

    export_std_vector<MineralSpecies>("MineralSpeciesVector");
}

} // namespace Reaktoro
