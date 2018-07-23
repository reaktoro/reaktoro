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

#include "PyMineralCatalyst.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Thermodynamics/Reactions/MineralCatalyst.hpp>

// PyReaktoro includes
#include <PyReaktoro/Common/PyConverters.hpp>

namespace Reaktoro {

auto export_MineralCatalyst() -> void
{
    py::class_<MineralCatalyst>("MineralCatalyst")
        .def(py::init<>())
        .def(py::init<std::string, std::string, double>())
        .def(py::init<std::string>())
        .def_readwrite("species", &MineralCatalyst::species)
        .def_readwrite("quantity", &MineralCatalyst::quantity)
        .def_readwrite("power", &MineralCatalyst::power)
        ;

    export_std_vector<MineralCatalyst>("MineralCatalystVector");
}

} // namespace Reaktoro
