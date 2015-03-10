// Reaktor is a C++ library for computational reaction modelling.
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

#include "PyPartition.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Partition.hpp>

namespace Reaktor {

auto export_Partition() -> void
{
    auto setEquilibriumSpecies1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setEquilibriumSpecies);
    auto setEquilibriumSpecies2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setEquilibriumSpecies);
    auto setKineticSpecies1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setKineticSpecies);
    auto setKineticSpecies2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setKineticSpecies);
    auto setInertSpecies1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setInertSpecies);
    auto setInertSpecies2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setInertSpecies);

    py::class_<Partition>("Partition", py::no_init)
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const Partition&>())
        .def("setEquilibriumSpecies", setEquilibriumSpecies1)
        .def("setEquilibriumSpecies", setEquilibriumSpecies2)
        .def("setEquilibriumPhases", &Partition::setEquilibriumPhases)
        .def("setKineticSpecies", setKineticSpecies1)
        .def("setKineticSpecies", setKineticSpecies2)
        .def("setKineticPhases", &Partition::setKineticPhases)
        .def("setInertSpecies", setInertSpecies1)
        .def("setInertSpecies", setInertSpecies2)
        .def("setInertPhases", &Partition::setInertPhases)
        .def("indicesEquilibriumSpecies", &Partition::indicesEquilibriumSpecies, py::return_value_policy<py::copy_const_reference>())
        .def("indicesKineticSpecies", &Partition::indicesKineticSpecies, py::return_value_policy<py::copy_const_reference>())
        .def("indicesInertSpecies", &Partition::indicesInertSpecies, py::return_value_policy<py::copy_const_reference>())
        ;
}

} // namespace Reaktor
