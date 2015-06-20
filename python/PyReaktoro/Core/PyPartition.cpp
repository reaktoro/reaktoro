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

#include "PyPartition.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>

namespace Reaktoro {

auto export_Partition() -> void
{
    using return_const_ref = py::return_value_policy<py::copy_const_reference>;

    auto setEquilibriumSpecies1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setEquilibriumSpecies);
    auto setEquilibriumSpecies2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setEquilibriumSpecies);
    auto setKineticSpecies1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setKineticSpecies);
    auto setKineticSpecies2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setKineticSpecies);
    auto setInertSpecies1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setInertSpecies);
    auto setInertSpecies2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setInertSpecies);

    py::class_<Partition>("Partition")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalSystem&, std::string>())
        .def("set", &Partition::set)
        .def("setEquilibriumSpecies", setEquilibriumSpecies1)
        .def("setEquilibriumSpecies", setEquilibriumSpecies2)
        .def("setEquilibriumPhases", &Partition::setEquilibriumPhases)
        .def("setKineticSpecies", setKineticSpecies1)
        .def("setKineticSpecies", setKineticSpecies2)
        .def("setKineticPhases", &Partition::setKineticPhases)
        .def("setInertSpecies", setInertSpecies1)
        .def("setInertSpecies", setInertSpecies2)
        .def("setInertPhases", &Partition::setInertPhases)
        .def("numEquilibriumSpecies", &Partition::numEquilibriumSpecies)
        .def("numKineticSpecies", &Partition::numKineticSpecies)
        .def("numInertSpecies", &Partition::numInertSpecies)
        .def("numEquilibriumElements", &Partition::numEquilibriumElements)
        .def("numKineticElements", &Partition::numKineticElements)
        .def("numInertElements", &Partition::numInertElements)
        .def("indicesEquilibriumSpecies", &Partition::indicesEquilibriumSpecies, return_const_ref())
        .def("indicesKineticSpecies", &Partition::indicesKineticSpecies, return_const_ref())
        .def("indicesInertSpecies", &Partition::indicesInertSpecies, return_const_ref())
        .def("indicesEquilibriumElements", &Partition::indicesEquilibriumElements, return_const_ref())
        .def("indicesKineticElements", &Partition::indicesKineticElements, return_const_ref())
        .def("indicesInertElements", &Partition::indicesInertElements, return_const_ref())
        .def("formulaMatrixEquilibriumSpecies", &Partition::formulaMatrixEquilibriumSpecies, return_const_ref())
        .def("formulaMatrixKineticSpecies", &Partition::formulaMatrixKineticSpecies, return_const_ref())
        .def("formulaMatrixInertSpecies", &Partition::formulaMatrixInertSpecies, return_const_ref())
        ;
}

} // namespace Reaktoro
