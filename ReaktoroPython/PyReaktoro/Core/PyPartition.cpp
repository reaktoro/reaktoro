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
    auto setEquilibriumSpecies1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setEquilibriumSpecies);
    auto setEquilibriumSpecies2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setEquilibriumSpecies);

    auto setEquilibriumPhases1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setEquilibriumPhases);
    auto setEquilibriumPhases2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setEquilibriumPhases);

    auto setKineticSpecies1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setKineticSpecies);
    auto setKineticSpecies2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setKineticSpecies);

    auto setKineticPhases1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setKineticPhases);
    auto setKineticPhases2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setKineticPhases);

    auto setInertSpecies1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setInertSpecies);
    auto setInertSpecies2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setInertSpecies);

    auto setInertPhases1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setInertPhases);
    auto setInertPhases2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setInertPhases);

    py::class_<Partition>("Partition")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalSystem&, std::string>())
        .def("set", &Partition::set)
        .def("setEquilibriumSpecies", setEquilibriumSpecies1)
        .def("setEquilibriumSpecies", setEquilibriumSpecies2)
        .def("setEquilibriumPhases", setEquilibriumPhases1)
        .def("setEquilibriumPhases", setEquilibriumPhases2)
        .def("setKineticSpecies", setKineticSpecies1)
        .def("setKineticSpecies", setKineticSpecies2)
        .def("setKineticPhases", setKineticPhases1)
        .def("setKineticPhases", setKineticPhases2)
        .def("setInertSpecies", setInertSpecies1)
        .def("setInertSpecies", setInertSpecies2)
        .def("setInertPhases", setInertPhases1)
        .def("setInertPhases", setInertPhases2)
        .def("numEquilibriumSpecies", &Partition::numEquilibriumSpecies)
        .def("numKineticSpecies", &Partition::numKineticSpecies)
        .def("numInertSpecies", &Partition::numInertSpecies)
        .def("numEquilibriumElements", &Partition::numEquilibriumElements)
        .def("numKineticElements", &Partition::numKineticElements)
        .def("numInertElements", &Partition::numInertElements)
        .def("indicesEquilibriumSpecies", &Partition::indicesEquilibriumSpecies, py::return_internal_reference<>())
        .def("indicesKineticSpecies", &Partition::indicesKineticSpecies, py::return_internal_reference<>())
        .def("indicesInertSpecies", &Partition::indicesInertSpecies, py::return_internal_reference<>())
        .def("indicesEquilibriumElements", &Partition::indicesEquilibriumElements, py::return_internal_reference<>())
        .def("indicesKineticElements", &Partition::indicesKineticElements, py::return_internal_reference<>())
        .def("indicesInertElements", &Partition::indicesInertElements, py::return_internal_reference<>())
        .def("formulaMatrixEquilibriumSpecies", &Partition::formulaMatrixEquilibriumSpecies, py::return_internal_reference<>())
        .def("formulaMatrixKineticSpecies", &Partition::formulaMatrixKineticSpecies, py::return_internal_reference<>())
        .def("formulaMatrixInertSpecies", &Partition::formulaMatrixInertSpecies, py::return_internal_reference<>())
        ;
}

} // namespace Reaktoro
