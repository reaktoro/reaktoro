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

    auto setFluidPhases1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setFluidPhases);
    auto setFluidPhases2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setFluidPhases);

    auto setSolidPhases1 = static_cast<void(Partition::*)(const Indices&)>(&Partition::setSolidPhases);
    auto setSolidPhases2 = static_cast<void(Partition::*)(const std::vector<std::string>&)>(&Partition::setSolidPhases);

    py::class_<Partition>("Partition")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
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
        .def("setFluidPhases", setFluidPhases1)
        .def("setFluidPhases", setFluidPhases2)
        .def("setSolidPhases", setSolidPhases1)
        .def("setSolidPhases", setSolidPhases2)
        .def("system", &Partition::system, py::return_internal_reference<>())
        .def("numFluidPhases", &Partition::numFluidPhases)
        .def("numFluidSpecies", &Partition::numFluidSpecies)
        .def("numSolidPhases", &Partition::numSolidPhases)
        .def("numSolidSpecies", &Partition::numSolidSpecies)
        .def("numEquilibriumSpecies", &Partition::numEquilibriumSpecies)
        .def("numEquilibriumFluidSpecies", &Partition::numEquilibriumFluidSpecies)
        .def("numEquilibriumSolidSpecies", &Partition::numEquilibriumSolidSpecies)
        .def("numKineticSpecies", &Partition::numKineticSpecies)
        .def("numKineticFluidSpecies", &Partition::numKineticFluidSpecies)
        .def("numKineticSolidSpecies", &Partition::numKineticSolidSpecies)
        .def("numInertSpecies", &Partition::numInertSpecies)
        .def("numInertFluidSpecies", &Partition::numInertFluidSpecies)
        .def("numInertSolidSpecies", &Partition::numInertSolidSpecies)
        .def("numEquilibriumElements", &Partition::numEquilibriumElements)
        .def("numEquilibriumFluidElements", &Partition::numEquilibriumFluidElements)
        .def("numEquilibriumSolidElements", &Partition::numEquilibriumSolidElements)
        .def("numKineticElements", &Partition::numKineticElements)
        .def("numKineticFluidElements", &Partition::numKineticFluidElements)
        .def("numKineticSolidElements", &Partition::numKineticSolidElements)
        .def("numInertElements", &Partition::numInertElements)
        .def("numInertFluidElements", &Partition::numInertFluidElements)
        .def("numInertSolidElements", &Partition::numInertSolidElements)
        .def("indicesFluidPhases", &Partition::indicesFluidPhases, py::return_internal_reference<>())
        .def("indicesFluidSpecies", &Partition::indicesFluidSpecies, py::return_internal_reference<>())
        .def("indicesSolidPhases", &Partition::indicesSolidPhases, py::return_internal_reference<>())
        .def("indicesSolidSpecies", &Partition::indicesSolidSpecies, py::return_internal_reference<>())
        .def("indicesEquilibriumSpecies", &Partition::indicesEquilibriumSpecies, py::return_internal_reference<>())
        .def("indicesEquilibriumFluidSpecies", &Partition::indicesEquilibriumFluidSpecies, py::return_internal_reference<>())
        .def("indicesEquilibriumSolidSpecies", &Partition::indicesEquilibriumSolidSpecies, py::return_internal_reference<>())
        .def("indicesKineticSpecies", &Partition::indicesKineticSpecies, py::return_internal_reference<>())
        .def("indicesKineticFluidSpecies", &Partition::indicesKineticFluidSpecies, py::return_internal_reference<>())
        .def("indicesKineticSolidSpecies", &Partition::indicesKineticSolidSpecies, py::return_internal_reference<>())
        .def("indicesInertSpecies", &Partition::indicesInertSpecies, py::return_internal_reference<>())
        .def("indicesInertFluidSpecies", &Partition::indicesInertFluidSpecies, py::return_internal_reference<>())
        .def("indicesInertSolidSpecies", &Partition::indicesInertSolidSpecies, py::return_internal_reference<>())
        .def("indicesEquilibriumElements", &Partition::indicesEquilibriumElements, py::return_internal_reference<>())
        .def("indicesEquilibriumFluidElements", &Partition::indicesEquilibriumFluidElements, py::return_internal_reference<>())
        .def("indicesEquilibriumSolidElements", &Partition::indicesEquilibriumSolidElements, py::return_internal_reference<>())
        .def("indicesKineticElements", &Partition::indicesKineticElements, py::return_internal_reference<>())
        .def("indicesKineticFluidElements", &Partition::indicesKineticFluidElements, py::return_internal_reference<>())
        .def("indicesKineticSolidElements", &Partition::indicesKineticSolidElements, py::return_internal_reference<>())
        .def("indicesInertElements", &Partition::indicesInertElements, py::return_internal_reference<>())
        .def("indicesInertFluidElements", &Partition::indicesInertFluidElements, py::return_internal_reference<>())
        .def("indicesInertSolidElements", &Partition::indicesInertSolidElements, py::return_internal_reference<>())
        .def("formulaMatrixEquilibriumPartition", &Partition::formulaMatrixEquilibriumPartition, py::return_internal_reference<>())
        .def("formulaMatrixEquilibriumFluidPartition", &Partition::formulaMatrixEquilibriumFluidPartition, py::return_internal_reference<>())
        .def("formulaMatrixEquilibriumSolidPartition", &Partition::formulaMatrixEquilibriumSolidPartition, py::return_internal_reference<>())
        .def("formulaMatrixKineticPartition", &Partition::formulaMatrixKineticPartition, py::return_internal_reference<>())
        .def("formulaMatrixKineticFluidPartition", &Partition::formulaMatrixKineticFluidPartition, py::return_internal_reference<>())
        .def("formulaMatrixKineticSolidPartition", &Partition::formulaMatrixKineticSolidPartition, py::return_internal_reference<>())
        .def("formulaMatrixInertPartition", &Partition::formulaMatrixInertPartition, py::return_internal_reference<>())
        .def("formulaMatrixInertFluidPartition", &Partition::formulaMatrixInertFluidPartition, py::return_internal_reference<>())
        .def("formulaMatrixInertSolidPartition", &Partition::formulaMatrixInertSolidPartition, py::return_internal_reference<>())
        ;
}

} // namespace Reaktoro
