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

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Partition.hpp>

namespace Reaktoro {

void exportPartition(py::module& m)
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

    py::class_<Partition>(m, "Partition")
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
        .def("system", &Partition::system, py::return_value_policy::reference_internal)
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
        .def("indicesFluidPhases", &Partition::indicesFluidPhases, py::return_value_policy::reference_internal)
        .def("indicesFluidSpecies", &Partition::indicesFluidSpecies, py::return_value_policy::reference_internal)
        .def("indicesSolidPhases", &Partition::indicesSolidPhases, py::return_value_policy::reference_internal)
        .def("indicesSolidSpecies", &Partition::indicesSolidSpecies, py::return_value_policy::reference_internal)
        .def("indicesEquilibriumSpecies", &Partition::indicesEquilibriumSpecies, py::return_value_policy::reference_internal)
        .def("indicesEquilibriumFluidSpecies", &Partition::indicesEquilibriumFluidSpecies, py::return_value_policy::reference_internal)
        .def("indicesEquilibriumSolidSpecies", &Partition::indicesEquilibriumSolidSpecies, py::return_value_policy::reference_internal)
        .def("indicesKineticSpecies", &Partition::indicesKineticSpecies, py::return_value_policy::reference_internal)
        .def("indicesKineticFluidSpecies", &Partition::indicesKineticFluidSpecies, py::return_value_policy::reference_internal)
        .def("indicesKineticSolidSpecies", &Partition::indicesKineticSolidSpecies, py::return_value_policy::reference_internal)
        .def("indicesInertSpecies", &Partition::indicesInertSpecies, py::return_value_policy::reference_internal)
        .def("indicesInertFluidSpecies", &Partition::indicesInertFluidSpecies, py::return_value_policy::reference_internal)
        .def("indicesInertSolidSpecies", &Partition::indicesInertSolidSpecies, py::return_value_policy::reference_internal)
        .def("indicesEquilibriumElements", &Partition::indicesEquilibriumElements, py::return_value_policy::reference_internal)
        .def("indicesEquilibriumFluidElements", &Partition::indicesEquilibriumFluidElements, py::return_value_policy::reference_internal)
        .def("indicesEquilibriumSolidElements", &Partition::indicesEquilibriumSolidElements, py::return_value_policy::reference_internal)
        .def("indicesKineticElements", &Partition::indicesKineticElements, py::return_value_policy::reference_internal)
        .def("indicesKineticFluidElements", &Partition::indicesKineticFluidElements, py::return_value_policy::reference_internal)
        .def("indicesKineticSolidElements", &Partition::indicesKineticSolidElements, py::return_value_policy::reference_internal)
        .def("indicesInertElements", &Partition::indicesInertElements, py::return_value_policy::reference_internal)
        .def("indicesInertFluidElements", &Partition::indicesInertFluidElements, py::return_value_policy::reference_internal)
        .def("indicesInertSolidElements", &Partition::indicesInertSolidElements, py::return_value_policy::reference_internal)
        .def("formulaMatrixEquilibriumPartition", &Partition::formulaMatrixEquilibriumPartition, py::return_value_policy::reference_internal)
        .def("formulaMatrixEquilibriumFluidPartition", &Partition::formulaMatrixEquilibriumFluidPartition, py::return_value_policy::reference_internal)
        .def("formulaMatrixEquilibriumSolidPartition", &Partition::formulaMatrixEquilibriumSolidPartition, py::return_value_policy::reference_internal)
        .def("formulaMatrixKineticPartition", &Partition::formulaMatrixKineticPartition, py::return_value_policy::reference_internal)
        .def("formulaMatrixKineticFluidPartition", &Partition::formulaMatrixKineticFluidPartition, py::return_value_policy::reference_internal)
        .def("formulaMatrixKineticSolidPartition", &Partition::formulaMatrixKineticSolidPartition, py::return_value_policy::reference_internal)
        .def("formulaMatrixInertPartition", &Partition::formulaMatrixInertPartition, py::return_value_policy::reference_internal)
        .def("formulaMatrixInertFluidPartition", &Partition::formulaMatrixInertFluidPartition, py::return_value_policy::reference_internal)
        .def("formulaMatrixInertSolidPartition", &Partition::formulaMatrixInertSolidPartition, py::return_value_policy::reference_internal)
        ;
}

} // namespace Reaktoro
