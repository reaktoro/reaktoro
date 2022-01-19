// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
using namespace Reaktoro;

void exportChemicalSystem(py::module& m)
{
    auto createChemicalSystem = [](const Database& db, py::args gphases)
    {
        Phases phases(db);
        for(auto phase : gphases)
        {
            try { phases.add(phase.cast<const GenericPhase&>()); }
            catch(...)
            {
                try { phases.add(phase.cast<const GenericPhasesGenerator&>()); }
                catch(...)
                {
                    errorif(true, "Could not create ChemicalSystem with phase object:\n", py::str(phase));
                }
            }
        }

        return ChemicalSystem(phases);
    };

    py::class_<ChemicalSystem>(m, "ChemicalSystem")
        .def(py::init<>())
        .def(py::init<const Phases&>())
        .def(py::init<const Database&, const Vec<Phase>&>())
        .def(py::init(createChemicalSystem))
        .def("database", &ChemicalSystem::database)
        .def("element", &ChemicalSystem::element)
        .def("elements", &ChemicalSystem::elements)
        .def("species", py::overload_cast<>(&ChemicalSystem::species, py::const_), py::return_value_policy::reference_internal)
        .def("species", py::overload_cast<Index>(&ChemicalSystem::species, py::const_), py::return_value_policy::reference_internal)
        .def("phase", &ChemicalSystem::phase)
        .def("phases", &ChemicalSystem::phases)
        .def("formulaMatrix", &ChemicalSystem::formulaMatrix, py::return_value_policy::reference_internal)
        .def("formulaMatrixElements", &ChemicalSystem::formulaMatrixElements, py::return_value_policy::reference_internal)
        .def("formulaMatrixCharge", &ChemicalSystem::formulaMatrixCharge, py::return_value_policy::reference_internal)
        ;
}
