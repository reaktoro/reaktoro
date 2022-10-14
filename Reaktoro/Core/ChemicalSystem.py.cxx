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

namespace rkt4py {

/// Creates a Phases object with given database and list of arguments expected to be Phase-like objects.
extern auto createPhases(Database const& db, py::args gphases) -> Phases;

/// Creates a ChemicalSystem object with given database and list of arguments expected to be Phase-like objects.
auto createChemicalSystem(Database const& db, py::args gphases) -> ChemicalSystem
{
    return ChemicalSystem(createPhases(db, gphases));
};

} // namespace rkt4py

void exportChemicalSystem(py::module& m)
{

    py::class_<ChemicalSystem>(m, "ChemicalSystem")
        .def(py::init<>())
        .def(py::init<Database const&, PhaseList const&>())
        .def(py::init<Database const&, PhaseList const&, SurfaceList const&>())
        .def(py::init<Database const&, PhaseList const&, ReactionList const&>())
        .def(py::init<Database const&, PhaseList const&, ReactionList const&, SurfaceList const&>())
        .def(py::init<Phases const&>())
        .def(py::init<Phases const&, Surfaces const&>())
        .def(py::init<Phases const&, Reactions const&>())
        .def(py::init<Phases const&, Reactions const&, Surfaces const&>())
        .def(py::init(&rkt4py::createChemicalSystem))
        .def("id", &ChemicalSystem::id)
        .def("database", &ChemicalSystem::database, return_internal_ref)
        .def("element", &ChemicalSystem::element, return_internal_ref)
        .def("elements", &ChemicalSystem::elements, return_internal_ref)
        .def("species", py::overload_cast<>(&ChemicalSystem::species, py::const_), return_internal_ref)
        .def("species", py::overload_cast<Index>(&ChemicalSystem::species, py::const_), return_internal_ref)
        .def("phase", &ChemicalSystem::phase, return_internal_ref)
        .def("phases", &ChemicalSystem::phases, return_internal_ref)
        .def("reaction", &ChemicalSystem::reaction, return_internal_ref)
        .def("reactions", &ChemicalSystem::reactions, return_internal_ref)
        .def("surface", &ChemicalSystem::surface, return_internal_ref)
        .def("surfaces", &ChemicalSystem::surfaces, return_internal_ref)
        .def("formulaMatrix", &ChemicalSystem::formulaMatrix, return_internal_ref)
        .def("formulaMatrixElements", &ChemicalSystem::formulaMatrixElements, return_internal_ref)
        .def("formulaMatrixCharge", &ChemicalSystem::formulaMatrixCharge, return_internal_ref)
        .def("stoichiometricMatrix", &ChemicalSystem::stoichiometricMatrix, return_internal_ref)
        ;
}
