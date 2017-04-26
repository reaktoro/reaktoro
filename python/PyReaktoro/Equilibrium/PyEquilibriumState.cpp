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

#include "PyEquilibriumState.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumState.hpp>

namespace Reaktoro {
namespace {

auto assignEquilibriumState(EquilibriumState& state, const EquilibriumState& other) -> void
{
    state = other;
}

auto cloneEquilibriumState(EquilibriumState& state) -> EquilibriumState
{
    return state;
}

}  // namespace

auto export_EquilibriumState() -> void
{
    py::class_<EquilibriumState, py::bases<ChemicalState>>("EquilibriumState")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalState&>())
        .def("assign", assignEquilibriumState)
        .def("clone", cloneEquilibriumState)
        .def("setElementDualPotentials", &EquilibriumState::setElementDualPotentials)
        .def("setSpeciesDualPotentials", &EquilibriumState::setSpeciesDualPotentials)
        .def("elementDualPotentials", &EquilibriumState::elementDualPotentials, py::return_internal_reference<>())
        .def("speciesDualPotentials", &EquilibriumState::speciesDualPotentials, py::return_internal_reference<>())
        .def("phaseStabilityIndices", &EquilibriumState::phaseStabilityIndices)
        .def("output", &EquilibriumState::output)
        .def(py::self_ns::str(py::self_ns::self))
        ;
}

} // namespace Reaktoro
