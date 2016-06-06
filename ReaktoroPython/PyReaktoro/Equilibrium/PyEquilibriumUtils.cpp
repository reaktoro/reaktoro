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

#include "PyEquilibriumSolver.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumUtils.hpp>

namespace Reaktoro {

auto export_EquilibriumUtils() -> void
{
    auto equilibrate1  = static_cast<EquilibriumResult (*)(EquilibriumState&)>(equilibrate);
    auto equilibrate2  = static_cast<EquilibriumResult (*)(EquilibriumState&, const Partition&)>(equilibrate);
    auto equilibrate3  = static_cast<EquilibriumResult (*)(EquilibriumState&, const EquilibriumOptions&)>(equilibrate);
    auto equilibrate4  = static_cast<EquilibriumResult (*)(EquilibriumState&, const Partition&, const EquilibriumOptions&)>(equilibrate);
    auto equilibrate5  = static_cast<EquilibriumResult (*)(EquilibriumState&, const EquilibriumProblem&)>(equilibrate);
    auto equilibrate6  = static_cast<EquilibriumResult (*)(EquilibriumState&, const EquilibriumProblem&, const EquilibriumOptions&)>(equilibrate);
    auto equilibrate7  = static_cast<EquilibriumState (*)(const EquilibriumProblem&)>(equilibrate);
    auto equilibrate8  = static_cast<EquilibriumState (*)(const EquilibriumProblem&, const EquilibriumOptions&)>(equilibrate);
    auto equilibrate9  = static_cast<EquilibriumResult (*)(EquilibriumState&, const EquilibriumInverseProblem&)>(equilibrate);
    auto equilibrate10 = static_cast<EquilibriumResult (*)(EquilibriumState&, const EquilibriumInverseProblem&, const EquilibriumOptions&)>(equilibrate);
    auto equilibrate11 = static_cast<EquilibriumState (*)(const EquilibriumInverseProblem&)>(equilibrate);
    auto equilibrate12 = static_cast<EquilibriumState (*)(const EquilibriumInverseProblem&, const EquilibriumOptions&)>(equilibrate);

    py::def("equilibrate", equilibrate1);
    py::def("equilibrate", equilibrate2);
    py::def("equilibrate", equilibrate3);
    py::def("equilibrate", equilibrate4);
    py::def("equilibrate", equilibrate5);
    py::def("equilibrate", equilibrate6);
    py::def("equilibrate", equilibrate7);
    py::def("equilibrate", equilibrate8);
    py::def("equilibrate", equilibrate9);
    py::def("equilibrate", equilibrate10);
    py::def("equilibrate", equilibrate11);
    py::def("equilibrate", equilibrate12);
}

} // namespace Reaktoro
