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

#include "PyEquilibriumSolver.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/Partition.hpp>
#include <Reaktoro/Equilibrium/EquilibriumOptions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumInverseProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
#include <Reaktoro/Equilibrium/EquilibriumUtils.hpp>

namespace Reaktoro {

auto export_EquilibriumUtils() -> void
{
    auto equilibrate1  = static_cast<EquilibriumResult (*)(ChemicalState&)>(equilibrate);
    auto equilibrate2  = static_cast<EquilibriumResult (*)(ChemicalState&, const Partition&)>(equilibrate);
    auto equilibrate3  = static_cast<EquilibriumResult (*)(ChemicalState&, const EquilibriumOptions&)>(equilibrate);
    auto equilibrate4  = static_cast<EquilibriumResult (*)(ChemicalState&, const Partition&, const EquilibriumOptions&)>(equilibrate);
    auto equilibrate5  = static_cast<EquilibriumResult (*)(ChemicalState&, const EquilibriumProblem&)>(equilibrate);
    auto equilibrate6  = static_cast<EquilibriumResult (*)(ChemicalState&, const EquilibriumProblem&, const EquilibriumOptions&)>(equilibrate);
    auto equilibrate7  = static_cast<ChemicalState (*)(const EquilibriumProblem&)>(equilibrate);
    auto equilibrate8  = static_cast<ChemicalState (*)(const EquilibriumProblem&, const EquilibriumOptions&)>(equilibrate);
    auto equilibrate9  = static_cast<EquilibriumResult (*)(ChemicalState&, const EquilibriumInverseProblem&)>(equilibrate);
    auto equilibrate10 = static_cast<EquilibriumResult (*)(ChemicalState&, const EquilibriumInverseProblem&, const EquilibriumOptions&)>(equilibrate);
    auto equilibrate11 = static_cast<ChemicalState (*)(const EquilibriumInverseProblem&)>(equilibrate);
    auto equilibrate12 = static_cast<ChemicalState (*)(const EquilibriumInverseProblem&, const EquilibriumOptions&)>(equilibrate);

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
