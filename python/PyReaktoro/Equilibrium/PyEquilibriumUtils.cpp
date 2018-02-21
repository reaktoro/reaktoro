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
#include <pybind11/eigen.h>
namespace py = pybind11;

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

void exportEquilibriumUtils(py::module& m)
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

    m.def("equilibrate", equilibrate1);
    m.def("equilibrate", equilibrate2);
    m.def("equilibrate", equilibrate3);
    m.def("equilibrate", equilibrate4);
    m.def("equilibrate", equilibrate5);
    m.def("equilibrate", equilibrate6);
    m.def("equilibrate", equilibrate7);
    m.def("equilibrate", equilibrate8);
    m.def("equilibrate", equilibrate9);
    m.def("equilibrate", equilibrate10);
    m.def("equilibrate", equilibrate11);
    m.def("equilibrate", equilibrate12);
}

} // namespace Reaktoro
