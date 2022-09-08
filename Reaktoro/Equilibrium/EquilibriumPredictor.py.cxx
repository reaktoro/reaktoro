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
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumPredictor.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSensitivity.hpp>
using namespace Reaktoro;

void exportEquilibriumProjector(py::module& m)
{
    py::class_<EquilibriumPredictor>(m, "EquilibriumPredictor")
        .def(py::init<ChemicalState const&, EquilibriumSensitivity const&>())
        .def("predict", py::overload_cast<ChemicalState&, EquilibriumConditions const&>(&EquilibriumPredictor::predict, py::const_), "Perform a first-order Taylor prediction of the chemical state at given conditions.")
        .def("predict", py::overload_cast<ChemicalState&, EquilibriumConditions const&, VectorXrConstRef const&>(&EquilibriumPredictor::predict, py::const_), "Perform a first-order Taylor prediction of the chemical state at given conditions.")
        .def("predictSpeciesChemicalPotential", &EquilibriumPredictor::predictSpeciesChemicalPotential, "Perform a first-order Taylor prediction of the chemical potential of a species at given conditions.")
        .def("referenceSpeciesChemicalPotential", &EquilibriumPredictor::referenceSpeciesChemicalPotential, "Return the chemical potential of a species at given reference conditions.")
        ;
}
