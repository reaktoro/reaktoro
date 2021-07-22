// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

void exportEquilibriumConditions(py::module& m);
void exportEquilibriumDims(py::module& m);
void exportEquilibriumOptions(py::module& m);
void exportEquilibriumProblem(py::module& m);
void exportEquilibriumRestrictions(py::module& m);
void exportEquilibriumResult(py::module& m);
void exportEquilibriumSolver(py::module& m);
void exportEquilibriumSpecs(py::module& m);
void exportEquilibriumUtils(py::module& m);

void exportEquilibrium(py::module& m)
{
    exportEquilibriumConditions(m);
    exportEquilibriumDims(m);
    exportEquilibriumOptions(m);
    exportEquilibriumRestrictions(m);
    exportEquilibriumProblem(m); // Ensure exportEquilibriumProblem is executed after exportEquilibriumConditions and exportEquilibriumRestrictions!
    exportEquilibriumResult(m);
    exportEquilibriumSolver(m);
    exportEquilibriumSpecs(m);
    exportEquilibriumUtils(m);
}
