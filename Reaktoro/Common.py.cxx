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

void exportConstants(py::module& m);
void exportInterpolationUtils(py::module& m);
void exportMemoization(py::module& m);
void exportParseUtils(py::module& m);
void exportStringList(py::module& m);
void exportStringUtils(py::module& m);
void exportTable(py::module& m);
void exportTimeUtils(py::module& m);
void exportTypes(py::module& m);
void exportUnits(py::module& m);
void exportWarnings(py::module& m);

void exportCommon(py::module& m)
{
    exportConstants(m);
    exportInterpolationUtils(m);
    exportMemoization(m);
    exportParseUtils(m);
    exportStringList(m);
    exportStringUtils(m);
    exportTable(m);
    exportTimeUtils(m);
    exportTypes(m);
    exportUnits(m);
    exportWarnings(m);
}
