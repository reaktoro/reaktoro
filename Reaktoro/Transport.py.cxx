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

void exportChemicalField(py::module& m);
void exportMesh(py::module& m);
void exportReactiveTransportProfiler(py::module& m);
void exportReactiveTransportResult(py::module& m);
void exportReactiveTransportAnalysis(py::module& m);
void exportReactiveTransportOptions(py::module& m);
void exportTransportSolver(py::module& m);

void exportTransport(py::module& m)
{
    exportChemicalField(m);
    exportMesh(m);
    exportReactiveTransportProfiler(m);
    exportReactiveTransportResult(m);
    exportReactiveTransportAnalysis(m);
    exportReactiveTransportOptions(m);
    exportTransportSolver(m);
}
