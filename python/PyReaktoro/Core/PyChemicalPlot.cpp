// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

namespace Reaktoro {

void exportChemicalPlot(py::module& m)
{
    auto y1 = static_cast<void(ChemicalPlot::*)(std::string)>(&ChemicalPlot::y);
    auto y2 = static_cast<void(ChemicalPlot::*)(std::string,std::string)>(&ChemicalPlot::y);
    auto points1 = static_cast<void(ChemicalPlot::*)(std::string, std::vector<double>, std::vector<double>)>(&ChemicalPlot::points);
    auto points2 = static_cast<void(ChemicalPlot::*)(std::string, std::string, std::string)>(&ChemicalPlot::points);
    auto showlegend1 = static_cast<void(ChemicalPlot::*)(bool)>(&ChemicalPlot::showlegend);
    auto showlegend2 = static_cast<bool(ChemicalPlot::*)() const>(&ChemicalPlot::showlegend);
    auto lshift = static_cast<ChemicalPlot&(ChemicalPlot::*)(std::string)>(&ChemicalPlot::operator<<);

    py::class_<ChemicalPlot>(m, "ChemicalPlot")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ReactionSystem&>())
        .def("name", &ChemicalPlot::name)
        .def("x", &ChemicalPlot::x)
        .def("y", y1)
        .def("y", y2)
        .def("points", points1)
        .def("points", points2)
        .def("legend", &ChemicalPlot::legend)
        .def("showlegend", showlegend1)
        .def("showlegeng", showlegend2)
        .def("title", &ChemicalPlot::title)
        .def("xlabel", &ChemicalPlot::xlabel)
        .def("ylabel", &ChemicalPlot::ylabel)
        .def("xtics", &ChemicalPlot::xtics)
        .def("ytics", &ChemicalPlot::ytics)
        .def("xformat", &ChemicalPlot::xformat)
        .def("yformat", &ChemicalPlot::yformat)
        .def("xlogscale", &ChemicalPlot::xlogscale, py::arg("base")=10)
        .def("ylogscale", &ChemicalPlot::ylogscale, py::arg("base")=10)
        .def("frequency", &ChemicalPlot::frequency)
        .def("__lshift__", lshift, py::return_value_policy::reference_internal)
        .def("open", &ChemicalPlot::open)
        .def("update", &ChemicalPlot::update)
        ;

//    exportstd_vector<ChemicalPlot>("ChemicalPlotVector");
}

} // namespace Reaktoro
