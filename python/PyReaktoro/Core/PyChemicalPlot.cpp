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

#include "PyChemicalOutput.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Core/ChemicalPlot.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ReactionSystem.hpp>

// PyReator includes
#include <PyReaktoro/Common/PyConverters.hpp>

namespace Reaktoro {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(xlogscale_overloads, xlogscale, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ylogscale_overloads, ylogscale, 0, 1)

auto export_ChemicalPlot() -> void
{
    auto y1 = static_cast<void(ChemicalPlot::*)(std::string)>(&ChemicalPlot::y);
    auto y2 = static_cast<void(ChemicalPlot::*)(std::string,std::string)>(&ChemicalPlot::y);
    auto points1 = static_cast<void(ChemicalPlot::*)(std::string, std::vector<double>, std::vector<double>)>(&ChemicalPlot::points);
    auto points2 = static_cast<void(ChemicalPlot::*)(std::string, std::string, std::string)>(&ChemicalPlot::points);
    auto showlegend1 = static_cast<void(ChemicalPlot::*)(bool)>(&ChemicalPlot::showlegend);
    auto showlegend2 = static_cast<bool(ChemicalPlot::*)() const>(&ChemicalPlot::showlegend);
    auto lshift = static_cast<ChemicalPlot&(ChemicalPlot::*)(std::string)>(&ChemicalPlot::operator<<);

    py::class_<ChemicalPlot>("ChemicalPlot")
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
        .def("xlogscale", &ChemicalPlot::xlogscale)
        .def("ylogscale", &ChemicalPlot::ylogscale)
        .def("frequency", &ChemicalPlot::frequency)
        .def("__lshift__", lshift, py::return_internal_reference<>())
        .def("open", &ChemicalPlot::open)
        .def("update", &ChemicalPlot::update)
        ;

    export_std_vector<ChemicalPlot>("ChemicalPlotVector");
}

} // namespace Reaktoro
