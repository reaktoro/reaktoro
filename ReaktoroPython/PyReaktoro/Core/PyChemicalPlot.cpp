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
    auto legend1 = static_cast<void(ChemicalPlot::*)(bool)>(&ChemicalPlot::legend);
    auto legend2 = static_cast<bool(ChemicalPlot::*)() const>(&ChemicalPlot::legend);
    auto points1 = static_cast<void(ChemicalPlot::*)(std::string, std::vector<double>, std::vector<double>)>(&ChemicalPlot::points);
    auto points2 = static_cast<void(ChemicalPlot::*)(std::string, std::string, std::string)>(&ChemicalPlot::points);
    auto lshift = static_cast<ChemicalPlot&(ChemicalPlot::*)(std::string)>(&ChemicalPlot::operator<<);

    py::class_<ChemicalPlot>("ChemicalPlot")
        .def(py::init<>())
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ReactionSystem&>())
        .def("name", &ChemicalPlot::name)
        .def("x", &ChemicalPlot::x)
        .def("y", &ChemicalPlot::y)
        .def("points", points1)
        .def("points", points2)
        .def("legend", legend1)
        .def("legeng", legend2)
        .def("title", &ChemicalPlot::title)
        .def("xlabel", &ChemicalPlot::xlabel)
        .def("ylabel", &ChemicalPlot::ylabel)
        .def("xtics", &ChemicalPlot::xtics)
        .def("ytics", &ChemicalPlot::ytics)
        .def("xformat", &ChemicalPlot::xformat)
        .def("yformat", &ChemicalPlot::yformat)
        .def("xlogscale", &ChemicalPlot::xlogscale)
        .def("ylogscale", &ChemicalPlot::ylogscale)
        .def("key", &ChemicalPlot::key)
        .def("frequency", &ChemicalPlot::frequency)
        .def("__lshift__", lshift, py::return_internal_reference<>())
        .def("open", &ChemicalPlot::open)
        .def("update", &ChemicalPlot::update)
        ;

    export_std_vector<ChemicalPlot>("ChemicalPlotVector");
}

} // namespace Reaktoro
