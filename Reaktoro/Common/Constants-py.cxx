// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Constants.hpp>
using namespace Reaktoro;

void exportConstants(py::module& m)
{
    m.attr("universalGasConstant") = universalGasConstant;
    m.attr("faradayConstant") = faradayConstant;
    m.attr("jouleToCalorie") = jouleToCalorie;
    m.attr("calorieToJoule") = calorieToJoule;
    m.attr("barToPascal") = barToPascal;
    m.attr("cubicCentimeterToCubicMeter") = cubicCentimeterToCubicMeter;
    m.attr("cubicMeterToCubicCentimeter") = cubicMeterToCubicCentimeter;
    m.attr("ln10") = ln10;
}
