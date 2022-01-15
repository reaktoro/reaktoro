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

// Reaktoro includes
#include <Reaktoro/Core/StateOfMatter.hpp>
using namespace Reaktoro;

void exportStateOfMatter(py::module& m)
{
    py::enum_<StateOfMatter>(m, "StateOfMatter")
        .value("Unspecified"   , StateOfMatter::Unspecified   , "When the state of matter of a phase is unspecified")
        .value("Solid"         , StateOfMatter::Solid         , "When the state of matter of a phase is solid")
        .value("Liquid"        , StateOfMatter::Liquid        , "When the state of matter of a phase is liquid")
        .value("Gas"           , StateOfMatter::Gas           , "When the state of matter of a phase is gas")
        .value("Supercritical" , StateOfMatter::Supercritical , "When the state of matter of a phase is supercritical")
        .value("Plasma"        , StateOfMatter::Plasma        , "When the state of matter of a phase is plama")
        .value("Fluid"         , StateOfMatter::Fluid         , "When the state of matter of a phase is either liquid, gas, or plasma")
        .value("Condensed"     , StateOfMatter::Condensed     , "When the state of matter of a phase is either liquid or solid")
        ;
}
