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

// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Common/ReactionEquation.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalProperty.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

class ChemicalPropertyDummy {};

void exportChemicalProperty(py::module& m)
{
    auto pE1 = static_cast<ChemicalPropertyFunction(*)(const ChemicalSystem&)>(ChemicalProperty::pE);
    auto pE2 = static_cast<ChemicalPropertyFunction(*)(const ChemicalSystem&, const ReactionEquation&)>(ChemicalProperty::pE);

    auto Eh1 = static_cast<ChemicalPropertyFunction(*)(const ChemicalSystem&)>(ChemicalProperty::Eh);
    auto Eh2 = static_cast<ChemicalPropertyFunction(*)(const ChemicalSystem&, const ReactionEquation&)>(ChemicalProperty::Eh);

    py::class_<ChemicalPropertyDummy>(m, "ChemicalProperty")
        .def_static("ionicStrength", &ChemicalProperty::ionicStrength)
        .def_static("pH", &ChemicalProperty::pH)
        .def_static("pE", pE1)
        .def_static("pE", pE2)
        .def_static("Eh", Eh1)
        .def_static("Eh", Eh2)
        .def_static("alkalinity", &ChemicalProperty::alkalinity)
        ;
}

} // namespace Reaktoro
