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

#include "PyChemicalProperties.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalPropertiesAqueousPhase.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {

auto export_ChemicalPropertiesAqueousPhase() -> void
{
    auto pE1 = static_cast<ChemicalScalar (ChemicalPropertiesAqueousPhase::*)() const>(&ChemicalPropertiesAqueousPhase::pE);
    auto pE2 = static_cast<ChemicalScalar (ChemicalPropertiesAqueousPhase::*)(std::string) const>(&ChemicalPropertiesAqueousPhase::pE);

    auto Eh1 = static_cast<ChemicalScalar (ChemicalPropertiesAqueousPhase::*)() const>(&ChemicalPropertiesAqueousPhase::Eh);
    auto Eh2 = static_cast<ChemicalScalar (ChemicalPropertiesAqueousPhase::*)(std::string) const>(&ChemicalPropertiesAqueousPhase::Eh);

    py::class_<ChemicalPropertiesAqueousPhase>("ChemicalPropertiesAqueousPhase", py::no_init)
        .def(py::init<const ChemicalProperties&>())
        .def("ionicStrength", &ChemicalPropertiesAqueousPhase::ionicStrength)
        .def("pH", &ChemicalPropertiesAqueousPhase::pH)
        .def("pe", pE1)
        .def("pe", pE2)
        .def("Eh", Eh1)
        .def("Eh", Eh2)
        ;
}

} // namespace Reaktoro
