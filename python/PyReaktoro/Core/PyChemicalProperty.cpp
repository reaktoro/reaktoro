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

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProperty.hpp>
#include <Reaktoro/Core/ChemicalProperties.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Common/ReactionEquation.hpp>

namespace Reaktoro {


void exportChemicalProperty(py::module& m)
{
	py::class_<double>(m, "ChemicalProperty")
		.def_static("ionicStrength", &ChemicalProperty::ionicStrength)
		.def_static("pH", &ChemicalProperty::pH)
		.def_static("pE", py::overload_cast<const ChemicalSystem&>(&ChemicalProperty::pE))
		.def_static("Eh", py::overload_cast<const ChemicalSystem&>(&ChemicalProperty::Eh))
		.def_static("alkalinity", &ChemicalProperty::alkalinity);

//    BOOST_PYTHON_FUNCTION_OVERLOADS(pE_overloads, pE, 1, 2);
//    BOOST_PYTHON_FUNCTION_OVERLOADS(Eh_overloads, Eh, 1, 2);
//
//    py::def("ionicStrength", &ionicStrength);
//    py::def("pH", &pH);
//    py::def("pE", ChemicalPropertyFunction(*)(const ChemicalSystem&, const ReactionEquation&), pE_overloads);
//    py::def("Eh", ChemicalPropertyFunction(*)(const ChemicalSystem&, const ReactionEquation&), Eh_overloads);
//    py::def("alkalinity", &alkalinity);
}

} // namespace Reaktoro
