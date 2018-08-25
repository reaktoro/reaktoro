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
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProperty.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Common/ReactionEquation.hpp>

namespace Reaktoro {

void exportChemicalProperty(py::module& m)
{
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
