// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalPropsPhase.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousProps.hpp>
using namespace Reaktoro;

void exportAqueousProps(py::module& m)
{
    py::class_<AqueousProps>(m, "AqueousProps")
        .def(py::init<const ChemicalSystem&>())
        .def(py::init<const ChemicalState&>())
        .def(py::init<const ChemicalProps&>())
        .def("setActivityModel", &AqueousProps::setActivityModel)
        .def("update", py::overload_cast<const ChemicalState&>(&AqueousProps::update))
        .def("update", py::overload_cast<const ChemicalProps&>(&AqueousProps::update))
        .def("temperature", &AqueousProps::temperature)
        .def("pressure", &AqueousProps::pressure)
        .def("elementMolality", &AqueousProps::elementMolality)
        .def("elementMolalities", &AqueousProps::elementMolalities)
        .def("speciesMolality", &AqueousProps::speciesMolality)
        .def("speciesMolalities", &AqueousProps::speciesMolalities)
        .def("ionicStrength", &AqueousProps::ionicStrength)
        .def("ionicStrengthEffective", &AqueousProps::ionicStrengthEffective)
        .def("ionicStrengthStoichiometric", &AqueousProps::ionicStrengthStoichiometric)
        .def("pH", &AqueousProps::pH)
        .def("pE", &AqueousProps::pE)
        .def("Eh", &AqueousProps::Eh)
        .def("alkalinity", &AqueousProps::alkalinity)
        .def("saturationSpecies", &AqueousProps::saturationSpecies)
        .def("saturationIndex", &AqueousProps::saturationIndex)
        .def("saturationIndexLn", &AqueousProps::saturationIndexLn)
        .def("saturationIndexLg", &AqueousProps::saturationIndexLg)
        .def("saturationIndices", &AqueousProps::saturationIndices)
        .def("saturationIndicesLn", &AqueousProps::saturationIndicesLn)
        .def("saturationIndicesLg", &AqueousProps::saturationIndicesLg)
        .def("phase", &AqueousProps::phase, py::return_value_policy::reference_internal)
        .def("output", py::overload_cast<std::ostream&>(&AqueousProps::output, py::const_))
        .def("output", py::overload_cast<const String&>(&AqueousProps::output, py::const_))
        .def("__repr__", [](const AqueousProps& self) { std::stringstream ss; ss << self; return ss.str(); })
        ;
}
