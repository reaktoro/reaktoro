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
        .def("phase", &AqueousProps::phase, py::return_value_policy::reference_internal)
        ;
}
