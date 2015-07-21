// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "PyGems.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Interfaces/Phreeqc.hpp>

namespace Reaktoro {

auto export_Phreeqc() -> void
{
    py::class_<Phreeqc>("Phreeqc")
        .def(py::init<>())
        .def(py::init<std::string, std::string>())
        .def("setTemperature", &Phreeqc::setTemperature)
        .def("setPressure", &Phreeqc::setPressure)
        .def("setSpeciesAmounts", &Phreeqc::setSpeciesAmounts)
        .def("numElements", &Phreeqc::numElements)
        .def("numSpecies", &Phreeqc::numSpecies)
        .def("numPhases", &Phreeqc::numPhases)
        .def("numSpeciesInPhase", &Phreeqc::numSpeciesInPhase)
        .def("elementName", &Phreeqc::elementName)
        .def("speciesName", &Phreeqc::speciesName)
        .def("phaseName", &Phreeqc::phaseName)
        .def("indexElement", &Phreeqc::indexElement)
        .def("indexSpecies", &Phreeqc::indexSpecies)
        .def("indexPhase", &Phreeqc::indexPhase)
        .def("indexPhaseWithSpecies", &Phreeqc::indexPhaseWithSpecies)
        .def("elementCoefficientInSpecies", &Phreeqc::elementCoefficientInSpecies)
        .def("speciesCharge", &Phreeqc::speciesCharge)
        .def("elementsInSpecies", &Phreeqc::elementsInSpecies)
        .def("elementMolarMass", &Phreeqc::elementMolarMass)
        .def("speciesMolarMass", &Phreeqc::speciesMolarMass)
        .def("temperature", &Phreeqc::temperature)
        .def("pressure", &Phreeqc::pressure)
        .def("speciesAmounts", &Phreeqc::speciesAmounts)
        .def("speciesAmount", &Phreeqc::speciesAmount)
        .def("speciesAmountsInPhase", &Phreeqc::speciesAmountsInPhase)
        .def("formulaMatrix", &Phreeqc::formulaMatrix)
        .def("standardGibbsEnergies", &Phreeqc::standardGibbsEnergies)
        .def("standardVolumes", &Phreeqc::standardVolumes)
        .def("activities", &Phreeqc::activities)
        .def("chemicalPotentials", &Phreeqc::chemicalPotentials)
        .def("phaseMolarVolumes", &Phreeqc::phaseMolarVolumes)
        ;
}

} // namespace Reaktoro
