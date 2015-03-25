// Reaktor is a C++ library for computational reaction modelling.
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

// Reaktor includes
#include <Reaktor/Interfaces/Phreeqx.hpp>

namespace Reaktor {

auto export_Phreeqx() -> void
{
    py::class_<Phreeqx>("Phreeqx")
        .def(py::init<>())
        .def(py::init<std::string, std::string>())
        .def("setTemperature", &Phreeqx::setTemperature)
        .def("setPressure", &Phreeqx::setPressure)
        .def("setSpeciesAmounts", &Phreeqx::setSpeciesAmounts)
        .def("numElements", &Phreeqx::numElements)
        .def("numSpecies", &Phreeqx::numSpecies)
        .def("numPhases", &Phreeqx::numPhases)
        .def("numSpeciesInPhase", &Phreeqx::numSpeciesInPhase)
        .def("elementName", &Phreeqx::elementName)
        .def("speciesName", &Phreeqx::speciesName)
        .def("phaseName", &Phreeqx::phaseName)
        .def("indexElement", &Phreeqx::indexElement)
        .def("indexSpecies", &Phreeqx::indexSpecies)
        .def("indexPhase", &Phreeqx::indexPhase)
        .def("indexPhaseWithSpecies", &Phreeqx::indexPhaseWithSpecies)
        .def("elementAtomsInSpecies", &Phreeqx::elementAtomsInSpecies)
        .def("speciesCharge", &Phreeqx::speciesCharge)
        .def("elementsInSpecies", &Phreeqx::elementsInSpecies)
        .def("elementMolarMass", &Phreeqx::elementMolarMass)
        .def("speciesMolarMass", &Phreeqx::speciesMolarMass)
        .def("temperature", &Phreeqx::temperature)
        .def("pressure", &Phreeqx::pressure)
        .def("speciesAmounts", &Phreeqx::speciesAmounts)
        .def("speciesAmount", &Phreeqx::speciesAmount)
        .def("speciesAmountsInPhase", &Phreeqx::speciesAmountsInPhase)
        .def("formulaMatrix", &Phreeqx::formulaMatrix)
        .def("standardGibbsEnergies", &Phreeqx::standardGibbsEnergies)
        .def("standardVolumes", &Phreeqx::standardVolumes)
        .def("activities", &Phreeqx::activities)
        .def("chemicalPotentials", &Phreeqx::chemicalPotentials)
        .def("phaseMolarVolumes", &Phreeqx::phaseMolarVolumes)
        ;
}

} // namespace Reaktor
