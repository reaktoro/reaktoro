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
#include <Reaktor/Interfaces/Gems.hpp>

namespace Reaktor {

auto export_Gems() -> void
{
    py::class_<GemsOptions>("GemsOptions")
        .def_readwrite("warmstart", &GemsOptions::warmstart)
        ;

    py::class_<Gems>("Gems")
        .def(py::init<>())
        .def(py::init<std::string>())
        .def("setTemperature", &Gems::setTemperature)
        .def("setPressure", &Gems::setPressure)
        .def("setSpeciesAmounts", &Gems::setSpeciesAmounts)
        .def("setElementAmounts", &Gems::setElementAmounts)
        .def("setOptions", &Gems::setOptions)
        .def("numElements", &Gems::numElements)
        .def("numSpecies", &Gems::numSpecies)
        .def("numPhases", &Gems::numPhases)
        .def("numSpeciesInPhase", &Gems::numSpeciesInPhase)
        .def("elementName", &Gems::elementName)
        .def("speciesName", &Gems::speciesName)
        .def("phaseName", &Gems::phaseName)
        .def("indexElement", &Gems::indexElement)
        .def("indexSpecies", &Gems::indexSpecies)
        .def("indexPhase", &Gems::indexPhase)
        .def("elementCoefficientInSpecies", &Gems::elementCoefficientInSpecies)
        .def("speciesCharge", &Gems::speciesCharge)
        .def("elementsInSpecies", &Gems::elementsInSpecies)
        .def("elementMolarMass", &Gems::elementMolarMass)
        .def("speciesMolarMass", &Gems::speciesMolarMass)
        .def("temperature", &Gems::temperature)
        .def("pressure", &Gems::pressure)
        .def("elementAmounts", &Gems::elementAmounts)
        .def("speciesAmounts", &Gems::speciesAmounts)
        .def("formulaMatrix", &Gems::formulaMatrix)
        .def("standardGibbsEnergies", &Gems::standardGibbsEnergies)
        .def("chemicalPotentials", &Gems::chemicalPotentials)
        .def("equilibrate", &Gems::equilibrate)
        .def("converged", &Gems::converged)
        .def("numIterations", &Gems::numIterations)
        .def("elapsedTime", &Gems::elapsedTime)
        ;
}

} // namespace Reaktor
