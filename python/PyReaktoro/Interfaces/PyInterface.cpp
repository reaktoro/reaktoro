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

#include <PyReaktoro/PyReaktoro.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Interfaces/Interface.hpp>

namespace Reaktoro {

class PyInterface : public Interface
{
public:
    auto temperature() const -> double
    {
        PYBIND11_OVERLOAD_PURE(double, Interface, temperature);
    }

    auto pressure() const -> double
    {
        PYBIND11_OVERLOAD_PURE(double, Interface, pressure);
    }

    auto speciesAmounts() const -> Vector
    {
        PYBIND11_OVERLOAD_PURE(Vector, Interface, speciesAmounts);
    }

    auto numElements() const -> unsigned
    {
        PYBIND11_OVERLOAD_PURE(unsigned, Interface, numElements);
    }

    auto numSpecies() const -> unsigned
    {
        PYBIND11_OVERLOAD_PURE(unsigned, Interface, numSpecies);
    }

    auto numPhases() const -> unsigned
    {
        PYBIND11_OVERLOAD_PURE(unsigned, Interface, numPhases);
    }

    auto numSpeciesInPhase(Index iphase) const -> unsigned
    {
        PYBIND11_OVERLOAD_PURE(unsigned, Interface, numSpeciesInPhase, iphase);
    }

    auto elementName(Index ielement) const -> std::string
    {
        PYBIND11_OVERLOAD_PURE(std::string, Interface, elementName, ielement);
    }

    auto elementMolarMass(Index ielement) const -> double
    {
        PYBIND11_OVERLOAD_PURE(double, Interface, elementMolarMass, ielement);
    }

    auto elementStoichiometry(Index ispecies, Index ielement) const -> double
    {
        PYBIND11_OVERLOAD_PURE(double, Interface, elementStoichiometry, ispecies, ielement);
    }

    auto speciesName(Index ispecies) const -> std::string
    {
        PYBIND11_OVERLOAD_PURE(std::string, Interface, speciesName, ispecies);
    }

    auto phaseName(Index iphase) const -> std::string
    {
        PYBIND11_OVERLOAD_PURE(std::string, Interface, phaseName, iphase);
    }

    auto properties(PhaseThermoModelResult& res, Index iphase, double T, double P) -> void
    {
        PYBIND11_OVERLOAD_PURE(void, Interface, properties, res, iphase, T, P);
    }

    auto properties(PhaseChemicalModelResult& res, Index iphase, double T, double P, VectorConstRef n) -> void
    {
        PYBIND11_OVERLOAD_PURE(void, Interface, properties, res, iphase, T, P, n);
    }

    auto clone() const -> std::shared_ptr<Interface>
    {
        PYBIND11_OVERLOAD_PURE(std::shared_ptr<Interface>, Interface, clone);
    }
};

void exportInterface(py::module& m)
{
    auto properties1 = static_cast<void (Interface::*)(ThermoModelResult&, double, double)>(&Interface::properties);
    auto properties2 = static_cast<void (Interface::*)(ChemicalModelResult&, double, double, VectorConstRef)>(&Interface::properties);

    py::class_<Interface, PyInterface>(m, "Interface")
        .def("temperature", &Interface::temperature)
        .def("pressure", &Interface::pressure)
        .def("speciesAmounts", &Interface::speciesAmounts)
        .def("numElements", &Interface::numElements)
        .def("numSpecies", &Interface::numSpecies)
        .def("numPhases", &Interface::numPhases)
        .def("numSpeciesInPhase", &Interface::numSpeciesInPhase)
        .def("elementName", &Interface::elementName)
        .def("elementMolarMass", &Interface::elementMolarMass)
        .def("elementStoichiometry", &Interface::elementStoichiometry)
        .def("speciesName", &Interface::speciesName)
        .def("phaseName", &Interface::phaseName)
        .def("properties", properties1)
        .def("properties", properties2)
        .def("clone", &Interface::clone)
        .def("formulaMatrix", &Interface::formulaMatrix)
        .def("indexElement", &Interface::indexElement)
        .def("indexSpecies", &Interface::indexSpecies)
        .def("indexPhase", &Interface::indexPhase)
        .def("indexPhaseWithSpecies", &Interface::indexPhaseWithSpecies)
        .def("indexFirstSpeciesInPhase", &Interface::indexFirstSpeciesInPhase)
        .def("system", &Interface::system)
        .def("state", &Interface::state)
        ;
}

} // namespace Reaktoro
