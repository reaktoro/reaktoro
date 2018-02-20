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
// but WITHOUT ANY WARRANTY without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

// Reaktoro includes
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Interfaces/Interface.hpp>

namespace Reaktoro {

//struct InterfaceWrapper : Interface, py::wrapper<Interface>
//{
//    auto temperature() const -> double
//    {
//        return this->get_override("temperature")();
//    }
//
//    auto pressure() const -> double
//    {
//        return this->get_override("pressure")();
//    }
//
//    auto speciesAmounts() const -> Vector
//    {
//        return this->get_override("speciesAmounts")();
//    }
//
//    auto numElements() const -> unsigned
//    {
//        return this->get_override("numElements")();
//    }
//
//    auto numSpecies() const -> unsigned
//    {
//        return this->get_override("numSpecies")();
//    }
//
//    auto numPhases() const -> unsigned
//    {
//        return this->get_override("numPhases")();
//    }
//
//    auto numSpeciesInPhase(Index iphase) const -> unsigned
//    {
//        return this->get_override("numSpeciesInPhase")(iphase);
//    }
//
//    auto elementName(Index ielement) const -> std::string
//    {
//        return this->get_override("elementName")(ielement);
//    }
//
//    auto elementMolarMass(Index ielement) const -> double
//    {
//        return this->get_override("elementMolarMass")(ielement);
//    }
//
//    auto elementStoichiometry(Index ispecies, Index ielement) const -> double
//    {
//        return this->get_override("elementStoichiometry")(ispecies, ielement);
//    }
//
//    auto speciesName(Index ispecies) const -> std::string
//    {
//        return this->get_override("speciesName")(ispecies);
//    }
//
//    auto phaseName(Index iphase) const -> std::string
//    {
//        return this->get_override("phaseName")(iphase);
//    }
//
//    auto properties(PhaseThermoModelResult& res, Index iphase, double T, double P) -> void
//    {
//        this->get_override("properties")(res, iphase, T, P);
//    }
//
//    auto properties(PhaseChemicalModelResult& res, Index iphase, double T, double P, VectorConstRef n) -> void
//    {
//        this->get_override("properties")(res, iphase, T, P, n);
//    }
//
//    auto clone() const -> std::shared_ptr<Interface>
//    {
//        return this->get_override("clone")();
//    }
//};

void exportInterface(py::module& m)
{
//    auto properties1 = static_cast<void (InterfaceWrapper::*)(PhaseThermoModelResult&, Index, double, double)>(&InterfaceWrapper::properties);
//    auto properties2 = static_cast<void (InterfaceWrapper::*)(PhaseChemicalModelResult&, Index, double, double, VectorConstRef)>(&InterfaceWrapper::properties);
//
//    py::class_<InterfaceWrapper, boost::noncopyable>(m, "Interface")
//        .def("temperature", py::pure_virtual(&Interface::temperature))
//        .def("pressure", py::pure_virtual(&Interface::pressure))
//        .def("speciesAmounts", py::pure_virtual(&Interface::speciesAmounts))
//        .def("numElements", py::pure_virtual(&Interface::numElements))
//        .def("numSpecies", py::pure_virtual(&Interface::numSpecies))
//        .def("numPhases", py::pure_virtual(&Interface::numPhases))
//        .def("numSpeciesInPhase", py::pure_virtual(&Interface::numSpeciesInPhase))
//        .def("elementName", py::pure_virtual(&Interface::elementName))
//        .def("elementMolarMass", py::pure_virtual(&Interface::elementMolarMass))
//        .def("elementStoichiometry", py::pure_virtual(&Interface::elementStoichiometry))
//        .def("speciesName", py::pure_virtual(&Interface::speciesName))
//        .def("phaseName", py::pure_virtual(&Interface::phaseName))
//        .def("properties", py::pure_virtual(properties1))
//        .def("properties", py::pure_virtual(properties2))
//        .def("clone", py::pure_virtual(&Interface::clone))
//        .def("formulaMatrix", &Interface::formulaMatrix)
//        .def("indexElement", &Interface::indexElement)
//        .def("indexSpecies", &Interface::indexSpecies)
//        .def("indexPhase", &Interface::indexPhase)
//        .def("indexPhaseWithSpecies", &Interface::indexPhaseWithSpecies)
//        .def("indexFirstSpeciesInPhase", &Interface::indexFirstSpeciesInPhase)
//        .def("system", &Interface::system)
//        .def("state", &Interface::state)
//        ;
}

} // namespace Reaktoro
