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

#include "PyInterface.hpp"

// Boost includes
#include <boost/python.hpp>
namespace py = boost::python;

// Reaktoro includes
#include <Reaktoro/Interfaces/Interface.hpp>

namespace Reaktoro {

struct InterfaceWrapper : Interface, py::wrapper<Interface>
{
    auto set(double T, double P) -> void
    {
        this->get_override("set")(T, P);
    }

    auto set(double T, double P, const Vector& n) -> void
    {
        this->get_override("set")(T, P, n);
    }

    auto temperature() const -> double
    {
        return this->get_override("temperature")();
    }

    auto pressure() const -> double
    {
        return this->get_override("pressure")();
    }

    auto speciesAmounts() const -> Vector
    {
        return this->get_override("speciesAmounts")();
    }

    auto numElements() const -> unsigned
    {
        return this->get_override("numElements")();
    }

    auto numSpecies() const -> unsigned
    {
        return this->get_override("numSpecies")();
    }

    auto numPhases() const -> unsigned
    {
        return this->get_override("numPhases")();
    }

    auto numSpeciesInPhase(Index iphase) const -> unsigned
    {
        return this->get_override("numSpeciesInPhase")(iphase);
    }

    auto elementName(Index ielement) const -> std::string
    {
        return this->get_override("elementName")(ielement);
    }

    auto elementMolarMass(Index ielement) const -> double
    {
        return this->get_override("elementMolarMass")(ielement);
    }

    auto elementStoichiometry(Index ispecies, Index ielement) const -> double
    {
        return this->get_override("elementStoichiometry")(ispecies, ielement);
    }

    auto speciesName(Index ispecies) const -> std::string
    {
        return this->get_override("speciesName")(ispecies);
    }

    auto phaseName(Index iphase) const -> std::string
    {
        return this->get_override("phaseName")(iphase);
    }

    auto phaseReferenceState(Index iphase) const -> std::string
    {
        return this->get_override("phaseReferenceState")(iphase);
    }

    auto standardMolarGibbsEnergies() const -> Vector
    {
        return this->get_override("standardMolarGibbsEnergies")();
    }

    auto standardMolarEnthalpies() const -> Vector
    {
        return this->get_override("standardMolarEnthalpies")();
    }

    auto standardMolarVolumes() const -> Vector
    {
        return this->get_override("standardMolarVolumes")();
    }

    auto standardMolarHeatCapacitiesConstP() const -> Vector
    {
        return this->get_override("standardMolarHeatCapacitiesConstP")();
    }

    auto standardMolarHeatCapacitiesConstV() const -> Vector
    {
        return this->get_override("standardMolarHeatCapacitiesConstV")();
    }

    auto lnActivityCoefficients() const -> Vector
    {
        return this->get_override("lnActivityCoefficients")();
    }

    auto lnActivities() const -> Vector
    {
        return this->get_override("lnActivities")();
    }

    auto phaseMolarVolumes() const -> Vector
    {
        if(py::override f = this->get_override("phaseMolarVolumes"))
            return f();
        return Interface::phaseMolarVolumes();
    }

    auto phaseResidualMolarGibbsEnergies() const -> Vector
    {
        if(py::override f = this->get_override("phaseResidualMolarGibbsEnergies"))
            return f();
        return Interface::phaseResidualMolarGibbsEnergies();
    }

    auto phaseResidualMolarEnthalpies() const -> Vector
    {
        if(py::override f = this->get_override("phaseResidualMolarEnthalpies"))
            return f();
        return Interface::phaseResidualMolarEnthalpies();
    }

    auto phaseResidualMolarHeatCapacitiesConstP() const -> Vector
    {
        if(py::override f = this->get_override("phaseResidualMolarHeatCapacitiesConstP"))
            return f();
        return Interface::phaseResidualMolarHeatCapacitiesConstP();
    }

    auto phaseResidualMolarHeatCapacitiesConstV() const -> Vector
    {
        if(py::override f = this->get_override("phaseResidualMolarHeatCapacitiesConstV"))
            return f();
        return Interface::phaseResidualMolarHeatCapacitiesConstV();
    }

    auto phaseMolarVolumesDefault() const -> Vector
    {
        return Interface::phaseMolarVolumes();
    }

    auto phaseResidualMolarGibbsEnergiesDefault() const -> Vector
    {
        return Interface::phaseResidualMolarGibbsEnergies();
    }

    auto phaseResidualMolarEnthalpiesDefault() const -> Vector
    {
        return Interface::phaseResidualMolarEnthalpies();
    }

    auto phaseResidualMolarHeatCapacitiesConstPDefault() const -> Vector
    {
        return Interface::phaseResidualMolarHeatCapacitiesConstP();
    }

    auto phaseResidualMolarHeatCapacitiesConstVDefault() const -> Vector
    {
        return Interface::phaseResidualMolarHeatCapacitiesConstV();
    }
};

auto export_Interface() -> void
{
    auto set1 = static_cast<void (InterfaceWrapper::*)(double,double)>(&InterfaceWrapper::set);
    auto set2 = static_cast<void (InterfaceWrapper::*)(double,double,const Vector&)>(&InterfaceWrapper::set);

    py::class_<InterfaceWrapper, boost::noncopyable>("Interface")
        .def("set", py::pure_virtual(set1))
        .def("set", py::pure_virtual(set2))
        .def("temperature", py::pure_virtual(&Interface::temperature))
        .def("pressure", py::pure_virtual(&Interface::pressure))
        .def("speciesAmounts", py::pure_virtual(&Interface::speciesAmounts))
        .def("numElements", py::pure_virtual(&Interface::numElements))
        .def("numSpecies", py::pure_virtual(&Interface::numSpecies))
        .def("numPhases", py::pure_virtual(&Interface::numPhases))
        .def("numSpeciesInPhase", py::pure_virtual(&Interface::numSpeciesInPhase))
        .def("elementName", py::pure_virtual(&Interface::elementName))
        .def("elementMolarMass", py::pure_virtual(&Interface::elementMolarMass))
        .def("elementStoichiometry", py::pure_virtual(&Interface::elementStoichiometry))
        .def("speciesName", py::pure_virtual(&Interface::speciesName))
        .def("phaseName", py::pure_virtual(&Interface::phaseName))
        .def("phaseReferenceState", py::pure_virtual(&Interface::phaseReferenceState))
        .def("standardMolarGibbsEnergies", py::pure_virtual(&Interface::standardMolarGibbsEnergies))
        .def("standardMolarEnthalpies", py::pure_virtual(&Interface::standardMolarEnthalpies))
        .def("standardMolarVolumes", py::pure_virtual(&Interface::standardMolarVolumes))
        .def("standardMolarHeatCapacitiesConstP", py::pure_virtual(&Interface::standardMolarHeatCapacitiesConstP))
        .def("standardMolarHeatCapacitiesConstV", py::pure_virtual(&Interface::standardMolarHeatCapacitiesConstV))
        .def("lnActivityCoefficients", py::pure_virtual(&Interface::lnActivityCoefficients))
        .def("lnActivities", py::pure_virtual(&Interface::lnActivities))
        .def("phaseMolarVolumes", &Interface::phaseMolarVolumes, &InterfaceWrapper::phaseMolarVolumesDefault)
        .def("phaseResidualMolarGibbsEnergies", &Interface::phaseResidualMolarGibbsEnergies, &InterfaceWrapper::phaseResidualMolarGibbsEnergiesDefault)
        .def("phaseResidualMolarEnthalpies", &Interface::phaseResidualMolarEnthalpies, &InterfaceWrapper::phaseResidualMolarEnthalpiesDefault)
        .def("phaseResidualMolarHeatCapacitiesConstP", &Interface::phaseResidualMolarHeatCapacitiesConstP, &InterfaceWrapper::phaseResidualMolarHeatCapacitiesConstPDefault)
        .def("phaseResidualMolarHeatCapacitiesConstV", &Interface::phaseResidualMolarHeatCapacitiesConstV, &InterfaceWrapper::phaseResidualMolarHeatCapacitiesConstVDefault)
        .def("formulaMatrix", &Interface::formulaMatrix)
        .def("indexElement", &Interface::indexElement)
        .def("indexSpecies", &Interface::indexSpecies)
        .def("indexPhase", &Interface::indexPhase)
        .def("indexPhaseWithSpecies", &Interface::indexPhaseWithSpecies)
        .def("indexFirstSpeciesInPhase", &Interface::indexFirstSpeciesInPhase)
        ;
}

} // namespace Reaktoro
