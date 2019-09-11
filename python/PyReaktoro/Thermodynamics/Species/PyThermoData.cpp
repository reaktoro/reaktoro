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

#include <Reaktoro/Thermodynamics/Species/ThermoData.hpp>

namespace Reaktoro {

//TODO: Fix for pybind bug, where it moves optional attributes, losing its dynamic content.
//https://github.com/pybind/pybind11/issues/1919
template <typename _Property>
_Property* get_optional(std::optional<_Property>& prop)
{
    if (prop.has_value())
        return &prop.value();
        
    return nullptr;
}

void exportThermoData(py::module& m)
{
    py::class_<SpeciesThermoData>(m, "SpeciesThermoData")
        .def(py::init<>())
        .def_property("properties", [] (SpeciesThermoData& self) { return get_optional(self.properties);}, 
                                    [] (SpeciesThermoData& self, SpeciesThermoInterpolatedProperties& property) {self.properties = property;} )
        .def_property("reaction", [] (SpeciesThermoData& self) { return get_optional(self.reaction);},
                                  [] (SpeciesThermoData& self, ReactionThermoInterpolatedProperties& property) {self.reaction = property;})
        .def_property("phreeqc", [] (SpeciesThermoData& self) { return get_optional(self.phreeqc);},
                                 [] (SpeciesThermoData& self, SpeciesThermoParamsPhreeqc& property) {self.phreeqc = property;})
        ;
    
    py::class_<AqueousSpeciesThermoData, SpeciesThermoData>(m, "AqueousSpeciesThermoData")
        .def(py::init<>())
        .def_property("hkf", [] (AqueousSpeciesThermoData& self) { return get_optional(self.hkf);},
                             [] (AqueousSpeciesThermoData& self, AqueousSpeciesThermoParamsHKF& hkf) {self.hkf = hkf;})
        ;

    py::class_<FluidSpeciesThermoData, SpeciesThermoData>(m, "FluidSpeciesThermoData")
        .def(py::init<>())
        .def_property("hkf", [] (FluidSpeciesThermoData& self) { return get_optional(self.hkf);},
                             [] (FluidSpeciesThermoData& self, FluidSpeciesThermoParamsHKF& hkf) {self.hkf = hkf;})
        ;

    m.attr("LiquidSpeciesThermoData") = m.attr("FluidSpeciesThermoData");
    m.attr("GaseousSpeciesThermoData") = m.attr("FluidSpeciesThermoData");
    
    {
        py::class_<MineralSpeciesThermoData, SpeciesThermoData> (m, "MineralSpeciesThermoData")
        .def(py::init<>())
        .def_property("hkf", [&] (MineralSpeciesThermoData& self) { return get_optional(self.hkf);},
                             [] (MineralSpeciesThermoData& self, MineralSpeciesThermoParamsHKF& hkf) {self.hkf = hkf;})
        ;
    }
    
}

void exportThermoDataProperties(py::module& m) 
{
    py::class_<SpeciesThermoInterpolatedProperties>(m, "SpeciesThermoInterpolatedProperties")
        .def(py::init<>())
        .def_readwrite("gibbs_energy", &SpeciesThermoInterpolatedProperties::gibbs_energy)
        .def_readwrite("helmholtz_energy", &SpeciesThermoInterpolatedProperties::helmholtz_energy)
        .def_readwrite("internal_energy", &SpeciesThermoInterpolatedProperties::internal_energy)
        .def_readwrite("enthalpy", &SpeciesThermoInterpolatedProperties::enthalpy)
        .def_readwrite("entropy", &SpeciesThermoInterpolatedProperties::entropy)
        .def_readwrite("volume", &SpeciesThermoInterpolatedProperties::volume)
        .def_readwrite("heat_capacity_cp", &SpeciesThermoInterpolatedProperties::heat_capacity_cp)
        .def_readwrite("heat_capacity_cv", &SpeciesThermoInterpolatedProperties::heat_capacity_cv)
        ;

    py::class_<ReactionThermoInterpolatedProperties, SpeciesThermoInterpolatedProperties>(m, "ReactionThermoInterpolatedProperties")
        .def(py::init<>())
        .def_readwrite("equation", &ReactionThermoInterpolatedProperties::equation)
        .def_readwrite("lnk", &ReactionThermoInterpolatedProperties::lnk)
        ;

    py::class_<SpeciesThermoParamsPhreeqc>(m, "SpeciesThermoParamsPhreeqc")
        .def(py::init<>())
        .def_readwrite("reaction", &SpeciesThermoParamsPhreeqc::reaction)
        ;

    py::class_<SpeciesThermoParamsPhreeqc::ReactionParams>(m, "ReactionParams")
        .def(py::init<>())
        .def_readwrite("equation", &SpeciesThermoParamsPhreeqc::ReactionParams::equation)
        .def_readwrite("log_k", &SpeciesThermoParamsPhreeqc::ReactionParams::log_k)
        .def_readwrite("delta_h", &SpeciesThermoParamsPhreeqc::ReactionParams::delta_h)
        .def_readwrite("analytic", &SpeciesThermoParamsPhreeqc::ReactionParams::analytic)
        ;

    py::class_<AqueousSpeciesThermoParamsHKF>(m, "AqueousSpeciesThermoParamsHKF")
        .def(py::init<>())
        .def_readwrite("Gf", &AqueousSpeciesThermoParamsHKF::Gf)
        .def_readwrite("Hf", &AqueousSpeciesThermoParamsHKF::Hf)
        .def_readwrite("Sr", &AqueousSpeciesThermoParamsHKF::Sr)
        .def_readwrite("a1", &AqueousSpeciesThermoParamsHKF::a1)
        .def_readwrite("a2", &AqueousSpeciesThermoParamsHKF::a2)
        .def_readwrite("a3", &AqueousSpeciesThermoParamsHKF::a3)
        .def_readwrite("a4", &AqueousSpeciesThermoParamsHKF::a4)
        .def_readwrite("c1", &AqueousSpeciesThermoParamsHKF::c1)
        .def_readwrite("c2", &AqueousSpeciesThermoParamsHKF::c2)
        .def_readwrite("wref", &AqueousSpeciesThermoParamsHKF::wref)
        ;

    py::class_<FluidSpeciesThermoParamsHKF>(m, "FluidSpeciesThermoParamsHKF")
        .def(py::init<>())
        .def_readwrite("Gf", &FluidSpeciesThermoParamsHKF::Gf)
        .def_readwrite("Hf", &FluidSpeciesThermoParamsHKF::Hf)
        .def_readwrite("Sr", &FluidSpeciesThermoParamsHKF::Sr)
        .def_readwrite("a", &FluidSpeciesThermoParamsHKF::a)
        .def_readwrite("b", &FluidSpeciesThermoParamsHKF::b)
        .def_readwrite("c", &FluidSpeciesThermoParamsHKF::c)
        .def_readwrite("Tmax", &FluidSpeciesThermoParamsHKF::Tmax)
        ;

    m.attr("LiquidSpeciesThermoParamsHKF") = m.attr("FluidSpeciesThermoParamsHKF");
    m.attr("GaseousSpeciesThermoParamsHKF") = m.attr("FluidSpeciesThermoParamsHKF");

    py::class_<MineralSpeciesThermoParamsHKF>(m, "MineralSpeciesThermoParamsHKF")
        .def(py::init<>())
        .def_readwrite("Gf", &MineralSpeciesThermoParamsHKF::Gf)
        .def_readwrite("Hf", &MineralSpeciesThermoParamsHKF::Hf)
        .def_readwrite("Sr", &MineralSpeciesThermoParamsHKF::Sr)
        .def_readwrite("Vr", &MineralSpeciesThermoParamsHKF::Vr)
        .def_readwrite("nptrans", &MineralSpeciesThermoParamsHKF::nptrans)
        .def_readwrite("a", &MineralSpeciesThermoParamsHKF::a)
        .def_readwrite("b", &MineralSpeciesThermoParamsHKF::b)
        .def_readwrite("c", &MineralSpeciesThermoParamsHKF::c)
        .def_readwrite("Ttr", &MineralSpeciesThermoParamsHKF::Ttr)
        .def_readwrite("Htr", &MineralSpeciesThermoParamsHKF::Htr)
        .def_readwrite("Vtr", &MineralSpeciesThermoParamsHKF::Vtr)
        .def_readwrite("dPdTtr", &MineralSpeciesThermoParamsHKF::dPdTtr)
        .def_readwrite("Tmax", &MineralSpeciesThermoParamsHKF::Tmax)
        ;
}

} // namespace Reaktoro
