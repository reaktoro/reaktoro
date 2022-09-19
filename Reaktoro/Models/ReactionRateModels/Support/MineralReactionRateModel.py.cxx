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
#include <Reaktoro/Core/Model.py.hxx>
#include <Reaktoro/Core/ReactionStandardThermoProps.hpp>
#include <Reaktoro/Models/ReactionRateModels/Support/MineralReactionRateModel.hpp>
using namespace Reaktoro;

void exportMineralReactionRateModel(py::module& m)
{
    py::class_<MineralReactionRateModelArgs>(m, "MineralReactionRateModelArgs")
        .def_property_readonly("props", [](MineralReactionRateModelArgs const& self) { return self.props; }, "The properties of the chemical system.")
        .def_property_readonly("aprops", [](MineralReactionRateModelArgs const& self) { return self.aprops; }, "The properties of the aqueous solution.")
        .def_property_readonly("T", [](MineralReactionRateModelArgs const& self) { return self.T; }, "The temperature of the system (in K).")
        .def_property_readonly("P", [](MineralReactionRateModelArgs const& self) { return self.P; }, "The pressure of the system (in Pa).")
        .def_property_readonly("pH", [](MineralReactionRateModelArgs const& self) { return self.pH; }, "The pH of the aqueous solution.")
        .def_property_readonly("Omega", [](MineralReactionRateModelArgs const& self) { return self.Omega; }, "The saturation index Omega = IAP}/K of the mineral reaction.")
        .def_property_readonly("area", [](MineralReactionRateModelArgs const& self) { return self.area; }, "The surface area between the mineral and the aqueous solution.")
        ;

    exportModel<ReactionRate, MineralReactionRateModelArgs>(m, "MineralReactionRateModel");
}
