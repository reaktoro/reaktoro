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
#include <Reaktoro/Core/Reactions.hpp>
#include <Reaktoro/Models/ReactionRateModels/Support/MineralReactions.hpp>
using namespace Reaktoro;

void exportMineralReactions(py::module& m)
{
    py::class_<MineralReactions>(m, "MineralReactions")
        .def(py::init<StringList const&>(), "Construct a MineralReactions object with given mineral names.")
        .def("setRateModel", py::overload_cast<MineralReactionRateModelGenerator const&>(&MineralReactions::setRateModel), "Set a common mineral reaction rate model generator for all minerals.")
        .def("setRateModel", py::overload_cast<String const&, MineralReactionRateModelGenerator const&>(&MineralReactions::setRateModel), "Set a mineral reaction rate model generator for a specific mineral.")
        .def("__call__", &MineralReactions::operator(), "Convert this MineralReactions object into a vector of Reaction objects.")
        ;

    py::implicitly_convertible<MineralReactions, Reactions>();
}
