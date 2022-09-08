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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Reaction.hpp>
#include <Reaktoro/Core/Reactions.hpp>
using namespace Reaktoro;

void exportReactions(py::module& m)
{
    auto createReactions = [](py::args reaction_generators)
    {
        Reactions reactions;
        for(auto generator : reaction_generators)
        {
            try { reactions.add(generator.cast<ReactionGenerator>()); }
            catch(...)
            {
                try { reactions.add(generator.cast<const Reaction&>()); }
                catch(...)
                {
                    errorif(true, "Could not create Reactions with reaction generator object:\n", py::str(generator));
                }
            }
        }

        return reactions;
    };

    auto add = [](Reactions& self, py::object generator)
    {
        try { self.add(generator.cast<ReactionGenerator>()); }
        catch(...)
        {
            try { self.add(generator.cast<const Reaction&>()); }
            catch(...)
            {
                errorif(true, "Could not add reaction generator object to Reactions object:\n", py::str(generator));
            }
        }
    };

    py::class_<Reactions>(m, "Reactions")
        .def(py::init<>())
        .def(py::init(createReactions))
        .def("add", add)
        .def("convert", &Reactions::convert)
        ;
}
