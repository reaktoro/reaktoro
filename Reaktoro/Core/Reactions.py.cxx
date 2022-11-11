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
    py::class_<GeneralReaction>(m, "GeneralReaction")
        .def(py::init<>())
        .def(py::init<String const&>())
        .def("setName", &GeneralReaction::setName, "Set the unique name of the reaction.")
        .def("setEquation", &GeneralReaction::setEquation, "Set the equation of the reaction as a formatted string.")
        .def("setRateModel", &GeneralReaction::setRateModel, "Set the reaction rate model of the reaction.")
        .def("set", &GeneralReaction::set, "Set the reaction rate model of the reaction.")
        .def("name", &GeneralReaction::name, "Return the name of the reaction.")
        .def("equation", &GeneralReaction::equation, "Return the reaction equation of the reaction.")
        .def("rateModel", &GeneralReaction::rateModel, "Return the reaction rate model of the reaction.")
        .def("convert", &GeneralReaction::operator(), "Convert this GeneralReaction object into a Reaction object.") // NOTE: Do not use __call__ here because pybind11 will gladly cast a Python GeneralReaction object to a std::function of any type without any runtime errors! When checking if an argument in a ChemicalSystem constructor is of type ReactionGenerator or SurfaceGenerator (both objects of class std::function), the Python GeneralReaction object will be sucessfully converted, which is not expected.
        ;

    auto add = [](Reactions& self, py::object generator)
    {
        try { self.add(generator.cast<Reaction const&>()); }
        catch(...) {
            try { self.add(generator.cast<GeneralReaction const&>()); }
            catch(...) {
                try { self.add(generator.cast<ReactionGenerator>()); }
                catch(...) {
                    errorif(true, "Could not add reaction generator object to Reactions object:\n", py::str(generator));
                }
            }
        }
    };

    auto createReactions = [](py::args reaction_generators)
    {
        Reactions reactions;
        for(auto generator : reaction_generators)
        {
            try { reactions.add(generator.cast<Reaction const&>()); }
            catch(...) {
                try { reactions.add(generator.cast<GeneralReaction const&>()); }
                catch(...) {
                    try { reactions.add(generator.cast<ReactionGenerator>()); }
                    catch(...) {
                        errorif(true, "Could not create Reactions with reaction generator object:\n", py::str(generator));
                    }
                }
            }
        }

        return reactions;
    };

    py::class_<Reactions>(m, "Reactions")
        .def(py::init<>())
        .def(py::init(createReactions))
        .def("add", add)
        .def("convert", &Reactions::convert)
        ;
}
