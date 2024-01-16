// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

auto isReactionRateModel(py::object obj) -> bool;
auto isReactionRateModelConvertible(py::object obj) -> bool;
auto isReactionRateModelGenerator(py::object obj) -> bool;
auto isReactionRateModelGeneratorConvertible(py::object obj) -> bool;
auto createReactionRateModel(py::object obj) -> ReactionRateModel;
auto createReactionRateModelGenerator(py::object obj) -> ReactionRateModelGenerator;

void exportReactions(py::module& m)
{
    py::class_<ReactionGeneratorArgs>(m, "ReactionGeneratorArgs")
        .def(py::init([]() { return ReactionGeneratorArgs{ Database(), SpeciesList(), PhaseList(), SurfaceList() }; }))
        .def(py::init<Database const&, SpeciesList const&, PhaseList const&, SurfaceList const&>())
        .def_property_readonly("database", [](ReactionGeneratorArgs const& self) { return self.database; }, "The thermodynamic database used to construct the chemical system where the reaction belongs to.")
        .def_property_readonly("species", [](ReactionGeneratorArgs const& self) { return self.species; }, "The species in the chemical system where the reaction belongs to.")
        .def_property_readonly("phases", [](ReactionGeneratorArgs const& self) { return self.phases; }, "The phases in the chemical system where the reaction belongs to.")
        .def_property_readonly("surfaces", [](ReactionGeneratorArgs const& self) { return self.surfaces; }, "The surfaces in the chemical system where the reaction belongs to.")
        ;

    auto setRateModel = [](GeneralReaction& self, py::object obj) -> GeneralReaction&
    {
        if(isReactionRateModel(obj))
            self.setRateModel(obj.cast<ReactionRateModel>());

        else if(isReactionRateModelConvertible(obj))
            self.setRateModel(createReactionRateModel(obj));

        else if(isReactionRateModelGenerator(obj))
            self.setRateModel(obj.cast<ReactionRateModelGenerator>());

        else if(isReactionRateModelGeneratorConvertible(obj))
            self.setRateModel(createReactionRateModelGenerator(obj));

        else errorif(true, "Expecting either a ReactionRateModel or ReactionRateModelGenerator object. Also possible is a Python callable object with a single argument of type ChemicalProps or ReactionRateModelGeneratorArgs.");

        return self;
    };

    py::class_<GeneralReaction>(m, "GeneralReaction")
        .def(py::init<>())
        .def(py::init<String const&>())
        .def("setName", &GeneralReaction::setName, return_internal_ref, "Set the unique name of the reaction.")
        .def("setEquation", &GeneralReaction::setEquation, return_internal_ref, "Set the equation of the reaction as a formatted string.")
        .def("setRateModel", setRateModel, return_internal_ref, "Set a reaction rate model or a reaction rate model generator for the reaction using a Python function.")
        .def("set", setRateModel, return_internal_ref, "Set a reaction rate model or a reaction rate model generator for the reaction using a Python function (equvalent to GeneralReaction.setRateModel).")
        .def("name", &GeneralReaction::name, return_internal_ref, "Return the name of the reaction.")
        .def("equation", &GeneralReaction::equation, return_internal_ref, "Return the reaction equation of the reaction.")
        .def("rateModel", &GeneralReaction::rateModel, return_internal_ref, "Return the reaction rate model of the reaction.")
        .def("rateModelGenerator", &GeneralReaction::rateModelGenerator, return_internal_ref, "Return the reaction rate model generator of the reaction.")
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
