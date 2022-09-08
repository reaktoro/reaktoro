// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright © 2014-2022 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// #pragma once

// // Reaktoro includes
// #include <Reaktoro/Common/StringList.hpp>
// #include <Reaktoro/Common/Types.hpp>
// #include <Reaktoro/Core/Database.hpp>
// #include <Reaktoro/Core/Reaction.hpp>
// #include <Reaktoro/Core/ReactionRateModel.hpp>

// namespace Reaktoro {

// /// The base type for all other classes defining more specific reactions.
// /// @ingroup Core
// class GenericReaction
// {
// public:
//     /// Construct a default GenericReaction object.
//     GenericReaction();

//     /// Construct a GenericReaction object with given species names.
//     explicit GenericReaction(const StringList& species);

//     /// Destroy this GenericReaction object.
//     virtual ~GenericReaction();

//     /// Set the unique name of the reaction.
//     auto setName(String name) -> GenericReaction&;

//     /// Set the rate model of the reaction.
//     auto setRateModel(const ReactionRateModelGenerator& model) -> GenericReaction&;

//     /// Set a unique name of the reaction (equivalent to GenericReaction::setName).
//     auto named(String name) -> GenericReaction&;

//     /// Set the rate model of the reaction (equivalent to GenericReaction::setRateModel).
//     auto set(const ReactionRateModelGenerator& model) -> GenericReaction&;

//     /// Return the name of the reaction.
//     auto name() const -> String;

//     /// Return the specified rate model of the reaction.
//     auto rateModel() const -> const ReactionRateModelGenerator&;

//     /// Convert this GenericReaction object into a Reaction object.
//     auto convert(const Database& db, const Strings& elements) const -> Reaction;

// private:
//     /// The name of the reaction.
//     String reaction_name;

//     /// The rate model of the reaction.
//     ReactionRateModelGenerator rate_model;
// };

// /// The base type for a generator of generic reactions.
// /// @ingroup Core
// class GenericReactionsGenerator
// {
// public:
//     /// Construct a default GenericReactionsGenerator object.
//     GenericReactionsGenerator();

//     /// Construct a GenericReactionsGenerator object with given species names.
//     explicit GenericReactionsGenerator(const StringList& species);

//     /// Destroy this GenericReactionsGenerator object.
//     virtual ~GenericReactionsGenerator();

//     /// Set the common rate model of the generated reactions.
//     auto setRateModel(const ReactionRateModelGenerator& model) -> GenericReactionsGenerator&;

//     /// Set the common rate model of the generated reactions (equivalent to GenericReactionsGenerator::setRateModel).
//     auto set(const ReactionRateModelGenerator& model) -> GenericReactionsGenerator&;

//     /// Return the specified common rate model of the generated reactions.
//     auto rateModel() const -> const ReactionRateModelGenerator&;

//     /// Convert this GenericReactionsGenerator object into a vector of GenericReaction objects.
//     auto convert(const Database& db, const Strings& elements) const -> Vec<GenericReaction>;

// private:
//     /// The names of the selected species to compose each generated phase.
//     Strings names;

//     /// The common rate model of the generated reactions.
//     ReactionRateModelGenerator rate_model;
// };

// /// The class used to define the reactions that will exist in the chemical system of interest.
// /// @ingroup Core
// class Reactions
// {
// public:
//     /// Construct a Reactions object.
//     /// @param db The database used to construct the species and elements in the reactions.
//     Reactions(const Database& db);

//     /// Construct a Reactions object with given database and generic reactions.
//     /// @param db The database used to construct the species and elements in the reactions.
//     /// @param gphases The generic reactions that will be converted into Reaction objects.
//     template<typename... GenericReaction>
//     Reactions(const Database& db, const GenericReaction&... gphases)
//     : Reactions(db)
//     {
//         static_assert(sizeof...(gphases) > 0);
//         addAux(gphases...);
//     }

//     /// Add a GenericReaction object into the Reactions container.
//     auto add(const GenericReaction& phase) -> void;

//     /// Add a GenericReactionsGenerator object into the Reactions container.
//     auto add(const GenericReactionsGenerator& generator) -> void;

//     /// Return the database object used to construct the species and elements in the reactions.
//     auto database() const -> const Database&;

//     /// Convert this Reactions object into a vector of Reaction objects.
//     auto convert() const -> Vec<Reaction>;

//     /// Convert this Reactions object into a vector of Reaction objects.
//     operator Vec<Reaction>() const;

// private:
//     /// The thermodynamic database used to deploy the Reaction objects from the GenericReaction ones.
//     Database db;

//     /// The GenericReaction objects collected so far with each call to Reactions::add method.
//     Vec<GenericReaction> generic_reactions;

//     /// The GenericReactionsGenerator objects collected so far with each call to Reactions::add method.
//     Vec<GenericReactionsGenerator> generators;

//     /// Add one or more GenericReaction or GenericReactionsGenerator objects into the Reactions container.
//     template<typename Arg, typename... Args>
//     auto addAux(const Arg& arg, const Args&... args) -> void
//     {
//         add(arg);
//         if constexpr (sizeof...(Args) > 0)
//             addAux(args...);
//     }
// };

// /// The class used to configure a mineral reaction.
// class MineralReaction : public GenericReaction
// {
// public:
//     /// Construct a default MineralReaction object.
//     MineralReaction() : GenericReaction() { initialize(); }

//     /// Construct an MineralReaction object with given species names.
//     explicit MineralReaction(const StringList& minerals) : GenericReaction(minerals) { initialize(); }

//     /// Initialize the default attributes of this MineralReaction object.
//     auto initialize() -> void
//     {
//         // setName("MineralReaction");
//         // setRateModel(rateModelIdealAqueous());
//     }
// };

// /// The class used to configure multiple mineral reactions simultaneously.
// class MineralReactions : public GenericPhasesGenerator
// {
// public:
//     /// Construct a default MineralReactions object.
//     MineralReactions() : GenericPhasesGenerator() { initialize(); }

//     /// Construct a MineralReactions object with given mineral names.
//     explicit MineralReactions(const StringList& minerals) : GenericPhasesGenerator(minerals) { initialize(); }

//     /// Initialize the default attributes of this MineralReactions object.
//     auto initialize() -> void
//     {
//         // setStateOfMatter(StateOfMatter::Solid);
//         // setAggregateState(AggregateState::Solid);
//         // setAdditionalAggregateStates({
//         //     AggregateState::CrystallineSolid,
//         //     AggregateState::AmorphousSolid
//         // });
//         // setActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid));
//         // setIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid));
//     }
// };
// } // namespace Reaktoro
