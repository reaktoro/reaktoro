// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealAqueous.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealGas.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealSolution.hpp>
#include <Reaktoro/Thermodynamics/Ideal/ActivityModelIdealIonExchange.hpp>

namespace Reaktoro {

/// The auxiliary type used to specify phase species to be determined from element symbols.
struct Speciate
{
    /// The symbols of the elements composing the species in a phase.
    Strings symbols;

    /// Add other element symbols into the speciation list.
    auto operator+=(const Strings& othersymbols) -> Speciate& { symbols = merge(symbols, othersymbols); return *this; }
};

/// The auxiliary function used to specify phase species to be determined from element symbols.
inline auto speciate(const StringList& symbols) { return Speciate{symbols}; }

/// The auxiliary type used to specify species that should be filtered out when contructing a phase.
struct Exclude
{
    /// The tags which species cannot have when populating the species in a phase.
    Strings tags;

    /// Add other tags symbols into the exclude list.
    auto operator+=(const Strings& othertags) -> Exclude& { tags = merge(tags, othertags); return *this; }
};

/// The auxiliary function used to specify species that should be filtered out when contructing a phase.
inline auto exclude(const StringList& tags) { return Exclude{tags}; }

/// The base type for all other classes defining more specific phases.
/// @ingroup Core
class GenericPhase
{
public:
    /// Construct a default GenericPhase object.
    GenericPhase();

    /// Construct a GenericPhase object with given species names.
    explicit GenericPhase(const StringList& species);

    /// Construct a GenericPhase object with given element symbols.
    explicit GenericPhase(const Speciate& elements);

    /// Construct a GenericPhase object with given element symbols excluding the species with provided tags.
    explicit GenericPhase(const Speciate& elements, const Exclude& withtags);

    /// Construct a GenericPhase object excluding the species with provided tags.
    explicit GenericPhase(const Exclude& withtags);

    /// Destroy this GenericPhase object.
    virtual ~GenericPhase();

    /// Set the unique name of the phase.
    auto setName(String name) -> GenericPhase&;

    /// Set the state of matter of the phase.
    auto setStateOfMatter(StateOfMatter option) -> GenericPhase&;

    /// Set the aggregate state of the species in the phase.
    auto setAggregateState(AggregateState option) -> GenericPhase&;

    /// Set additional aggregate states to be considered when searching for species in a database.
    auto setAdditionalAggregateStates(const Vec<AggregateState>& options) -> GenericPhase&;

    /// Set the activity model of the phase.
    auto setActivityModel(const ActivityModelGenerator& model) -> GenericPhase&;

    /// Set the ideal activity model of the phase.
    auto setIdealActivityModel(const ActivityModelGenerator& model) -> GenericPhase&;

    /// Set a unique name of the phase (equivalent to GenericPhase::setName).
    auto named(String name) -> GenericPhase&;

    /// Set the state of matter of the phase (equivalent to GenericPhase::setStateOfMatter).
    auto set(StateOfMatter option) -> GenericPhase&;

    /// Set the aggregate state of the species in the phase (equivalent to GenericPhase::setAggregateState).
    auto set(AggregateState option) -> GenericPhase&;

    /// Set the activity model of the phase (equivalent to GenericPhase::setActivityModel).
    auto set(const ActivityModelGenerator& model) -> GenericPhase&;

    /// Return the name of the phase.
    auto name() const -> String;

    /// Return the state of matter of the phase.
    auto stateOfMatter() const -> StateOfMatter;

    /// Return the common aggregate state of the species composing the phase.
    auto aggregateState() const -> AggregateState;

    /// Return the additional aggregate states to be considered when searching for species in a database.
    auto additionalAggregateStates() const -> const Vec<AggregateState>&;

    /// Return the names of the selected species to compose the phase (empty if not given).
    auto species() const -> const Strings&;

    /// Return the element symbols for automatic species selection (empty if not given).
    auto elements() const -> const Strings&;

    /// Return the specified activity model of the phase.
    auto activityModel() const -> const ActivityModelGenerator&;

    /// Return the specified ideal activity model of the phase.
    auto idealActivityModel() const -> const ActivityModelGenerator&;

    /// Convert this GenericPhase object into a Phase object.
    auto convert(const Database& db, const Strings& elements) const -> Phase;

private:
    /// The name of the phase.
    String phasename;

    /// The state of matter of the phase.
    StateOfMatter stateofmatter = StateOfMatter::Solid;

    /// The aggregate state of the species in the phase.
    AggregateState aggregatestate = AggregateState::Undefined;

    /// The additional aggregate states used for searching species in the database.
    Vec<AggregateState> other_aggregate_states;

    /// The names of the selected species to compose the phase.
    Strings names;

    /// The element symbols for automatic selection of the species composing the phase.
    Strings symbols;

    /// The tags that indicate species to be excluded from the phase.
    Strings excludetags;

    /// The activity model of the phase.
    ActivityModelGenerator activity_model;

    /// The ideal activity model of the phase.
    ActivityModelGenerator ideal_activity_model;
};

/// The base type for a generator of generic phases with a single species.
/// @ingroup Core
class GenericPhasesGenerator
{
public:
    /// Construct a default GenericPhasesGenerator object.
    GenericPhasesGenerator();

    /// Construct a GenericPhasesGenerator object with given species names.
    explicit GenericPhasesGenerator(const StringList& species);

    /// Construct a GenericPhasesGenerator object with given element symbols.
    explicit GenericPhasesGenerator(const Speciate& elements);

    /// Construct a GenericPhasesGenerator object with given element symbols excluding the species with provided tags.
    explicit GenericPhasesGenerator(const Speciate& elements, const Exclude& withtags);

    /// Construct a GenericPhasesGenerator object excluding the species with provided tags.
    explicit GenericPhasesGenerator(const Exclude& withtags);

    /// Destroy this GenericPhasesGenerator object.
    virtual ~GenericPhasesGenerator();

    /// Set the common state of matter of the generated phases.
    auto setStateOfMatter(StateOfMatter option) -> GenericPhasesGenerator&;

    /// Set the common aggregate state of the species in the generated phases.
    auto setAggregateState(AggregateState option) -> GenericPhasesGenerator&;

    /// Set additional aggregate states to be considered when searching for species in a database.
    auto setAdditionalAggregateStates(const Vec<AggregateState>& options) -> GenericPhasesGenerator&;

    /// Set the common activity model of the generated phases.
    auto setActivityModel(const ActivityModelGenerator& model) -> GenericPhasesGenerator&;

    /// Set the common ideal activity model of the generated phases.
    auto setIdealActivityModel(const ActivityModelGenerator& model) -> GenericPhasesGenerator&;

    /// Set the common state of matter of the generated phases (equivalent to GenericPhasesGenerator::setStateOfMatter).
    auto set(StateOfMatter option) -> GenericPhasesGenerator&;

    /// Set the common aggregate state of the species in the generated phases (equivalent to GenericPhasesGenerator::setAggregateState).
    auto set(AggregateState option) -> GenericPhasesGenerator&;

    /// Set the common activity model of the generated phases (equivalent to GenericPhasesGenerator::setActivityModel).
    auto set(const ActivityModelGenerator& model) -> GenericPhasesGenerator&;

    /// Return the common state of matter of the generated phases.
    auto stateOfMatter() const -> StateOfMatter;

    /// Return the common aggregate state of the species composing the generated phases.
    auto aggregateState() const -> AggregateState;

    /// Return the additional aggregate states to be considered when searching for species in a database.
    auto additionalAggregateStates() const -> const Vec<AggregateState>&;

    /// Return the names of the selected species to compose the generated phase (empty if not given).
    auto species() const -> const Strings&;

    /// Return the element symbols for automatic species selection that will compose the generated phases (empty if not given).
    auto elements() const -> const Strings&;

    /// Return the specified common activity model of the generated phases.
    auto activityModel() const -> const ActivityModelGenerator&;

    /// Return the specified common ideal activity model of the generated phases.
    auto idealActivityModel() const -> const ActivityModelGenerator&;

    /// Convert this GenericPhasesGenerator object into a vector of GenericPhase objects.
    auto convert(const Database& db, const Strings& elements) const -> Vec<GenericPhase>;

private:
    /// The common state of matter of the generated phases.
    StateOfMatter stateofmatter = StateOfMatter::Solid;

    /// The common aggregate state of the species composing the generated phases.
    AggregateState aggregatestate = AggregateState::Undefined;

    /// The additional aggregate states used for searching species in the database.
    Vec<AggregateState> other_aggregate_states;

    /// The names of the selected species to compose each generated phase.
    Strings names;

    /// The element symbols for automatic selection of the species composing the generated phases.
    Strings symbols;

    /// The tags that indicate species to be excluded from the phase.
    Strings excludetags;

    /// The common activity model of the generated phases.
    ActivityModelGenerator activity_model;

    /// The common ideal activity model of the generated phases.
    ActivityModelGenerator ideal_activity_model;
};

/// The class used to define the phases that will constitute the chemical system of interest.
/// @ingroup Core
class Phases
{
public:
    /// Construct a Phases object.
    /// @param db The database used to construct the species and elements in the phases.
    Phases(const Database& db);

    /// Construct a Phases object with given database and generic phases.
    /// @param db The database used to construct the species and elements in the phases.
    /// @param gphases The generic phases that will be converted into Phase objects.
    template<typename... GenericPhases>
    Phases(const Database& db, const GenericPhases&... gphases)
    : Phases(db)
    {
        static_assert(sizeof...(gphases) > 0);
        addAux(gphases...);
    }

    /// Add a GenericPhase object into the Phases container.
    auto add(const GenericPhase& phase) -> void;

    /// Add a GenericPhasesGenerator object into the Phases container.
    auto add(const GenericPhasesGenerator& generator) -> void;

    /// Return the database object used to construct the species and elements in the phases.
    auto database() const -> const Database&;

    /// Convert this Phases object into a vector of Phase objects.
    auto convert() const -> Vec<Phase>;

    /// Convert this Phases object into a vector of Phase objects.
    operator Vec<Phase>() const;

private:
    /// The thermodynamic database used to deploy the Phase objects from the GenericPhase ones.
    Database db;

    /// The GenericPhase objects collected so far with each call to Phases::add method.
    Vec<GenericPhase> genericphases;

    /// The GenericPhaseGenerator objects collected so far with each call to Phases::add method.
    Vec<GenericPhasesGenerator> generators;

    /// Add one or more GenericPhase or GenericPhasesGenerator objects into the Phases container.
    template<typename Arg, typename... Args>
    auto addAux(const Arg& arg, const Args&... args) -> void
    {
        add(arg);
        if constexpr (sizeof...(Args) > 0)
            addAux(args...);
    }
};

/// The class used to configure an aqueous solution phase.
class AqueousPhase : public GenericPhase
{
public:
    /// Construct a default AqueousPhase object.
    AqueousPhase() : GenericPhase() { initialize(); }

    /// Construct an AqueousPhase object with given species names.
    explicit AqueousPhase(const StringList& species) : GenericPhase(species) { initialize(); }

    /// Construct an AqueousPhase object with given element symbols.
    explicit AqueousPhase(Speciate elements) : GenericPhase(elements += {"H", "O"}) { initialize(); }

    /// Construct an AqueousPhase object with given element symbols and tags indicating which species must be excluded from the final list.
    explicit AqueousPhase(Speciate elements, const Exclude& withtags) : GenericPhase(elements += {"H", "O"}, withtags) { initialize(); }

    /// Construct an AqueousPhase object with tags indicating which species must be excluded from the final list.
    explicit AqueousPhase(const Exclude& withtags) : GenericPhase(speciate("H O"), withtags) { initialize(); }

    /// Initialize the default attributes of this AqueousPhase object.
    auto initialize() -> void
    {
        setName("AqueousPhase");
        setStateOfMatter(StateOfMatter::Liquid);
        setAggregateState(AggregateState::Aqueous);
        setActivityModel(ActivityModelIdealAqueous());
        setIdealActivityModel(ActivityModelIdealAqueous());
    }
};

/// The class used to configure a gaseous solution phase.
class GaseousPhase : public GenericPhase
{
public:
    /// Construct a default GaseousPhase object.
    GaseousPhase() : GenericPhase() { initialize(); }

    /// Construct a GaseousPhase object with given species names.
    explicit GaseousPhase(const StringList& species) : GenericPhase(species) { initialize(); }

    /// Construct a GaseousPhase object with given element symbols.
    explicit GaseousPhase(const Speciate& elements) : GenericPhase(elements) { initialize(); }

    /// Construct a GaseousPhase object with given element symbols and tags indicating which species must be excluded from the final list.
    explicit GaseousPhase(const Speciate& elements, const Exclude& withtags) : GenericPhase(elements, withtags) { initialize(); }

    /// Construct a GaseousPhase object with tags indicating which species must be excluded from the final list.
    explicit GaseousPhase(const Exclude& withtags) : GenericPhase(withtags) { initialize(); }

    /// Initialize the default attributes of this GaseousPhase object.
    auto initialize() -> void
    {
        setName("GaseousPhase");
        setStateOfMatter(StateOfMatter::Gas);
        setAggregateState(AggregateState::Gas);
        setActivityModel(ActivityModelIdealGas());
        setIdealActivityModel(ActivityModelIdealGas());
    }
};

/// The class used to configure a liquid solution phase.
class LiquidPhase : public GenericPhase
{
public:
    /// Construct a default LiquidPhase object.
    LiquidPhase() : GenericPhase() { initialize(); }

    /// Construct a LiquidPhase object with given species names.
    explicit LiquidPhase(const StringList& species) : GenericPhase(species) { initialize(); }

    /// Construct a LiquidPhase object with given element symbols.
    explicit LiquidPhase(const Speciate& elements) : GenericPhase(elements) { initialize(); }

    /// Construct a LiquidPhase object with given element symbols excluding the species with provided tags.
    explicit LiquidPhase(const Speciate& elements, const Exclude& withtags) : GenericPhase(elements, withtags) { initialize(); };

    /// Construct a LiquidPhase object excluding the species with provided tags.
    explicit LiquidPhase(const Exclude& withtags) : GenericPhase(withtags) { initialize(); };

    /// Initialize the default attributes of this LiquidPhase object.
    auto initialize() -> void
    {
        setName("LiquidPhase");
        setStateOfMatter(StateOfMatter::Liquid);
        setAggregateState(AggregateState::Liquid);
        setActivityModel(ActivityModelIdealSolution());
        setIdealActivityModel(ActivityModelIdealSolution());
    }
};

/// The class used to configure a solid solution phase.
class SolidPhase : public GenericPhase
{
public:
    /// Construct a SolidPhase object with given species names.
    explicit SolidPhase(const StringList& species) : GenericPhase(species) { initialize(); }

    /// Initialize the default attributes of this SolidPhase object.
    auto initialize() -> void
    {
        String phasename;
        auto i = 0;
        for(auto&& name : species())
            phasename += (i++ == 0) ? name : "-" + name;

        setName("SolidPhase");
        setStateOfMatter(StateOfMatter::Solid);
        setAggregateState(AggregateState::Solid);
        setActivityModel(ActivityModelIdealSolution());
        setIdealActivityModel(ActivityModelIdealSolution());
    }
};

/// The class used to configure a pure mineral phase.
class MineralPhase : public GenericPhase
{
public:
    /// Construct a default MineralPhase object.
    explicit MineralPhase(String mineral) : GenericPhase(mineral) { initialize(); }

    /// Initialize the default attributes of this MineralPhase object.
    auto initialize() -> void
    {
        setName(species().front());
        setStateOfMatter(StateOfMatter::Solid);
        setAggregateState(AggregateState::Solid);
        setAdditionalAggregateStates({
            AggregateState::CrystallineSolid,
            AggregateState::AmorphousSolid
        });
        setActivityModel(ActivityModelIdealSolution());
        setIdealActivityModel(ActivityModelIdealSolution());
    }
};

/// The class used to configure automatic selection of pure mineral phases.
class MineralPhases : public GenericPhasesGenerator
{
public:
    /// Construct a default MineralPhases object.
    MineralPhases() : GenericPhasesGenerator() { initialize(); }

    /// Construct a MineralPhases object with given species names.
    explicit MineralPhases(const StringList& species) : GenericPhasesGenerator(species) { initialize(); }

    /// Construct a MineralPhases object with given element symbols.
    explicit MineralPhases(const Speciate& elements) : GenericPhasesGenerator(elements) { initialize(); }

    /// Construct a MineralPhases object with given element symbols excluding the species with provided tags
    explicit MineralPhases(const Speciate& elements, const Exclude& withtags) : GenericPhasesGenerator(elements, withtags) { initialize(); };

    /// Construct a MineralPhases object excluding the species with provided tags
    explicit MineralPhases(const Exclude& withtags) : GenericPhasesGenerator(withtags) { initialize(); };

    /// Initialize the default attributes of this MineralPhases object.
    auto initialize() -> void
    {
        setStateOfMatter(StateOfMatter::Solid);
        setAggregateState(AggregateState::Solid);
        setAdditionalAggregateStates({
            AggregateState::CrystallineSolid,
            AggregateState::AmorphousSolid
        });
        setActivityModel(ActivityModelIdealSolution());
        setIdealActivityModel(ActivityModelIdealSolution());
    }
};

/// The class used to configure a pure condensed phase.
class CondensedPhase : public GenericPhase
{
public:
    /// Construct a default CondensedPhase object.
    explicit CondensedPhase(String species) : GenericPhase(species) { initialize(); }

    /// Initialize the default attributes of this CondensedPhase object.
    auto initialize() -> void
    {
        setName(species().front());
        setStateOfMatter(StateOfMatter::Condensed);
        setAggregateState(AggregateState::CondensedPhase);
        setAdditionalAggregateStates({
            AggregateState::Liquid,
            AggregateState::LiquidCrystal,
            AggregateState::Solid,
            AggregateState::CrystallineSolid,
            AggregateState::AmorphousSolid
        });
        setActivityModel(ActivityModelIdealSolution());
        setIdealActivityModel(ActivityModelIdealSolution());
    }
};

/// The class used to configure automatic selection of pure condensed phases.
class CondensedPhases : public GenericPhasesGenerator
{
public:
    /// Construct a default CondensedPhases object.
    CondensedPhases() : GenericPhasesGenerator() { initialize(); }

    /// Construct a CondensedPhases object with given species names.
    explicit CondensedPhases(const StringList& species) : GenericPhasesGenerator(species) { initialize(); }

    /// Construct a CondensedPhases object with given element symbols.
    explicit CondensedPhases(const Speciate& elements) : GenericPhasesGenerator(elements) { initialize(); }

    /// Construct a CondensedPhases object with given element symbols excluding the species with provided tags
    explicit CondensedPhases(const Speciate& elements, const Exclude& withtags) : GenericPhasesGenerator(elements, withtags) { initialize(); };

    /// Construct a CondensedPhases object excluding the species with provided tags
    explicit CondensedPhases(const Exclude& withtags) : GenericPhasesGenerator(withtags) { initialize(); };

    /// Initialize the default attributes of this CondensedPhases object.
    auto initialize() -> void
    {
        setStateOfMatter(StateOfMatter::Condensed);
        setAggregateState(AggregateState::CondensedPhase);
        setAdditionalAggregateStates({
            AggregateState::Liquid,
            AggregateState::LiquidCrystal,
            AggregateState::Solid,
            AggregateState::CrystallineSolid,
            AggregateState::AmorphousSolid
        });
        setActivityModel(ActivityModelIdealSolution());
        setIdealActivityModel(ActivityModelIdealSolution());
    }
};

/// The class used to configure an ion exchange phase.
class IonExchangePhase : public GenericPhase
{
public:

    /// Construct an IonExchangePhase object with given species names.
    explicit IonExchangePhase(const StringList& species) : GenericPhase(species) { initialize(); }

    /// Initialize the default attributes of this IonExchangePhase object.
    auto initialize() -> void
    {
        setName("IonExchangePhase");
        setStateOfMatter(StateOfMatter::Solid);
        setAggregateState(AggregateState::IonExchange);
        setActivityModel(ActivityModelIdealIonExchange());
        setIdealActivityModel(ActivityModelIdealIonExchange());

    }
};

} // namespace Reaktoro
