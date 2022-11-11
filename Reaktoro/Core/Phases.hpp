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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/StringList.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealAqueous.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealGas.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealSolution.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelIdealIonExchange.hpp>

namespace Reaktoro {

/// The auxiliary type used to specify phase species to be determined from element symbols.
struct Speciate
{
    /// The symbols of the elements composing the species in a phase.
    Strings symbols;

    /// Add other element symbols into the speciation list.
    auto operator+=(Strings const& othersymbols) -> Speciate& { symbols = merge(symbols, othersymbols); return *this; }
};

/// The auxiliary function used to specify phase species to be determined from element symbols.
auto speciate(StringList const& symbols) -> Speciate;

/// The auxiliary type used to specify species that should be filtered out when contructing a phase.
struct Exclude
{
    /// The tags which species cannot have when populating the species in a phase.
    Strings tags;

    /// Add other tags symbols into the exclude list.
    auto operator+=(Strings const& othertags) -> Exclude& { tags = merge(tags, othertags); return *this; }
};

/// The auxiliary function used to specify species that should be filtered out when contructing a phase.
auto exclude(StringList const& tags) -> Exclude;

/// The base type for all other classes defining more specific phases.
/// @ingroup Core
class GeneralPhase
{
public:
    /// Construct a default GeneralPhase object.
    GeneralPhase();

    /// Construct a GeneralPhase object with given species names.
    explicit GeneralPhase(StringList const& species);

    /// Construct a GeneralPhase object with given element symbols.
    explicit GeneralPhase(Speciate const& elements);

    /// Construct a GeneralPhase object with given element symbols excluding the species with provided tags.
    explicit GeneralPhase(Speciate const& elements, Exclude const& withtags);

    /// Construct a GeneralPhase object excluding the species with provided tags.
    explicit GeneralPhase(Exclude const& withtags);

    /// Destroy this GeneralPhase object.
    virtual ~GeneralPhase();

    /// Set the unique name of the phase.
    auto setName(String name) -> GeneralPhase&;

    /// Set the state of matter of the phase.
    auto setStateOfMatter(StateOfMatter option) -> GeneralPhase&;

    /// Set the aggregate state of the species in the phase.
    auto setAggregateState(AggregateState option) -> GeneralPhase&;

    /// Set additional aggregate states to be considered when searching for species in a database.
    auto setAdditionalAggregateStates(Vec<AggregateState> const& options) -> GeneralPhase&;

    /// Set the activity model of the phase.
    auto setActivityModel(ActivityModelGenerator const& model) -> GeneralPhase&;

    /// Set the ideal activity model of the phase.
    auto setIdealActivityModel(ActivityModelGenerator const& model) -> GeneralPhase&;

    /// Set a unique name of the phase (equivalent to GeneralPhase::setName).
    auto named(String name) -> GeneralPhase&;

    /// Set the state of matter of the phase (equivalent to GeneralPhase::setStateOfMatter).
    auto set(StateOfMatter option) -> GeneralPhase&;

    /// Set the aggregate state of the species in the phase (equivalent to GeneralPhase::setAggregateState).
    auto set(AggregateState option) -> GeneralPhase&;

    /// Set the activity model of the phase (equivalent to GeneralPhase::setActivityModel).
    auto set(ActivityModelGenerator const& model) -> GeneralPhase&;

    /// Return the name of the phase.
    auto name() const -> String;

    /// Return the state of matter of the phase.
    auto stateOfMatter() const -> StateOfMatter;

    /// Return the common aggregate state of the species composing the phase.
    auto aggregateState() const -> AggregateState;

    /// Return the additional aggregate states to be considered when searching for species in a database.
    auto additionalAggregateStates() const -> Vec<AggregateState> const&;

    /// Return the names of the selected species to compose the phase (empty if not given).
    auto species() const -> Strings const&;

    /// Return the element symbols for automatic species selection (empty if not given).
    auto elements() const -> Strings const&;

    /// Return the specified activity model of the phase.
    auto activityModel() const -> ActivityModelGenerator const&;

    /// Return the specified ideal activity model of the phase.
    auto idealActivityModel() const -> ActivityModelGenerator const&;

    /// Convert this GeneralPhase object into a Phase object.
    auto convert(Database const& db, Strings const& elements) const -> Phase;

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

/// The base type for a generator of general phases with a single species.
/// @ingroup Core
class GeneralPhasesGenerator
{
public:
    /// Construct a default GeneralPhasesGenerator object.
    GeneralPhasesGenerator();

    /// Construct a GeneralPhasesGenerator object with given species names.
    explicit GeneralPhasesGenerator(StringList const& species);

    /// Construct a GeneralPhasesGenerator object with given element symbols.
    explicit GeneralPhasesGenerator(Speciate const& elements);

    /// Construct a GeneralPhasesGenerator object with given element symbols excluding the species with provided tags.
    explicit GeneralPhasesGenerator(Speciate const& elements, Exclude const& withtags);

    /// Construct a GeneralPhasesGenerator object excluding the species with provided tags.
    explicit GeneralPhasesGenerator(Exclude const& withtags);

    /// Destroy this GeneralPhasesGenerator object.
    virtual ~GeneralPhasesGenerator();

    /// Set the common state of matter of the generated phases.
    auto setStateOfMatter(StateOfMatter option) -> GeneralPhasesGenerator&;

    /// Set the common aggregate state of the species in the generated phases.
    auto setAggregateState(AggregateState option) -> GeneralPhasesGenerator&;

    /// Set additional aggregate states to be considered when searching for species in a database.
    auto setAdditionalAggregateStates(Vec<AggregateState> const& options) -> GeneralPhasesGenerator&;

    /// Set the common activity model of the generated phases.
    auto setActivityModel(ActivityModelGenerator const& model) -> GeneralPhasesGenerator&;

    /// Set the common ideal activity model of the generated phases.
    auto setIdealActivityModel(ActivityModelGenerator const& model) -> GeneralPhasesGenerator&;

    /// Set the common state of matter of the generated phases (equivalent to GeneralPhasesGenerator::setStateOfMatter).
    auto set(StateOfMatter option) -> GeneralPhasesGenerator&;

    /// Set the common aggregate state of the species in the generated phases (equivalent to GeneralPhasesGenerator::setAggregateState).
    auto set(AggregateState option) -> GeneralPhasesGenerator&;

    /// Set the common activity model of the generated phases (equivalent to GeneralPhasesGenerator::setActivityModel).
    auto set(ActivityModelGenerator const& model) -> GeneralPhasesGenerator&;

    /// Return the common state of matter of the generated phases.
    auto stateOfMatter() const -> StateOfMatter;

    /// Return the common aggregate state of the species composing the generated phases.
    auto aggregateState() const -> AggregateState;

    /// Return the additional aggregate states to be considered when searching for species in a database.
    auto additionalAggregateStates() const -> Vec<AggregateState> const&;

    /// Return the names of the selected species to compose the generated phase (empty if not given).
    auto species() const -> Strings const&;

    /// Return the element symbols for automatic species selection that will compose the generated phases (empty if not given).
    auto elements() const -> Strings const&;

    /// Return the specified common activity model of the generated phases.
    auto activityModel() const -> ActivityModelGenerator const&;

    /// Return the specified common ideal activity model of the generated phases.
    auto idealActivityModel() const -> ActivityModelGenerator const&;

    /// Convert this GeneralPhasesGenerator object into a vector of GeneralPhase objects.
    auto convert(Database const& db, Strings const& elements) const -> Vec<GeneralPhase>;

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

template <typename T, typename... Ts>
constexpr auto _areGeneralPhasesImpl()
{
    constexpr auto aux = isBaseOf<GeneralPhase, T> || isBaseOf<GeneralPhasesGenerator, T>;
    if constexpr (sizeof...(Ts))
        return aux && _areGeneralPhasesImpl<Ts...>();
    else return aux;
}

/// Used to determine if `T` and all types in `Ts` are either GeneralPhase or GeneralPhaseGenerator.
template<typename T, typename... Ts>
constexpr auto areGeneralPhases = _areGeneralPhasesImpl<T, Ts...>();

/// The class used to define the phases that will constitute the chemical system of interest.
/// @ingroup Core
class Phases
{
public:
    /// Construct a Phases object.
    /// @param db The database used to construct the species and elements in the phases.
    Phases(Database const& db);

    /// Construct a Phases object with given database and general phases.
    /// @param db The database used to construct the species and elements in the phases.
    /// @param gphases The general phases that will be converted into Phase objects.
    template<typename... GeneralPhases, EnableIf<areGeneralPhases<GeneralPhases...>>...>
    explicit Phases(Database const& db, GeneralPhases const&... gphases)
    : Phases(db)
    {
        static_assert(sizeof...(gphases) > 0);
        addAux(gphases...);
    }

    // TODO: Implement `auto add(Phase const& phase) -> void` as well, in case the user provides a
    // Phase object. This will need a new data member `phases` of type `Vec<Phase>`.

    /// Add a GeneralPhase object into the Phases container.
    auto add(GeneralPhase const& phase) -> void;

    /// Add a GeneralPhasesGenerator object into the Phases container.
    auto add(GeneralPhasesGenerator const& generator) -> void;

    /// Return the database object used to construct the species and elements in the phases.
    auto database() const -> Database const&;

    /// Return the GeneralPhase objects collected so far with each call to Phases::add method.
    auto generalPhases() const -> Vec<GeneralPhase> const&;

    /// Return the GeneralPhaseGenerator objects collected so far with each call to Phases::add method.
    auto generalPhasesGenerators() const -> Vec<GeneralPhasesGenerator> const&;

    /// Convert this Phases object into a vector of Phase objects.
    auto convert() const -> Vec<Phase>;

    /// Convert this Phases object into a vector of Phase objects.
    operator Vec<Phase>() const;

private:
    /// The thermodynamic database used to deploy the Phase objects from the GeneralPhase ones.
    Database db;

    /// The GeneralPhase objects collected so far with each call to Phases::add method.
    Vec<GeneralPhase> generalphases;

    /// The GeneralPhaseGenerator objects collected so far with each call to Phases::add method.
    Vec<GeneralPhasesGenerator> generators;

    /// Add one or more GeneralPhase or GeneralPhasesGenerator objects into the Phases container.
    template<typename Arg, typename... Args>
    auto addAux(Arg const& arg, Args const&... args) -> void
    {
        add(arg);
        if constexpr (sizeof...(Args) > 0)
            addAux(args...);
    }
};

/// The class used to configure an aqueous solution phase.
class AqueousPhase : public GeneralPhase
{
public:
    /// Construct a default AqueousPhase object.
    AqueousPhase() : GeneralPhase() { initialize(); }

    /// Construct an AqueousPhase object with given species names.
    explicit AqueousPhase(StringList const& species) : GeneralPhase(species) { initialize(); }

    /// Construct an AqueousPhase object with given element symbols.
    explicit AqueousPhase(Speciate elements) : GeneralPhase(elements += {"H", "O"}) { initialize(); }

    /// Construct an AqueousPhase object with given element symbols and tags indicating which species must be excluded from the final list.
    explicit AqueousPhase(Speciate elements, Exclude const& withtags) : GeneralPhase(elements += {"H", "O"}, withtags) { initialize(); }

    /// Construct an AqueousPhase object with tags indicating which species must be excluded from the final list.
    explicit AqueousPhase(Exclude const& withtags) : GeneralPhase(speciate("H O"), withtags) { initialize(); }

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
class GaseousPhase : public GeneralPhase
{
public:
    /// Construct a default GaseousPhase object.
    GaseousPhase() : GeneralPhase() { initialize(); }

    /// Construct a GaseousPhase object with given species names.
    explicit GaseousPhase(StringList const& species) : GeneralPhase(species) { initialize(); }

    /// Construct a GaseousPhase object with given element symbols.
    explicit GaseousPhase(Speciate const& elements) : GeneralPhase(elements) { initialize(); }

    /// Construct a GaseousPhase object with given element symbols and tags indicating which species must be excluded from the final list.
    explicit GaseousPhase(Speciate const& elements, Exclude const& withtags) : GeneralPhase(elements, withtags) { initialize(); }

    /// Construct a GaseousPhase object with tags indicating which species must be excluded from the final list.
    explicit GaseousPhase(Exclude const& withtags) : GeneralPhase(withtags) { initialize(); }

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
class LiquidPhase : public GeneralPhase
{
public:
    /// Construct a default LiquidPhase object.
    LiquidPhase() : GeneralPhase() { initialize(); }

    /// Construct a LiquidPhase object with given species names.
    explicit LiquidPhase(StringList const& species) : GeneralPhase(species) { initialize(); }

    /// Construct a LiquidPhase object with given element symbols.
    explicit LiquidPhase(Speciate const& elements) : GeneralPhase(elements) { initialize(); }

    /// Construct a LiquidPhase object with given element symbols excluding the species with provided tags.
    explicit LiquidPhase(Speciate const& elements, Exclude const& withtags) : GeneralPhase(elements, withtags) { initialize(); };

    /// Construct a LiquidPhase object excluding the species with provided tags.
    explicit LiquidPhase(Exclude const& withtags) : GeneralPhase(withtags) { initialize(); };

    /// Initialize the default attributes of this LiquidPhase object.
    auto initialize() -> void
    {
        setName("LiquidPhase");
        setStateOfMatter(StateOfMatter::Liquid);
        setAggregateState(AggregateState::Liquid);
        setActivityModel(ActivityModelIdealSolution(StateOfMatter::Liquid));
        setIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Liquid));
    }
};

/// The class used to configure a solid solution phase.
class SolidPhase : public GeneralPhase
{
public:
    /// Construct a SolidPhase object with given species names.
    explicit SolidPhase(StringList const& species) : GeneralPhase(species) { initialize(); }

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
        setAdditionalAggregateStates({AggregateState::CrystallineSolid});
        setActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid));
        setIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid));
    }
};

/// The class used to configure a pure mineral phase.
class MineralPhase : public GeneralPhase
{
public:
    /// Construct a default MineralPhase object.
    explicit MineralPhase(String mineral) : GeneralPhase(mineral) { initialize(); }

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
        setActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid));
        setIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid));
    }
};

/// The class used to configure automatic selection of pure mineral phases.
class MineralPhases : public GeneralPhasesGenerator
{
public:
    /// Construct a default MineralPhases object.
    MineralPhases() : GeneralPhasesGenerator() { initialize(); }

    /// Construct a MineralPhases object with given species names.
    explicit MineralPhases(StringList const& species) : GeneralPhasesGenerator(species) { initialize(); }

    /// Construct a MineralPhases object with given element symbols.
    explicit MineralPhases(Speciate const& elements) : GeneralPhasesGenerator(elements) { initialize(); }

    /// Construct a MineralPhases object with given element symbols excluding the species with provided tags
    explicit MineralPhases(Speciate const& elements, Exclude const& withtags) : GeneralPhasesGenerator(elements, withtags) { initialize(); };

    /// Construct a MineralPhases object excluding the species with provided tags
    explicit MineralPhases(Exclude const& withtags) : GeneralPhasesGenerator(withtags) { initialize(); };

    /// Initialize the default attributes of this MineralPhases object.
    auto initialize() -> void
    {
        setStateOfMatter(StateOfMatter::Solid);
        setAggregateState(AggregateState::Solid);
        setAdditionalAggregateStates({
            AggregateState::CrystallineSolid,
            AggregateState::AmorphousSolid
        });
        setActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid));
        setIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid));
    }
};

/// The class used to configure a pure condensed phase.
class CondensedPhase : public GeneralPhase
{
public:
    /// Construct a default CondensedPhase object.
    explicit CondensedPhase(String species) : GeneralPhase(species) { initialize(); }

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
        setActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)); // TODO: Create ActivityModelIdealCondensedPhase(melting_temperature) or rely on a MeltingTemperature attribute in the species block of the database
        setIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid));
    }
};

/// The class used to configure automatic selection of pure condensed phases.
class CondensedPhases : public GeneralPhasesGenerator
{
public:
    /// Construct a default CondensedPhases object.
    CondensedPhases() : GeneralPhasesGenerator() { initialize(); }

    /// Construct a CondensedPhases object with given species names.
    explicit CondensedPhases(StringList const& species) : GeneralPhasesGenerator(species) { initialize(); }

    /// Construct a CondensedPhases object with given element symbols.
    explicit CondensedPhases(Speciate const& elements) : GeneralPhasesGenerator(elements) { initialize(); }

    /// Construct a CondensedPhases object with given element symbols excluding the species with provided tags
    explicit CondensedPhases(Speciate const& elements, Exclude const& withtags) : GeneralPhasesGenerator(elements, withtags) { initialize(); };

    /// Construct a CondensedPhases object excluding the species with provided tags
    explicit CondensedPhases(Exclude const& withtags) : GeneralPhasesGenerator(withtags) { initialize(); };

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
        setActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid)); // TODO: Create ActivityModelIdealCondensedPhase(melting_temperature)
        setIdealActivityModel(ActivityModelIdealSolution(StateOfMatter::Solid));
    }
};

/// The class used to configure an ion exchange phase.
class IonExchangePhase : public GeneralPhase
{
public:

    /// Construct an IonExchangePhase object with given species names.
    explicit IonExchangePhase(StringList const& species) : GeneralPhase(species) { initialize(); }

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
