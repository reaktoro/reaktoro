// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <Reaktoro/Core/Phase.hpp>
#include <Reaktoro/Core/ThermoEngine.hpp>
#include <Reaktoro/Models/ActivityModelIdeal.hpp>

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
inline auto speciate(const StringList& symbols) { return Speciate{symbols}; };

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

    /// Destroy this GenericPhase object.
    virtual ~GenericPhase();

    /// Set a unique name for the phase.
    auto setName(String name) -> GenericPhase&;

    /// Set the state of matter for the phase.
    auto setStateOfMatter(StateOfMatter option) -> GenericPhase&;

    /// Set the aggregate state of the species in the phase.
    auto setAggregateState(AggregateState option) -> GenericPhase&;

    /// Set the activity model for the phase.
    auto setActivityModel(const ActivityModel& model) -> GenericPhase&;

    /// Set a unique name for the phase (equivalent to GenericPhase::setName).
    auto named(String name) -> GenericPhase&;

    /// Set the state of matter for the phase (equivalent to GenericPhase::setStateOfMatter).
    auto set(StateOfMatter option) -> GenericPhase&;

    /// Set the aggregate state of the species in the phase (equivalent to GenericPhase::setAggregateState).
    auto set(AggregateState option) -> GenericPhase&;

    /// Set the activity model for the phase (equivalent to GenericPhase::setActivityModel).
    auto set(const ActivityModel& model) -> GenericPhase&;

    /// Return the name of the phase.
    auto name() const -> String;

    /// Return the state of matter of the phase.
    auto stateOfMatter() const -> StateOfMatter;

    /// Return the common aggregate state of the species composing the phase.
    auto aggregateState() const -> AggregateState;

    /// Return the names of the selected species to compose the phase (empty if not given).
    auto species() const -> const Strings&;

    /// Return the element symbols for automatic species selection (empty if not given).
    auto elements() const -> const Strings&;

    /// Return the specified activity model of the phase.
    auto activityModel() const -> const ActivityModel&;

    /// Convert this GenericPhase object into a Phase object.
    auto convert(const ThermoEngine& engine, const Strings& elements) const -> Phase;

private:
    /// The name of the phase.
    String phasename;

    /// The state of matter of the phase.
    StateOfMatter stateofmatter = StateOfMatter::Solid;

    /// The aggregate state of the species in the phase.
    AggregateState aggregatestate = AggregateState::Undefined;

    /// The names of the selected species to compose the phase.
    Strings names;

    /// The element symbols for automatic selection of the species composing the phase.
    Strings symbols;

    /// The activity model of the phase.
    ActivityModel activitymodel;
};

/// The base type for all other classes defining generic pure phases at once.
/// @ingroup Core
class GenericPhases
{
public:
    /// Construct a default GenericPhases object.
    GenericPhases();

    /// Construct a GenericPhases object with given species names.
    explicit GenericPhases(const StringList& species);

    /// Construct a GenericPhases object with given element symbols.
    explicit GenericPhases(const Speciate& elements);

    /// Destroy this GenericPhases object.
    virtual ~GenericPhases();

    /// Set the common state of matter for the pure phases.
    auto setStateOfMatter(StateOfMatter option) -> GenericPhases&;

    /// Set the common aggregate state of the species in the pure phases.
    auto setAggregateState(AggregateState option) -> GenericPhases&;

    /// Set the common activity model for the pure phases.
    auto setActivityModel(const ActivityModel& model) -> GenericPhases&;

    /// Set the common state of matter for the pure phases (equivalent to GenericPhases::setStateOfMatter).
    auto set(StateOfMatter option) -> GenericPhases&;

    /// Set the common aggregate state of the species in the pure phases (equivalent to GenericPhases::setAggregateState).
    auto set(AggregateState option) -> GenericPhases&;

    /// Set the common activity model for the pure phases (equivalent to GenericPhases::setActivityModel).
    auto set(const ActivityModel& model) -> GenericPhases&;

    /// Return the common state of matter of the phase.
    auto stateOfMatter() const -> StateOfMatter;

    /// Return the common aggregate state of the species composing the pure phases.
    auto aggregateState() const -> AggregateState;

    /// Return the names of the selected species to compose the pure phase (empty if not given).
    auto species() const -> const Strings&;

    /// Return the element symbols for automatic species selection that will compose the pure phases (empty if not given).
    auto elements() const -> const Strings&;

    /// Return the specified common activity model of the pure phases.
    auto activityModel() const -> const ActivityModel&;

    /// Convert this GenericPhases object into a vector of GenericPhase objects.
    auto convert(const ThermoEngine& engine, const Strings& elements) const -> Vec<GenericPhase>;

private:
    /// The common state of matter of the pure phases.
    StateOfMatter stateofmatter = StateOfMatter::Solid;

    /// The common aggregate state of the species composing the pure phases.
    AggregateState aggregatestate = AggregateState::Undefined;

    /// The names of the selected species to compose the pure phases.
    Strings names;

    /// The element symbols for automatic selection of the species composing the pure phases.
    Strings symbols;

    /// The common activity model of the pure phases.
    ActivityModel activitymodel;
};

/// The class used to define the phases that will constitute the chemical system of interest.
/// @ingroup Core
class Phases
{
public:
    /// Construct a Phases object with given GenericPhase and GenericPhases objects.
    template<typename... Args>
    Phases(const ThermoEngine& engine, const Args&... phases)
    : engine(engine)
    {
        collectElements(phases...);
        addPhases(phases...);
        fixDuplicatePhaseNames();
    }

private:
    /// The thermodynamic engine used to deploy the Phase objects from the GenericPhase ones.
    ThermoEngine engine;

    /// The GenericPhase objects collected so far with each call to Phases::add methods.
    std::deque<GenericPhase> genericphases;

    /// The element symbols collected so far with each call to Phases::add methods.
    Strings elements;

private:
    /// Collect the element symbols in a GenericPhase or GenericPhases object.
    template<typename Arg>
    auto collectElements(const Arg& phase) -> void
    {
        static_assert(std::is_same_v<Arg, GenericPhase>);
        static_assert(std::is_same_v<Arg, GenericPhases>);

        if(phase.elements().size())
            elements = merge(elements, phase.elements());

        if(phase.species().size())
            for(auto&& s : engine.database().species().withNames(phase.species()))
                elements = merge(elements, s.elements().symbols());
    }

    /// Collect the element symbols in one or more GenericPhase or GenericPhases objects.
    template<typename Arg, typename... Args>
    auto collectElements(const Arg& phase, const Args&... phases) -> void
    {
        collectElements(phase);
        collectElements(phases...);
    }

    /// Add a GenericPhase object into the Phases container.
    auto addPhases(const GenericPhase& phase) -> void
    {
        genericphases.push_back(phase);
    }

    /// Add all GenericPhase objects into the Phases container.
    auto addPhases(const Vec<GenericPhase>& phases) -> void
    {
        genericphases.insert(genericphases.end(), phases.begin(), phases.end());
    }

    /// Add all underlying GenericPhase objects, in given GenericPhases object, into the Phases container.
    auto addPhases(const GenericPhases& phases) -> void
    {
        addPhases(phases.convert(engine, elements));
    }

    /// Add one or more GenericPhase or GenericPhases objects into the Phases container.
    template<typename Arg, typename... Args>
    auto addPhases(const Arg& phase, const Args&... phases) -> void
    {
        static_assert(std::is_same_v<Arg, GenericPhase>);
        static_assert(std::is_same_v<Arg, GenericPhases>);

        addPhases(phase);
        addPhases(phases...);
    }

    /// Replace duplicate phase names with unique names.
    auto fixDuplicatePhaseNames() -> void
    {
        Strings phasenames = vectorize(genericphases, lambda(x, x.name()));
        phasenames = makeunique(phasenames, "!");
        auto i = 0;
        for(auto& phase : genericphases)
            phase.setName(phasenames[i++]);
    }
};

/// The class used to configure an aqueous solution phase.
class AqueousSolution : public GenericPhase
{
public:
    /// Construct a default AqueousSolution object.
    AqueousSolution() : GenericPhase() { initialize(); }

    /// Construct a AqueousSolution object with given species names.
    explicit AqueousSolution(const StringList& species) : GenericPhase(species) { initialize(); }

    /// Construct a AqueousSolution object with given element symbols.
    explicit AqueousSolution(Speciate elements) : GenericPhase(elements += {"H", "O"}) { initialize(); }

    /// Initialize the default attributes of this AqueousSolution object.
    auto initialize() -> void
    {
        setName("AqueousSolution");
        setStateOfMatter(StateOfMatter::Liquid);
        setAggregateState(AggregateState::Aqueous);
        setActivityModel(ActivityModelIdealAqueousSolution());
    }
};

/// The class used to configure a gaseous solution phase.
class GaseousSolution : public GenericPhase
{
public:
    /// Construct a default GaseousSolution object.
    GaseousSolution() : GenericPhase() { initialize(); }

    /// Construct a GaseousSolution object with given species names.
    explicit GaseousSolution(const StringList& species) : GenericPhase(species) { initialize(); }

    /// Construct a GaseousSolution object with given element symbols.
    explicit GaseousSolution(const Speciate& elements) : GenericPhase(elements) { initialize(); }

    /// Initialize the default attributes of this GaseousSolution object.
    auto initialize() -> void
    {
        setName("GaseousSolution");
        setStateOfMatter(StateOfMatter::Gas);
        setAggregateState(AggregateState::Gas);
        setActivityModel(ActivityModelIdealGaseousSolution());
    }
};

/// The class used to configure a liquid solution phase.
class LiquidSolution : public GenericPhase
{
public:
    /// Construct a default LiquidSolution object.
    LiquidSolution() : GenericPhase() { initialize(); }

    /// Construct a LiquidSolution object with given species names.
    explicit LiquidSolution(const StringList& species) : GenericPhase(species) { initialize(); }

    /// Construct a LiquidSolution object with given element symbols.
    explicit LiquidSolution(const Speciate& elements) : GenericPhase(elements) { initialize(); }

    /// Initialize the default attributes of this LiquidSolution object.
    auto initialize() -> void
    {
        setName("LiquidSolution");
        setStateOfMatter(StateOfMatter::Liquid);
        setAggregateState(AggregateState::Liquid);
        setActivityModel(ActivityModelIdealSolution());
    }
};

/// The class used to configure a solid solution phase.
class SolidSolution : public GenericPhase
{
public:
    /// Construct a SolidSolution object with given species names.
    explicit SolidSolution(const StringList& species) : GenericPhase(species) { initialize(); }

    /// Initialize the default attributes of this SolidSolution object.
    auto initialize() -> void
    {
        String phasename;
        auto i = 0;
        for(auto&& name : species())
            phasename += (i++ == 0) ? name : "-" + name;

        setName("SolidSolution");
        setStateOfMatter(StateOfMatter::Solid);
        setAggregateState(AggregateState::Solid);
        setActivityModel(ActivityModelIdealSolution());
    }
};

/// The class used to configure a pure mineral phase.
class Mineral : public GenericPhase
{
public:
    /// Construct a default Mineral object.
    explicit Mineral(String mineral) : GenericPhase(mineral) { initialize(); }

    /// Initialize the default attributes of this Mineral object.
    auto initialize() -> void
    {
        setName(species().front());
        setStateOfMatter(StateOfMatter::Solid);
        setAggregateState(AggregateState::Solid);
        setActivityModel(ActivityModelIdealSolution());
    }
};

/// The class used to configure automatic selection of pure mineral phases.
class MineralPhases : public GenericPhases
{
public:
    /// Construct a default MineralPhases object.
    MineralPhases() : GenericPhases() { initialize(); }

    /// Construct a MineralPhases object with given species names.
    explicit MineralPhases(const StringList& species) : GenericPhases(species) { initialize(); }

    /// Construct a MineralPhases object with given element symbols.
    explicit MineralPhases(const Speciate& elements) : GenericPhases(elements) { initialize(); }

    /// Initialize the default attributes of this MineralPhases object.
    auto initialize() -> void
    {
        setStateOfMatter(StateOfMatter::Solid);
        setAggregateState(AggregateState::Solid);
        setActivityModel(ActivityModelIdealSolution());
    }
};

} // namespace Reaktoro
