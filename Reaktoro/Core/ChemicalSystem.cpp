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

#include "ChemicalSystem.hpp"

// C++ includes
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {
namespace detail {

auto computeChemicalSystemID() -> Index
{
    static thread_local Index counter = 0;
    return counter++;
}

template<typename NamedObjects>
auto fixDuplicateNames(NamedObjects& objects)
{
    Strings names = vectorize(objects, RKT_LAMBDA(x, x.name()));
    names = makeunique(names, "!");
    for(auto i = 0; i < objects.size(); ++i)
        if(objects[i].name() != names[i])
            objects[i] = objects[i].withName(names[i]);
}

} // namespace detail

struct ChemicalSystem::Impl
{
    /// The unique identification number of the chemical system
    Index id = detail::computeChemicalSystemID();

    /// The database used to construct the chemical system.
    Database database;

    /// The list of phases in the system.
    PhaseList phases;

    /// The list of species in the system.
    SpeciesList species;

    /// The list of elements in the system.
    ElementList elements;

    /// The list of reactions in the system.
    ReactionList reactions;

    /// The list of surfaces in the system.
    SurfaceList surfaces;

    /// The formula matrix of the species in the system with respect to its elements.
    MatrixXd formula_matrix;

    /// The stoichiometric matrix of the reactions in the system with respect to its species.
    MatrixXd stoichiometric_matrix;

    /// Construct a default ChemicalSystem::Impl object.
    Impl()
    {}

    /// Construct a ChemicalSystem::Impl object with given database and phases.
    Impl(Database const& database, PhaseList const& phases)
    : Impl(database, phases, ReactionList{}, SurfaceList{})
    {
    }

    /// Construct a ChemicalSystem::Impl object with given database, phases, and surfaces.
    Impl(Database const& database, PhaseList const& phases, SurfaceList const& surfaces)
    : Impl(database, phases, ReactionList{}, surfaces)
    {
    }

    /// Construct a ChemicalSystem::Impl object with given database, phases, and reactions (determine surfaces from reactions and phases).
    Impl(Database const& database, PhaseList const& phases, ReactionList const& reactions)
    : Impl(database, phases, reactions, SurfaceList{})
    {
    }

    /// Construct a ChemicalSystem::Impl object with given database, phases, reactions, and surfaces.
    Impl(Database const& database0, PhaseList const& phases0, ReactionList const& reactions0, SurfaceList const& surfaces0)
    : database(database0), phases(phases0), reactions(reactions0), surfaces(surfaces0)
    {
        errorif(database.species().empty(), "Expecting at least one species in the Database object provided when creating a ChemicalSystem object.");
        errorif(phases.empty(), "Expecting at least one phase when creating a ChemicalSystem object, but none was provided.");

        species = phases.species();
        elements = species.elements();
        formula_matrix = detail::assembleFormulaMatrix(species, elements);
        stoichiometric_matrix = detail::assembleStoichiometricMatrix(reactions, species);

        detail::fixDuplicateNames(phases);
        detail::fixDuplicateNames(species);
        detail::fixDuplicateNames(reactions);
        detail::fixDuplicateNames(surfaces);
    }

    /// Construct a ChemicalSystem::Impl object with given database, phases, reactions, and surfaces.
    Impl(Database const& database, PhaseList const& phases, Reactions const& reactions)
    : Impl(database, phases, ReactionList(reactions.convert({database, phases.species(), phases, SurfaceList{}})))
    {
    }

    /// Construct a ChemicalSystem::Impl object with given database, phases, reactions, and surfaces.
    Impl(Database const& database, PhaseList const& phases, Surfaces const& surfaces)
    : Impl(database, phases, SurfaceList(surfaces.convert(phases)))
    {
    }

    /// Construct a ChemicalSystem::Impl object with given database, phases, reactions, and surfaces.
    Impl(Database const& database, PhaseList const& phases, Reactions const& reactions, Surfaces const& surfaces)
    : Impl(database, phases, reactions, surfaces.convert(phases))
    {
    }

    /// Construct a ChemicalSystem::Impl object with given database, phases, reactions, and surfaces.
    Impl(Database const& database, PhaseList const& phases, Reactions const& reactions, SurfaceList const& surfaces)
    : Impl(database, phases, reactions.convert({database, phases.species(), phases, surfaces}), surfaces)
    {
    }
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(Database const& database, PhaseList const& phases)
: pimpl(new Impl(database, phases))
{}

ChemicalSystem::ChemicalSystem(Database const& database, PhaseList const& phases, SurfaceList const& surfaces)
: pimpl(new Impl(database, phases, surfaces))
{}

ChemicalSystem::ChemicalSystem(Database const& database, PhaseList const& phases, ReactionList const& reactions)
: pimpl(new Impl(database, phases, reactions))
{}

ChemicalSystem::ChemicalSystem(Database const& database, PhaseList const& phases, ReactionList const& reactions, SurfaceList const& surfaces)
: pimpl(new Impl(database, phases, reactions, surfaces))
{}

ChemicalSystem::ChemicalSystem(Phases const& phases)
: pimpl(new Impl(phases.database(), phases.convert()))
{}

ChemicalSystem::ChemicalSystem(Phases const& phases, Surfaces const& surfaces)
: pimpl(new Impl(phases.database(), phases.convert(), surfaces))
{}

ChemicalSystem::ChemicalSystem(Phases const& phases, Reactions const& reactions)
: pimpl(new Impl(phases.database(), phases.convert(), reactions))
{}

ChemicalSystem::ChemicalSystem(Phases const& phases, Reactions const& reactions, Surfaces const& surfaces)
: pimpl(new Impl(phases.database(), phases.convert(), reactions, surfaces))
{}

auto ChemicalSystem::id() const -> Index
{
    return pimpl->id;
}

auto ChemicalSystem::database() const -> Database const&
{
    return pimpl->database;
}

auto ChemicalSystem::element(Index index) const -> Element const&
{
    return pimpl->elements[index];
}

auto ChemicalSystem::elements() const -> ElementList const&
{
    return pimpl->elements;
}

auto ChemicalSystem::species(Index index) const -> Species const&
{
    return pimpl->species[index];
}

auto ChemicalSystem::species() const -> SpeciesList const&
{
    return pimpl->species;
}

auto ChemicalSystem::phase(Index index) const -> Phase const&
{
    return pimpl->phases[index];
}

auto ChemicalSystem::phases() const -> PhaseList const&
{
    return pimpl->phases;
}

auto ChemicalSystem::reaction(Index index) const -> Reaction const&
{
    return pimpl->reactions[index];
}

auto ChemicalSystem::reactions() const -> ReactionList const&
{
    return pimpl->reactions;
}

auto ChemicalSystem::surface(Index index) const -> Surface const&
{
    return pimpl->surfaces[index];
}

auto ChemicalSystem::surfaces() const -> SurfaceList const&
{
    return pimpl->surfaces;
}

auto ChemicalSystem::formulaMatrix() const -> MatrixXdConstRef
{
    return pimpl->formula_matrix;
}

auto ChemicalSystem::formulaMatrixElements() const -> MatrixXdConstRef
{
    return pimpl->formula_matrix.topRows(pimpl->elements.size());
}

auto ChemicalSystem::formulaMatrixCharge() const -> MatrixXdConstRef
{
    return pimpl->formula_matrix.bottomRows(1);
}

auto ChemicalSystem::stoichiometricMatrix() const -> MatrixXdConstRef
{
    return pimpl->stoichiometric_matrix;
}

auto operator<<(std::ostream& out, ChemicalSystem const& system) -> std::ostream&
{
    // auto const& phases = system.phases();
    // auto const& species = system.species();
    // auto const& elements = system.elements();

    // const unsigned num_phases = phases.size();
    // const unsigned bar_size = std::max(unsigned(4), num_phases) * 25;
    // const String bar1(bar_size, '=');
    // const String bar2(bar_size, '-');

    // std::size_t max_size = 0;
    // for(auto const& phase : phases)
    //     max_size = std::max(max_size, phase.species().size());

    // out << bar1 << std::endl;
    // for(auto const& phase : phases)
    //     out << std::setw(25) << std::left << phase.name();
    // out << std::endl;
    // out << bar2 << std::endl;
    // for(unsigned i = 0; ; ++i)
    // {
    //     if(max_size <= i)
    //         break;

    //     for(auto const& phase : phases)
    //     {
    //         if(i < phase.species().size())
    //             out << std::setw(25) << std::left << phase.species(i).name();
    //         else
    //             out << std::setw(25) << std::left << "";
    //     }

    //     out << std::endl;
    // }

    // out << bar1 << std::endl;
    // out << std::setw(25) << std::left << "Index";
    // out << std::setw(25) << std::left << "Species";
    // out << std::setw(25) << std::left << "Element";
    // out << std::setw(25) << std::left << "Phase";
    // out << std::endl;
    // out << bar2 << std::endl;

    // for(unsigned i = 0; ; ++i)
    // {
    //     if(elements.size() <= i && species.size() <= i && phases.size() <= i)
    //         break;

    //     out << std::setw(25) << std::left << i;
    //     out << std::setw(25) << std::left << (i < species.size() ? species[i].name() : "");
    //     out << std::setw(25) << std::left << (i < elements.size() ? elements[i].name() : "");
    //     out << std::setw(25) << std::left << (i < phases.size() ? phases[i].name() : "");
    //     out << std::endl;
    // }

    // out << bar1 << std::endl;

    return out;
}

} // namespace Reaktoro
