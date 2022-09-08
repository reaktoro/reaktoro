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

#include "ChemicalSystem.hpp"

// C++ includes
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {
namespace detail {

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

    /// The formula matrix of the species in the system with respect to its elements.
    MatrixXd formula_matrix;

    /// The stoichiometric matrix of the reactions in the system with respect to its species.
    MatrixXd stoichiometric_matrix;

    /// The phase index pairs defining the kinetically controlled reacting interfaces in the system.
    Pairs<Index, Index> reacting_phase_interfaces;

    /// Construct a default ChemicalSystem::Impl object.
    Impl()
    {}

    /// Construct a ChemicalSystem::Impl object with given database and phases.
    Impl(const Database& database, const Vec<Phase>& phaselist)
    : Impl(database, phaselist, {})
    {
    }

    /// Construct a ChemicalSystem::Impl object with given database, phases, and reactions.
    Impl(const Database& database, const Vec<Phase>& phaselist, const Vec<Reaction>& reactionlist)
    : database(database), phases(phaselist), reactions(reactionlist)
    {
        errorif(database.species().empty(), "Expecting at least one species in the Database object provided when creating a ChemicalSystem object.");
        errorif(phases.empty(), "Expecting at least one phase when creating a ChemicalSystem object, but none was provided.");

        species = phases.species();
        elements = species.elements();
        formula_matrix = detail::assembleFormulaMatrix(species, elements);
        stoichiometric_matrix = detail::assembleStoichiometricMatrix(reactions, species);
        reacting_phase_interfaces = detail::determineReactingPhaseInterfaces(reactions, phases);

        detail::fixDuplicateNames(phases);
        detail::fixDuplicateNames(species);
        detail::fixDuplicateNames(reactions);
    }
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const Database& database, const Vec<Phase>& phases)
: pimpl(new Impl(database, phases))
{}

ChemicalSystem::ChemicalSystem(const Database& database, const Vec<Phase>& phases, const Vec<Reaction>& reactions)
: pimpl(new Impl(database, phases, reactions))
{}

ChemicalSystem::ChemicalSystem(const Phases& phases)
: pimpl(new Impl(phases.database(), phases))
{}

ChemicalSystem::ChemicalSystem(const Phases& phases, const Reactions& reactions)
: pimpl(new Impl(phases.database(), phases, reactions.convert(ChemicalSystem(phases))))
{}

auto ChemicalSystem::database() const -> const Database&
{
    return pimpl->database;
}

auto ChemicalSystem::element(Index index) const -> const Element&
{
    return elements()[index];
}

auto ChemicalSystem::elements() const -> const ElementList&
{
    return pimpl->elements;
}

auto ChemicalSystem::species(Index index) const -> const Species&
{
    return species()[index];
}

auto ChemicalSystem::species() const -> const SpeciesList&
{
    return pimpl->species;
}

auto ChemicalSystem::phase(Index index) const -> const Phase&
{
    return phases()[index];
}

auto ChemicalSystem::phases() const -> const PhaseList&
{
    return pimpl->phases;
}

auto ChemicalSystem::reaction(Index index) const -> const Reaction&
{
    return reactions()[index];
}

auto ChemicalSystem::reactions() const -> const ReactionList&
{
    return pimpl->reactions;
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

auto ChemicalSystem::reactingPhaseInterfaces() const -> const Pairs<Index, Index>&
{
    return pimpl->reacting_phase_interfaces;
}

auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&
{
    // const auto& phases = system.phases();
    // const auto& species = system.species();
    // const auto& elements = system.elements();

    // const unsigned num_phases = phases.size();
    // const unsigned bar_size = std::max(unsigned(4), num_phases) * 25;
    // const String bar1(bar_size, '=');
    // const String bar2(bar_size, '-');

    // std::size_t max_size = 0;
    // for(const auto& phase : phases)
    //     max_size = std::max(max_size, phase.species().size());

    // out << bar1 << std::endl;
    // for(const auto& phase : phases)
    //     out << std::setw(25) << std::left << phase.name();
    // out << std::endl;
    // out << bar2 << std::endl;
    // for(unsigned i = 0; ; ++i)
    // {
    //     if(max_size <= i)
    //         break;

    //     for(const auto& phase : phases)
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
