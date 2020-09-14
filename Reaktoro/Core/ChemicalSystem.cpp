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

#include "ChemicalSystem.hpp"

// C++ includes
#include <cassert>
#include <iostream>
#include <iomanip>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {
namespace detail {

/// Replace duplicate phase names with unique names.
auto fixDuplicatePhaseNames(Vec<Phase>& phases)
{
    Strings phasenames = vectorize(phases, RKT_LAMBDA(x, x.name()));
    phasenames = makeunique(phasenames, "!");
    for(auto i = 0; i < phases.size(); ++i)
        if(phases[i].name() != phasenames[i])
            phases[i] = phases[i].withName(phasenames[i]);
}

/// Replace duplicate species names with unique names.
auto fixDuplicateSpeciesNames(Vec<Species>& species)
{
    Strings speciesnames = vectorize(species, RKT_LAMBDA(x, x.name()));
    speciesnames = makeunique(speciesnames, "!");
    for(auto i = 0; i < species.size(); ++i)
        if(species[i].name() != speciesnames[i])
            species[i] = species[i].withName(speciesnames[i]);
}

/// Collect all Element objects in the given Species objects.
auto collectElements(const Vec<Species>& species) -> Vec<Element>
{
    Map<String, Element> collected;
    for(const auto& s : species)
        for(const auto& [element, coeff] : s.elements())
            collected.emplace(element.symbol(), element);
    auto elements = vectorize(collected, RKT_LAMBDA(x, x.second));
    std::sort(elements.begin(), elements.end(),
        [](auto l, auto r) { return l.molarMass() < r.molarMass(); }); // sort in ascending order of molar mass (atomic weight)
    return elements;
}

/// Collect all Species objects in the given Phase objects.
auto collectSpecies(const Vec<Phase>& phases) -> Vec<Species>
{
    auto num_species = 0;
    for(const auto& phase : phases)
        num_species += phase.species().size();

    Vec<Species> species;
    species.reserve(num_species);
    for(const auto& phase : phases)
        species = concatenate(species, phase.species().data());
    return species;
}

/// Return the formula matrix of the species with respect to given elements.
auto formulaMatrix(const Vec<Species>& species, const Vec<Element>& elements) -> MatrixXd
{
    const auto num_elements = elements.size();
    const auto num_components = num_elements + 1;
    const auto num_species = species.size();
    MatrixXd A(num_components, num_species);
    for(auto i = 0; i < num_species; ++i)
        for(auto j = 0; j < num_elements; ++j)
            A(j, i) = species[i].elements().coefficient(elements[j].symbol());
    for(auto i = 0; i < num_species; ++i)
        A(num_elements, i) = species[i].charge();
    return A;
}

} // namespace detail

struct ChemicalSystem::Impl
{
    /// The database used to construct the chemical system.
    Database database;

    /// The list of phases in the system.
    Vec<Phase> phases;

    /// The list of species in the system.
    Vec<Species> species;

    /// The list of elements in the system.
    Vec<Element> elements;

    /// The formula matrix of the system.
    MatrixXd formula_matrix;

    /// Construct a default ChemicalSystem::Impl object.
    Impl()
    {}

    /// Construct a ChemicalSystem::Impl object with given phases.
    Impl(const Database& database, const Vec<Phase>& phaselist)
    : database(database), phases(phaselist)
    {
        species = detail::collectSpecies(phases);
        elements = detail::collectElements(species);
        formula_matrix = detail::formulaMatrix(species, elements);

        detail::fixDuplicatePhaseNames(phases);
        detail::fixDuplicateSpeciesNames(species);
    }
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const Database& database, const Vec<Phase>& phases)
: pimpl(new Impl(database, phases))
{}

auto ChemicalSystem::database() const -> const Database&
{
    return pimpl->database;
}

auto ChemicalSystem::element(Index index) const -> const Element&
{
    return elements()[index];
}

auto ChemicalSystem::elements() const -> ElementListConstRef
{
    return pimpl->elements;
}

auto ChemicalSystem::species(Index index) const -> const Species&
{
    return species()[index];
}

auto ChemicalSystem::species() const -> SpeciesListConstRef
{
    return pimpl->species;
}

auto ChemicalSystem::phase(Index index) const -> const Phase&
{
    return phases()[index];
}

auto ChemicalSystem::phases() const -> PhaseListConstRef
{
    return pimpl->phases;
}

auto ChemicalSystem::formulaMatrix() const -> MatrixXdConstRef
{
    return pimpl->formula_matrix;
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
