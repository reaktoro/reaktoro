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
#include <iostream>
#include <iomanip>
#include <set>

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
#include <Reaktoro/Core/Utils.hpp>

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

auto collectElements(const std::vector<Species>& species) -> std::vector<Element>
{
    std::set<Element> elements;
    for(const Species& iter : species)
        for(const auto& pair : iter.elements())
            elements.insert(pair.first);
    return std::vector<Element>(elements.begin(), elements.end());
}

auto collectSpecies(const std::vector<Phase>& phases) -> std::vector<Species>
{
    unsigned num_species = 0;
    for(const Phase& phase : phases)
        num_species += phase.species().size();

    std::vector<Species> list;
    list.reserve(num_species);
    for(const Phase& phase : phases)
        for(const Species& iter : phase.species())
            list.push_back(iter);
    return list;
}

} // namespace detail

struct ChemicalSystem::Impl
{
    /// The list of phases in the system
    std::vector<Phase> phases;

    /// The list of species in the system
    std::vector<Species> species;

    /// The list of elements in the system
    std::vector<Element> elements;

    // /// The thermodynamic model of the system
    // ThermoModel thermo_model;

    // /// The chemical model of the system
    // ChemicalModel chemical_model;

    /// The formula matrix of the system
    MatrixXd formula_matrix;

    Impl()
    {}

    Impl(const std::vector<Phase>& phaselist)
    {
        initializePhasesSpeciesElements(phaselist);
        initializeFormulaMatrix();
    }

    auto initializePhasesSpeciesElements(const std::vector<Phase>& phaselist) -> void
    {
        // phases = fixDuplicatedSpeciesNames(phaselist);
        species = detail::collectSpecies(phases);
        elements = detail::collectElements(species);
    }

    auto initializeFormulaMatrix() -> void
    {
        const auto num_elements = elements.size();
        const auto num_species = species.size();
        formula_matrix.resize(num_elements, num_species);
        for(auto i = 0; i < num_species; ++i)
            for(auto j = 0; j < num_elements; ++j)
                formula_matrix(j, i) = species[i].elements().coefficient(elements[j].symbol());
    }
};

ChemicalSystem::ChemicalSystem()
: pimpl(new Impl())
{}

ChemicalSystem::ChemicalSystem(const std::vector<Phase>& phases)
: pimpl(new Impl(phases))
{}

auto ChemicalSystem::numElements() const -> unsigned
{
    return elements().size();
}

auto ChemicalSystem::numSpecies() const -> unsigned
{
    return species().size();
}

auto ChemicalSystem::numSpeciesInPhase(Index iphase) const -> unsigned
{
    return phase(iphase).species().size();
}

auto ChemicalSystem::numPhases() const -> unsigned
{
    return phases().size();
}

auto ChemicalSystem::element(Index index) const -> const Element&
{
    return elements()[index];
}

auto ChemicalSystem::element(std::string name) const -> const Element&
{
    return element(indexElementWithError(name));
}

auto ChemicalSystem::elements() const -> const std::vector<Element>&
{
    return pimpl->elements;
}

auto ChemicalSystem::species(Index index) const -> const Species&
{
    return species()[index];
}

auto ChemicalSystem::species(std::string name) const -> const Species&
{
    return species(indexSpeciesWithError(name));
}

auto ChemicalSystem::species() const -> const std::vector<Species>&
{
    return pimpl->species;
}

auto ChemicalSystem::phase(Index index) const -> const Phase&
{
    Assert(index < numPhases(),
        "Could not get a reference to a Phase instance with given index.",
        "The given index " + std::to_string(index) + " is out of bounds.")
    return phases()[index];
}

auto ChemicalSystem::phase(std::string name) const -> const Phase&
{
    return phase(indexPhaseWithError(name));
}

auto ChemicalSystem::phases() const -> const std::vector<Phase>&
{
    return pimpl->phases;
}

auto ChemicalSystem::formulaMatrix() const -> MatrixXdConstRef
{
    return pimpl->formula_matrix;
}

auto ChemicalSystem::indexElement(std::string name) const -> Index
{
    return indexfn(elements(), RKT_LAMBDA(x, x.name() == name));
}

auto ChemicalSystem::indexElementWithError(std::string name) const -> Index
{
    const Index index = indexElement(name);

    Assert(index < numElements(),
        "Could not get the index of element `" + name + "`.",
        "There is no element called `" + name + "` in the system.");

    return index;
}

auto ChemicalSystem::indexSpecies(std::string name) const -> Index
{
    return indexfn(species(), RKT_LAMBDA(x, x.name() == name));
}

auto ChemicalSystem::indexSpeciesWithError(std::string name) const -> Index
{
    const Index index = indexSpecies(name);
    Assert(index < numSpecies(),
        "Could not get the index of species `" + name + "`.",
        "There is no species called `" + name + "` in the system.");
    return index;
}

auto ChemicalSystem::indexSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    return indexfn(species(), RKT_LAMBDA(s, contains(names, s.name())));
}

auto ChemicalSystem::indexSpeciesAnyWithError(const std::vector<std::string>& names) const -> Index
{
    const Index index = indexSpeciesAny(names);
    Assert(index < numSpecies(),
        "Could not get the index of the species with "
        "any of the following names `" + join(names, ", ") + "`.",
        "There is no species in the system with any of these names.");
    return index;
}

auto ChemicalSystem::indexPhase(std::string name) const -> Index
{
    return indexfn(phases(), RKT_LAMBDA(x, x.name() == name));
}

auto ChemicalSystem::indexPhaseWithError(std::string name) const -> Index
{
    const Index index = indexPhase(name);

    Assert(index < numPhases(),
        "Could not get the index of phase `" + name + "`.",
        "There is no phase called `" + name + "` in the system.");

    return index;
}

auto ChemicalSystem::indexPhaseWithSpecies(Index index) const -> Index
{
    unsigned counter = 0;
    for(unsigned i = 0; i < numPhases(); ++i)
    {
        counter += numSpeciesInPhase(i);
        if(counter > index) return i;
    }
    return numPhases();
}

auto ChemicalSystem::indexFirstSpeciesInPhase(Index iphase) const -> unsigned
{
    unsigned counter = 0;
    for(unsigned i = 0; i < iphase; ++i)
        counter += phase(i).species().size();
    return counter;
}

auto ChemicalSystem::indicesElements(const std::vector<std::string>& names) const -> Indices
{
    Indices indices;
    indices.reserve(names.size());
    for(const auto& name : names)
        indices.push_back(indexElementWithError(name));
    return indices;
}

auto ChemicalSystem::indicesElementsInSpecies(Index index) const -> Indices
{
    Indices indices;
    for(const auto& pair : species(index).elements())
        indices.push_back(indexElement(pair.first.name()));
    return indices;
}

auto ChemicalSystem::indicesElementsInSpecies(const Indices& ispecies) const -> Indices
{
    std::set<Index> ielements;
    for(const Index& i : ispecies)
    {
        const Indices& indices = indicesElementsInSpecies(i);
        ielements.insert(indices.begin(), indices.end());
    }
    return Indices(ielements.begin(), ielements.end());
}

auto ChemicalSystem::indicesSpecies(const std::vector<std::string>& names) const -> Indices
{
    Indices indices;
    indices.reserve(names.size());
    for(const auto& name : names)
        indices.push_back(indexSpeciesWithError(name));
    return indices;
}

auto ChemicalSystem::indicesSpeciesInPhases(const Indices& indices) const -> Indices
{
    Indices res;
    res.reserve(numSpecies());
    for(Index iphase : indices)
    {
        const Index ifirst = indexFirstSpeciesInPhase(iphase);
        const Index nspecies = numSpeciesInPhase(iphase);
        for(unsigned i = 0; i < nspecies; ++i)
            res.push_back(ifirst + i);
    }
    return res;
}

auto ChemicalSystem::indicesPhases(const std::vector<std::string>& names) const -> Indices
{
    Indices indices;
    indices.reserve(names.size());
    for(const auto& name : names)
        indices.push_back(indexPhaseWithError(name));
    return indices;
}

auto ChemicalSystem::indicesPhasesWithSpecies(const Indices& ispecies) const -> Indices
{
    std::set<Index> iphases;
    for(Index i : ispecies)
        iphases.insert(indexPhaseWithSpecies(i));
    return Indices(iphases.begin(), iphases.end());
}

auto ChemicalSystem::indicesFluidPhases() const -> Indices
{
    Indices indices;
    indices.reserve(numPhases());
    for(Index i = 0; i < numPhases(); ++i)
        if(phase(i).stateOfMatter() != StateOfMatter::Solid)
            indices.push_back(i);
    return indices;
}

auto ChemicalSystem::indicesFluidSpecies() const -> Indices
{
    return indicesSpeciesInPhases(indicesFluidPhases());
}

auto ChemicalSystem::indicesSolidPhases() const -> Indices
{
    Indices indices;
    indices.reserve(numPhases());
    for(Index i = 0; i < numPhases(); ++i)
        if(phase(i).stateOfMatter() == StateOfMatter::Solid)
            indices.push_back(i);
    return indices;
}

auto ChemicalSystem::indicesSolidSpecies() const -> Indices
{
    return indicesSpeciesInPhases(indicesSolidPhases());
}

auto operator<<(std::ostream& out, const ChemicalSystem& system) -> std::ostream&
{
    const auto& phases = system.phases();
    const auto& species = system.species();
    const auto& elements = system.elements();

    const unsigned num_phases = phases.size();
    const unsigned bar_size = std::max(unsigned(4), num_phases) * 25;
    const std::string bar1(bar_size, '=');
    const std::string bar2(bar_size, '-');

    std::size_t max_size = 0;
    for(const auto& phase : phases)
        max_size = std::max(max_size, phase.species().size());

    out << bar1 << std::endl;
    for(const auto& phase : phases)
        out << std::setw(25) << std::left << phase.name();
    out << std::endl;
    out << bar2 << std::endl;
    for(unsigned i = 0; ; ++i)
    {
        if(max_size <= i)
            break;

        for(const auto& phase : phases)
        {
            if(i < phase.species().size())
                out << std::setw(25) << std::left << phase.species(i).name();
            else
                out << std::setw(25) << std::left << "";
        }

        out << std::endl;
    }

    out << bar1 << std::endl;
    out << std::setw(25) << std::left << "Index";
    out << std::setw(25) << std::left << "Species";
    out << std::setw(25) << std::left << "Element";
    out << std::setw(25) << std::left << "Phase";
    out << std::endl;
    out << bar2 << std::endl;

    for(unsigned i = 0; ; ++i)
    {
        if(elements.size() <= i && species.size() <= i && phases.size() <= i)
            break;

        out << std::setw(25) << std::left << i;
        out << std::setw(25) << std::left << (i < species.size() ? species[i].name() : "");
        out << std::setw(25) << std::left << (i < elements.size() ? elements[i].name() : "");
        out << std::setw(25) << std::left << (i < phases.size() ? phases[i].name() : "");
        out << std::endl;
    }

    out << bar1 << std::endl;

    return out;
}

} // namespace Reaktoro
