// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "Partition.hpp"

// C++ includes
#include <set>

// Reaktor includes
#include <Reaktor/Common/SetUtils.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/Utils.hpp>

namespace Reaktor {
namespace {

// Collect the names of the species in a collection of phases
auto speciesNamesInPhases(const ChemicalSystem& system, const std::vector<std::string>& phases) -> std::vector<std::string>
{
    std::vector<std::string> names;
    names.reserve(system.numSpecies());
    for(std::string phase : phases)
        for(Species species : system.phase(phase).species())
            names.push_back(species.name());
    return names;
}

} // namespace

struct Partition::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The universe of indices
    Indices universe;

    /// The indices of the equilibrium species
    Indices indices_equilibrium_species;

    /// The indices of the kinetic species
    Indices indices_kinetic_species;

    /// The indices of the inert species
    Indices indices_inert_species;

    /// The indices of the equilibrium elements
    Indices indices_equilibrium_elements;

    /// The indices of the kinetic elements
    Indices indices_kinetic_elements;

    /// The indices of the inert elements
    Indices indices_inert_elements;

    /// The formula matrix of the equilibrium partition
    Matrix formula_matrix_equilibrium;

    /// The formula matrix of the kinetic partition
    Matrix formula_matrix_kinetic;

    /// The formula matrix of the inert partition
    Matrix formula_matrix_inert;

    Impl()
    {}

    Impl(const ChemicalSystem& system)
    : system(system)
    {
        universe = range(system.species().size());

        setEquilibriumSpecies(universe);
    }

    auto setEquilibriumSpecies(const Indices& ispecies) -> void
    {
        indices_equilibrium_species = ispecies;
        indices_inert_species       = difference(indices_inert_species, ispecies);
        indices_kinetic_species     = difference(universe, unify(ispecies, indices_inert_species));
        finalise();
    }

    auto setEquilibriumSpecies(const std::vector<std::string>& species) -> void
    {
        setEquilibriumSpecies(system.indicesSpecies(species));
    }

    auto setEquilibriumPhases(const std::vector<std::string>& phases) -> void
    {
        setEquilibriumSpecies(speciesNamesInPhases(system, phases));
    }

    auto setKineticSpecies(const Indices& ispecies) -> void
    {
        indices_kinetic_species     = ispecies;
        indices_inert_species       = difference(indices_inert_species, ispecies);
        indices_equilibrium_species = difference(universe, unify(ispecies, indices_inert_species));
        finalise();
    }

    auto setKineticSpecies(const std::vector<std::string>& species) -> void
    {
        setKineticSpecies(system.indicesSpecies(species));
    }

    auto setKineticPhases(const std::vector<std::string>& phases) -> void
    {
        setKineticSpecies(speciesNamesInPhases(system, phases));
    }

    auto setInertSpecies(const Indices& ispecies) -> void
    {
        indices_inert_species       = ispecies;
        indices_equilibrium_species = difference(indices_equilibrium_species, ispecies);
        indices_kinetic_species     = difference(indices_kinetic_species, ispecies);
        finalise();
    }

    auto setInertSpecies(const std::vector<std::string>& species) -> void
    {
        setInertSpecies(system.indicesSpecies(species));
    }

    auto setInertPhases(const std::vector<std::string>& phases) -> void
    {
        setInertSpecies(speciesNamesInPhases(system, phases));
    }

    auto finalise() -> void
    {
        indices_equilibrium_elements = system.indicesElementsInSpecies(indices_equilibrium_species);
        indices_kinetic_elements = system.indicesElementsInSpecies(indices_kinetic_species);
        indices_inert_elements = system.indicesElementsInSpecies(indices_inert_species);

        formula_matrix_equilibrium = submatrix(system.formulaMatrix(), indices_equilibrium_elements, indices_equilibrium_species);
        formula_matrix_kinetic = submatrix(system.formulaMatrix(), indices_kinetic_elements, indices_kinetic_species);
        formula_matrix_inert = submatrix(system.formulaMatrix(), indices_inert_elements, indices_inert_species);
    }
};

Partition::Partition()
: pimpl(new Impl())
{}

Partition::Partition(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

Partition::Partition(const Partition& other)
: pimpl(new Impl(*other.pimpl))
{}

Partition::~Partition()
{}

auto Partition::operator=(Partition other) -> Partition&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Partition::setEquilibriumSpecies(const Indices& ispecies) -> void
{
    pimpl->setEquilibriumSpecies(ispecies);
}

auto Partition::setEquilibriumSpecies(const std::vector<std::string>& species) -> void
{
    pimpl->setEquilibriumSpecies(species);
}

auto Partition::setEquilibriumPhases(const std::vector<std::string>& phases) -> void
{
    pimpl->setEquilibriumPhases(phases);
}

auto Partition::setKineticSpecies(const Indices& ispecies) -> void
{
    pimpl->setKineticSpecies(ispecies);
}

auto Partition::setKineticSpecies(const std::vector<std::string>& species) -> void
{
    pimpl->setKineticSpecies(species);
}

auto Partition::setKineticPhases(const std::vector<std::string>& phases) -> void
{
    pimpl->setKineticPhases(phases);
}

auto Partition::setInertSpecies(const Indices& ispecies) -> void
{
    pimpl->setInertSpecies(ispecies);
}

auto Partition::setInertSpecies(const std::vector<std::string>& species) -> void
{
    pimpl->setInertSpecies(species);
}

auto Partition::setInertPhases(const std::vector<std::string>& phases) -> void
{
    pimpl->setInertPhases(phases);
}

auto Partition::numEquilibriumSpecies() const -> unsigned
{
    return pimpl->indices_equilibrium_species.size();
}

auto Partition::numKineticSpecies() const -> unsigned
{
    return pimpl->indices_kinetic_species.size();
}

auto Partition::numInertSpecies() const -> unsigned
{
    return pimpl->indices_inert_species.size();
}

auto Partition::numEquilibriumElements() const -> unsigned
{
    return pimpl->indices_equilibrium_elements.size();
}

auto Partition::numKineticElements() const -> unsigned
{
    return pimpl->indices_kinetic_elements.size();
}

auto Partition::numInertElements() const -> unsigned
{
    return pimpl->indices_inert_elements.size();
}

auto Partition::indicesEquilibriumSpecies() const -> const Indices&
{
    return pimpl->indices_equilibrium_species;
}

auto Partition::indicesKineticSpecies() const -> const Indices&
{
    return pimpl->indices_kinetic_species;
}

auto Partition::indicesInertSpecies() const -> const Indices&
{
    return pimpl->indices_inert_species;
}

auto Partition::indicesEquilibriumElements() const -> const Indices&
{
    return pimpl->indices_equilibrium_elements;
}

auto Partition::indicesKineticElements() const -> const Indices&
{
    return pimpl->indices_kinetic_elements;
}

auto Partition::indicesInertElements() const -> const Indices&
{
    return pimpl->indices_inert_elements;
}

auto Partition::formulaMatrixEquilibriumSpecies() const -> const Matrix&
{
    return pimpl->formula_matrix_equilibrium;
}

auto Partition::formulaMatrixKineticSpecies() const -> const Matrix&
{
    return pimpl->formula_matrix_kinetic;
}

auto Partition::formulaMatrixInertSpecies() const -> const Matrix&
{
    return pimpl->formula_matrix_inert;
}

} // namespace Reaktor
