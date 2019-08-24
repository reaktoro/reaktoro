// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

#include "Partition.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>

namespace Reaktoro {
namespace {

// Collect the indices of the species in a collection of phases
auto indicesSpeciesInPhases(const ChemicalSystem& system, const Indices& iphases) -> Indices
{
    Indices indices;
    indices.reserve(system.numSpecies());
    for(Index iphase : iphases)
    {
        const Index ifirst = system.indexFirstSpeciesInPhase(iphase);
        const Index nspecies = system.numSpeciesInPhase(iphase);
        for(unsigned i = 0; i < nspecies; ++i)
            indices.push_back(ifirst + i);
    }
    return indices;
}

// Collect the indices of the species in a collection of phases
auto indicesSpeciesInPhases(const ChemicalSystem& system, const std::vector<std::string>& phases) -> Indices
{
    return indicesSpeciesInPhases(system, system.indicesPhases(phases));
}

} // namespace

struct Partition::Impl
{
    /// The chemical system instance
    ChemicalSystem system;

    /// The universe for the indices of all species
    Indices universe_species;

    /// The universe for the indices of all phases
    Indices universe_phases;


    /// The indices of the phases in the fluid partition
    Indices indices_fluid_phases;

    /// The indices of the species in the fluid partition
    Indices indices_fluid_species;

    /// The indices of the phases in the solid partition
    Indices indices_solid_phases;

    /// The indices of the species in the solid partition
    Indices indices_solid_species;


    /// The indices of the species in the equilibrium partition
    Indices indices_equilibrium_species;

    /// The indices of the elements in the equilibrium partition
    Indices indices_equilibrium_elements;

    /// The indices of the species in the equilibrium-fluid partition
    Indices indices_equilibrium_fluid_species;

    /// The indices of the elements in the equilibrium-fluid partition
    Indices indices_equilibrium_fluid_elements;

    /// The indices of the species in the equilibrium-solid partition
    Indices indices_equilibrium_solid_species;

    /// The indices of the elements in the equilibrium-solid partition
    Indices indices_equilibrium_solid_elements;


    /// The indices of the species in the kinetic partition
    Indices indices_kinetic_species;

    /// The indices of the elements in the kinetic partition
    Indices indices_kinetic_elements;

    /// The indices of the species in the kinetic-fluid partition
    Indices indices_kinetic_fluid_species;

    /// The indices of the elements in the kinetic-fluid partition
    Indices indices_kinetic_fluid_elements;

    /// The indices of the species in the kinetic-solid partition
    Indices indices_kinetic_solid_species;

    /// The indices of the elements in the kinetic-solid partition
    Indices indices_kinetic_solid_elements;


    /// The indices of the species in the inert partition
    Indices indices_inert_species;

    /// The indices of the elements in the inert partition
    Indices indices_inert_elements;

    /// The indices of the species in the inert-fluid partition
    Indices indices_inert_fluid_species;

    /// The indices of the elements in the inert-fluid partition
    Indices indices_inert_fluid_elements;

    /// The indices of the species in the inert-solid partition
    Indices indices_inert_solid_species;

    /// The indices of the elements in the inert-solid partition
    Indices indices_inert_solid_elements;


    /// The formula matrix of the equilibrium partition
    Matrix formula_matrix_equilibrium;

    /// The formula matrix of the equilibrium-fluid partition
    Matrix formula_matrix_equilibrium_fluid;

    /// The formula matrix of the equilibrium-solid partition
    Matrix formula_matrix_equilibrium_solid;


    /// The formula matrix of the kinetic partition
    Matrix formula_matrix_kinetic;

    /// The formula matrix of the kinetic-fluid partition
    Matrix formula_matrix_kinetic_fluid;

    /// The formula matrix of the kinetic-solid partition
    Matrix formula_matrix_kinetic_solid;


    /// The formula matrix of the inert partition
    Matrix formula_matrix_inert;

    /// The formula matrix of the inert-fluid partition
    Matrix formula_matrix_inert_fluid;

    /// The formula matrix of the inert-solid partition
    Matrix formula_matrix_inert_solid;


    Impl()
    {}

    Impl(const ChemicalSystem& system)
    : system(system)
    {
        universe_species = range(system.species().size());
        universe_phases = range(system.phases().size());

        indices_equilibrium_species = universe_species;
        indices_solid_phases = system.indicesSolidPhases();
        indices_fluid_phases = system.indicesFluidPhases();

        finalise();
    }

    auto setEquilibriumSpecies(const Indices& ispecies) -> void
    {
        indices_equilibrium_species = ispecies;
        indices_inert_species       = difference(indices_inert_species, ispecies);
        indices_kinetic_species     = difference(universe_species, unify(ispecies, indices_inert_species));
        finalise();
    }

    auto setKineticSpecies(const Indices& ispecies) -> void
    {
        indices_kinetic_species     = ispecies;
        indices_inert_species       = difference(indices_inert_species, ispecies);
        indices_equilibrium_species = difference(universe_species, unify(ispecies, indices_inert_species));
        finalise();
    }

    auto setInertSpecies(const Indices& ispecies) -> void
    {
        indices_inert_species       = ispecies;
        indices_equilibrium_species = difference(indices_equilibrium_species, ispecies);
        indices_kinetic_species     = difference(indices_kinetic_species, ispecies);
        finalise();
    }

    auto setFluidPhases(const Indices& indices) -> void
    {
        indices_fluid_phases = indices;
        indices_solid_phases = difference(universe_phases, indices);
        finalise();
    }

    auto setSolidPhases(const Indices& indices) -> void
    {
        indices_solid_phases = indices;
        indices_fluid_phases = difference(universe_phases, indices);
        finalise();
    }

    auto finalise() -> void
    {
        indices_fluid_species = indicesSpeciesInPhases(system, indices_fluid_phases);
        indices_solid_species = indicesSpeciesInPhases(system, indices_solid_phases);

        indices_equilibrium_fluid_species = intersect(indices_equilibrium_species, indices_fluid_species);
        indices_equilibrium_solid_species = intersect(indices_equilibrium_species, indices_solid_species);

        indices_kinetic_fluid_species = intersect(indices_kinetic_species, indices_fluid_species);
        indices_kinetic_solid_species = intersect(indices_kinetic_species, indices_solid_species);

        indices_inert_fluid_species = intersect(indices_inert_species, indices_fluid_species);
        indices_inert_solid_species = intersect(indices_inert_species, indices_solid_species);

        indices_equilibrium_elements       = system.indicesElementsInSpecies(indices_equilibrium_species);
        indices_equilibrium_fluid_elements = system.indicesElementsInSpecies(indices_equilibrium_fluid_species);
        indices_equilibrium_solid_elements = system.indicesElementsInSpecies(indices_equilibrium_solid_species);

        indices_kinetic_elements       = system.indicesElementsInSpecies(indices_kinetic_species);
        indices_kinetic_fluid_elements = system.indicesElementsInSpecies(indices_kinetic_fluid_species);
        indices_kinetic_solid_elements = system.indicesElementsInSpecies(indices_kinetic_solid_species);

        indices_inert_elements       = system.indicesElementsInSpecies(indices_inert_species);
        indices_inert_fluid_elements = system.indicesElementsInSpecies(indices_inert_fluid_species);
        indices_inert_solid_elements = system.indicesElementsInSpecies(indices_inert_solid_species);

        formula_matrix_equilibrium       = submatrix(system.formulaMatrix(), indices_equilibrium_elements, indices_equilibrium_species);
        formula_matrix_equilibrium_fluid = submatrix(system.formulaMatrix(), indices_equilibrium_fluid_elements, indices_equilibrium_fluid_species);
        formula_matrix_equilibrium_solid = submatrix(system.formulaMatrix(), indices_equilibrium_solid_elements, indices_equilibrium_solid_species);

        formula_matrix_kinetic       = submatrix(system.formulaMatrix(), indices_kinetic_elements, indices_kinetic_species);
        formula_matrix_kinetic_fluid = submatrix(system.formulaMatrix(), indices_kinetic_fluid_elements, indices_kinetic_fluid_species);
        formula_matrix_kinetic_solid = submatrix(system.formulaMatrix(), indices_kinetic_solid_elements, indices_kinetic_solid_species);

        formula_matrix_inert       = submatrix(system.formulaMatrix(), indices_inert_elements, indices_inert_species);
        formula_matrix_inert_fluid = submatrix(system.formulaMatrix(), indices_inert_fluid_elements, indices_inert_fluid_species);
        formula_matrix_inert_solid = submatrix(system.formulaMatrix(), indices_inert_solid_elements, indices_inert_solid_species);
    }
};

Partition::Partition()
: pimpl(new Impl())
{}

Partition::Partition(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

auto Partition::setEquilibriumSpecies(const Indices& ispecies) -> void
{
    pimpl->setEquilibriumSpecies(ispecies);
}

auto Partition::setEquilibriumSpecies(const std::vector<std::string>& species) -> void
{
    pimpl->setEquilibriumSpecies(pimpl->system.indicesSpecies(species));
}

auto Partition::setEquilibriumPhases(const Indices& iphases) -> void
{
    pimpl->setEquilibriumSpecies(indicesSpeciesInPhases(pimpl->system, iphases));
}

auto Partition::setEquilibriumPhases(const std::vector<std::string>& phases) -> void
{
    pimpl->setEquilibriumSpecies(indicesSpeciesInPhases(pimpl->system, phases));
}

auto Partition::setKineticSpecies(const Indices& ispecies) -> void
{
    pimpl->setKineticSpecies(ispecies);
}

auto Partition::setKineticSpecies(const std::vector<std::string>& species) -> void
{
    pimpl->setKineticSpecies(pimpl->system.indicesSpecies(species));
}

auto Partition::setKineticPhases(const Indices& iphases) -> void
{
    pimpl->setKineticSpecies(indicesSpeciesInPhases(pimpl->system, iphases));
}

auto Partition::setKineticPhases(const std::vector<std::string>& phases) -> void
{
    pimpl->setKineticSpecies(indicesSpeciesInPhases(pimpl->system, phases));
}

auto Partition::setInertSpecies(const Indices& ispecies) -> void
{
    pimpl->setInertSpecies(ispecies);
}

auto Partition::setInertSpecies(const std::vector<std::string>& species) -> void
{
    pimpl->setInertSpecies(pimpl->system.indicesSpecies(species));
}

auto Partition::setInertPhases(const Indices& iphases) -> void
{
    pimpl->setInertSpecies(indicesSpeciesInPhases(pimpl->system, iphases));
}

auto Partition::setInertPhases(const std::vector<std::string>& phases) -> void
{
    pimpl->setInertSpecies(indicesSpeciesInPhases(pimpl->system, phases));
}

auto Partition::setFluidPhases(const Indices& indices) -> void
{
    pimpl->setFluidPhases(indices);
}

auto Partition::setFluidPhases(const std::vector<std::string>& names) -> void
{
    pimpl->setFluidPhases(pimpl->system.indicesPhases(names));
}

auto Partition::setSolidPhases(const Indices& indices) -> void
{
    pimpl->setSolidPhases(indices);
}

auto Partition::setSolidPhases(const std::vector<std::string>& names) -> void
{
    pimpl->setSolidPhases(pimpl->system.indicesPhases(names));
}

auto Partition::system() const -> const ChemicalSystem&
{
    return pimpl->system;
}

auto Partition::numFluidPhases() const -> unsigned
{
    return pimpl->indices_fluid_phases.size();
}

auto Partition::numFluidSpecies() const -> unsigned
{
    return pimpl->indices_fluid_species.size();
}

auto Partition::numSolidPhases() const -> unsigned
{
    return pimpl->indices_solid_phases.size();
}

auto Partition::numSolidSpecies() const -> unsigned
{
    return pimpl->indices_solid_species.size();
}

auto Partition::numEquilibriumSpecies() const -> unsigned
{
    return pimpl->indices_equilibrium_species.size();
}

auto Partition::numEquilibriumFluidSpecies() const -> unsigned
{
    return pimpl->indices_equilibrium_fluid_species.size();
}

auto Partition::numEquilibriumSolidSpecies() const -> unsigned
{
    return pimpl->indices_equilibrium_solid_species.size();
}

auto Partition::numKineticSpecies() const -> unsigned
{
    return pimpl->indices_kinetic_species.size();
}

auto Partition::numKineticFluidSpecies() const -> unsigned
{
    return pimpl->indices_kinetic_fluid_species.size();
}

auto Partition::numKineticSolidSpecies() const -> unsigned
{
    return pimpl->indices_kinetic_solid_species.size();
}

auto Partition::numInertSpecies() const -> unsigned
{
    return pimpl->indices_inert_species.size();
}

auto Partition::numInertFluidSpecies() const -> unsigned
{
    return pimpl->indices_inert_fluid_species.size();
}

auto Partition::numInertSolidSpecies() const -> unsigned
{
    return pimpl->indices_inert_solid_species.size();
}

auto Partition::numEquilibriumElements() const -> unsigned
{
    return pimpl->indices_equilibrium_elements.size();
}

auto Partition::numEquilibriumFluidElements() const -> unsigned
{
    return pimpl->indices_equilibrium_fluid_elements.size();
}

auto Partition::numEquilibriumSolidElements() const -> unsigned
{
    return pimpl->indices_equilibrium_solid_elements.size();
}

auto Partition::numKineticElements() const -> unsigned
{
    return pimpl->indices_kinetic_elements.size();
}

auto Partition::numKineticFluidElements() const -> unsigned
{
    return pimpl->indices_kinetic_fluid_elements.size();
}

auto Partition::numKineticSolidElements() const -> unsigned
{
    return pimpl->indices_kinetic_solid_elements.size();
}

auto Partition::numInertElements() const -> unsigned
{
    return pimpl->indices_inert_elements.size();
}

auto Partition::numInertFluidElements() const -> unsigned
{
    return pimpl->indices_inert_fluid_elements.size();
}

auto Partition::numInertSolidElements() const -> unsigned
{
    return pimpl->indices_inert_solid_elements.size();
}

auto Partition::indicesFluidPhases() const -> const Indices&
{
    return pimpl->indices_fluid_phases;
}

auto Partition::indicesFluidSpecies() const -> const Indices&
{
    return pimpl->indices_fluid_species;
}

auto Partition::indicesSolidPhases() const -> const Indices&
{
    return pimpl->indices_solid_phases;
}

auto Partition::indicesSolidSpecies() const -> const Indices&
{
    return pimpl->indices_solid_species;
}

auto Partition::indicesEquilibriumSpecies() const -> const Indices&
{
    return pimpl->indices_equilibrium_species;
}

auto Partition::indicesEquilibriumFluidSpecies() const -> const Indices&
{
    return pimpl->indices_equilibrium_fluid_species;
}

auto Partition::indicesEquilibriumSolidSpecies() const -> const Indices&
{
    return pimpl->indices_equilibrium_solid_species;
}

auto Partition::indicesKineticSpecies() const -> const Indices&
{
    return pimpl->indices_kinetic_species;
}

auto Partition::indicesKineticFluidSpecies() const -> const Indices&
{
    return pimpl->indices_kinetic_fluid_species;
}

auto Partition::indicesKineticSolidSpecies() const -> const Indices&
{
    return pimpl->indices_kinetic_solid_species;
}

auto Partition::indicesInertSpecies() const -> const Indices&
{
    return pimpl->indices_inert_species;
}

auto Partition::indicesInertFluidSpecies() const -> const Indices&
{
    return pimpl->indices_inert_fluid_species;
}

auto Partition::indicesInertSolidSpecies() const -> const Indices&
{
    return pimpl->indices_inert_solid_species;
}

auto Partition::indicesEquilibriumElements() const -> const Indices&
{
    return pimpl->indices_equilibrium_elements;
}

auto Partition::indicesEquilibriumFluidElements() const -> const Indices&
{
    return pimpl->indices_equilibrium_fluid_elements;
}

auto Partition::indicesEquilibriumSolidElements() const -> const Indices&
{
    return pimpl->indices_equilibrium_solid_elements;
}

auto Partition::indicesKineticElements() const -> const Indices&
{
    return pimpl->indices_kinetic_elements;
}

auto Partition::indicesKineticFluidElements() const -> const Indices&
{
    return pimpl->indices_kinetic_fluid_elements;
}

auto Partition::indicesKineticSolidElements() const -> const Indices&
{
    return pimpl->indices_kinetic_solid_elements;
}

auto Partition::indicesInertElements() const -> const Indices&
{
    return pimpl->indices_inert_elements;
}

auto Partition::indicesInertFluidElements() const -> const Indices&
{
    return pimpl->indices_inert_fluid_elements;
}

auto Partition::indicesInertSolidElements() const -> const Indices&
{
    return pimpl->indices_inert_solid_elements;
}

auto Partition::formulaMatrixEquilibriumPartition() const -> MatrixConstRef
{
    return pimpl->formula_matrix_equilibrium;
}

auto Partition::formulaMatrixEquilibriumFluidPartition() const -> MatrixConstRef
{
    return pimpl->formula_matrix_equilibrium_fluid;
}

auto Partition::formulaMatrixEquilibriumSolidPartition() const -> MatrixConstRef
{
    return pimpl->formula_matrix_equilibrium_solid;
}

auto Partition::formulaMatrixKineticPartition() const -> MatrixConstRef
{
    return pimpl->formula_matrix_kinetic;
}

auto Partition::formulaMatrixKineticFluidPartition() const -> MatrixConstRef
{
    return pimpl->formula_matrix_kinetic_fluid;
}

auto Partition::formulaMatrixKineticSolidPartition() const -> MatrixConstRef
{
    return pimpl->formula_matrix_kinetic_solid;
}

auto Partition::formulaMatrixInertPartition() const -> MatrixConstRef
{
    return pimpl->formula_matrix_inert;
}

auto Partition::formulaMatrixInertFluidPartition() const -> MatrixConstRef
{
    return pimpl->formula_matrix_inert_fluid;
}

auto Partition::formulaMatrixInertSolidPartition() const -> MatrixConstRef
{
    return pimpl->formula_matrix_inert_solid;
}

bool Partition::operator==(const Partition& partition) const
{
    return partition.pimpl == pimpl;
}

bool Partition::operator!=(const Partition& partition) const
{
    return partition.pimpl != pimpl;
}

} // namespace Reaktoro
