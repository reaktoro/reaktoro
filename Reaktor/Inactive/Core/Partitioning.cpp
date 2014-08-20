/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "Partitioning.hpp"

// C++ includes
#include <set>

// Reaktor includes
#include <Reaktor/Common/Exception.hpp>
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Core/ReactionSystem.hpp>
#include <Reaktor/Utils/MatrixUtils.hpp>
#include <Reaktor/Utils/SetUtils.hpp>
#include <Reaktor/Utils/StringUtils.hpp>

namespace Reaktor {
namespace internal {

void uninitialisedPartitionError()
{
    Exception exception;
    exception.error << "Cannot set the partitioning instance.";
    exception.reason << "The partitioning instance has not been initialized with a chemical system instance.";
    raise(exception);
}

void unknownSpeciesError(const std::string& species)
{
    Exception exception;
    exception.error << "Cannot set the partitioning instance.";
    exception.reason << "The provided species " << species << " is unknown to the chemical system. Did you forget to include it?";
    raise(exception);
}

void checkExistenceOfSpecies(const std::vector<std::string>& subset, const std::vector<std::string>& set)
{
    for(const auto& species : subset)
    {
        if(not contained(species, set))
        {
            Exception exception;
            exception.error << "Cannot set the partitioning instance.";
            exception.reason << "The provided species " << species << " is unknown to the chemical system. Did you forget to include it?";
            raise(exception);
        }
    }
}

} /* namespace internal */

class Partitioning::Impl
{
public:
    /// The chemical system instance
    ChemicalSystem system$;

    /// The species in the chemical system
    std::vector<std::string> species$;

    /// The equilibrium species in the chemical system
    std::vector<std::string> equilibrium_species$;

    /// The kinetic species in the chemical system
    std::vector<std::string> kinetic_species$;

    /// The inert species in the chemical system
    std::vector<std::string> inert_species$;

    /// The equilibrium elements in the chemical system
    std::vector<std::string> equilibrium_elements$;

    /// The kinetic elements in the chemical system
    std::vector<std::string> kinetic_elements$;

    /// The inert elements in the chemical system
    std::vector<std::string> inert_elements$;

    /// The indices of the equilibrium species
    Indices idx_equilibrium_species$;

    /// The indices of the kinetic species
    Indices idx_kinetic_species$;

    /// The indices of the inert species
    Indices idx_inert_species$;

    /// The indices of the equilibrium elements
    Indices idx_equilibrium_elements$;

    /// The indices of the kinetic elements
    Indices idx_kinetic_elements$;

    /// The indices of the inert elements
    Indices idx_inert_elements$;

    /// The indices of the phases that contains equilibrium species
    Indices idx_phases_with_equilibrium_species$;

    /// The indices of the phases that contains kinetic species
    Indices idx_phases_with_kinetic_species$;

    /// The indices of the phases that contains inert species
    Indices idx_phases_with_inert_species$;

public:
    Impl() = delete;

    Impl(const ChemicalSystem& system)
    : system$(system), species$(system.speciesNames())
    {
        setEquilibriumSpecies(species$);
    }

    ~Impl()
    {}

    auto idxPhasesWithSpecies(const Indices& indices) -> Indices
    {
        std::set<Index> idx_phases;
        for(const Index& idx_species : indices)
            idx_phases.insert(system$.idxPhaseWithSpecies(idx_species));
        return Indices(idx_phases.begin(), idx_phases.end());
    }

    auto idxElementsInPartition(const Indices& idx_species) -> Indices
    {
        std::set<Index> indices;

        for(const Index& i : idx_species)
        {
            const Indices& idx_elements = system$.idxElementsInSpecies(i);
            indices.insert(idx_elements.begin(), idx_elements.end());
        }

        return Indices(indices.begin(), indices.end());
    }

    auto finalise() -> void
    {
        // Initialise the indices of the equilibrium, kinetic and inert elements
        idx_equilibrium_elements$ = idxElementsInPartition(idx_equilibrium_species$);
        idx_kinetic_elements$     = idxElementsInPartition(idx_kinetic_species$);
        idx_inert_elements$       = idxElementsInPartition(idx_inert_species$);

        // Initialise the names of the equilibrium, kinetic and inert elements
        equilibrium_elements$ = extract(system$.elements(), idx_equilibrium_elements$);
        kinetic_elements$     = extract(system$.elements(), idx_kinetic_elements$);
        inert_elements$       = extract(system$.elements(), idx_inert_elements$);

        // Initialise the indices of the phases that contains equilibrium, kinetic and inert elements
        idx_phases_with_equilibrium_species$ = idxPhasesWithSpecies(idx_equilibrium_species$);
        idx_phases_with_kinetic_species$     = idxPhasesWithSpecies(idx_kinetic_species$);
        idx_phases_with_inert_species$       = idxPhasesWithSpecies(idx_inert_species$);
    }

    auto setEquilibriumSpecies(const std::vector<std::string>& species) -> void
    {
        if(species$.empty()) internal::uninitialisedPartitionError();

        internal::checkExistenceOfSpecies(species, species$);

        equilibrium_species$     = species;
        inert_species$           = difference(inert_species$, equilibrium_species$);
        kinetic_species$         = difference(species$, unify(equilibrium_species$, inert_species$));
        idx_equilibrium_species$ = find(equilibrium_species$, species$);
        idx_kinetic_species$     = find(kinetic_species$, species$);
        idx_inert_species$       = find(inert_species$, species$);

        finalise();
    }

    auto setKineticSpecies(const std::vector<std::string>& species) -> void
    {
        if(species$.empty()) internal::uninitialisedPartitionError();

        internal::checkExistenceOfSpecies(species, species$);

        kinetic_species$         = species;
        inert_species$           = difference(inert_species$, kinetic_species$);
        equilibrium_species$     = difference(species$, unify(kinetic_species$, inert_species$));
        idx_equilibrium_species$ = find(equilibrium_species$, species$);
        idx_kinetic_species$     = find(kinetic_species$, species$);
        idx_inert_species$       = find(inert_species$, species$);

        finalise();
    }

    auto setInertSpecies(const std::vector<std::string>& species) -> void
    {
        if(species$.empty()) internal::uninitialisedPartitionError();

        internal::checkExistenceOfSpecies(species, species$);

        inert_species$           = species;
        equilibrium_species$     = difference(equilibrium_species$, inert_species$);
        kinetic_species$         = difference(kinetic_species$, inert_species$);
        idx_equilibrium_species$ = find(equilibrium_species$, species$);
        idx_kinetic_species$     = find(kinetic_species$, species$);
        idx_inert_species$       = find(inert_species$, species$);

        finalise();
    }

    auto numSpecies() const -> unsigned
    {
        return species$.size();
    }

    auto numEquilibriumSpecies() const -> unsigned
    {
        return equilibrium_species$.size();
    }

    auto numKineticSpecies() const -> unsigned
    {
        return kinetic_species$.size();
    }

    auto numInertSpecies() const -> unsigned
    {
        return inert_species$.size();
    }

    auto numEquilibriumElements() const -> unsigned
    {
        return equilibrium_elements$.size();
    }

    auto numKineticElements() const -> unsigned
    {
        return kinetic_elements$.size();
    }

    auto numInertElements() const -> unsigned
    {
        return inert_elements$.size();
    }

    auto numPhasesWithEquilibriumSpecies() const -> unsigned
    {
        return idx_phases_with_equilibrium_species$.size();
    }

    auto numPhasesWithKineticSpecies() const -> unsigned
    {
        return idx_phases_with_kinetic_species$.size();
    }

    auto numPhasesWithInertSpecies() const -> unsigned
    {
        return idx_phases_with_inert_species$.size();
    }

    auto setEquilibriumSpecies(const std::string& names) -> void
    {
        setEquilibriumSpecies(split(names, " "));
    }

    auto setKineticSpecies(const std::string& names) -> void
    {
        setKineticSpecies(split(names, " "));
    }

    auto setInertSpecies(const std::string& names) -> void
    {
        setInertSpecies(split(names, " "));
    }

    auto idxInertSpecies() const -> const Indices&
    {
        return idx_inert_species$;
    }

    auto idxSpecies(const std::string& name) const -> Index
    {
        return find(name, species$);
    }

    auto idxSpecies(const std::vector<std::string>& names) const -> Indices
    {
        return find(names, species$);
    }

    auto idxEquilibriumSpecies(const std::string& name) const -> Index
    {
        return find(name, equilibrium_species$);
    }

    auto idxEquilibriumSpecies(const std::vector<std::string>& names) const -> Indices
    {
        return find(names, equilibrium_species$);
    }

    auto idxKineticSpecies(const std::string& name) const -> Index
    {
        return find(name, kinetic_species$);
    }

    auto idxKineticSpecies(const std::vector<std::string>& names) const -> Indices
    {
        return find(names, kinetic_species$);
    }

    auto idxInertSpecies(const std::string& name) const -> Index
    {
        return find(name, inert_species$);
    }

    auto idxInertSpecies(const std::vector<std::string>& names) const -> Indices
    {
        return find(names, inert_species$);
    }

    auto idxEquilibriumElement(const std::string& element) const -> Index
    {
        return find(element, equilibrium_elements$);
    }

    auto idxEquilibriumElements(const std::vector<std::string>& elements) const -> Indices
    {
        return find(elements, equilibrium_elements$);
    }

    auto idxKineticElement(const std::string& element) const -> Index
    {
        return find(element, kinetic_elements$);
    }

    auto idxKineticElements(const std::vector<std::string>& elements) const -> Indices
    {
        return find(elements, kinetic_elements$);
    }

    auto idxInertElement(const std::string& element) const -> Index
    {
        return find(element, inert_elements$);
    }

    auto idxInertElements(const std::vector<std::string>& elements) const -> Indices
    {
        return find(elements, inert_elements$);
    }

    auto idxPhasesWithEquilibriumSpecies() const -> const Indices&
    {
        return idx_phases_with_equilibrium_species$;
    }

    auto idxPhasesWithKineticSpecies() const -> const Indices&
    {
        return idx_phases_with_kinetic_species$;
    }

    auto idxPhasesWithInertSpecies() const -> const Indices&
    {
        return idx_phases_with_inert_species$;
    }

    auto setEquilibriumRows(const Vector& values, Vector& vec) const -> void
    {
        setRows(idx_equilibrium_species$, values, vec);
    }

    auto setKineticRows(const Vector& values, Vector& vec) const -> void
    {
        setRows(idx_kinetic_species$, values, vec);
    }

    auto setInertRows(const Vector& values, Vector& vec) const -> void
    {
        setRows(idx_inert_species$, values, vec);
    }

    auto equilibriumRows(const Vector& vec) const -> Vector
    {
        return rows(idx_equilibrium_species$, vec);
    }

    auto kineticRows(const Vector& vec) const -> Vector
    {
        return rows(idx_kinetic_species$, vec);
    }

    auto inertRows(const Vector& vec) const -> Vector
    {
        return rows(idx_inert_species$, vec);
    }

    auto equilibriumCols(const Matrix& mat) const -> Matrix
    {
        return cols(idx_equilibrium_species$, mat);
    }

    auto kineticCols(const Matrix& mat) const -> Matrix
    {
        return cols(idx_kinetic_species$, mat);
    }

    auto inertCols(const Matrix& mat) const -> Matrix
    {
        return cols(idx_inert_species$, mat);
    }

    auto equilibriumRowsCols(const Matrix& mat) const -> Matrix
    {
        return submatrix(idx_equilibrium_species$, idx_equilibrium_species$, mat);
    }

    auto kineticRowsCols(const Matrix& mat) const -> Matrix
    {
        return submatrix(idx_kinetic_species$, idx_kinetic_species$, mat);
    }

    auto inertRowsCols(const Matrix& mat) const -> Matrix
    {
        return submatrix(idx_inert_species$, idx_inert_species$, mat);
    }

    auto equilibriumFormulaMatrix(const ChemicalSystem& system) const -> Matrix
    {
        return submatrix(idx_equilibrium_elements$, idx_equilibrium_species$, system.formulaMatrix());
    }

    auto kineticFormulaMatrix(const ChemicalSystem& system) const -> Matrix
    {
        return submatrix(idx_kinetic_elements$, idx_kinetic_species$, system.formulaMatrix());
    }

    auto inertFormulaMatrix(const ChemicalSystem& system) const -> Matrix
    {
        return submatrix(idx_inert_elements$, idx_inert_species$, system.formulaMatrix());
    }

    auto equilibriumStoichiometricMatrix(const ReactionSystem& reactions) const -> Matrix
    {
        return cols(idx_equilibrium_species$, reactions.stoichiometricMatrix());
    }

    auto kineticStoichiometricMatrix(const ReactionSystem& reactions) const -> Matrix
    {
        return cols(idx_kinetic_species$, reactions.stoichiometricMatrix());
    }

    auto inertStoichiometricMatrix(const ReactionSystem& reactions) const -> Matrix
    {
        return cols(idx_inert_species$, reactions.stoichiometricMatrix());
    }
};

Partitioning::Partitioning(const ChemicalSystem& system)
: pimpl(new Impl(system))
{}

Partitioning::Partitioning(const Partitioning& other)
: pimpl(new Impl(*other.pimpl))
{}

Partitioning::~Partitioning()
{}

auto Partitioning::operator=(Partitioning other) -> Partitioning&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Partitioning::setEquilibriumSpecies(const std::vector<std::string>& equilibrium_species) -> void
{
    pimpl->setEquilibriumSpecies(equilibrium_species);
}

auto Partitioning::setEquilibriumSpecies(const std::string& equilibrium_species) -> void
{
    pimpl->setEquilibriumSpecies(equilibrium_species);
}

auto Partitioning::setKineticSpecies(const std::vector<std::string>& kinetic_species) -> void
{
    pimpl->setKineticSpecies(kinetic_species);
}

auto Partitioning::setKineticSpecies(const std::string& kinetic_species) -> void
{
    pimpl->setKineticSpecies(kinetic_species);
}

auto Partitioning::setInertSpecies(const std::vector<std::string>& inert_species) -> void
{
    pimpl->setInertSpecies(inert_species);
}

auto Partitioning::setInertSpecies(const std::string& inert_species) -> void
{
    pimpl->setInertSpecies(inert_species);
}

auto Partitioning::numSpecies() const -> unsigned
{
    return pimpl->numSpecies();
}

auto Partitioning::numEquilibriumSpecies() const -> unsigned
{
    return pimpl->numEquilibriumSpecies();
}

auto Partitioning::numKineticSpecies() const -> unsigned
{
    return pimpl->numKineticSpecies();
}

auto Partitioning::numInertSpecies() const -> unsigned
{
    return pimpl->numInertSpecies();
}

auto Partitioning::numEquilibriumElements() const -> unsigned
{
    return pimpl->numEquilibriumElements();
}

auto Partitioning::numKineticElements() const -> unsigned
{
    return pimpl->numKineticElements();
}

auto Partitioning::numInertElements() const -> unsigned
{
    return pimpl->numInertElements();
}

auto Partitioning::numPhasesWithEquilibriumSpecies() const -> unsigned
{
    return pimpl->numPhasesWithEquilibriumSpecies();
}

auto Partitioning::numPhasesWithKineticSpecies() const -> unsigned
{
    return pimpl->numPhasesWithKineticSpecies();
}

auto Partitioning::numPhasesWithInertSpecies() const -> unsigned
{
    return pimpl->numPhasesWithInertSpecies();
}

auto Partitioning::species() const -> const std::vector<std::string>&
{
    return pimpl->species$;
}

auto Partitioning::equilibriumSpecies() const -> const std::vector<std::string>&
{
    return pimpl->equilibrium_species$;
}

auto Partitioning::kineticSpecies() const -> const std::vector<std::string>&
{
    return pimpl->kinetic_species$;
}

auto Partitioning::inertSpecies() const -> const std::vector<std::string>&
{
    return pimpl->inert_species$;
}

auto Partitioning::equilibriumElements() const -> const std::vector<std::string>&
{
    return pimpl->equilibrium_elements$;
}

auto Partitioning::kineticElements() const -> const std::vector<std::string>&
{
    return pimpl->kinetic_elements$;
}

auto Partitioning::inertElements() const -> const std::vector<std::string>&
{
    return pimpl->inert_elements$;
}

auto Partitioning::idxEquilibriumSpecies() const -> const Indices&
{
    return pimpl->idx_equilibrium_species$;
}

auto Partitioning::idxKineticSpecies() const -> const Indices&
{
    return pimpl->idx_kinetic_species$;
}

auto Partitioning::idxInertSpecies() const -> const Indices&
{
    return pimpl->idx_inert_species$;
}

auto Partitioning::idxEquilibriumElements() const -> const Indices&
{
    return pimpl->idx_equilibrium_elements$;
}

auto Partitioning::idxKineticElements() const -> const Indices&
{
    return pimpl->idx_kinetic_elements$;
}

auto Partitioning::idxInertElements() const -> const Indices&
{
    return pimpl->idx_inert_elements$;
}

auto Partitioning::idxSpecies(const std::string& species) const -> Index
{
    return pimpl->idxSpecies(species);
}

auto Partitioning::idxSpecies(const std::vector<std::string>& species) const -> Indices
{
    return pimpl->idxSpecies(species);
}

auto Partitioning::idxEquilibriumSpecies(const std::string& species) const -> Index
{
    return pimpl->idxEquilibriumSpecies(species);
}

auto Partitioning::idxEquilibriumSpecies(const std::vector<std::string>& species) const -> Indices
{
    return pimpl->idxEquilibriumSpecies(species);
}

auto Partitioning::idxKineticSpecies(const std::string& species) const -> Index
{
    return pimpl->idxKineticSpecies(species);
}

auto Partitioning::idxKineticSpecies(const std::vector<std::string>& species) const -> Indices
{
    return pimpl->idxKineticSpecies(species);
}

auto Partitioning::idxInertSpecies(const std::string& species) const -> Index
{
    return pimpl->idxInertSpecies(species);
}

auto Partitioning::idxInertSpecies(const std::vector<std::string>& species) const -> Indices
{
    return pimpl->idxInertSpecies(species);
}

auto Partitioning::idxEquilibriumElement(const std::string& element) const -> Index
{
    return pimpl->idxEquilibriumElement(element);
}

auto Partitioning::idxEquilibriumElements(const std::vector<std::string>& elements) const -> Indices
{
    return pimpl->idxEquilibriumElements(elements);
}

auto Partitioning::idxKineticElement(const std::string& element) const -> Index
{
    return pimpl->idxKineticElement(element);
}

auto Partitioning::idxKineticElements(const std::vector<std::string>& elements) const -> Indices
{
    return pimpl->idxKineticElements(elements);
}

auto Partitioning::idxInertElement(const std::string& element) const -> Index
{
    return pimpl->idxInertElement(element);
}

auto Partitioning::idxInertElements(const std::vector<std::string>& elements) const -> Indices
{
    return pimpl->idxInertElements(elements);
}

auto Partitioning::idxPhasesWithEquilibriumSpecies() const -> const Indices&
{
    return pimpl->idxPhasesWithEquilibriumSpecies();
}

auto Partitioning::idxPhasesWithKineticSpecies() const -> const Indices&
{
    return pimpl->idxPhasesWithKineticSpecies();
}

auto Partitioning::idxPhasesWithInertSpecies() const -> const Indices&
{
    return pimpl->idxPhasesWithInertSpecies();
}

auto Partitioning::setEquilibriumRows(const Vector& values, Vector& vec) const -> void
{
    pimpl->setEquilibriumRows(values, vec);
}

auto Partitioning::setKineticRows(const Vector& values, Vector& vec) const -> void
{
    pimpl->setKineticRows(values, vec);
}

auto Partitioning::setInertRows(const Vector& values, Vector& vec) const -> void
{
    pimpl->setInertRows(values, vec);
}

auto Partitioning::equilibriumRows(const Vector& vec) const -> Vector
{
    return pimpl->equilibriumRows(vec);
}

auto Partitioning::kineticRows(const Vector& vec) const -> Vector
{
    return pimpl->kineticRows(vec);
}

auto Partitioning::inertRows(const Vector& vec) const -> Vector
{
    return pimpl->inertRows(vec);
}

auto Partitioning::equilibriumCols(const Matrix& mat) const -> Matrix
{
    return pimpl->equilibriumCols(mat);
}

auto Partitioning::kineticCols(const Matrix& mat) const -> Matrix
{
    return pimpl->kineticCols(mat);
}

auto Partitioning::inertCols(const Matrix& mat) const -> Matrix
{
    return pimpl->inertCols(mat);
}

auto Partitioning::equilibriumRowsCols(const Matrix& mat) const -> Matrix
{
    return pimpl->equilibriumRowsCols(mat);
}

auto Partitioning::kineticRowsCols(const Matrix& mat) const -> Matrix
{
    return pimpl->kineticRowsCols(mat);
}

auto Partitioning::inertRowsCols(const Matrix& mat) const -> Matrix
{
    return pimpl->inertRowsCols(mat);
}

auto Partitioning::equilibriumFormulaMatrix(const ChemicalSystem& system) const -> Matrix
{
    return pimpl->equilibriumFormulaMatrix(system);
}

auto Partitioning::kineticFormulaMatrix(const ChemicalSystem& system) const -> Matrix
{
    return pimpl->kineticFormulaMatrix(system);
}

auto Partitioning::inertFormulaMatrix(const ChemicalSystem& system) const -> Matrix
{
    return pimpl->inertFormulaMatrix(system);
}

auto Partitioning::equilibriumStoichiometricMatrix(const ReactionSystem& reactions) const -> Matrix
{
    return pimpl->equilibriumStoichiometricMatrix(reactions);
}

auto Partitioning::kineticStoichiometricMatrix(const ReactionSystem& reactions) const -> Matrix
{
    return pimpl->kineticStoichiometricMatrix(reactions);
}

auto Partitioning::inertStoichiometricMatrix(const ReactionSystem& reactions) const -> Matrix
{
    return pimpl->inertStoichiometricMatrix(reactions);
}

} /* namespace Reaktor */
