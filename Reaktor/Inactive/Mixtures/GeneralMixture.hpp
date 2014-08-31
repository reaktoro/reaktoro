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

#pragma once

// C++ includes
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/PartialVector.hpp>

namespace Reaktor {

/**
 * Provides a base of implementation for the mixture classes
 *
 * @ingroup Mixtures
 */
template<class SpeciesType>
class GeneralMixture
{
public:
    /**
     * Constructs a default GeneralMixture instance
     */
    GeneralMixture();

    /**
     * Constructs a GeneralMixture instance with given species
     * @param species The names of the species in the mixture
     */
    GeneralMixture(const std::vector<SpeciesType>& species);

    /**
     * Destroys the instance
     */
    virtual ~GeneralMixture();

    /**
     * Gets the number of species in the mixture
     */
    auto numSpecies() const -> unsigned;

    /**
     * Gets the species that compose the mixture
     * @return The species that compose the mixture
     */
    auto species() const -> const std::vector<SpeciesType>&;

    /**
     * Gets a species in the mixture
     * @param idx_species The index of the species
     * @return The species with given index
     */
    auto species(const Index& idx_species) const -> const SpeciesType&;

    /**
     * Gets the names of the species in the mixture
     */
    auto speciesNames() const -> std::vector<std::string>;

    /**
     * Gets the index of a species in the mixture
     * @param species The name of the species
     * @return The index of the given species if found. The number of species in otherwise.
     */
    auto idxSpecies(const std::string& species) const -> Index;

    /**
     * Calculates the molar fractions of the species and its molar derivatives
     * @param n The molar abundance of the species (in units of mol)
     * @return The molar fractions and its molar derivatives
     */
    auto molarFractions(const Vector& n) const -> PartialVector;

private:
    /// The name of the species in the mixture
    std::vector<SpeciesType> m_species;
};

template<class SpeciesType>
inline GeneralMixture<SpeciesType>::GeneralMixture(const std::vector<SpeciesType>& species)
: m_species(species)
{}

template<class SpeciesType>
inline GeneralMixture<SpeciesType>::GeneralMixture()
{}

template<class SpeciesType>
inline GeneralMixture<SpeciesType>::~GeneralMixture()
{}

template<class SpeciesType>
inline auto GeneralMixture<SpeciesType>::numSpecies() const -> unsigned
{
    return m_species.size();
}

template<class SpeciesType>
inline auto GeneralMixture<SpeciesType>::species(const Index& idx_species) const -> const SpeciesType&
{
    return m_species[idx_species];
}

template<class SpeciesType>
inline auto GeneralMixture<SpeciesType>::speciesNames() const -> std::vector<std::string>
{
    std::vector<std::string> names(m_species.size());

    for(unsigned i = 0; i < names.size(); ++i)
        names[i] = m_species[i].name();

    return names;
}

template<class SpeciesType>
inline auto GeneralMixture<SpeciesType>::species() const -> const std::vector<SpeciesType>&
{
    return m_species;
}

template<class SpeciesType>
inline auto GeneralMixture<SpeciesType>::idxSpecies(const std::string& name) const -> Index
{
    for(Index i = 0; i < m_species.size(); ++i)
        if(m_species[i].name() == name) return i;
    return numSpecies();
}

template<class SpeciesType>
inline auto GeneralMixture<SpeciesType>::molarFractions(const Vector& n) const -> PartialVector
{
    const unsigned num_species = numSpecies();
    const double nt = n.sum();
    const Vector x = n / nt;
    Matrix dxdn(num_species, num_species);
    for(unsigned i = 0; i < num_species; ++i)
        dxdn.row(i).setConstant(-x[i]/nt);
    for(unsigned i = 0; i < num_species; ++i)
        dxdn(i, i) += 1.0/nt;
    return partialVector(x, dxdn);
}

} // namespace Reaktor
