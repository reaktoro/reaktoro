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

#pragma once

// C++ includes
#include <string>
#include <vector>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>

namespace Reaktor {

/// Provide a base of implementation for the mixture classes.
/// @ingroup Mixtures
template<class SpeciesType>
class GeneralMixture
{
public:
    /// Construct a default GeneralMixture instance
    GeneralMixture();

    /// Construct a GeneralMixture instance with given species
    /// @param species The names of the species in the mixture
    GeneralMixture(const std::vector<SpeciesType>& species);

    /// Destroy the instance
    virtual ~GeneralMixture();

    /// Return the number of species in the mixture
    auto numSpecies() const -> unsigned;

    /// Return the species that compose the mixture
    /// @return The species that compose the mixture
    auto species() const -> const std::vector<SpeciesType>&;

    /// Return a species in the mixture
    /// @param index The index of the species
    /// @return The species with given index
    auto species(const Index& index) const -> const SpeciesType&;

    /// Return the index of a species in the mixture
    /// @param name The name of the species in the solution
    /// @return The index of the species if found, or the number of species otherwise
    auto indexSpecies(const std::string& name) const -> Index;

    /// Return the names of the species in the mixture
    auto namesSpecies() const -> std::vector<std::string>;

    /// Return the charges of the species in the mixture
    auto chargesSpecies() const -> Vector;

    /// Calculates the molar fractions of the species and its molar derivatives
    /// @param n The molar abundance of the species (in units of mol)
    /// @return The molar fractions and its molar derivatives
    auto molarFractions(const Vector& n) const -> ChemicalVector;

private:
    /// The name of the species in the mixture
    std::vector<SpeciesType> m_species;
};

template<class SpeciesType>
GeneralMixture<SpeciesType>::GeneralMixture(const std::vector<SpeciesType>& species)
: m_species(species)
{}

template<class SpeciesType>
GeneralMixture<SpeciesType>::GeneralMixture()
{}

template<class SpeciesType>
GeneralMixture<SpeciesType>::~GeneralMixture()
{}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::numSpecies() const -> unsigned
{
    return m_species.size();
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::species() const -> const std::vector<SpeciesType>&
{
    return m_species;
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::species(const Index& index) const -> const SpeciesType&
{
    return m_species[index];
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::indexSpecies(const std::string& name) const -> Index
{
    for(Index i = 0; i < m_species.size(); ++i)
        if(m_species[i].name == name) return i;
    return numSpecies();
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::namesSpecies() const -> std::vector<std::string>
{
    std::vector<std::string> names(m_species.size());
    for(unsigned i = 0; i < names.size(); ++i)
        names[i] = m_species[i].name;
    return names;
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::chargesSpecies() const -> Vector
{
    const unsigned nspecies = numSpecies();
    Vector charges(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        charges[i] = m_species[i].charge;
    return charges;
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::molarFractions(const Vector& n) const -> ChemicalVector
{
    const unsigned nspecies = n.size();
    const double nt = n.sum();
    Vector x = zeros(nspecies);
    Matrix dxdt = zeros(nspecies);
    Matrix dxdp = zeros(nspecies);
    Matrix dxdn = zeros(nspecies, nspecies);

    if(nt == 0.0)
        return {x, dxdt, dxdp, dxdn};

    x = n/nt;

    for(unsigned i = 0; i < nspecies; ++i)
    {
        dxdn.row(i).fill(-x[i]/nt);
        dxdn(i, i) += 1.0/nt;
    }

    return {x, dxdt, dxdp, dxdn};
}

} // namespace Reaktor
