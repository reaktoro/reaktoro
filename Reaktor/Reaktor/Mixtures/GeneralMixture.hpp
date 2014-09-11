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
#include <Reaktor/Common/ThermoVector.hpp>

namespace Reaktor {

/// A type that describes a general mixture of species, used as the base class of real mixtures
/// @ingroup Mixtures
template<class SpeciesType>
class GeneralMixture
{
public:
	/// Construct a default GeneralMixture instance
    GeneralMixture();

	/// Construct a GeneralMixture instance with given species
	/// @param species The species that compose the mixture
    GeneralMixture(const std::vector<SpeciesType>& species);

	/// Destroy the GeneralMixture instance
    virtual ~GeneralMixture();

	/// Get the species that compose the mixture
    auto species() const -> const std::vector<SpeciesType>&;

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
auto GeneralMixture<SpeciesType>::species() const -> const std::vector<SpeciesType>&
{
    return m_species;
}

/// Get the number of species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
template<class SpeciesType>
auto numSpecies(const GeneralMixture<SpeciesType>& mixture) -> unsigned
{
    return mixture.species().size();
}

/// Get the names of the species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
template<class SpeciesType>
auto namesSpecies(const GeneralMixture<SpeciesType>& mixture) -> std::vector<std::string>
{
	const unsigned num_species = numSpecies(mixture);
    std::vector<std::string> names(num_species);
    for(unsigned i = 0; i < num_species; ++i)
        names[i] = mixture.species()[i].name;
    return names;
}

/// Get the index of a species in a mixture
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
/// @param name The name of the species in the mixture
/// @return The index of the species if found, or the number of species otherwise
template<class SpeciesType>
auto indexSpecies(const GeneralMixture<SpeciesType>& mixture, const std::string& name) -> Index
{
	const auto& species = mixture.species();
    for(Index i = 0; i < species.size(); ++i)
        if(species[i].name == name) return i;
    return species.size();
}

/// Calculate the molar fractions of the species and its molar derivatives
/// @param mixture The mixture (e.g., aqueous, gaseous, mineral, etc.)
/// @param n The molar abundance of the species (in units of mol)
/// @return The molar fractions and its molar derivatives
template<class SpeciesType>
auto molarFractions(const GeneralMixture<SpeciesType>& mixture, const Vector& n) -> ThermoVector
{
    const unsigned num_species = numSpecies(mixture);
    const double nt = arma::sum(n);
    ThermoVector x(num_species);
    x.val = n/nt;
    for(unsigned i = 0; i < num_species; ++i)
        x.ddn.row(i).fill(-x.val[i]/nt);
    for(unsigned i = 0; i < num_species; ++i)
    	x.ddn(i, i) += 1.0/nt;
    return x;
}

} // namespace Reaktor
