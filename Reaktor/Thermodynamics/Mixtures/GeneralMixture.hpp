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
#include <Reaktor/Common/ChemicalScalar.hpp>
#include <Reaktor/Common/ChemicalVector.hpp>

namespace Reaktor {

/// A type used to describe the state of a mixture
struct MixtureState
{
    /// The temperature of the mixture (in units of K)
    double T;

    /// The pressure of the mixture (in units of Pa)
    double P;

    /// The amounts of the species in the mixture (in units of mol)
    Vector n;

    /// The molar fractions of the species in the mixture and their partial derivatives
    ChemicalVector x;
};

/// Compare two MixtureState instances for equality
auto operator==(const MixtureState& l, const MixtureState& r) -> bool;

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

    /// Set the name of the mixture.
    auto setName(std::string name) -> void;

    /// Destroy the instance
    virtual ~GeneralMixture();

    /// Return the number of species in the mixture
    auto numSpecies() const -> unsigned;

    /// Return the name of the mixture.
    auto name() const -> std::string;

    /// Return the species that compose the mixture
    /// @return The species that compose the mixture
    auto species() const -> const std::vector<SpeciesType>&;

    /// Return a species in the mixture
    /// @param index The index of the species
    /// @return The species with given index
    auto species(const Index& index) const -> const SpeciesType&;

    /// Return the index of a species in the mixture
    /// @param name The name of the species in the mixture
    /// @return The index of the species if found, or the number of species otherwise
    auto indexSpecies(const std::string& name) const -> Index;

    /// Return the names of the species in the mixture
    auto namesSpecies() const -> std::vector<std::string>;

    /// Return the charges of the species in the mixture
    auto chargesSpecies() const -> Vector;

    /// Calculates the molar fractions of the species and their partial derivatives
    /// @param n The molar abundance of the species (in units of mol)
    /// @return The molar fractions and their partial derivatives
    auto molarFractions(const Vector& n) const -> ChemicalVector;

    /// Calculate the state of the mixture.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of bar)
    /// @param n The molar amounts of the species in the mixture (in units of mol)
    auto state(double T, double P, const Vector& n) const -> MixtureState;

private:
    /// The name of mixture
    std::string _name;

    /// The species in the mixture
    std::vector<SpeciesType> _species;
};

template<class SpeciesType>
GeneralMixture<SpeciesType>::GeneralMixture(const std::vector<SpeciesType>& species)
: _species(species)
{}

template<class SpeciesType>
GeneralMixture<SpeciesType>::GeneralMixture()
{}

template<class SpeciesType>
GeneralMixture<SpeciesType>::~GeneralMixture()
{}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::setName(std::string name) -> void
{
    return _name = name;
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::numSpecies() const -> unsigned
{
    return _species.size();
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::name() const -> std::string
{
    return _name;
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::species() const -> const std::vector<SpeciesType>&
{
    return _species;
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::species(const Index& index) const -> const SpeciesType&
{
    return _species[index];
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::indexSpecies(const std::string& name) const -> Index
{
    for(Index i = 0; i < _species.size(); ++i)
        if(_species[i].name() == name) return i;
    return numSpecies();
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::namesSpecies() const -> std::vector<std::string>
{
    std::vector<std::string> names(_species.size());
    for(unsigned i = 0; i < names.size(); ++i)
        names[i] = _species[i].name();
    return names;
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::chargesSpecies() const -> Vector
{
    const unsigned nspecies = numSpecies();
    Vector charges(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        charges[i] = _species[i].charge();
    return charges;
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::molarFractions(const Vector& n) const -> ChemicalVector
{
    const unsigned nspecies = n.size();
    ChemicalVector x(nspecies, nspecies);
    const double nt = n.sum();
    if(nt == 0.0) return x;
    x.val = n/nt;
    for(unsigned i = 0; i < nspecies; ++i)
    {
        x.ddn.row(i).fill(-x.val[i]/nt);
        x.ddn(i, i) += 1.0/nt;
    }
    return x;
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::state(double T, double P, const Vector& n) const -> MixtureState
{
    MixtureState res;
    res.T = T;
    res.P = P;
    res.n = n;
    res.x = molarFractions(n);
    return res;
}

} // namespace Reaktor
