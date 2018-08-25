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

#pragma once

// C++ includes
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ChemicalVector.hpp>
#include <Reaktoro/Common/ScalarTypes.hpp>
#include <Reaktoro/Common/SetUtils.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Core/Utils.hpp>

namespace Reaktoro {

/// A type used to describe the state of a mixture
struct MixtureState
{
    /// The temperature of the mixture (in units of K)
    Temperature T;

    /// The pressure of the mixture (in units of Pa)
    Pressure P;

    /// The mole fractions of the species in the mixture and their partial derivatives
    ChemicalVector x;
};

/// Compare two MixtureState instances for equality
inline auto operator==(const MixtureState& l, const MixtureState& r) -> bool
{
    return l.T == r.T && r.P == r.P && l.x == r.x;
}

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
    explicit GeneralMixture(const std::vector<SpeciesType>& species);

    /// Destroy the instance
    virtual ~GeneralMixture();

    /// Set the name of the mixture.
    auto setName(std::string name) -> void;

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

    /// Return the index of the first species in the mixture with any of the given names.
    /// @param names The tentative names of the species in the mixture.
    /// @return The index of the species if found, or the number of species otherwise.
    auto indexSpeciesAny(const std::vector<std::string>& names) const -> Index;

    /// Return the names of the species in the mixture
    auto namesSpecies() const -> std::vector<std::string>;

    /// Return the charges of the species in the mixture
    auto chargesSpecies() const -> Vector;

    /// Calculates the mole fractions of the species and their partial derivatives
    /// @param n The molar abundance of the species (in units of mol)
    /// @return The mole fractions and their partial derivatives
    auto moleFractions(VectorConstRef n) const -> ChemicalVector;

    /// Calculate the state of the mixture.
    /// @param T The temperature (in units of K)
    /// @param P The pressure (in units of Pa)
    /// @param n The molar amounts of the species in the mixture (in units of mol)
    auto state(Temperature T, Pressure P, VectorConstRef n) const -> MixtureState;

private:
    /// The name of mixture
    std::string _name;

    /// The species in the mixture
    std::vector<SpeciesType> _species;
};

template<class SpeciesType>
GeneralMixture<SpeciesType>::GeneralMixture()
{}

template<class SpeciesType>
GeneralMixture<SpeciesType>::GeneralMixture(const std::vector<SpeciesType>& species)
: _species(species)
{}

template<class SpeciesType>
GeneralMixture<SpeciesType>::~GeneralMixture()
{}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::setName(std::string name) -> void
{
    _name = name;
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
    return index(name, _species);
}

template<class SpeciesType>
auto GeneralMixture<SpeciesType>::indexSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    return indexAny(names, _species);
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
auto GeneralMixture<SpeciesType>::moleFractions(VectorConstRef n) const -> ChemicalVector
{
    const unsigned nspecies = numSpecies();
    if(nspecies == 1)
    {
        ChemicalVector x(1);
        x.val[0] = 1.0;
        return x;
    }
    ChemicalVector x(nspecies);
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
auto GeneralMixture<SpeciesType>::state(Temperature T, Pressure P, VectorConstRef n) const -> MixtureState
{
    MixtureState res;
    res.T = T;
    res.P = P;
    res.x = moleFractions(n);
    return res;
}

} // namespace Reaktoro
