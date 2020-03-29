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

#include "GeneralMixture.hpp"

namespace Reaktoro {

/// Compare two MixtureState instances for equality
auto operator==(const MixtureState& l, const MixtureState& r) -> bool
{
    return l.T == r.T && r.P == r.P && l.x == r.x;
}

GeneralMixture::GeneralMixture()
{}

GeneralMixture::GeneralMixture(const std::vector<Species>& species)
: _species(species)
{}

GeneralMixture::~GeneralMixture()
{}

auto GeneralMixture::setName(std::string name) -> void
{
    _name = name;
}

auto GeneralMixture::numSpecies() const -> unsigned
{
    return _species.size();
}

auto GeneralMixture::name() const -> std::string
{
    return _name;
}

auto GeneralMixture::species() const -> const std::vector<Species>&
{
    return _species;
}

auto GeneralMixture::species(const Index& index) const -> const Species&
{
    return _species[index];
}

auto GeneralMixture::indexSpecies(const std::string& name) const -> Index
{
    return index(name, _species);
}

auto GeneralMixture::indexSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    return indexAny(names, _species);
}

auto GeneralMixture::namesSpecies() const -> std::vector<std::string>
{
    std::vector<std::string> names(_species.size());
    for(unsigned i = 0; i < names.size(); ++i)
        names[i] = _species[i].name();
    return names;
}

auto GeneralMixture::chargesSpecies() const -> VectorXr
{
    const unsigned nspecies = numSpecies();
    VectorXr charges(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        charges[i] = _species[i].charge();
    return charges;
}

auto GeneralMixture::moleFractions(VectorXrConstRef n) const -> VectorXr
{
    const unsigned nspecies = numSpecies();
    if(nspecies == 1)
        return VectorXr::Ones(1);
    const real nsum = n.sum();
    if(nsum == 0.0) return VectorXr::Zero(nspecies);
    return n/nsum;
}

auto GeneralMixture::state(real T, real P, VectorXrConstRef n) const -> MixtureState
{
    MixtureState res;
    res.T = T;
    res.P = P;
    res.x = moleFractions(n);
    return res;
}

} // namespace Reaktoro
