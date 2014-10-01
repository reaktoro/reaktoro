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

namespace Reaktor {

template<typename Mixture>
auto numSpecies(const Mixture& mixture) -> unsigned
{
    return mixture.size();
}

template<typename Mixture>
auto speciesNames(const Mixture& mixture) -> std::vector<std::string>
{
    const unsigned nspecies = numSpecies(mixture);
    std::vector<std::string> names(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        names[i] = mixture[i].name;
    return names;
}

template<typename Mixture>
auto speciesCharges(const Mixture& mixture) -> Vector
{
    const unsigned nspecies = numSpecies(mixture);
    Vector charges(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        charges[i] = mixture[i].charge;
    return charges;
}

template<class Mixture>
auto speciesIndex(const Mixture& mixture, const std::string& name) -> Index
{
    for(Index i = 0; i < numSpecies(mixture); ++i)
        if(mixture[i].name == name) return i;
    return mixture.size();
}

} // namespace Reaktor
