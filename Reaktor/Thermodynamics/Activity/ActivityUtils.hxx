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

template<typename Solution>
auto numSpecies(const Solution& solution) -> unsigned
{
    return solution.size();
}

template<typename Solution>
auto speciesNames(const Solution& solution) -> std::vector<std::string>
{
    const unsigned nspecies = numSpecies(solution);
    std::vector<std::string> names(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        names[i] = solution[i].name;
    return names;
}

template<typename Solution>
auto speciesCharges(const Solution& solution) -> Vector
{
    const unsigned nspecies = numSpecies(solution);
    Vector charges(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        charges[i] = solution[i].charge;
    return charges;
}

template<class Solution>
auto speciesIndex(const Solution& solution, const std::string& name) -> Index
{
    for(Index i = 0; i < numSpecies(solution); ++i)
        if(solution[i].name == name) return i;
    return solution.size();
}

} // namespace Reaktor
