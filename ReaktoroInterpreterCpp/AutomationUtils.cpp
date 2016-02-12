// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

#include "Interpreter.hpp"

#include <ReaktoroInterpreterCpp/ParserUtils.hpp>

namespace Reaktoro {
namespace {

/// Return all entity names it can find in list of EntityValueUnits objects.
template<typename Triplet>
auto collectEntities(std::vector<std::string>& list, const std::vector<Triplet>& triplets) -> void
{
    for(const Triplet& t : triplets)
        list.push_back(t.entity);
}

/// Return all titrant names it can find in a list of EquilibriumConstraintNode objects.
template<typename Constraint>
auto collectTitrants(std::vector<std::string>& list, const std::vector<Constraint>& constraints) -> void
{
    for(const Constraint& c : constraints)
    {
        list.push_back(c.entity);
        list.push_back(c.titrant1);
        list.push_back(c.titrant2);
    }
}

} // namespace

auto collectCompounds(const EquilibriumNode& e) -> std::vector<std::string>
{
    std::vector<std::string> list;
    collectEntities(list, e.mixture);
    collectTitrants(list, e.pH);
    collectTitrants(list, e.species_amounts);
    collectTitrants(list, e.species_activities);
    collectEntities(list, e.inert_species);
    remove(list, [](std::string x) { return x.empty(); });
    return list;
}

} // namespace Reaktoro
