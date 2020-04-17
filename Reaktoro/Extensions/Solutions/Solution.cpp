// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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

#include "Solution.hpp"

namespace Reaktoro {

/// Compare two SolutionState instances for equality
auto operator==(const SolutionState& l, const SolutionState& r) -> bool
{
    return l.T == r.T && r.P == r.P && l.x == r.x;
}

Solution::Solution()
{}

Solution::Solution(std::vector<Species> species)
: m_species(std::move(species))
{}

auto Solution::setName(std::string name) -> void
{
    m_name = name;
}

auto Solution::numSpecies() const -> unsigned
{
    return m_species.size();
}

auto Solution::name() const -> std::string
{
    return m_name;
}

auto Solution::species() const -> const std::vector<Species>&
{
    return m_species;
}

auto Solution::species(const Index& index) const -> const Species&
{
    return m_species[index];
}

auto Solution::indexSpecies(const std::string& name) const -> Index
{
    return indexfn(m_species, RKT_LAMBDA(s, s.name() == name));
}

auto Solution::indexSpeciesAny(const std::vector<std::string>& names) const -> Index
{
    return indexAny(names, m_species);
}

auto Solution::namesSpecies() const -> std::vector<std::string>
{
    std::vector<std::string> names(m_species.size());
    for(unsigned i = 0; i < names.size(); ++i)
        names[i] = m_species[i].name();
    return names;
}

auto Solution::charges() const -> VectorXr
{
    const unsigned nspecies = numSpecies();
    VectorXr charges(nspecies);
    for(unsigned i = 0; i < nspecies; ++i)
        charges[i] = m_species[i].charge();
    return charges;
}

auto Solution::moleFractions(VectorXrConstRef n) const -> VectorXr
{
    const unsigned nspecies = numSpecies();
    if(nspecies == 1)
        return VectorXr::Ones(1);
    const real nsum = n.sum();
    if(nsum == 0.0) return VectorXr::Zero(nspecies);
    return n/nsum;
}

auto Solution::state(real T, real P, VectorXrConstRef n) const -> SolutionState
{
    SolutionState res;
    res.T = T;
    res.P = P;
    res.x = moleFractions(n);
    return res;
}

} // namespace Reaktoro
