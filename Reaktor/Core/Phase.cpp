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

#include "Phase.hpp"

// C++ includes
#include <algorithm>

namespace Reaktor {

struct Phase::Impl
{
    PhaseData data;
};

Phase::Phase()
: pimpl(new Impl())
{}

Phase::Phase(const PhaseData& data)
: pimpl(new Impl{data})
{}

Phase::Phase(std::string name, std::vector<Species> species)
: Phase(PhaseData{name, species})
{}

auto Phase::name() const -> const std::string&
{
    return pimpl->data.name;
}

auto Phase::species() const -> const std::vector<Species>&
{
    return pimpl->data.species;
}

auto operator<(const Phase& lhs, const Phase& rhs) -> bool
{
	return lhs.name() < rhs.name();
}

auto operator==(const Phase& lhs, const Phase& rhs) -> bool
{
	return lhs.name() == rhs.name();
}

auto collectSpecies(const std::vector<Phase>& phases) -> std::vector<Species>
{
    unsigned num_species = 0;
    for(const Phase& phase : phases)
        num_species += phase.species().size();

    std::vector<Species> list;
    list.reserve(num_species);
    for(const Phase& phase : phases)
        for(const Species& iter : phase.species())
            list.push_back(iter);
    return list;
}

} // namespace Reaktor
