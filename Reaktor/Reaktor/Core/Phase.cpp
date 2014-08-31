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

// Reaktor includes
#include <Reaktor/Core/Species.hpp>
#include <Reaktor/Core/SpeciesUtils.hpp>

namespace Reaktor {

struct Phase::Impl
{
    Impl()
    {}

    /// The name of the phase
    std::string name;

    /// The chemical species that compose the phase
    std::vector<Species> species;

    /// The concentration function of the phase
    Concentration concentration;
};

Phase::Phase()
: pimpl(new Impl())
{}

Phase::Phase(const Phase& other)
: pimpl(new Impl(*other.pimpl))
{}

Phase::~Phase()
{}

auto Phase::operator=(Phase other) -> Phase&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Phase::setName(const std::string& name) -> Phase&
{
    pimpl->name = name;
    return *this;
}

auto Phase::setSpecies(const std::vector<Species>& species) -> Phase&
{
    pimpl->species = species;
    return *this;
}

auto Phase::setConcentration(const Concentration& concentration) -> Phase&
{
    pimpl->concentration = concentration;
    return *this;
}

auto Phase::name() const -> const std::string&
{
    return pimpl->name;
}

auto Phase::species() const -> const std::vector<Species>&
{
    return pimpl->species;
}

auto Phase::concentration() const -> const Concentration&
{
    return pimpl->concentration;
}

auto Phase::operator==(const Phase& phase) const -> bool
{
    return name() == phase.name();
}

auto operator<<(std::ostream& out, const Phase& phase) -> std::ostream&
{
    out << phase.name() << std::endl;
    out << "  ";
    const auto& species = phase.species();
    for(unsigned i = 0; i < species.size(); ++i)
        out << (i > 0 ? ", " : "") << species[i].name();
    return out;
}

} // namespace Reaktor
