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

#include "Reaction.hpp"

// C++ includes
#include <cmath>
#include <iomanip>
#include <sstream>

// Reaktor includes
#include <Reaktor/Common/Macros.hpp>
#include <Reaktor/Common/ThermoScalar.hpp>
#include <Reaktor/Common/ThermoVector.hpp>

namespace Reaktor {

struct Reaction::Impl
{
    /// The names of the reacting species of the reaction
    std::vector<std::string> species;

    /// The indices of the reacting species of the reaction
    Indices indices;

    /// The stoichiometries of the reacting species of the reaction
    std::vector<double> stoichiometries;

    /// The thermodynamic model of the reaction
    ReactionThermoModel thermo_model;

    /// The kinetics model of the reaction
    ReactionKineticsModel kinetics_model;
};

Reaction::Reaction()
: pimpl(new Impl())
{}

Reaction::Reaction(const Reaction& other)
: pimpl(new Impl(*other.pimpl))
{}

Reaction::~Reaction()
{}

auto Reaction::operator=(Reaction other) -> Reaction&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto Reaction::setSpecies(const std::vector<std::string>& species) -> Reaction&
{
	pimpl->species = species;
    return *this;
}

auto Reaction::setIndices(const Indices& indices) -> Reaction&
{
	pimpl->indices = indices;
    return *this;
}

auto Reaction::setStoichiometries(const std::vector<double>& stoichiometries) -> Reaction&
{
	pimpl->stoichiometries = stoichiometries;
    return *this;
}

auto Reaction::setThermoModel(const ReactionThermoModel& thermo_model) -> Reaction&
{
	pimpl->thermo_model = thermo_model;
    return *this;
}

auto Reaction::setKineticsModel(const ReactionKineticsModel& kinetics_model) -> Reaction&
{
	pimpl->kinetics_model = kinetics_model;
    return *this;
}

auto Reaction::thermoModel() const -> const ReactionThermoModel&
{
    return pimpl->thermo_model;
}

auto Reaction::kineticsModel() const -> const ReactionKineticsModel&
{
    return pimpl->kinetics_model;
}

} // namespace Reaktor
