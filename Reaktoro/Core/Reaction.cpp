// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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

#include "Reaction.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>

namespace Reaktoro {

struct Reaction::Impl
{
    /// The name that uniquely identifies this reaction.
    String name;

    /// The equation of the reaction with its species and stoichiometric coefficients.
    ReactionEquation equation;

    /// The function that computes the equilibrium constant of the reaction (in natural log).
    EquilibriumConstantFn lnKfn;

    /// The function that computes the rate of the reaction (in mol/s).
    ReactionRateFn ratefn;
};

Reaction::Reaction()
: pimpl(new Impl())
{}

auto Reaction::clone() const -> Reaction
{
    Reaction reaction;
    *reaction.pimpl = *pimpl;
    return reaction;
}

auto Reaction::withName(String name) const -> Reaction
{
    Reaction copy = clone();
    copy.pimpl->name = name;
    return copy;
}

auto Reaction::withEquation(const ReactionEquation& equation) const -> Reaction
{
    Reaction copy = clone();
    copy.pimpl->equation = equation;
    return copy;
}

auto Reaction::withEquilibriumConstantFn(const EquilibriumConstantFn& fn) const -> Reaction
{
    Reaction copy = clone();
    copy.pimpl->lnKfn = fn;
    return copy;
}

auto Reaction::withRateFn(const ReactionRateFn& fn) const -> Reaction
{
    Reaction copy = clone();
    copy.pimpl->ratefn = fn;
    return copy;
}

auto Reaction::name() const -> String
{
    return pimpl->name;
}

auto Reaction::equation() const -> const ReactionEquation&
{
    return pimpl->equation;
}

auto Reaction::equilibriumConstantFn() const -> const EquilibriumConstantFn&
{
    return pimpl->lnKfn;
}

auto Reaction::rateFn() const -> const ReactionRateFn&
{
    return pimpl->ratefn;
}

auto operator<(const Reaction& lhs, const Reaction& rhs) -> bool
{
    return lhs.name() < rhs.name();
}

auto operator==(const Reaction& lhs, const Reaction& rhs) -> bool
{

    return lhs.name() == rhs.name();
}

} // namespace Reaktoro
