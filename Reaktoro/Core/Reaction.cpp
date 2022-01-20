// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Units.hpp>

namespace Reaktoro {

const auto R = universalGasConstant;

struct Reaction::Impl
{
    /// The name that uniquely identifies this reaction.
    String name;

    /// The equation of the reaction with its species and stoichiometric coefficients.
    ReactionEquation equation;

    /// The function that computes the rate of the reaction (in mol/s).
    ReactionRateFn ratefn;

    /// Construct a default Reaction::Impl object
    Impl()
    : ratefn([](const ChemicalProps&) { return 0.0; })
    {
    }

    /// Calculate the complete set of thermodynamic properties of the reaction.
    auto props(real T, real P) const -> ReactionProps
    {
        ReactionProps rprops;
        rprops.T = T;
        rprops.P = P;

        for(const auto& [species, coeff] : equation)
        {
            const auto sprops = species.props(T, P);
            rprops.dG0  += coeff * sprops.G0;
            rprops.dH0  += coeff * sprops.H0;
            rprops.dV0  += coeff * sprops.V0;
            rprops.dVT0 += coeff * sprops.VT0;
            rprops.dVP0 += coeff * sprops.VP0;
            rprops.dCp0 += coeff * sprops.Cp0;
            rprops.dCv0 += coeff * sprops.Cv0;
            rprops.dU0  += coeff * sprops.U0;
            rprops.dS0  += coeff * sprops.S0;
            rprops.dA0  += coeff * sprops.A0;
        }

        rprops.lgK = -rprops.dG0/(R*T*ln10);

        return rprops;
    }
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

auto Reaction::rateFn() const -> const ReactionRateFn&
{
    return pimpl->ratefn;
}

auto Reaction::props(real T, real P) const -> ReactionProps
{
    return pimpl->props(T, P);
}

auto Reaction::props(real T, Chars unitT, real P, Chars unitP) const -> ReactionProps
{
    T = units::convert(T, unitT, "K");
    P = units::convert(P, unitP, "Pa");
    return props(T, P);
}

auto Reaction::rate(const ChemicalProps& props) const -> real
{
    return pimpl->ratefn(props);
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
