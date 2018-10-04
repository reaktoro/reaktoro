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

#pragma once

// C++ includes
#include <string>
#include <map>

namespace Reaktoro {

/// Define a type that describes the equation of a reaction.
/// The equation of a reaction is assumed as a sequence of pairs
/// (species, stoichiometry). It is shown below how the equation
/// of reaction @f$\mathrm{CO_{2}(g)+H_{2}O\rightleftharpoons H^{+}+HCO_{3}^{-}}@f$
/// can be defined by two equivalent ways:
/// ~~~~~~~~~~~~~~
/// ReactionEquation equation1 = {{"CO2(g)", -1}, {"H2O(l)", -1}, {"H+", 1}, {"HCO3-", 1}};
/// ReactionEquation equation2 = "CO2(g) + H2O(l) = H+ + HCO3-";
/// ~~~~~~~~~~~~~~
class ReactionEquation
{
public:
    /// Construct a default ReactionEquation instance
    ReactionEquation();

    /// Construct a ReactionEquation instance by parsing a string.
    /// Below are examples of how to set a reaction equation via a formatted string.
    /// ~~~
    /// ReactionEquation equation1("Calcite + H+ = Ca++ + HCO3-");
    /// ReactionEquation equation2("CO2(g) + H2O(l) = H+ + HCO3-");
    /// ReactionEquation equation3("Dolomite + 2*H+ = Ca++ + Mg++ + 2*HCO3-");
    /// ~~~
    /// Note that unity stoichiometry coefficients can be ommited from the equation. 
    /// The operator `*` must be used when this is not the case.
    /// @param equation The string representing the rection equation
    ReactionEquation(std::string equation);

    /// Construct a ReactionEquation instance from a list of species names and a list of stoichiometries
    /// @param species The names of the participating chemical species
    /// @param coeffs The stoichiometries of the participating chemical species
    ReactionEquation(const std::map<std::string, double>& equation);

    /// Return true if the rection equation is empty.
    auto empty() const -> bool;

    /// Return the number of species in the reaction equation.
    auto numSpecies() const -> unsigned;

    /// Return the stoichiometry of a species in the reaction equation.
    /// @param species The name of the species.
    auto stoichiometry(std::string species) const -> double;

    /// Return the reaction equation as a map of species names and stoichiometries.
    auto equation() const -> const std::map<std::string, double>&;

    /// Convert the ReactionEquation instance into a string
    operator std::string() const;

private:
    /// The string representation of the reaction equation
    std::string equation_str;

    /// The reaction equation represented as a map of species names and their stoichiometries
    std::map<std::string, double> equation_map;
};

/// Output a ReactionEquation instance
auto operator<<(std::ostream& out, const ReactionEquation& equation) -> std::ostream&;

/// Return begin const iterator of a ReactionEquation instance
inline auto begin(const Reaktoro::ReactionEquation& equation) -> decltype(equation.equation().begin())
{
    return equation.equation().begin();
}

/// Return begin iterator of a ReactionEquation instance
inline auto begin(Reaktoro::ReactionEquation& equation) -> decltype(equation.equation().begin())
{
    return equation.equation().begin();
}

/// Return end const iterator of a ReactionEquation instance
inline auto end(const Reaktoro::ReactionEquation& equation) -> decltype(equation.equation().end())
{
    return equation.equation().end();
}

/// Return end iterator of a ReactionEquation instance
inline auto end(Reaktoro::ReactionEquation& equation) -> decltype(equation.equation().end())
{
    return equation.equation().end();
}

} // namespace Reaktoro

