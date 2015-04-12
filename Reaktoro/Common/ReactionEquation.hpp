// Reaktoro is a C++ library for computational reaction modelling.
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

// C++ includes
#include <string>
#include <vector>

namespace Reaktoro {

/// Define a type that describes the equation of a reaction.
/// The equation of a reaction is assumed as a sequence of pairs
/// (species, stoichiometry). It is shown below how the equation
/// of reaction @f$\mathrm{CO_{2}(g)+H_{2}O\rightleftharpoons H^{+}+HCO_{3}^{-}}@f$
/// can be defined by two equivalent ways:
/// ~~~~~~~~~~~~~~
/// ReactionEquation equation1 = {{"CO2(g)", -1}, {"H2O(l)", -1}, {"H+", 1}, {"HCO3-", 1}};
/// ReactionEquation equation2 = "-1:CO2(g) -1:H2O(l) 1:H+ 1:HCO3-";
/// ~~~~~~~~~~~~~~
/// Note that the stoichiometry of a species in a reaction follows the following sign
/// convention: *positive* for products, *negative* for reactants.
class ReactionEquation : public std::vector<std::pair<std::string, double>>
{
public:
    /// Define an alias for the base class
    using Base = std::vector<std::pair<std::string, double>>;

    /// Inherit all constructors of the base class
    using Base::Base;

    /// Construct a default ReactionEquation instance
    ReactionEquation();

    /// Construct a ReactionEquation instance by parsing a string.
    /// The string representing a reaction equation must have the format
    /// @c "N1:S1 N2:S2 ... Ni:Si", where @c Si represents a chemical
    /// species (e.g., H2O(l), CO2(aq), Calcite) and @c Ni the stoichiometry
    /// of species @c Si.
    /// @param formula The string representing the elemental formula of a species
    ReactionEquation(const std::string& equation);

    /// Construct a ReactionEquation instance from a list of species names and a list of stoichiometries
    /// @param species The names of the participating chemical species
    /// @param coeffs The stoichiometries of the participating chemical species
    ReactionEquation(const std::vector<std::string>& species, const std::vector<double>& stoichiometries);

    /// Return the stoichiometry of a species in the reaction equation.
    /// @param species The name of the species.
    auto stoichiometry(std::string species) const -> double;

    /// Convert the ReactionEquation instance into a string
    operator std::string() const;
};

/// Output a ReactionEquation instance
auto operator<<(std::ostream& out, const ReactionEquation& equation) -> std::ostream&;

} // namespace Reaktoro
