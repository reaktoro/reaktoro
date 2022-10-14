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

#pragma once

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Core/SpeciesList.hpp>

namespace Reaktoro {

/// A type used to represent the equation of a reaction.
/// The equation of a reaction is assumed as a sequence of pairs (species,
/// stoichiometry). It is shown below how the equation of reaction
/// @f$\mathrm{CO_{2}(g)+H_{2}O\rightleftharpoons H^{+}+HCO_{3}^{-}}@f$ can be
/// defined by two equivalent ways:
/// ~~~~~~~~~~~~~~
/// ReactionEquation equation1 = {{"CO2(g)", -1}, {"H2O(l)", -1}, {"H+", 1}, {"HCO3-", 1}};
/// ReactionEquation equation2 = "CO2(g) + H2O(l) = H+ + HCO3-";
/// ~~~~~~~~~~~~~~
class ReactionEquation
{
public:
    /// Construct a default ReactionEquation object.
    ReactionEquation();

    /// Construct an ReactionEquation object with given species and respective stoichiometric coefficients.
    ReactionEquation(Pairs<Species, double> const& species);

    /// Construct a ReactionEquation object by parsing a formatted string.
    /// Below are examples of how to create a ReactionEquation object via a
    /// formatted string.
    /// ~~~
    /// ReactionEquation equation1("Calcite + H+ = Ca++ + HCO3-");
    /// ReactionEquation equation2("CO2(g) + H2O(l) = H+ + HCO3-");
    /// ReactionEquation equation3("Dolomite + 2*H+ = Ca++ + Mg++ + 2*HCO3-");
    /// ~~~
    /// Note that unity stoichiometric coefficients can be ommited from the
    /// equation. The operator `*` must be used when this is not the case.
    ///@{
    ReactionEquation(String const& equation);
    ReactionEquation(Chars equation);
    ///@}

    /// Construct a ReactionEquation object by parsing a formatted string.
    /// This method is similar to ReactionEquation(String const&), but the Species objects in the
    /// constructed ReactionEquation object are fetched from a given list of species.
    ReactionEquation(String const& equation, SpeciesList const& species);

    /// Return true if the rection equation is empty.
    auto empty() const -> bool;

    /// Return the number of species in the reaction equation.
    auto size() const -> Index;

    /// Return the species in the reaction equation.
    auto species() const -> Vec<Species>;

    /// Return the stoichiometric coefficients of the species in the reaction equation.
    auto coefficients() const -> Vec<double>;

    /// Return the stoichiometric coefficient of a species in the reaction equation.
    auto coefficient(const String& name) const -> double;

    /// Convert this ReactionEquation object into a string.
    operator String() const;

private:
    /// The species and their stoichiometric coefficients in the reaction equation.
    Pairs<Species, double> m_species;

public:
    /// Return begin const iterator of this ReactionEquation object (for STL compatibility reasons).
    inline auto begin() const { return m_species.begin(); }

    /// Return begin iterator of this ReactionEquation object (for STL compatibility reasons).
    inline auto begin() { return m_species.begin(); }

    /// Return end const iterator of this ReactionEquation object (for STL compatibility reasons).
    inline auto end() const { return m_species.end(); }

    /// Return end iterator of this ReactionEquation object (for STL compatibility reasons).
    inline auto end() { return m_species.end(); }
};

/// Return true if a Species object is less than another for sorting reasons.
auto operator<(const ReactionEquation& lhs, const ReactionEquation& rhs) -> bool;

/// Return true if two ReactionEquation objects are equal.
auto operator==(const ReactionEquation& lhs, const ReactionEquation& rhs) -> bool;

/// Output a ReactionEquation object
auto operator<<(std::ostream& out, const ReactionEquation& equation) -> std::ostream&;

} // namespace Reaktoro

