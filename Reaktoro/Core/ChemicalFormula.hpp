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

namespace Reaktoro {

/// A type used to represent the chemical formula of a chemical species.
class ChemicalFormula
{
public:
    /// Construct a default ChemicalFormula object.
    ChemicalFormula();

    /// Construct a ChemicalFormula object with given formula.
    /// @param formula The chemical formula of the species (e.g., `H2O`, `CaCO3`, `CO3--`, `CO3-2`).
    ChemicalFormula(const char* formula);

    /// Construct a ChemicalFormula object with given formula.
    /// @param formula The chemical formula of the species (e.g., `H2O`, `CaCO3`, `CO3--`, `CO3-2`).
    ChemicalFormula(String formula);

    /// Construct a ChemicalFormula object with given formula, element symbols, and charge.
    /// @param formula The chemical formula of the species (e.g., `HCO3-`).
    /// @param symbols The element symbols and their coefficients (e.g., `{{"H", 1}, {"C", 1}, {"O", 3}}` for `HCO3-`).
    /// @param charge The electric charge in the chemical formula (e.g., `-1` for `HCO3-`).
    ChemicalFormula(String formula, Pairs<String, double> symbols, double charge);

    /// Return the chemical formula of the substance as a string.
    auto str() const -> const String&;

    /// Return the element symbols and their coefficients in the chemical formula.
    auto elements() const -> const Pairs<String, double>&;

    /// Return the element symbols in the chemical formula.
    auto symbols() const -> Strings;

    /// Return the coefficients of the elements in the chemical formula.
    auto coefficients() const -> Vec<double>;

    /// Return the coefficient of an element symbol in the chemical formula.
    auto coefficient(const String& symbol) const -> double;

    /// Return the electric charge of the chemical formula.
    auto charge() const -> double;

    /// Return the molar mass of the chemical formula (in kg/mol).
    /// @warning An error is thrown if the formula contains symbols that do not
    /// exist in the periodic table.
    auto molarMass() const -> double;

    /// Return true if another chemical formula is equivalent to this one.
    /// Two chemical formulas are equivalent if they have the same elemental
    /// composition. Consider `Ca++` and `Ca+2` for example. These two formulas
    /// are equivalent, despite the different convention for indicating charge.
    /// Another example is `CaCO3` and `Ca(CO3)`. All these formulas have
    /// identical symbols and coefficients, and are thus equivalent.
    auto equivalent(const ChemicalFormula& other) const -> bool;

    /// Return true if two chemical formulas are equivalent.
    static auto equivalent(const ChemicalFormula& f1, const ChemicalFormula& f2) -> bool;

    /// Convert this ChemicalFormula object into a string.
    operator String() const;

    /// Convert this ChemicalFormula object into a Pairs<String, double>.
    operator Pairs<String, double>() const;

private:
    struct Impl;

    SharedPtr<Impl> pimpl;
};

/// Output a ChemicalFormula object
auto operator<<(std::ostream& out, const ChemicalFormula& formula) -> std::ostream&;

/// Compare two ChemicalFormula objects for less than
auto operator<(const ChemicalFormula& lhs, const ChemicalFormula& rhs) -> bool;

/// Compare two ChemicalFormula objects for equality
auto operator==(const ChemicalFormula& lhs, const ChemicalFormula& rhs) -> bool;

} // namespace Reaktoro
