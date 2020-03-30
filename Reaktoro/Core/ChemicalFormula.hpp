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
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace Reaktoro {

// Forward declarations
class Element;

/// A type used to represent the chemical formula of a chemical species.
class ChemicalFormula
{
public:
    /// Alias type for a vector of element symbols and their coefficients in a chemical formula.
    using ElementSymbols = std::vector<std::pair<std::string, double>>;

    /// Alias type for a vector of elements and their coefficients in a chemical formula.
    using Elements = std::vector<std::pair<Element, double>>;

    /// Construct a default ChemicalFormula object.
    ChemicalFormula();

    /// Construct a ChemicalFormula object with given formula.
    /// @param formula The chemical formula of the species (e.g., `H2O`, `CaCO3`, `CO3--`, `CO3-2`).
    explicit ChemicalFormula(const char* formula);

    /// Construct a ChemicalFormula object with given formula.
    /// @param formula The chemical formula of the species (e.g., `H2O`, `CaCO3`, `CO3--`, `CO3-2`).
    explicit ChemicalFormula(std::string formula);

    /// Construct a ChemicalFormula object with given formula and parsed formula.
    /// @param formula The chemical formula of the species (e.g., `H2O`, `CaCO3`, `CO3--`, `CO3-2`).
    /// @param elements The element symbols and their coefficients in the formula, e.g., {{"H", 2}, {"O", 1}} for `H2O`.
    explicit ChemicalFormula(std::string formula, ElementSymbols symbols);

    /// Return the chemical formula of the substance as a string.
    auto str() const -> const std::string&;

    /// Return the elements in the chemical formula and their coefficients.
    auto elements() const -> const Elements&;

    /// Return the electric charge of the chemical formula.
    auto charge() const -> double;

    /// Return the molar mass of the chemical formula (in kg/mol).
    auto molarMass() const -> double;

    /// Return the coefficient of an element symbol in the chemical formula.
    auto coefficient(const std::string& symbol) const -> double;

    /// Return true if another chemical formula is equivalent to this one.
    /// Equivalency is defined as two formulas with the same elemental composition.
    /// For example, `Ca++` and `Ca+2`; and `CaCO3` and `Ca(CO3)`.
    auto equivalent(const ChemicalFormula& other) const -> bool;

    /// Convert this ChemicalFormula object into a string.
    operator std::string() const;

private:
    struct Impl;

    std::shared_ptr<Impl> pimpl;
};

/// Compare two ChemicalFormula objects for less than
auto operator<(const ChemicalFormula& lhs, const ChemicalFormula& rhs) -> bool;

/// Compare two ChemicalFormula objects for equality
auto operator==(const ChemicalFormula& lhs, const ChemicalFormula& rhs) -> bool;

/// Return a map from element symbol to its coefficient in a chemical formula.
/// Successfully parsing a chemical formula requires that the first letter in
/// the formula is uppercase and all others in lowercase. Thus, even chemical
/// formulas such as `AaBbb` or `(Aa2Bbb4)Cc6` are supported.
/// There are two ways for specifying electric charge in a chemical formula:
///   1. as a suffix containing as many symbols `+` and `-` as there are charges (e.g., `Fe+++`, `Ca++`, `CO3--`); or
///   2. as a suffix containing the symbol `+` or `-` followed by the number of charges (e.g., `Fe+3`, `Ca+2`, `Na+`)
/// Note that number 1 is optional for the second format (e.g., `Na+` and `Na+1` are equivalent).
/// In both formats, (1) and (2), the symbol `+` is used for positively charged substances, and `-` for negatively charged ones.
/// See below several examples of parsing different chemical formulas:
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// using namespace Reaktoro;
/// auto formula01 = parseChemicalFormula("H2O");
/// auto formula02 = parseChemicalFormula("CaCl2");
/// auto formula03 = parseChemicalFormula("MgCO3");
/// auto formula04 = parseChemicalFormula("(CaMg)(CO3)2");
/// auto formula05 = parseChemicalFormula("Fe3Al2Si3O12");
/// auto formula06 = parseChemicalFormula("Na+");
/// auto formula07 = parseChemicalFormula("Ca++");
/// auto formula08 = parseChemicalFormula("Fe+++");
/// auto formula09 = parseChemicalFormula("Fe+3");
/// auto formula10 = parseChemicalFormula("CO3--");
/// auto formula11 = parseChemicalFormula("CO3-2");
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
auto parseChemicalFormula(const std::string& formula) -> std::vector<std::pair<std::string, double>>;

/// Return the electric charge in a chemical formula.
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// using namespace Reaktoro;
/// auto charge01 = parseElectricCharge("H2O");          \\ charge01 ===  0
/// auto charge02 = parseElectricCharge("CaCl2");        \\ charge02 ===  0
/// auto charge03 = parseElectricCharge("MgCO3");        \\ charge03 ===  0
/// auto charge04 = parseElectricCharge("(CaMg)(CO3)2"); \\ charge04 ===  0
/// auto charge05 = parseElectricCharge("Fe3Al2Si3O12"); \\ charge05 ===  0
/// auto charge06 = parseElectricCharge("Na+");          \\ charge06 ===  1
/// auto charge07 = parseElectricCharge("Ca++");         \\ charge07 ===  2
/// auto charge08 = parseElectricCharge("Fe+++");        \\ charge08 ===  3
/// auto charge09 = parseElectricCharge("Fe+3");         \\ charge09 ===  3
/// auto charge10 = parseElectricCharge("CO3--");        \\ charge10 === -2
/// auto charge11 = parseElectricCharge("CO3-2");        \\ charge11 === -2
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// @see parseChemicalFormula
auto parseElectricCharge(const std::string& formula) -> double;

/// Return true if the two given chemical formulas are equivalent.
auto equivalentChemicalFormulas(const std::string& formula1, const std::string& formula2) -> bool;

} // namespace Reaktoro
