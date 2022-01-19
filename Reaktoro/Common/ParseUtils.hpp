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

/// Parse a reaction equation
auto parseReaction(const String& equation) -> Pairs<String, double>;

/// Parse a formatted string containing pairs of numbers and strings.
/// See below several examples of parsing different formatted strings:
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// using namespace Reaktoro;
/// auto pairs0 = parseNumberStringPairs("2:H 1:O");
/// auto pairs1 = parseNumberStringPairs("1:Ca 2:Cl");
/// auto pairs2 = parseNumberStringPairs("1:Mg 1:C 3:O");
/// auto pairs3 = parseNumberStringPairs("1:Ca 1:Mg 2:C 6:O");
/// auto pairs4 = parseNumberStringPairs("3:Fe 2:Al 3:Si 12:O");
/// auto pairs5 = parseNumberStringPairs("1:Na+ 1:Cl-");
/// auto pairs6 = parseNumberStringPairs("1:Ca++ 1:CO3--");
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
auto parseNumberStringPairs(const String& str) -> Pairs<String, double>;

/// Return the element symbols and their coefficients in a chemical formula.
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
auto parseChemicalFormula(const String& formula) -> Pairs<String, double>;

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
auto parseElectricCharge(const String& formula) -> double;

} // namespace Reaktoro
