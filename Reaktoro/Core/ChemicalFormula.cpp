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

#include "ChemicalFormula.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Singletons/PeriodicTable.hpp>

namespace Reaktoro {
namespace detail {

using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

auto parseElementAtom(string::iterator begin, string::iterator end) -> pair<string, string::iterator>
{
    error(!isupper(*begin), "The first character in a chemical formula must be in uppercase.");
    if(begin == end) return {"", begin};
    auto endelement = find_if(begin + 1, end, [](char c){return isupper(c) || !isalpha(c);});
    string element = string(begin, endelement);
    return {element, endelement};
}

auto parseNumAtoms(string::iterator begin, string::iterator end) -> pair<double, string::iterator>
{
    if(begin == end) return {1.0, begin};
    if(!(isdigit(*begin) || *begin == '.')) return {1.0, begin};
    auto endnumber = find_if(begin, end, [](char c){return !(isdigit(c) || c == '.');});
    double number = atof(string(begin, endnumber).c_str());
    return {number, endnumber};
}

auto findMatchedParenthesis(string::iterator begin, string::iterator end) -> string::iterator
{
    if(begin == end) return end;
    int level = 0;
    for(auto iter = begin+1; iter < end; ++iter)
    {
        level = (*iter == '(') ? level + 1 : level;
        level = (*iter == ')') ? level - 1 : level;
        if(*iter == ')' && level == -1)
            return iter;
    }
    return end;
}

auto findMatchedParenthesisReverse(string::iterator begin, string::iterator end) -> string::iterator
{
    if(begin == end) return begin;
    int level = 0;
    for(auto iter = end-2; iter != begin; --iter)
    {
        level = (*iter == ')') ? level + 1 : level;
        level = (*iter == '(') ? level - 1 : level;
        if(*iter == '(' && level == -1)
            return iter;
    }
    return begin;
}

auto parseChemicalFormulaAux(string::iterator begin, string::iterator end, unordered_map<string, double>& result, double scalar) -> void
{
    if(begin == end) return;

    if(*begin != '(' && *begin != '.' && !isupper(*begin))
    {
        parseChemicalFormulaAux(begin + 1, end, result, scalar);
    }
    else if(*begin == '(')
    {
        string::iterator begin1 = begin + 1;
        string::iterator end1   = findMatchedParenthesis(begin, end);

        auto res = parseNumAtoms(end1 + 1, end);

        const double number = res.first;

        string::iterator begin2 = res.second;
        string::iterator end2   = end;

        parseChemicalFormulaAux(begin1, end1, result, scalar * number);
        parseChemicalFormulaAux(begin2, end2, result, scalar);
    }
    else if(*begin == '.')
    {
        auto res = parseNumAtoms(begin + 1, end);

        const double number = res.first;

        string::iterator begin1 = res.second;
        string::iterator end1   = end;

        parseChemicalFormulaAux(begin1, end1, result, scalar * number);
    }
    else
    {
        auto res1 = parseElementAtom(begin, end);
        auto res2 = parseNumAtoms(res1.second, end);

        string element = res1.first;
        double natoms  = res2.first;

        if(result.count(element)) result[element] += scalar * natoms;
        else result[element] = scalar * natoms;

        string::iterator begin1 = res2.second;
        string::iterator end1   = end;

        parseChemicalFormulaAux(begin1, end1, result, scalar);
    }
}

auto parseChemicalFormula(string formula) -> std::unordered_map<std::string, double>
{
    // Parse the formula for elements and their coefficients (without charge)
    unordered_map<string, double> result;
    parseChemicalFormulaAux(formula.begin(), formula.end(), result, 1.0);
    return result;
}

auto removeAggregateStateSuffix(string formula) -> string
{
    if(formula.back() != ')') return formula;

    string::iterator begin = formula.begin();
    string::iterator end = formula.end();
    begin = findMatchedParenthesisReverse(begin, end);

    for(auto iter = begin; iter < end; ++iter)
        if(isupper(*begin)) // no upper case char - otherwise, there are element symbols and thus not aggregate state
            return formula;

    const auto n = begin - formula.begin();  // H+(aq) => n = 2

    return formula.substr(0, n);
}

auto parseElectricChargeModeSignNumber(string formula) -> double
{
    size_t ipos = formula.find_last_of('+');
    size_t ineg = formula.find_last_of('-');
    size_t imin = std::min(ipos, ineg);

    if(imin == string::npos)
        return 0.0;

    int sign = (imin == ipos) ? +1 : -1;

    if(imin + 1 == formula.size())
        return sign;

    string digits = formula.substr(imin + 1);

    return sign * stod(digits);
}

auto parseElectricChargeModeMultipleSigns(string formula) -> double
{
    const auto sign = formula.back();
    const auto signval = sign == '+' ? 1 : (sign == '-' ? -1 : 0);
    const auto ilast = formula.size() - 1;
    if(sign != 0)
    {
        auto i = 0;
        while(formula[ilast - i] == sign && i <= ilast)
            ++i;
        const auto charge = i * signval;
        return charge;
    }
    else return 0;
}

auto parseElectricChargeModeNumberSign(string formula) -> double
{
    if(formula.back() != ')') return 0.0;

    size_t iparbegin = formula.rfind('(');

    if(iparbegin == string::npos) return 0.0;

    const auto isign = formula.size() - 2;
    const auto sign = formula[isign] == '+' ? +1.0 : formula[isign] == '-' ? -1.0 : 0.0;

    if(sign == 0.0) return 0.0;

    string digits(formula.begin() + iparbegin + 1, formula.end() - 2);

    if(digits.empty()) return sign;

    return sign * stod(digits);
}

auto parseElectricCharge(string formula) -> double
{
    double charge;

    formula = removeAggregateStateSuffix(formula);

    charge = parseElectricChargeModeMultipleSigns(formula); if(charge != 0.0) return charge;
    charge = parseElectricChargeModeNumberSign(formula); if(charge != 0.0) return charge;
    charge = parseElectricChargeModeSignNumber(formula); if(charge != 0.0) return charge;

    return 0.0;
}

} // namespace detail

auto parseChemicalFormula(const std::string& formula) -> std::unordered_map<std::string, double>
{
    return detail::parseChemicalFormula(formula);
}

auto parseElectricCharge(const std::string& formula) -> double
{
    return detail::parseElectricCharge(formula);
}

struct ChemicalFormula::Impl
{
    /// The chemical formula of the substance.
    std::string formula;

    /// The element symbols and their coefficients in the chemical formula, e.g., {{"H", 2}, {"O", 1}} for `H2O`.
    ElementSymbols symbols;

    /// The elements and their coefficients in the chemical formula.
    Elements elements;

    /// The electric charge in the chemical formula.
    double charge = {};

    /// The molar mass of the chemical formula.
    double molar_mass = {};

    /// Construct an object of type Impl.
    Impl()
    {
    }

    /// Construct an object of type Impl with given formula.
    Impl(std::string formula)
    : Impl(formula, parseChemicalFormula(formula))
    {
    }

    /// Construct an object of type Impl with given data.
    Impl(std::string formula, ElementSymbols symbols)
    : formula(formula), symbols(symbols)
    {
        elements.reserve(symbols.size());
        for(auto&& [symbol, coeff] : symbols)
        {
            const auto element = PeriodicTable::elementWithSymbol(symbol);
            error(!element.has_value(), "Cannot construct ChemicalFormula object with formula ",
                formula, ". PeriodicTable contains no element with symbol ", symbol, ". "
                "Use method PeriodicTable::append (in C++) or PeriodicTable.append (in Python) "
                "to add a new Element with this symbol.");
            elements.emplace_back(element.value(), coeff);
        }

        charge = parseElectricCharge(formula);
        molar_mass = 0.0;
        for(auto&& [element, coeff] : elements)
            molar_mass += element.molarMass() * coeff;
    }

    /// Return the coefficient of an element symbol in the chemical formula.
    auto coefficient(const std::string& symbol) const -> double
    {
        const auto iter = symbols.find(symbol);
        if(iter != symbols.end()) return iter->second;
        return 0.0;
    }
};

ChemicalFormula::ChemicalFormula()
: pimpl(new Impl())
{}

ChemicalFormula::ChemicalFormula(const char* formula)
: ChemicalFormula(std::string(formula))
{
}

ChemicalFormula::ChemicalFormula(std::string formula)
: pimpl(new Impl(formula))
{
}

ChemicalFormula::ChemicalFormula(std::string formula, ElementSymbols symbols)
: pimpl(new Impl(formula, symbols))
{
}

auto ChemicalFormula::str() const -> const std::string&
{
    return pimpl->formula;
}

auto ChemicalFormula::elements() const -> const Elements&
{
    return pimpl->elements;
}

auto ChemicalFormula::symbols() const -> const ElementSymbols&
{
    return pimpl->symbols;
}

auto ChemicalFormula::charge() const -> double
{
    return pimpl->charge;
}

auto ChemicalFormula::molarMass() const -> double
{
    return pimpl->molar_mass;
}

auto ChemicalFormula::coefficient(const std::string& symbol) const -> double
{
    return pimpl->coefficient(symbol);
}

auto ChemicalFormula::equivalent(const ChemicalFormula& other) const -> bool
{
    return symbols() == other.symbols() && charge() == other.charge();
}

auto ChemicalFormula::equivalent(const ChemicalFormula& f1, const ChemicalFormula& f2) -> bool
{
    return f1.equivalent(f2);
}

ChemicalFormula::operator std::string() const
{
    return str();
}

auto operator<(const ChemicalFormula& lhs, const ChemicalFormula& rhs) -> bool
{
    return lhs.str() < rhs.str();
}

auto operator==(const ChemicalFormula& lhs, const ChemicalFormula& rhs) -> bool
{
    return lhs.str() == rhs.str();
}

} // namespace Reaktoro
