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

// C++ includes
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
using std::endl;
using std::pow;
using std::find_if;
using std::map;
using std::ostream;
using std::string;
using std::vector;
using std::shared_ptr;
using std::stringstream;

namespace units {
namespace internal {

inline std::size_t find(const char* s, const std::string& str, std::size_t pos)
{
    auto res = str.find(s, pos);
    return (res == str.npos) ? str.size() : res;
}

inline std::string numerator(std::string str, std::size_t pos)
{
    if(pos >= str.size()) return "";
    auto next = find("*", str, pos);
    if(str[pos] == '/') return numerator(str, next + 1);
    return str.substr(pos, next - pos).append("*").append(numerator(str, next + 1));
}

inline std::string denominator(std::string str, std::size_t pos)
{
    if(pos >= str.size()) return "";
    auto next = find("*", str, pos);
    if(str[pos] != '/') return denominator(str, next + 1);
    return str.substr(pos + 1, next - pos - 1).append("*").append(denominator(str, next + 1));
}

struct StringUnit
{
    double factor;
    string symbol;
    int power;
};

struct TemperatureUnit
{
    double factor;
    double translate;
    string symbol;
};

using DerivedUnit = vector<StringUnit>;

const double atto  = 1.0e-18;
const double femto = 1.0e-15;
const double pico  = 1.0e-12;
const double nano  = 1.0e-9;
const double micro = 1.0e-6;
const double milli = 1.0e-3;
const double centi = 1.0e-2;
const double deci  = 1.0e-1;
const double deca  = 1.0e+1;
const double hecto = 1.0e+2;
const double kilo  = 1.0e+3;
const double mega  = 1.0e+6;
const double giga  = 1.0e+9;
const double tera  = 1.0e+12;
const double peta  = 1.0e+15;
const double exa   = 1.0e+18;
const double pi    = 3.14159265359;

map<string, DerivedUnit> derivedUnitsMap =
{
    // Nondimensional Unit
    {"1"            , {{1, "1", 1}}},
    {"percent"      , {{centi, "1", 1}}},
    {"permille"     , {{milli, "1", 1}}},
    {"ppm"          , {{micro, "1", 1}}},
    {"ppb"          , {{nano,  "1", 1}}},
    {"ppt"          , {{pico, "1", 1}}},

    // Mass Units
    {"g"            , {{1, "g", 1}}},
    {"kg"           , {{kilo , "g", 1}}},
    {"mg"           , {{milli, "g", 1}}},
    {"mcg"          , {{micro, "g", 1}}},
    {"ug"           , {{micro, "g", 1}}},
    {"ng"           , {{nano , "g", 1}}},
    {"tonne"        , {{kilo , "kg" , 1}}},
    {"lb"           , {{453.59237, "g", 1}}},
    {"amu"          , {{1.660538921e-27, "kg" , 1}}},
    {"ton"          , {{2000 , "lb" , 1}}},
    {"gram"         , {{1, "g", 1}}},
    {"kilogram"     , {{1, "kg" , 1}}},
    {"lbm"          , {{1, "lb" , 1}}},

    // Length Units
    {"m"            , {{1, "m", 1}}},
    {"km"           , {{kilo , "m", 1}}},
    {"dm"           , {{deci, "m", 1}}},
    {"cm"           , {{centi, "m", 1}}},
    {"mm"           , {{milli, "m", 1}}},
    {"um"           , {{micro, "m", 1}}},
    {"nm"           , {{nano , "m", 1}}},
    {"Angstrom"     , {{1.0e-10, "m", 1}}},
    {"in"           , {{25.4 , "mm" , 1}}},
    {"mil"          , {{milli, "in" , 1}}},
    {"ft"           , {{12 , "in" , 1}}},
    {"yd"           , {{3, "ft" , 1}}},
    {"mi"           , {{5280 , "ft" , 1}}},
    {"nauticalmile" , {{1852 , "m", 1}}},
    {"lightyear"    , {{9460730472580800 , "m", 1}}},
    {"au"           , {{149597870700 , "m", 1}}},
    {"feet"         , {{1, "ft" , 1}}},
    {"foot"         , {{1, "ft" , 1}}},
    {"inch"         , {{1, "in" , 1}}},
    {"meter"        , {{1, "m", 1}}},
    {"mile"         , {{1, "mi" , 1}}},
    {"yard"         , {{1, "yd" , 1}}},

    // Time Units
    {"s"            , {{1, "s", 1}}},
    {"ms"           , {{milli, "s", 1}}},
    {"us"           , {{micro, "s", 1}}},
    {"ns"           , {{nano , "s", 1}}},
    {"minute"       , {{60 , "s", 1}}},
    {"min"          , {{60 , "s", 1}}},
    {"hour"         , {{3600 , "s", 1}}},
    {"h"            , {{3600 , "s", 1}}},
    {"day"          , {{24 , "hour" , 1}}},
    {"week"         , {{7, "day", 1}}},
    {"month"        , {{30.43685 , "day", 1}}},
    {"year"         , {{365.2422 , "day", 1}}},
    {"yr"           , {{365.2422 , "day", 1}}},
    {"second"       , {{1, "s", 1}}},

    // Electrical Current Units
    {"A"            , {{1, "A", 1}}},
    {"MA"           , {{mega , "A", 1}}},
    {"kA"           , {{kilo , "A", 1}}},
    {"mA"           , {{milli, "A", 1}}},
    {"uA"           , {{micro, "A", 1}}},
    {"nA"           , {{nano , "A", 1}}},
    {"pA"           , {{pico , "A", 1}}},
    {"ampere"       , {{1, "A", 1}}},

    // Amount of Substance Units
    {"mol"          , {{1, "mol", 1}}},
    {"mmol"         , {{milli, "mol", 1}}},
    {"umol"         , {{micro, "mol", 1}}},
    {"mcmol"        , {{micro, "mol", 1}}},
    {"kmol"         , {{kilo , "mol", 1}}},
    {"mole"         , {{1, "mol", 1}}},

    // Luminous Intensity Units
    {"cd"           , {{1, "cd" , 1}}},
    {"candela"      , {{1, "cd" , 1}}},

    // Angle Units
    {"radian"       , {{1, "radian" , 1}}},
    {"mradian"      , {{milli, "radian" , 1}}},
    {"uradian"      , {{micro, "radian" , 1}}},
    {"deg"          , {{pi/180 , "radian" , 1}}},
    {"degree"       , {{1, "deg", 1}}},

    // Solid Angle Units
    {"sr"           , {{1, "sr" , 1}}},
    {"steradian"    , {{1, "sr" , 1}}},

    // Area Units
    {"m2"           , {{1, "m", 2}}},
    {"m^2"          , {{1, "m", 2}}},
    {"nm2"          , {{1, "nm" , 2}}},
    {"nm^2"         , {{1, "nm" , 2}}},
    {"um2"          , {{1, "um" , 2}}},
    {"um^2"         , {{1, "um" , 2}}},
    {"mm2"          , {{1, "mm" , 2}}},
    {"mm^2"         , {{1, "mm" , 2}}},
    {"cm2"          , {{1, "cm" , 2}}},
    {"cm^2"         , {{1, "cm" , 2}}},
    {"dm2"          , {{1, "dm" , 2}}},
    {"dm^2"         , {{1, "dm" , 2}}},
    {"in2"          , {{1, "in" , 2}}},
    {"in^2"         , {{1, "in" , 2}}},
    {"ft2"          , {{1, "ft" , 2}}},
    {"ft^2"         , {{1, "ft" , 2}}},
    {"hectare"      , {{10000, "m2" , 1}}},
    {"acre"         , {{4046.8726, "m2" , 1}}},
    {"barn"         , {{1.0e-28, "m2" , 1}}},

    // Volume Units
    {"m3"           , {{1, "m", 3}}},
    {"m^3"          , {{1, "m", 3}}},
    {"nm3"          , {{1, "nm" , 3}}},
    {"nm^3"         , {{1, "nm" , 3}}},
    {"um3"          , {{1, "um" , 3}}},
    {"um^3"         , {{1, "um" , 3}}},
    {"mm3"          , {{1, "mm" , 3}}},
    {"mm^3"         , {{1, "mm" , 3}}},
    {"cm3"          , {{1, "cm" , 3}}},
    {"cm^3"         , {{1, "cm" , 3}}},
    {"dm3"          , {{1, "dm" , 3}}},
    {"dm^3"         , {{1, "dm" , 3}}},
    {"km3"          , {{1, "km" , 3}}},
    {"km^3"         , {{1, "km" , 3}}},
    {"in3"          , {{1, "in" , 3}}},
    {"in^3"         , {{1, "in" , 3}}},
    {"ft3"          , {{1, "ft" , 3}}},
    {"ft^3"         , {{1, "ft" , 3}}},
    {"l"            , {{milli, "m3" , 1}}},
    {"ml"           , {{milli, "l", 1}}},
    {"ul"           , {{micro, "l", 1}}},
    {"mcl"          , {{micro, "l", 1}}},
    {"L"            , {{milli, "m3" , 1}}},
    {"mL"           , {{milli, "L", 1}}},
    {"uL"           , {{micro, "L", 1}}},
    {"mcL"          , {{micro, "l", 1}}},
    {"gal"          , {{4.54609, "l", 1}}},
    {"quart"        , {{0.25 , "gal", 1}}},
    {"cc"           , {{1, "cm3", 1}}},
    {"gallon"       , {{1, "gal", 1}}},
    {"liter"        , {{1, "l", 1}}},
    {"litre"        , {{1, "l", 1}}},

    // Velocity Units
    {"mps"          , {{1 , "m", 1}, {1 , "s", -1}}},
    {"cmps"         , {{1 , "cm" , 1}, {1 , "s", -1}}},
    {"mph"          , {{1 , "mile" , 1}, {1 , "hour" , -1}}},
    {"knot"         , {{1 , "nauticalmile" , 1}, {1 , "hour" , -1}}},
    {"fps"          , {{1 , "foot" , 1}, {1 , "s", -1}}},

    // Acceleration Units
    {"mpss"         , {{1 , "mps", 1}, {1 , "s", -1}}},
    {"fpss"         , {{1 , "fps", 1}, {1 , "s", -1}}},
    {"gravity"      , {{9.80665, "mpss" , 1}}},

    // Frequency Units
    {"Hz"           , {{1, "s", -1}}},
    {"kHz"          , {{kilo , "Hz" , 1}}},
    {"MHz"          , {{mega , "Hz" , 1}}},
    {"GHz"          , {{giga , "Hz" , 1}}},
    {"rpm"          , {{1, "minute" , -1}}},
    {"hertz"        , {{1, "Hz" , 1}}},

    // Flow Rate Units
    {"gps"          , {{1 , "gal", 1}, {1 , "s", -1}}},
    {"gpm"          , {{1 , "gal", 1}, {1 , "minute" , -1}}},
    {"gph"          , {{1 , "gal", 1}, {1 , "hour" , -1}}},
    {"cfs"          , {{1 , "ft3", 1}, {1 , "s", -1}}},
    {"cfm"          , {{1 , "ft3", 1}, {1 , "minute" , -1}}},
    {"lps"          , {{1 , "liter", 1}, {1 , "s", -1}}},
    {"lpm"          , {{1 , "liter", 1}, {1 , "minute" , -1}}},

    // Force Units
    {"N"            , {{1 , "kg" , 1}, {1 , "m", 1}, {1 , "s" , -2}}},
    {"kN"           , {{kilo, "N", 1}}},
    {"lbf"          , {{4.44822162 , "N", 1}}},
    {"dyne"         , {{1.0e-5 , "N", 1}}},
    {"newton"       , {{1, "N", 1}}},

    // Pressure Units
    {"Pa"           , {{1 , "N", 1}, {1 , "m2" , -1}}},
    {"kPa"          , {{kilo , "Pa" , 1}}},
    {"MPa"          , {{mega , "Pa" , 1}}},
    {"GPa"          , {{giga , "Pa" , 1}}},
    {"atm"          , {{101325 , "Pa" , 1}}},
    {"mmHg"         , {{133.32239, "Pa" , 1}}},
    {"inHg"         , {{3386.3886, "Pa" , 1}}},
    {"psi"          , {{1 , "lbf", 1}, {1 , "in2", -1}}},
    {"kpsi"         , {{kilo , "psi", 1}}},
    {"Mpsi"         , {{mega , "psi", 1}}},
    {"psf"          , {{1 , "lbf", 1}, {1 , "ft2", -1}}},
    {"bar"          , {{100000 , "Pa" , 1}}},
    {"mbar"         , {{100 , "Pa" , 1}}},
    {"kbar"         , {{1000 , "bar" , 1}}},
    {"Mbar"         , {{1000000 , "bar" , 1}}},
    {"torr"         , {{133.322368 , "Pa" , 1}}},
    {"inH2O"        , {{249.08891, "Pa" , 1}}},
    {"ftH2O"        , {{12 , "inH2O", 1}}},
    {"pascal"       , {{1, "Pa" , 1}}},

    // Energy Units
    {"J"            , {{1 , "N", 1}, {1 , "m", 1}}},
    {"kJ"           , {{kilo, "J", 1}}},
    {"MJ"           , {{mega, "J", 1}}},
    {"erg"          , {{1 , "dyne" , 1}, {1 , "cm" , 1}}},
    {"cal"          , {{4.184 , "J", 1}}},
    {"kcal"         , {{kilo , "cal", 1}}},
    {"eV"           , {{1.60217657e-19 , "J", 1}}},
    {"keV"          , {{kilo , "eV" , 1}}},
    {"MeV"          , {{mega , "eV" , 1}}},
    {"GeV"          , {{giga , "eV" , 1}}},
    {"btu"          , {{1055.05585 , "J", 1}}},
    {"calorie"      , {{1, "cal", 1}}},
    {"joule"        , {{1, "J", 1}}},

    // Power Units
    {"W"            , {{1 , "J", 1}, {1 , "s", -1}}},
    {"watt"        , {{1, "W", 1}}},

    // Charge Units
    {"C"            , {{1 , "A", 1}, {1 , "s", 1}}},
    {"coulomb"      , {{1, "C", 1}}},
    {"mC"           , {{milli, "C", 1}}},

    // Voltage Units
    {"V"            , {{1 , "J", 1}, {1 , "C", -1}}},
    {"volt"         , {{1, "V", 1}}},
    {"mV"           , {{milli, "V", 1}}},

    // Resistance Units
    {"ohm"          , {{1 , "V", 1}, {1 , "A", -1}}},

    // Conductance Units
    {"S"            , {{1, "ohm", -1}}},
    {"siemens"      , {{1, "S", 1}}},

    // Capacitance Units
    {"F"            , {{1 , "C", 1}, {1 , "V", -1}}},
    {"pF"           , {{pico , "F", 1}}},
    {"nF"           , {{nano , "F", 1}}},
    {"uF"           , {{micro, "F", 1}}},
    {"mF"           , {{milli, "F", 1}}},
    {"farad"        , {{1, "F", 1}}},

    // Magnetic Flux Units
    {"Wb"           , {{1 , "V", 1}, {1 , "s", 1}}},
    {"weber"        , {{1, "Wb" , 1}}},

    // Magnetic Flux Density Units
    {"T"            , {{1 , "Wb" , 1}, {1 , "m", 2}}},
    {"gauss"        , {{0.0001 , "T", 1}}},
    {"tesla"        , {{1, "T", 1}}},

    // Inductance Units
    {"H"            , {{1 , "Wb" , 1}, {1 , "A", -1}}},
    {"mH"           , {{milli, "H", 1}}},
    {"uH"           , {{micro, "H", 1}}},
    {"henry"        , {{1, "H", 1}}},

    // Molality Units
    {"molal"        , {{1 , "mol", 1}, {1 , "kg" , -1}}},
    {"mmolal"       , {{milli, "molal", 1}}},
    {"umolal"       , {{micro, "molal", 1}}},
    {"mcmolal"      , {{micro, "molal", 1}}},
    {"nmolal"       , {{nano, "molal", 1}}},

    // Molarity Units
    {"molar"        , {{1 , "mol", 1}, {1 , "liter" , -1}}},
    {"mmolar"       , {{milli, "molar", 1}}},
    {"umolar"       , {{micro, "molar", 1}}},
    {"mcmolar"      , {{micro, "molar", 1}}},
    {"nmolar"       , {{nano, "molar", 1}}},

    // Equivalent Units
    {"eq"           , {{1, "eq", 1}}},
    {"meq"          , {{milli, "eq", 1}}},
    {"ueq"          , {{micro, "eq", 1}}},
    {"mceq"         , {{micro, "eq", 1}}},
    {"neq"          , {{nano, "eq", 1}}}
};

map<string, TemperatureUnit> temperatureUnitsMap =
{
    {"K"          , {1, 0, "K"}},
    {"degC"       , {1, -273.15, "K"}},
    {"degF"       , {1.8, 32, "degC"}},
    {"degR"       , {1, 459.67, "degF"}},
    {"kelvin"     , {1, 0, "K"}},
    {"celsius"    , {1, 0, "degC"}},
    {"fahrenheit" , {1, 0, "degF"}},
    {"rankine"    , {1, 0, "degR"}}
};

double toKelvin(double value, const string& from)
{
    if(from == "K") return value;
    auto unit = temperatureUnitsMap[from];
    return toKelvin((value - unit.translate)/unit.factor, unit.symbol);
}

double fromKelvin(double value, const string& to)
{
    if(to == "K") return value;
    auto unit = temperatureUnitsMap[to];
    return unit.factor * fromKelvin(value, unit.symbol) + unit.translate;
}

inline void checkTemperatureUnit(const string& symbol)
{
    if(!temperatureUnitsMap.count(symbol))
    {
        stringstream error; error << "*** Error *** there is no such temperature unit named: " << symbol << ".";
        throw std::runtime_error(error.str());
    }
}

inline void checkDerivedUnit(const string& symbol)
{
    if(!derivedUnitsMap.count(symbol))
    {
        stringstream error; error << "*** Error *** there is no such unit named: " << symbol << ".";
        throw std::runtime_error(error.str());
    }
}

double convertTemperature(double value, const string& from, const string& to)
{
    checkTemperatureUnit(from);
    checkTemperatureUnit(to);

    return fromKelvin(toKelvin(value, from), to);
}

double factor(const string& symbol)
{
    if(temperatureUnitsMap.count(symbol)) return 1.0;
    checkDerivedUnit(symbol);
    const auto& units = derivedUnitsMap[symbol];
    if(units.front().symbol == symbol) return pow(units.front().factor, units.front().power);
    double res = 1.0;
    for(const auto& unit : units)
        res *= unit.factor * pow(factor(unit.symbol), unit.power);
    return res;
}

double factor(const DerivedUnit& derivedUnit, unsigned pos)
{
    if(pos == derivedUnit.size()) return 1.0;

    return derivedUnit[pos].factor * pow(factor(derivedUnit[pos].symbol), derivedUnit[pos].power) * factor(derivedUnit, pos+1);
}

double factor(const DerivedUnit& derivedUnit)
{
    return factor(derivedUnit, 0);
}

void removeZero(map<string, int>& res)
{
    for(auto iter = res.begin(); iter != res.end();)
        if(iter->second == 0) res.erase(iter++); else ++iter;
}

void dimension(const string& symbol, int multpower, map<string, int>& res)
{
    if(temperatureUnitsMap.count(symbol)) { res[symbol] += multpower; return; }
    checkDerivedUnit(symbol);
    const auto& units = derivedUnitsMap[symbol];
    if(units.front().symbol == symbol) res[units.front().symbol] += units.front().power * multpower;
    else for(const auto& unit : units)
        dimension(unit.symbol, unit.power * multpower, res);
}

void dimension(const DerivedUnit& derivedUnit, unsigned pos, map<string, int>& res)
{
    if(pos == derivedUnit.size()) return;
    dimension(derivedUnit[pos].symbol, derivedUnit[pos].power, res);
    dimension(derivedUnit, pos + 1, res);
}

map<string, int> dimension(const string& symbol)
{
    map<string, int> res;
    dimension(symbol, 1, res);
    removeZero(res);
    return res;
}

map<string, int> dimension(const DerivedUnit& derivedUnit)
{
    map<string, int> res;
    dimension(derivedUnit, 0, res);
    removeZero(res);
    return res;
}

struct Node
{
    Node() : left(nullptr), right(nullptr) {}
    Node(string str) : str(str), left(nullptr), right(nullptr) {}

    string str;
    std::shared_ptr<Node> left;
    std::shared_ptr<Node> right;
};

std::size_t findMatchedParenthesisFromRight(const string& str, unsigned pos)
{
    int level = 0;
    for(int iter = pos-1; iter >= 0; --iter)
    {
        level = (str[iter] == ')') ? level + 1 : level;
        level = (str[iter] == '(') ? level - 1 : level;
        if(str[iter] == '(' && level == -1)
            return iter;
    }
    return string::npos;
}

inline void unmatchedParenthesisError()
{
    throw std::runtime_error("*** Error *** there is unmatched parenthesis in the unit string.");
}

std::shared_ptr<Node> parseUnit(const string& symbol, unsigned pos)
{
    if(symbol.empty()) return nullptr;

    if(symbol[pos] == ')')
    {
        auto i = findMatchedParenthesisFromRight(symbol, pos);
        if(i == string::npos) unmatchedParenthesisError();
        if(i > 0)
        {
            string op_str    = symbol.substr(i-1, 1);
            string left_str  = symbol.substr(0, i-1);
            string right_str = symbol.substr(i+1, pos-i-1);

            std::shared_ptr<Node> parent(new Node(op_str));
            parent->left  = parseUnit(left_str, left_str.size()-1);
            parent->right = parseUnit(right_str, right_str.size()-1);
            return parent;
        }
        else
        {
            string right_str = symbol.substr(i+1, pos-i-1);
            return parseUnit(right_str, right_str.size()-1);
        }
    }
    else
    {
        auto found = symbol.find_last_of("*/", pos);

        if(found == string::npos)
            return std::shared_ptr<Node>(new Node(symbol.substr(0, pos+1)));

        std::shared_ptr<Node> right(new Node());
        right->str = symbol.substr(found+1, pos-found);

        std::shared_ptr<Node> parent(new Node());
        parent->str = symbol[found];
        parent->right = right;
        parent->left  = parseUnit(symbol, found-1);

        return parent;
    }
}

void parseUnit(const std::shared_ptr<Node>& root, DerivedUnit& derivedUnit, int sign)
{
    if(root->str == "*")
    {
        parseUnit(root->left,  derivedUnit, sign);
        parseUnit(root->right, derivedUnit, sign);
    }
    else if(root->str == "/")
    {
        parseUnit(root->left,  derivedUnit,  sign);
        parseUnit(root->right, derivedUnit, -sign);
    }
    else
    {
        derivedUnit.push_back({1, root->str, sign});
    }
}

DerivedUnit parseUnit(const string& symbol)
{
    std::shared_ptr<Node> root = parseUnit(symbol, symbol.size()-1);
    DerivedUnit derivedUnit;
    parseUnit(root, derivedUnit, 1);
    return derivedUnit;
}

inline void checkConvertibleUnits(const DerivedUnit& from, const DerivedUnit& to, const string& strfrom, const string& strto)
{
    if(dimension(from) != dimension(to))
    {
        stringstream error; error << "*** Error *** the dimensions of the units " << strfrom << " and " << strto << " do not match.";
        throw std::runtime_error(error.str());
    }
}

} // namespace internal

auto slope(const std::string& from, const std::string& to) -> double
{
    if(internal::temperatureUnitsMap.count(from) && internal::temperatureUnitsMap.count(to))
        return internal::convertTemperature(1.0, from, to) - internal::convertTemperature(0.0, from, to);
    auto parsed_from = internal::parseUnit(from);
    auto parsed_to   = internal::parseUnit(to);
    internal::checkConvertibleUnits(parsed_from, parsed_to, from, to);
    return internal::factor(parsed_from)/internal::factor(parsed_to);
}

auto intercept(const std::string& from, const std::string& to) -> double
{
    if(internal::temperatureUnitsMap.count(from) && internal::temperatureUnitsMap.count(to))
        return internal::convertTemperature(0.0, from, to);
    return 0.0;
}

bool convertible(const std::string& from, const std::string& to)
{
    if(internal::temperatureUnitsMap.count(from) && internal::temperatureUnitsMap.count(to))
        return true;
    auto parsed_from = internal::parseUnit(from);
    auto parsed_to   = internal::parseUnit(to);
    return dimension(parsed_from) == dimension(parsed_to);
}

} // namespace units
