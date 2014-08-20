/*
 * StringUtils.hpp
 *
 *  Created on: 28 Oct 2013
 *      Author: allan
 */

#pragma once

// C++ includes
#include <string>

// Units++ includes
#include "Unit.hpp"

namespace units {

/**
 * Converts a numeric value from a unit to another
 * @param value The value
 * @param from The string representing the unit from which the conversion is made
 * @param to The string representing the unit to which the conversion is made
 * @return The converted value
 */
double convert(double value, const std::string& from, const std::string& to);

/**
 * Checks if two units are convertible among each other
 * @return True if they are convertible, false otherwise
 */
bool convertible(const std::string& from, const std::string& to);

/**
 * Stringfy a strongly typed unit
 */
template<typename U>
inline std::string symbol(Unit<U> u);

namespace internal {

template<typename U>
inline std::string str(Unit<U> u)
{
    return symbol(U()); // use the symbol function defined in the macro DefineUnit
}

template<typename U>
inline std::string str(Inv<U> u)
{
    return std::string("/").append(str(U()));
}

template<typename U1, typename U2>
inline std::string str(Mult<U1, U2> c)
{
    return str(U1()).append("*").append(str(U2()));
}

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

} /* namespace internal */

template<typename U>
inline std::string symbol(Unit<U> u)
{
    auto str = internal::str(U());            // the intermediate string representation of unit U
    auto num = internal::numerator(str, 0);   // form the string containing only the numerator units
    auto den = internal::denominator(str, 0); // form the string containing only the denominator units

    num = num.substr(0, num.size()-1); // remove the leading dot '.'
    den = den.substr(0, den.size()-1); // remove the leading dot '.'

    if(den.empty()) return num;
    else return num.append("/(").append(den).append(")");
}

} /* namespace units */
