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

// C++ includes
#include <functional>
#include <sstream>
#include <string>
#include <vector>
#include <variant>

namespace Reaktoro {
namespace internal {

template<typename T>
auto operator<<(std::ostream& out, std::vector<T> const& values) -> std::ostream&
{
    for(auto i = 0; i < values.size(); ++i)
        out << ((i == 0) ? "" : ", ") << values[i];
    return out;
}

template <typename Arg>
auto stringfy(std::ostringstream& ss, std::string const& sep, Arg const& item) -> void
{
    ss << item;
}

template <typename Arg, typename... Args>
auto stringfy(std::ostringstream& ss, std::string const& sep, Arg const& item, Args... items) -> void
{
    ss << item << sep;
    stringfy(ss, sep, items...);
}

} // namespace internal

/// Concatenate the arguments into a string using a given separator string.
template <typename... Args>
auto stringfy(std::string const& sep, Args... items) -> std::string
{
    std::ostringstream ss;
    internal::stringfy(ss, sep, items...);
    return ss.str();
}

/// Stringfy a `std::variant` object.
template <typename... Args>
auto stringfy(std::variant<Args...> const& var) -> std::string
{
    std::ostringstream ss;
    std::visit([&](auto arg){ ss << arg; }, var);
    return ss.str();
}

/// Concatenate the arguments into a string without any separator string.
template <typename... Args>
auto str(Args... items) -> std::string
{
    return stringfy("", items...);
}

/// Set the global precision used when converting floating-point numbers to string.
auto precision(int precision) -> void;

/// Return the precision used when converting floating-point numbers to string.
auto precision() -> int;

/// Return a string representation for a number in fixed format.
/// @param num The number to be converted to a string
/// @param precision The precision in the conversion (negative value results in global precision using method @ref precision)
auto strfix(double num, int precision = -1) -> std::string;

/// Return a string representation for a number in scientific format.
/// @param num The number to be converted to a string
/// @param precision The precision in the conversion (negative value results in global precision using method @ref precision)
auto strsci(double num, int precision = -1) -> std::string;

/// Return a new string where `substr` occurrences are replaced by `newsubstr`.
auto replace(std::string original, std::string substr, std::string newsubstr) -> std::string;

/// Return a string with lower case characters.
auto lowercase(std::string str) -> std::string;

/// Return a string with upper case characters.
auto uppercase(std::string str) -> std::string;

/// Trim the string from start (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
auto trimleft(std::string str) -> std::string;

/// Trim the string from end (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
auto trimright(std::string str) -> std::string;

/// Trim the string from both ends (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
auto trim(std::string str) -> std::string;

/// Split the string on every occurrence of the specified delimiters and apply a transform function.
auto split(std::string const& str, std::string const& delims, std::function<std::string(std::string)> transform) -> std::vector<std::string>;

/// Split the string on every occurrence of the specified delimiters
auto split(std::string const& str, std::string const& delims = " ") -> std::vector<std::string>;

/// Join several strings into one.
auto join(std::vector<std::string> const& strs, std::string sep = " ") -> std::string;

/// Convert the string into a floating point number
auto tofloat(std::string const& str) -> double;

/// Return a list of words with duplicate names converted to unique ones.
auto makeunique(std::vector<std::string> words, std::string suffix) -> std::vector<std::string>;

/// Returns true if string `str` starts with `substr`, or any other given sub string in `substrs`.
template<typename SubStr, typename... SubStrs>
auto startswith(std::string const& str, SubStr substr, SubStrs... substrs)
{
    if constexpr(sizeof...(SubStrs) > 0)
        return str.find(substr) == 0 || startswith(str, substrs...);
    return str.find(substr) == 0;
}

} // namespace Reaktoro
