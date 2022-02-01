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

namespace Reaktoro {
namespace internal {

template<typename T>
auto operator<<(std::ostream& out, const std::vector<T>& values) -> std::ostream&
{
    for(auto i = 0; i < values.size(); ++i)
        out << ((i == 0) ? "" : ", ") << values[i];
    return out;
}

template <typename Arg>
auto stringfy(std::ostringstream& ss, const std::string& sep, const Arg& item) -> void
{
    ss << item;
}

template <typename Arg, typename... Args>
auto stringfy(std::ostringstream& ss, const std::string& sep, const Arg& item, Args... items) -> void
{
    ss << item << sep;
    stringfy(ss, sep, items...);
}

} // namespace internal

/// Concatenate the arguments into a string using a given separator string.
template <typename... Args>
auto stringfy(const std::string& sep, Args... items) -> std::string
{
    std::ostringstream ss;
    internal::stringfy(ss, sep, items...);
    return ss.str();
}

/// Concatenate the arguments into a string without any separator string.
template <typename... Args>
auto str(Args... items) -> std::string
{
    return stringfy("", items...);
}

/// Return a string representation for a number in fixed format.
auto strfix(double num, int precision = 6) -> std::string;

/// Return a string representation for a number in scientific format.
auto strsci(double num, int precision = 6) -> std::string;

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
auto split(const std::string& str, const std::string& delims, std::function<std::string(std::string)> transform) -> std::vector<std::string>;

/// Split the string on every occurrence of the specified delimiters
auto split(const std::string& str, const std::string& delims = " ") -> std::vector<std::string>;

/// Join several strings into one.
auto join(const std::vector<std::string>& strs, std::string sep = " ") -> std::string;

/// Convert the string into a floating point number
auto tofloat(const std::string& str) -> double;

/// Return a list of words with duplicate names converted to unique ones.
auto makeunique(std::vector<std::string> words, std::string suffix) -> std::vector<std::string>;

} // namespace Reaktoro


