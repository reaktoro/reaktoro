// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {
namespace internal {

template<typename T>
auto operator<<(std::ostream& out, const Vec<T>& values) -> std::ostream&
{
    for(auto i = 0; i < values.size(); ++i)
        out << ((i == 0) ? "" : ", ") << values[i];
    return out;
}

template <typename Arg>
auto stringfy(std::stringstream& ss, const String& sep, const Arg& item) -> void
{
    ss << item;
}

template <typename Arg, typename... Args>
auto stringfy(std::stringstream& ss, const String& sep, const Arg& item, Args... items) -> void
{
    ss << item << sep;
    stringfy(ss, sep, items...);
}

} // namespace internal

/// Concatenate the arguments into a string using a given separator string.
template <typename... Args>
auto stringfy(const String& sep, Args... items) -> String
{
    std::stringstream ss;
    internal::stringfy(ss, sep, items...);
    return ss.str();
}

/// Concatenate the arguments into a string without any separator string.
template <typename... Args>
auto str(Args... items) -> String
{
    return stringfy("", items...);
}

/// Return a new string where `substr` occurrences are replaced by `newsubstr`.
inline auto replace(String original, String substr, String newsubstr) -> String
{
    if(substr.empty()) return original;
    auto pos = original.find(substr);
    while(pos != String::npos)
    {
        original.replace(pos, substr.size(), newsubstr);
        pos = original.find(substr, pos + newsubstr.size());
    }
    return original;
}

/// Return a string with lower case characters.
inline auto lowercase(String str) -> String
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

/// Return a string with upper case characters.
inline auto uppercase(String str) -> String
{
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

/// Trim the string from start (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
inline auto trimleft(String str) -> String
{
    str.erase(str.begin(), std::find_if(str.begin(), str.end(),
        [](unsigned char ch) { return !std::isspace(ch); }));
    return str;
}

/// Trim the string from end (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
inline auto trimright(String str) -> String
{
    str.erase(std::find_if(str.rbegin(), str.rend(),
        [](unsigned char ch) { return !std::isspace(ch); }).base(), str.end());
    return str;
}

/// Trim the string from both ends (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
inline auto trim(String str) -> String
{
    return trimleft(trimright(str));
}

/// Split the string on every occurrence of the specified delimiters and apply a transform function.
inline auto split(const String& str, const String& delims, std::function<String(String)> transform) -> Strings
{
    Strings words;
    std::size_t start = 0, end = 0;
    while(end != String::npos)
    {
        end = str.find_first_of(delims, start);
        String word = str.substr(start, end - start);
        if(word != "") words.push_back(transform ? transform(word) : word);
        start = end + 1;
    }
    return words;
}

/// Split the string on every occurrence of the specified delimiters
inline auto split(const String& str, const String& delims = " ") -> Strings
{
    return split(str, delims, {});
}

/// Join several strings into one.
inline auto join(const Strings& strs, String sep = " ") -> String
{
    String res;
    for(auto i = 0; i < strs.size(); ++i)
        res += (i == 0 ? "" : sep) + strs[i];
    return res;
}

/// Convert the string into a floating point number
inline auto tofloat(const String& str) -> double
{
    return atof(str.c_str());
}

/// Return a list of words with duplicate names converted to unique ones.
inline auto makeunique(Strings words, String suffix) -> Strings
{
    const auto size = words.size();
    for(auto i = 1; i < size; ++i)
    {
        auto uniqueword = words[i];
        for(auto j = 0; j < i; ++j)
            if(uniqueword == words[j])
                uniqueword += suffix;
        words[i] = uniqueword;
    }
    return words;
}

} // namespace Reaktoro


