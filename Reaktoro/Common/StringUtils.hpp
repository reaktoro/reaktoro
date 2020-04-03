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
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <functional>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

namespace Reaktoro {
namespace internal {

template <typename Arg>
auto stringfy(std::stringstream& ss, const std::string& sep, const Arg& item)
{
    ss << item;
}

template <typename Arg, typename... Args>
auto stringfy(std::stringstream& ss, const std::string& sep, const Arg& item, Args... items) -> void
{
    ss << item << sep;
    stringfy(ss, sep, items...);
}

} // namespace internal

/// Concatenate the arguments into a string using a given separator string.
template <typename... Args>
auto stringfy(const std::string& sep, Args... items) -> std::string
{
    std::stringstream ss;
    internal::stringfy(ss, sep, items...);
    return ss.str();
}

/// Concatenate the arguments into a string without any separator string.
template <typename... Args>
auto str(Args... items) -> std::string
{
    return stringfy("", items...);
}

/// Return a new string where `substr` occurrences are replaced by `newsubstr`.
inline auto replace(std::string original, std::string substr, std::string newsubstr) -> std::string
{
	auto pos = original.find(substr);
	while(pos != std::string::npos)
	{
		original.replace(pos, substr.size(), newsubstr);
		pos = original.find(substr, pos + newsubstr.size());
	}
    return original;
}

/// Return a string with lower case characters.
inline auto lowercase(std::string str) -> std::string
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

/// Return a string with upper case characters.
inline auto uppercase(std::string str) -> std::string
{
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

/// Trim the string from start (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
inline auto leftTrim(std::string& str) -> std::string&
{
    str.erase(str.begin(),
        std::find_if(str.begin(), str.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return str;
}

/// Trim the string from end (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
inline auto rightTrim(std::string& str) -> std::string&
{
    str.erase(std::find_if(str.rbegin(), str.rend(),
        std::not1(std::ptr_fun<int, int>(std::isspace))).base(), str.end());
    return str;
}

/// Trim the string from both ends (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
inline auto trim(std::string& str) -> std::string&
{
    return leftTrim(rightTrim(str));
}

/// Split the string on every occurrence of the specified delimiters
inline auto split(const std::string& str, const std::string& delims,
	std::function<std::string&(std::string&)> transform) -> std::vector<std::string>
{
	std::vector<std::string> words;
	std::size_t start = 0, end = 0;
	while(end != std::string::npos)
	{
		end = str.find_first_of(delims, start);
		std::string word = str.substr(start, end - start);
		if(word != "") words.push_back(transform ? transform(word) : word);
		start = end + 1;
	}
	return words;
}

/// Split the string on every occurrence of the specified delimiters
inline auto split(const std::string& str, const std::string& delims = " ") -> std::vector<std::string>
{
    return split(str, delims, {});
}

/// Split the string on every occurrence of the specified delimiters and trim each word
inline auto splitrim(const std::string& str, const std::string& delims = " ") -> std::vector<std::string>
{
    return split(str, delims, trim);
}

/// Join several strings into one.
inline auto join(const std::vector<std::string>& strs, std::string delim = " ") -> std::string
{
    std::string res;
    for(unsigned i = 0; i < strs.size(); ++i)
        res = res + (i > 0 ? delim : "") + strs[i];
    return res;
}

/// Convert the string into a floating point number
inline auto tofloat(const std::string& str) -> double
{
    return atof(str.c_str());
}

/// Convert the string into a list of floating point numbers
inline auto tofloats(const std::string& str, const std::string& delims = " ") -> std::vector<double>
{
    std::vector<double> values;
    for(const std::string& num : split(str, delims))
        values.push_back(tofloat(num));
    return values;
}

} // namespace Reaktoro


