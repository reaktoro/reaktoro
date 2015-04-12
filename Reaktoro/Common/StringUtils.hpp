// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// C++ includes
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <locale>
#include <string>
#include <vector>

namespace Reaktoro {

/// Trims the string from start (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
inline auto leftTrim(std::string& str) -> std::string&
{
    str.erase(str.begin(),
        std::find_if(str.begin(), str.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return str;
}

/// Trims the string from end (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
inline auto rightTrim(std::string& str) -> std::string&
{
    str.erase(std::find_if(str.rbegin(), str.rend(),
        std::not1(std::ptr_fun<int, int>(std::isspace))).base(), str.end());
    return str;
}

/// Trims the string from both ends (taken from http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring)
inline auto trim(std::string& str) -> std::string&
{
    return leftTrim(rightTrim(str));
}

/// Splits the string on every occurrence of the specified delimiters
inline auto split(const std::string& str, const std::string& delims = " ") -> std::vector<std::string>
{
    std::vector<std::string> words;
    std::size_t start = 0, end = 0;
    while(end != std::string::npos)
    {
        end = str.find_first_of(delims, start);
        std::string word = str.substr(start, end - start);
        if(word != "") words.push_back(word);
        start = end + 1;
    }
    return words;
}

/// Converts the string into a floating point number
inline auto tofloat(const std::string& str) -> double
{
    return atof(str.c_str());
}

/// Converts the string into a list of floating point numbers
inline auto tofloats(const std::string& str, const std::string& delims = " ") -> std::vector<double>
{
    std::vector<double> values;
    for(const std::string& num : split(str, delims))
        values.push_back(tofloat(num));
    return values;
}

} // namespace Reaktoro


