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

#include "StringUtils.hpp"

// C++ includes
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iomanip>

namespace Reaktoro {

static int __precision = 4;

auto precision(int precision) -> void
{
    __precision = precision;
}

auto precision() -> int
{
    return __precision;
}

auto strfix(double num, int precision) -> std::string
{
    precision = precision < 0 ? __precision : precision;
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << num;
    return ss.str();
}

auto strsci(double num, int precision) -> std::string
{
    precision = precision < 0 ? __precision : precision;
    std::ostringstream ss;
    ss << std::scientific << std::setprecision(precision) << num;
    return ss.str();
}

auto replace(std::string original, std::string substr, std::string newsubstr) -> std::string
{
    if(substr.empty()) return original;
    auto pos = original.find(substr);
    while(pos != std::string::npos)
    {
        original.replace(pos, substr.size(), newsubstr);
        pos = original.find(substr, pos + newsubstr.size());
    }
    return original;
}

auto lowercase(std::string str) -> std::string
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

auto uppercase(std::string str) -> std::string
{
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

auto trimleft(std::string str) -> std::string
{
    str.erase(str.begin(), std::find_if(str.begin(), str.end(),
        [](unsigned char ch) { return !std::isspace(ch); }));
    return str;
}

auto trimright(std::string str) -> std::string
{
    str.erase(std::find_if(str.rbegin(), str.rend(),
        [](unsigned char ch) { return !std::isspace(ch); }).base(), str.end());
    return str;
}

auto trim(std::string str) -> std::string
{
    return trimleft(trimright(str));
}

auto split(const std::string& str, const std::string& delims, std::function<std::string(std::string)> transform) -> std::vector<std::string>
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

auto split(const std::string& str, const std::string& delims) -> std::vector<std::string>
{
    return split(str, delims, {});
}

auto join(const std::vector<std::string>& strs, std::string sep) -> std::string
{
    std::string res;
    for(auto i = 0; i < strs.size(); ++i)
        res += (i == 0 ? "" : sep) + strs[i];
    return res;
}

auto tofloat(const std::string& str) -> double
{
    return atof(str.c_str());
}

auto makeunique(std::vector<std::string> words, std::string suffix) -> std::vector<std::string>
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


