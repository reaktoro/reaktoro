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

#include "StringList.hpp"

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {
namespace {

/// Return a string iterator to the position with a given token.
auto findNextToken(std::string::iterator begin, std::string::iterator end, char token) -> std::string::iterator
{
    if(*begin == token || begin == end) return begin;
    return findNextToken(begin + 1, end, token);
}

/// Return a string iterator to the next matching bracket.
/// This function assumes that `begin` points to an opening bracket `(`.
auto findMatchingBracket(std::string::iterator begin, std::string::iterator end) -> std::string::iterator
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

/// A recursive implementation that breaks a string into several strings at space.
/// This method treats brackets in a special case by ignoring spaces inside brackets.
/// Thus, the string `"H2O elementAmount(Ca units=mol)"` is broken into strings
/// `"H2O"` and `"elementAmount(Ca units=mol)"`.
auto convertStringToStringsHelper(
    std::string::iterator begin,
    std::string::iterator end,
    std::string::iterator current,
    char token,
    std::vector<std::string>& res) -> void
{
    if(current == end)
    {
        res.push_back(std::string(begin, end));
    }
    else if(*current == token)
    {
        res.push_back(std::string(begin, current));
        convertStringToStringsHelper(current + 1, end, current + 1, token, res);
    }
    else if(*current != '(')
    {
        convertStringToStringsHelper(begin, end, current + 1, token, res);
    }
    else if(*current == '(')
    {
        convertStringToStringsHelper(begin, end, findMatchingBracket(current, end), token, res);
    }
}

/// Convert a formatted string into a vector of strings
auto convertStringToStrings(std::string str, char token) -> std::vector<std::string>
{
    std::vector<std::string> strings;
    if(str.size())
        convertStringToStringsHelper(str.begin(), str.end(), str.begin(), token, strings);
    return strings;
}

} // namespace

StringList::StringList()
{}

StringList::StringList(std::initializer_list<std::string> strings)
: m_strings(strings.begin(), strings.end())
{}

StringList::StringList(std::vector<std::string> strings)
: m_strings(std::move(strings))
{}

StringList::StringList(const char* strings)
: m_strings(convertStringToStrings(trim(strings), ' '))
{}

StringList::StringList(const char* strings, char token)
: m_strings(convertStringToStrings(trim(strings), token))
{}

StringList::StringList(std::string strings)
: m_strings(convertStringToStrings(trim(strings), ' '))
{}

StringList::StringList(std::string strings, char token)
: m_strings(convertStringToStrings(trim(strings), token))
{}

auto StringList::empty() const -> bool
{
    return m_strings.empty();
}

auto StringList::size() const -> std::size_t
{
    return m_strings.size();
}

auto StringList::data() const -> const std::vector<std::string>&
{
    return m_strings;
}

auto StringList::operator[](std::size_t index) -> std::string&
{
    return m_strings[index];
}

auto StringList::operator[](std::size_t index) const -> const std::string&
{
    return m_strings[index];
}

StringList::operator std::vector<std::string>() const
{
    return data();
}

} // namespace Reaktoro

