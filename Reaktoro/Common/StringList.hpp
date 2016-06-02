// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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
#include <string>

namespace Reaktoro {

/// A class for representing a list of strings with special constructors.
class StringList
{
public:
    /// Construct a default StringList instance.
    StringList();

    /// Construct a StringList instance by breaking the words within separated by a token.
    /// This method converts a given string into a list of strings separated by a token.
    /// However, it ignores the token inside brackets. If the token is space, then the
    /// string `"hello (to you)"` are split into `"hello"` and `"(to you)"`.
    /// @param str The string containing words separated by `token`.
    /// @param token The token used to separate the words in `str`.
    StringList(std::string str, char token);

    /// Construct a StringList instance with a container of strings.
    template<typename StringContainer>
    StringList(const StringContainer& strings);

    /// Destroy this StringList instance.
    virtual ~StringList();

    /// Return a vector of strings.
    auto strings() const -> const std::vector<std::string>&;

private:
    std::vector<std::string> _strings;
};

template<typename StringContainer>
StringList::StringList(const StringContainer& strings)
: _strings(strings.begin(), strings.end())
{}

} // namespace Reaktoro
