// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
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
#include <string>
#include <vector>

namespace Reaktoro {

/// A class for representing a list of strings with special constructors.
class StringList
{
public:
    /// Construct a default StringList instance.
    StringList();

    /// Construct a StringList instance with a initializer list of strings.
    StringList(std::initializer_list<std::string> strings);

    /// Construct a StringList instance with a vector of strings.
    StringList(std::vector<std::string> strings);

    /// Construct a StringList instance by breaking the words within separated by space.
    /// This method converts a given string into a list of strings separated by a token.
    /// However, it ignores the token inside brackets. If the token is space, then the
    /// string `"hello (to you)"` are split into `"hello"` and `"(to you)"`.
    /// @param str The string containing words separated by space.
    StringList(const char* str);

    /// Construct a StringList instance by breaking the words within separated by a token.
    /// This method converts a given string into a list of strings separated by a token.
    /// However, it ignores the token inside brackets. If the token is space, then the
    /// string `"hello (to you)"` are split into `"hello"` and `"(to you)"`.
    /// @param str The string containing words separated by `token`.
    /// @param token The token used to separate the words in `str`.
    StringList(const char* str, char token);

    /// Construct a StringList instance by breaking the words within separated by space.
    /// This method converts a given string into a list of strings separated by a token.
    /// However, it ignores the token inside brackets. If the token is space, then the
    /// string `"hello (to you)"` are split into `"hello"` and `"(to you)"`.
    /// @param str The string containing words separated by space.
    StringList(std::string str);

    /// Construct a StringList instance by breaking the words within separated by a token.
    /// This method converts a given string into a list of strings separated by a token.
    /// However, it ignores the token inside brackets. If the token is space, then the
    /// string `"hello (to you)"` are split into `"hello"` and `"(to you)"`.
    /// @param str The string containing words separated by `token`.
    /// @param token The token used to separate the words in `str`.
    StringList(std::string str, char token);

    /// Return the number of strings.
    auto size() const -> std::size_t;

    /// Return the vector of strings.
    auto data() const -> const std::vector<std::string>&;

    /// Return the string at the given index position.
    auto operator[](std::size_t index) -> std::string&;

    /// Return the string at the given index position.
    auto operator[](std::size_t index) const -> const std::string&;

    /// Convert this StringList instance into a std::vector<std::string>.
    operator std::vector<std::string>() const;

private:
    /// The list of strings separated by given token
    std::vector<std::string> m_strings;

public:
    /// Return begin const iterator of this StringList instance (for STL compatibility reasons).
    inline auto begin() const { return data().begin(); }

    /// Return begin iterator of this StringList instance (for STL compatibility reasons).
    inline auto begin() { return data().begin(); }

    /// Return end const iterator of this StringList instance (for STL compatibility reasons).
    inline auto end() const { return data().end(); }

    /// Return end iterator of this StringList instance (for STL compatibility reasons).
    inline auto end() { return data().end(); }

    /// The type of the value stored in a StringList instance (for STL compatibility reasons).
    using value_type = std::string;
};

} // namespace Reaktoro
