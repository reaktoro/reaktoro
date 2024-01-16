// Reaktoro is a unified framework for modeling chemically reactive phases.
//
// Copyright © 2014-2024 Allan Leal
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

#include "Embedded.hpp"

// CMakeRC includes
#include <cmrc/cmrc.hpp>

CMRC_DECLARE(ReaktoroEmbedded);

namespace Reaktoro {

auto Embedded::get(String const& path) -> String
{
    return getAsString(path);
}

auto Embedded::getAsString(String const& path) -> String
{
    const auto [begin, end] = getAsStringView(path);
    return String(begin, end);
}

auto Embedded::getAsStringView(String const& path) -> Pair<Chars, Chars>
{
    auto fs = cmrc::ReaktoroEmbedded::get_filesystem();
    auto file = fs.open("embedded/" + path);
    return { file.begin(), file.end() };
}

} // namespace Reaktoro
