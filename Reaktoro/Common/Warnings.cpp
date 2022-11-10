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

#include "Warnings.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Algorithms.hpp>
#include <Reaktoro/Common/Types.hpp>

namespace Reaktoro {

auto getDisabledWarnings() -> Vec<int>&
{
    static thread_local Vec<int> disabled_warnings;
    return disabled_warnings;
}

auto Warnings::enable(int warningid) -> bool
{
    auto& disabled = getDisabledWarnings();
    auto pos = std::find(disabled.begin(), disabled.end(), warningid);
    if(pos < disabled.end())
    {
        disabled.erase(pos);
        return true;
    }
    return false;
}

auto Warnings::disable(int warningid) -> bool
{
    auto& disabled = getDisabledWarnings();
    if(!contains(disabled, warningid))
    {
        disabled.push_back(warningid);
        return true;
    }
    return false;
}

auto Warnings::isEnabled(int warningid) -> bool
{
    return !isDisabled(warningid);
}

auto Warnings::isDisabled(int warningid) -> bool
{
    auto const& disabled = getDisabledWarnings();
    return contains(disabled, warningid);
}

} // namespace Reaktoro
