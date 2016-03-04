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

#include "EquilibriumOptions.hpp"

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro{

EquilibriumOptions::EquilibriumOptions()
{}

EquilibriumOptions::EquilibriumOptions(const char* str)
{
    auto words = split(str);
    for(auto word : words)
    {
        auto pair = split(word, "=");
        if(pair.front() == "output")
            optimum.output.active = pair.back() == "true" ? true : false;
    }
}

} // namespace Reaktoro



