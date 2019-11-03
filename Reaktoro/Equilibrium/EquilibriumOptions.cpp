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
