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

#include "ThermoEngine.hpp"

// Reaktoro includes
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro {

ThermoEngine::ThermoEngine(const Database& database)
: db(database)
{
}

auto ThermoEngine::database() const -> const Database&
{
    return db;
}

auto ThermoEngine::standardThermoProps(Temperature T, Pressure P, const std::vector<Species>& species) const -> std::vector<StandardThermoProps>
{
    std::vector<StandardThermoProps> res;
    res.reserve(species.size());
    for(const auto& s : species)
        res.push_back(standardThermoProps(T, P, s));
    return res;
}

} // namespace Reaktoro
