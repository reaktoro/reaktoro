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

#include "ThermoEngine.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Core/Species.hpp>

namespace Reaktoro {

ThermoEngine::ThermoEngine(const Database& db, const StandardThermoPropsFn& model)
: db(db), model(model)
{
    // Assert given standard thermodynamic model function is not empty
    Assert(model, "Failure at ThermoEngine::ThermoEngine(const Database&, const StandardThermoPropsFn&).",
        "Given StandardThermoPropsFn object is empty.");
}

auto ThermoEngine::database() const -> const Database&
{
    return db;
}

auto ThermoEngine::standardThermoPropsFn() const -> const StandardThermoPropsFn&
{
    return model;
}

auto ThermoEngine::props(real T, real P, const Species& species) const -> StandardThermoProps
{
    return model(T, P, species);
}

} // namespace Reaktoro
