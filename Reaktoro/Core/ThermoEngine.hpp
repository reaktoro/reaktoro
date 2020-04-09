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

// Reaktoro includes
#include <Reaktoro/Common/Real.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/StandardThermoModel.hpp>

namespace Reaktoro {

/// The thermodynamic engine base class for standard thermodynamic properties calculation.
class ThermoEngine
{
public:
    /// Construct a ThermoEngine object.
    /// @param db The thermodynamic database
    /// @param model The standard thermodynamic model function
    ThermoEngine(const Database& db, const StandardThermoPropsFn& model);

    /// Return the database of this thermodynamic engine.
    auto database() const -> const Database&;

    /// Return the standard thermodynamic model function of this thermodynamic engine.
    auto standardThermoModelFn() const -> const StandardThermoPropsFn&;

    /// Return the standard thermodynamic properties of a species at given temperature and pressure.
    /// @param T The temperature for the calculation (in unit of K)
    /// @param P The pressure for the calculation (in unit of Pa)
    /// @param species The species object.
    auto standardThermoProps(real T, real P, const Species& species) const -> StandardThermoProps;

private:
    /// The thermodynamic database of this thermodynamic engine
    Database db;

    /// The standard thermodynamic model function of this thermodynamic engine
    StandardThermoPropsFn model;
};

} // namespace Reaktoro
