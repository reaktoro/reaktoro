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

#pragma once

// C++ includes
#include <memory>

// Reaktoro includes
#include <Reaktoro/Core/ThermoEngine.hpp>

namespace Reaktoro {

// Forward declarations
class DatabaseSupcrt;

/// A type defining the calculation of standard thermodynamic properties based on SUPCRT.
class ThermoEngineSupcrt : public ThermoEngine
{
public:
    /// Construct a ThermoEngineSupcrt object.
    explicit ThermoEngineSupcrt(const DatabaseSupcrt& database);

    /// Construct a copy of a ThermoEngineSupcrt object.
    ThermoEngineSupcrt(const ThermoEngineSupcrt& other);

    /// Destroy this ThermoEngineSupcrt object.
    ~ThermoEngineSupcrt();

    /// Assign a ThermoEngineSupcrt object to this.
    auto operator=(ThermoEngineSupcrt other) -> ThermoEngineSupcrt&;

    /// Return the standard thermodynamic properties of a species at given temperature and pressure.
    /// @param T The temperature for the calculation (in unit of K)
    /// @param P The pressure for the calculation (in unit of Pa)
    /// @param species The species object.
    auto standardThermoProps(Temperature T, Pressure P, const Species& species) const -> StandardThermoProps;

private:
    struct Impl;

    std::unique_ptr<Impl> pimpl;
};

} // namespace Reaktoro
