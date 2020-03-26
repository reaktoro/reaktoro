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

// Reaktoro includes
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Core/Database.hpp>

namespace Reaktoro {

// Forward declarations
class Species;

/// A type to store the standard thermodynamic properties of a single species.
struct StandardThermoProps
{
    /// The standard molar Gibbs energy @f$G^{\circ}@f$ of the species (in unit of J/mol)
    ThermoScalar G0;

    /// The standard molar Helmholtz energy @f$A^{\circ}@f$ of the species (in unit of J/mol)
    ThermoScalar A0;

    /// The standard molar internal energy @f$U^{\circ}@f$ of the species (in unit of J/mol)
    ThermoScalar U0;

    /// The standard molar enthalpy @f$H^{\circ}@f$ of the species (in unit of J/mol)
    ThermoScalar H0;

    /// The standard molar entropy @f$S^{\circ}@f$ of the species (in unit of J/K)
    ThermoScalar S0;

    /// The standard molar volume @f$V^{\circ}@f$ of the species (in unit of m3/mol)
    ThermoScalar V0;

    /// The standard molar isobaric heat capacity @f$C_{P}^{\circ}@f$ of the species (in unit of J/(mol·K))
    ThermoScalar Cp0;

    /// The standard molar isochoric heat capacity @f$C_{V}^{\circ}@f$ of the species (in unit of J/(mol·K))
    ThermoScalar Cv0;
};

/// A base type defining the interface for calculation of standard thermodynamic properties.
class ThermoEngine
{
public:
    /// Construct a ThermoEngine object.
    explicit ThermoEngine(const Database& database);

    /// Return the database in this thermodynamic engine.
    auto database() const -> const Database&;

    /// Return the standard thermodynamic properties of a species at given temperature and pressure.
    /// @param T The temperature for the calculation (in unit of K)
    /// @param P The pressure for the calculation (in unit of Pa)
    /// @param species The species object.
    virtual auto standardThermoProps(Temperature T, Pressure P, const Species& species) const -> StandardThermoProps = 0;

    /// Return the standard thermodynamic properties of multiple species at given temperature and pressure.
    /// @param T The temperature for the calculation (in unit of K)
    /// @param P The pressure for the calculation (in unit of Pa)
    /// @param species The species object.
    virtual auto standardThermoProps(Temperature T, Pressure P, const std::vector<Species>& species) const -> std::vector<StandardThermoProps>;

private:
    /// The thermodynamic database
    Database db;
};

} // namespace Reaktoro
