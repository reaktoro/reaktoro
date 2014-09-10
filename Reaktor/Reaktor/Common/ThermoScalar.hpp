// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

#pragma once

// Reaktor includes
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class ThermoVectorRow;

/// A type that defines a scalar thermodynamic quantity with its partial temperature, pressure and molar derivatives
class ThermoScalar
{
public:
    /// Construct a default ThermoScalar instance
    ThermoScalar();

    /// Construct a ThermoScalar instance from a row of a ThermoVector instance
    /// @param row The ThermoVectorRow instance from which the ThermoScalar instance is built
    ThermoScalar(const ThermoVectorRow& row);

    /// Set the scalar value of the thermodynamic quantity
    auto val(double val) -> ThermoScalar&;

    /// Set the partial derivative of the thermodynamic quantity with respect to temperature (temperature in units of K)
    auto ddt(double ddt) -> ThermoScalar&;

    /// Set the partial derivative of the thermodynamic quantity with respect to pressure (pressure in units of Pa)
    auto ddp(double ddp) -> ThermoScalar&;

    /// Set the partial derivative of the thermodynamic quantity with respect to composition (composition in units of mol)
    auto ddn(const Vector& ddn) -> ThermoScalar&;

    /// Get the scalar value of the thermodynamic quantity
    auto val() const -> const double&;

    /// Get the partial derivative of the thermodynamic quantity with respect to temperature (temperature in units of K)
    auto ddt() const -> const double&;

    /// Get the partial derivative of the thermodynamic quantity with respect to pressure (pressure in units of Pa)
    auto ddp() const -> const double&;

    /// Get the partial derivative of the thermodynamic quantity with respect to composition (composition in units of mol)
    auto ddn() const -> const Vector&;

private:
    /// The scalar value of the thermodynamic quantity
    double m_val;

    /// The partial derivative of the thermodynamic quantity with respect to temperature (temperature in units of K)
    double m_ddt;

    /// The partial derivative of the thermodynamic quantity with respect to pressure (pressure in units of Pa)
    double m_ddp;

    /// The partial derivative of the thermodynamic quantity with respect to composition (composition in units of mol)
    Vector m_ddn;
};

/// Get the scalar value of the thermodynamic quantity
inline auto val(const ThermoScalar& scalar) -> const double&
{
    return scalar.val();
}

/// Get the partial derivative of the thermodynamic quantity with respect to temperature (temperature in units of K)
inline auto ddt(const ThermoScalar& scalar) -> const double&
{
    return scalar.ddt();
}

/// Get the partial derivative of the thermodynamic quantity with respect to pressure (pressure in units of Pa)
inline auto ddp(const ThermoScalar& scalar) -> const double&
{
    return scalar.ddp();
}

/// Get the partial derivative of the thermodynamic quantity with respect to composition (composition in units of mol)
inline auto ddn(const ThermoScalar& scalar) -> const Vector&
{
    return scalar.ddn();
}

} // namespace Reaktor
