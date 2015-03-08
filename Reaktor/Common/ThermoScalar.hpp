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

// C++ includes
#include <functional>

namespace Reaktor {

// Forward declarations
struct ThermoVectorRow;
struct ThermoVectorConstRow;

/// Describe a thermodynamic property value and its partial derivatives w.r.t. temperature and pressure
class ThermoScalar
{
public:
    /// Construct a default ThermoScalar instance
    ThermoScalar();

    /// Construct a ThermoScalar instance
    /// @param val The value of the thermodynamic property
    /// @param ddt The partial temperature derivative of the thermodynamic property
    /// @param ddp The partial pressure derivative of the thermodynamic property
    ThermoScalar(double val, double ddt, double ddp);

    /// Get the value of the thermodynamic property
    auto val() const -> double;

    /// Get the partial temperature derivative of the thermodynamic property
    auto ddt() const -> double;

    /// Get the partial pressure derivative of the thermodynamic property
    auto ddp() const -> double;

    /// Assign a row of a ThermoVector instance to this ThermoScalar instance
    auto operator=(const ThermoVectorRow& row) -> ThermoScalar&;

    /// Assign a row of a ThermoVector instance to this ThermoScalar instance
    auto operator=(const ThermoVectorConstRow& row) -> ThermoScalar&;

private:
    /// The value of the thermodynamic property
    double m_val = 0.0;

    /// The partial temperature derivative of the thermodynamic property
    double m_ddt = 0.0;

    /// The partial pressure derivative of the thermodynamic property
    double m_ddp = 0.0;
};

/// Compares two ThermoScalar instances for equality
auto operator==(const ThermoScalar& l, const ThermoScalar& r) -> bool;

/// A type used to define the function signature for the calculation of a thermodynamic property.
/// @see ThermoScalar, ThermoVector, ThermoVectorFunction
using ThermoScalarFunction = std::function<ThermoScalar(double, double)>;

} // namespace Reaktor
