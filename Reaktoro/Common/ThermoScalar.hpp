// Reaktoro is a C++ library for computational reaction modelling.
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

namespace Reaktoro {

// Forward declarations
class ThermoVectorRow;
class ThermoVectorConstRow;

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

    /// Assign a row of a ThermoVector instance to this ThermoScalar instance
    auto operator=(const ThermoVectorRow& row) -> ThermoScalar&;

    /// Assign a row of a ThermoVector instance to this ThermoScalar instance
    auto operator=(const ThermoVectorConstRow& row) -> ThermoScalar&;

    /// Assign-addition of a ThermoScalar instance
    auto operator+=(const ThermoScalar& other) -> ThermoScalar&;

    /// Assign-subtraction of a ThermoScalar instance
    auto operator-=(const ThermoScalar& other) -> ThermoScalar&;

    /// Assign-multiplication of a ThermoScalar instance
    auto operator*=(double scalar) -> ThermoScalar&;

    /// Assign-division of a ThermoScalar instance
    auto operator/=(double scalar) -> ThermoScalar&;

    /// The value of the thermodynamic property
    double val = 0.0;

    /// The partial temperature derivative of the thermodynamic property
    double ddt = 0.0;

    /// The partial pressure derivative of the thermodynamic property
    double ddp = 0.0;
};

/// A type used to define the function signature for the calculation of a thermodynamic property.
/// @see ThermoScalar, ThermoVector, ThermoVectorFunction
using ThermoScalarFunction = std::function<ThermoScalar(double, double)>;

/// Compare two ThermoScalar instances for equality
auto operator==(const ThermoScalar& l, const ThermoScalar& r) -> bool;

/// Unary addition operator for a ThermoScalar instance
auto operator+(const ThermoScalar& l) -> ThermoScalar;

/// Unary subtraction operator for a ThermoScalar instance
auto operator-(const ThermoScalar& l) -> ThermoScalar;

/// Add two ThermoScalar instances
auto operator+(const ThermoScalar& l, const ThermoScalar& r) -> ThermoScalar;

/// Subtract two ThermoScalar instances
auto operator-(const ThermoScalar& l, const ThermoScalar& r) -> ThermoScalar;

/// Left-multiply a ThermoScalar instance by a scalar
auto operator*(double scalar, const ThermoScalar& r) -> ThermoScalar;

/// Right-multiply a ThermoScalar instance by a scalar
auto operator*(const ThermoScalar& l, double scalar) -> ThermoScalar;

/// Multiply two ThermoScalar instances
auto operator*(const ThermoScalar& l, const ThermoScalar& r) -> ThermoScalar;

/// Left-divide a ThermoScalar instance by a scalar
auto operator/(double scalar, const ThermoScalar& r) -> ThermoScalar;

/// Right-divide a ThermoScalar instance by a scalar
auto operator/(const ThermoScalar& l, double scalar) -> ThermoScalar;

/// Divide a ThermoScalar instance by another
auto operator/(const ThermoScalar& l, const ThermoScalar& r) -> ThermoScalar;

/// Return the square root of a ThermoScalar instance
auto sqrt(const ThermoScalar& l) -> ThermoScalar;

/// Return the power of a ThermoScalar instance
auto pow(const ThermoScalar& l, double power) -> ThermoScalar;

/// Return the natural exponential of a ThermoScalar instance
auto exp(const ThermoScalar& l) -> ThermoScalar;

/// Return the natural log of a ThermoScalar instance
auto log(const ThermoScalar& l) -> ThermoScalar;

/// Return the log10 of a ThermoScalar instance
auto log10(const ThermoScalar& l) -> ThermoScalar;

} // namespace Reaktoro
