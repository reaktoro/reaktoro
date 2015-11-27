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

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ThermoScalar;
class ChemicalVectorRow;
class ChemicalVectorRowConst;

/// A type that defines a scalar chemical property.
/// A chemical property means here any property that depends on
/// temperature, pressure and composition. A ChemicalScalar instance
/// not only holds the value of the chemical property, but also is partial
/// temperature, pressure and molar derivatives.
/// @see ChemicalVector
class ChemicalScalar
{
public:
    /// Construct a default ChemicalScalar instance.
    ChemicalScalar();

    /// Construct a ChemicalScalar instance with given number of species.
    explicit ChemicalScalar(unsigned nspecies);

    /// Construct a ChemicalScalar instance.
    /// @param val The scalar value of the thermodynamic quantity
    /// @param ddt The partial temperature derivative of the thermodynamic quantity
    /// @param ddp The partial pressure derivative of the thermodynamic quantity
    /// @param ddn The partial molar derivatives of the thermodynamic quantity
    ChemicalScalar(double val, double ddt, double ddp, const Vector& ddn);

    /// Construct a ChemicalScalar instance from a ChemicalVectorRow instance.
    /// @param row The row of a ChemicalVector instance
    ChemicalScalar(const ChemicalVectorRow& row);

    /// Construct a ChemicalScalar instance from a ChemicalVectorConstRow instance.
    /// @param row The row of a const ChemicalVector instance
    ChemicalScalar(const ChemicalVectorRowConst& row);

    /// Assign a ThermoScalar instance to this ChemicalScalar instance.
    auto operator=(const ThermoScalar& scalar) -> ChemicalScalar&;

    /// Assign a row of a ChemicalVector instance to this ChemicalScalar instance.
    auto operator=(const ChemicalVectorRow& row) -> ChemicalScalar&;

    /// Assign a row of a ChemicalVector instance to this ChemicalScalar instance.
    auto operator=(const ChemicalVectorRowConst& row) -> ChemicalScalar&;

    /// Assign-addition of a ChemicalScalar instance.
    auto operator+=(const ChemicalScalar& other) -> ChemicalScalar&;

    /// Assign-addition of a ThermoScalar instance.
    auto operator+=(const ThermoScalar& other) -> ChemicalScalar&;

    /// Assign-addition of a scalar.
    auto operator+=(double scalar) -> ChemicalScalar&;

    /// Assign-subtraction of a ChemicalScalar instance.
    auto operator-=(const ChemicalScalar& other) -> ChemicalScalar&;

    /// Assign-subtraction of a ThermoScalar instance.
    auto operator-=(const ThermoScalar& other) -> ChemicalScalar&;

    /// Assign-subtraction of a scalar.
    auto operator-=(double scalar) -> ChemicalScalar&;

    /// Assign-multiplication of a ChemicalScalar instance.
    auto operator*=(double scalar) -> ChemicalScalar&;

    /// Assign-division of a ChemicalScalar instance.
    auto operator/=(double scalar) -> ChemicalScalar&;

    /// The scalar value of the chemical property.
    double val;

    /// The partial temperature derivative of the chemical property.
    double ddt;

    /// The partial pressure derivative of the chemical property.
    double ddp;

    /// The partial molar derivatives of the chemical property.
    Vector ddn;
};

/// A type used to define the function signature for the calculation of a chemical property.
/// @see ChemicalScalar, ChemicalVectorFunction
using ChemicalScalarFunction = std::function<ChemicalScalar(double, double, const Vector&)>;

/// Compare two ChemicalScalar instances for equality
auto operator==(const ChemicalScalar& l, const ChemicalScalar& r) -> bool;

/// Unary addition operator for a ChemicalScalar instance
auto operator+(const ChemicalScalar& l) -> ChemicalScalar;

/// Unary subtraction operator for a ChemicalScalar instance
auto operator-(const ChemicalScalar& l) -> ChemicalScalar;

/// Add two ChemicalScalar instances
auto operator+(const ChemicalScalar& l, const ChemicalScalar& r) -> ChemicalScalar;

/// Add a ChemicalScalar instance and a ThermoScalar instance
auto operator+(const ChemicalScalar& l, const ThermoScalar& r) -> ChemicalScalar;

/// Add a ThermoScalar instance and a ChemicalScalar instance
auto operator+(const ThermoScalar& l, const ChemicalScalar& r) -> ChemicalScalar;

/// Add a ChemicalScalar instance and a scalar
auto operator+(const ChemicalScalar& l, double scalar) -> ChemicalScalar;

/// Add a scalar and a ChemicalScalar instance
auto operator+(double scalar, const ChemicalScalar& r) -> ChemicalScalar;

/// Subtract two ChemicalScalar instances
auto operator-(const ChemicalScalar& l, const ChemicalScalar& r) -> ChemicalScalar;

/// Subtract a ChemicalScalar instance and a ThermoScalar instance
auto operator-(const ChemicalScalar& l, const ThermoScalar& r) -> ChemicalScalar;

/// Subtract a ThermoScalar instance and a ChemicalScalar instance
auto operator-(const ThermoScalar& l, const ChemicalScalar& r) -> ChemicalScalar;

/// Subtract a ChemicalScalar instance and a scalar
auto operator-(const ChemicalScalar& l, double scalar) -> ChemicalScalar;

/// Subtract a scalar and a ChemicalScalar instance
auto operator-(double scalar, const ChemicalScalar& r) -> ChemicalScalar;

/// Left-multiply a ChemicalScalar instance by a scalar
auto operator*(double scalar, const ChemicalScalar& r) -> ChemicalScalar;

/// Right-multiply a ChemicalScalar instance by a scalar
auto operator*(const ChemicalScalar& l, double scalar) -> ChemicalScalar;

/// Left-multiply a ChemicalScalar instance by a ThermoScalar
auto operator*(const ThermoScalar& l, const ChemicalScalar& r) -> ChemicalScalar;

/// Right-multiply a ChemicalScalar instance by a ThermoScalar
auto operator*(const ChemicalScalar& l, const ThermoScalar& r) -> ChemicalScalar;

/// Multiply two ChemicalScalar instances
auto operator*(const ChemicalScalar& l, const ChemicalScalar& r) -> ChemicalScalar;

/// Left-divide a ChemicalScalar instance by a scalar
auto operator/(double scalar, const ChemicalScalar& r) -> ChemicalScalar;

/// Right-divide a ChemicalScalar instance by a scalar
auto operator/(const ChemicalScalar& l, double scalar) -> ChemicalScalar;

/// Left-divide a ChemicalScalar instance by a ThermoScalar
auto operator/(const ThermoScalar& l, const ChemicalScalar& r) -> ChemicalScalar;

/// Right-divide a ChemicalScalar instance by a ThermoScalar
auto operator/(const ChemicalScalar& l, const ThermoScalar& r) -> ChemicalScalar;

/// Divide a ChemicalScalar instance by another
auto operator/(const ChemicalScalar& l, const ChemicalScalar& r) -> ChemicalScalar;

/// Return the square root of a ChemicalScalar instance
auto sqrt(const ChemicalScalar& l) -> ChemicalScalar;

/// Return the power of a ChemicalScalar instance
auto pow(const ChemicalScalar& l, double power) -> ChemicalScalar;

/// Return the natural exponential of a ChemicalScalar instance
auto exp(const ChemicalScalar& l) -> ChemicalScalar;

/// Return the natural log of a ChemicalScalar instance
auto log(const ChemicalScalar& l) -> ChemicalScalar;

/// Return the log10 of a ChemicalScalar instance
auto log10(const ChemicalScalar& l) -> ChemicalScalar;

} // namespace Reaktoro
