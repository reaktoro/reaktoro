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
#include <Reaktor/Common/Index.hpp>
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

// Forward declarations
class ChemicalScalar;
class ChemicalVectorConstRow;
class ChemicalVectorRow;
class ThermoScalar;
class ThermoVector;

/// A type that defines a vector chemical property.
/// A chemical property means here any property that depends on
/// temperature, pressure and composition. A ChemicalVector instance
/// not only holds the values of chemical properties, but also their partial
/// temperature, pressure and molar derivatives.
/// @see ChemicalScalar
class ChemicalVector
{
public:
    /// Construct a default ChemicalVector instance
    ChemicalVector();

    /// Construct a ChemicalVector instance with given dimensions
    /// @param nrows The number of rows of the vector quantities
    /// @param nrows The number of columns of the matrix quantities
    ChemicalVector(unsigned nrows, unsigned ncols);

    /// Construct a ChemicalVector instance with given data members
    /// @param val The vector value of the chemical property
    /// @param ddt The partial temperature derivatives of the vector chemical property
    /// @param ddp The partial pressure derivative of the vector chemical property
    /// @param ddn The partial molar derivatives of the vector chemical property
    ChemicalVector(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn);

    /// Get a reference of a row of this ChemicalVector instance
    auto row(unsigned irow) -> ChemicalVectorRow;

    /// Get a const reference of a row of this ChemicalVector instance
    auto row(unsigned irow) const -> ChemicalVectorConstRow;

    /// Assign-addition of a ChemicalVector instance.
    auto operator+=(const ChemicalVector& other) -> ChemicalVector&;

    /// Assign-addition of a ThermoVector instance.
    auto operator+=(const ThermoVector& other) -> ChemicalVector&;

    /// Assign-subtraction of a ChemicalVector instance.
    auto operator-=(const ChemicalVector& other) -> ChemicalVector&;

    /// Assign-subtraction of a ThermoVector instance.
    auto operator-=(const ThermoVector& other) -> ChemicalVector&;

    /// Assign-multiplication of a ChemicalVector instance.
    auto operator*=(double scalar) -> ChemicalVector&;

    /// Assign-division of a ChemicalVector instance.
    auto operator/=(double scalar) -> ChemicalVector&;

    /// The vector value of the chemical property
    Vector val;

    /// The partial temperature derivatives of the vector chemical property
    Vector ddt;

    /// The partial pressure derivative of the vector chemical property
    Vector ddp;

    /// The partial molar derivatives of the vector chemical property
    Matrix ddn;
};

/// An auxiliary type for the representation of the view of a row of a ChemicalVector instance
class ChemicalVectorRow
{
public:
    ChemicalVectorRow(ChemicalVector& vector, unsigned irow);
    auto operator=(const ChemicalVectorRow& row) -> ChemicalVectorRow&;
    auto operator=(const ChemicalScalar& scalar) -> ChemicalVectorRow&;
    double& val;
    double& ddt;
    double& ddp;
    decltype(std::declval<Matrix>().row(0)) ddn;
};

/// An auxiliary type for the representation of the const view of a row of a ChemicalVector instance
class ChemicalVectorConstRow
{
public:
    ChemicalVectorConstRow(const ChemicalVector& vector, unsigned irow);
    const double& val;
    const double& ddt;
    const double& ddp;
    decltype(std::declval<const Matrix>().row(0)) ddn;
};

/// A type used to define the function signature for the calculation of a vector of chemical properties.
/// @see ChemicalVector, ChemicalScalarFunction
using ChemicalVectorFunction = std::function<ChemicalVector(double, double, const Vector&)>;

/// Compare two ChemicalVector instances for equality
auto operator==(const ChemicalVector& l, const ChemicalVector& r) -> bool;

/// Unary addition operator for a ChemicalVector instance
auto operator+(const ChemicalVector& l) -> ChemicalVector;

/// Unary subtraction operator for a ChemicalVector instance
auto operator-(const ChemicalVector& l) -> ChemicalVector;

/// Add two ChemicalVector instances
auto operator+(const ChemicalVector& l, const ChemicalVector& r) -> ChemicalVector;

/// Add a ChemicalVector instance and a ThermoVector instance
auto operator+(const ChemicalVector& l, const ThermoVector& r) -> ChemicalVector;

/// Add a ThermoVector instance and a ChemicalVector instance
auto operator+(const ThermoVector& l, const ChemicalVector& r) -> ChemicalVector;

/// Subtract two ChemicalVector instances
auto operator-(const ChemicalVector& l, const ChemicalVector& r) -> ChemicalVector;

/// Subtract a ChemicalVector instance and a ThermoVector instance
auto operator-(const ChemicalVector& l, const ThermoVector& r) -> ChemicalVector;

/// Subtract a ThermoVector instance and a ChemicalVector instance
auto operator-(const ThermoVector& l, const ChemicalVector& r) -> ChemicalVector;

/// Left-multiply a ChemicalVector instance by a scalar
auto operator*(double scalar, const ChemicalVector& r) -> ChemicalVector;

/// Right-multiply a ChemicalVector instance by a scalar
auto operator*(const ChemicalVector& l, double scalar) -> ChemicalVector;

/// Left-multiply a ChemicalVector instance by a ThermoScalar instance
auto operator*(const ThermoScalar& scalar, const ChemicalVector& r) -> ChemicalVector;

/// Right-multiply a ChemicalVector instance by a ThermoScalar instance
auto operator*(const ChemicalVector& l, const ThermoScalar& scalar) -> ChemicalVector;

/// Multiply two ChemicalVector instances
auto operator*(const ChemicalVector& l, const ChemicalVector& r) -> ChemicalVector;

/// Left-divide a ChemicalVector instance by a scalar
auto operator/(double scalar, const ChemicalVector& r) -> ChemicalVector;

/// Right-divide a ChemicalVector instance by a scalar
auto operator/(const ChemicalVector& l, double scalar) -> ChemicalVector;

/// Divide a ChemicalVector instance by another
auto operator/(const ChemicalVector& l, const ChemicalVector& r) -> ChemicalVector;

/// Return the natural exponential of a ChemicalVector instance
auto exp(const ChemicalVector& l) -> ChemicalVector;

/// Return the natural log of a ChemicalVector instance
auto log(const ChemicalVector& l) -> ChemicalVector;

} // namespace Reaktor
