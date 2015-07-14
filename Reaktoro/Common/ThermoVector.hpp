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
class ThermoVectorRows;
class ThermoVectorConstRows;
class ThermoVectorConstRow;
class ThermoVectorRow;

/// Describe the thermodynamic properties and their partial temperature and pressure
/// derivatives of a collection of species or reactions.
/// ingroup Common
class ThermoVector
{
public:
    /// Construct a default ThermoVector instance
    ThermoVector();

    /// Construct a ThermoVector instance with given dimension
    /// @param nrows The number of rows of the vector quantities
    explicit ThermoVector(unsigned nrows);

    /// Construct a ThermoVector instance
    /// @param val The values of the thermodynamic properties
    /// @param ddt The partial temperature derivatives of the thermodynamic properties
    /// @param ddp The partial pressure derivatives of the thermodynamic properties
    ThermoVector(const Vector& val, const Vector& ddt, const Vector& ddp);

    /// Resize the ThermoVector instance.
    auto resize(unsigned nrows) -> void;

    /// Return a reference of a row of this ThermoVector instance
    auto row(unsigned irow) -> ThermoVectorRow;

    /// Return a const reference of a row of this ThermoVector instance
    auto row(unsigned irow) const -> ThermoVectorConstRow;

    /// Return a reference of a block of this ThermoVector instance
    auto rows(unsigned irow, unsigned nrows) -> ThermoVectorRows;

    /// Return a const reference of a block of this ThermoVector instance
    auto rows(unsigned irow, unsigned nrows) const -> ThermoVectorConstRows;

    /// Assign-addition of a ThermoVector instance.
    auto operator+=(const ThermoVector& other) -> ThermoVector&;

    /// Assign-subtraction of a ThermoVector instance.
    auto operator-=(const ThermoVector& other) -> ThermoVector&;

    /// Assign-multiplication of a ChemicalVector instance.
    auto operator*=(double scalar) -> ThermoVector&;

    /// Assign-division of a ChemicalVector instance.
    auto operator/=(double scalar) -> ThermoVector&;

    /// Return a reference of a row of this ThermoVector instance
    auto operator[](unsigned irow) -> ThermoVectorRow;

    /// Return a const reference of a row of this ThermoVector instance
    auto operator[](unsigned irow) const -> ThermoVectorConstRow;

    /// The values of the thermodynamic properties
    Vector val;

    /// The partial temperature derivatives of the thermodynamic properties
    Vector ddt;

    /// The partial pressure derivatives of the thermodynamic properties
    Vector ddp;
};

/// An auxiliary type for the representation of the view of a row of a ThermoVector instance
class ThermoVectorRow
{
public:
    ThermoVectorRow(ThermoVector& vector, unsigned irow);
    auto operator=(const ThermoScalar& property) -> ThermoVectorRow&;
    double& val;
    double& ddt;
    double& ddp;
    operator ThermoScalar();
};

/// An auxiliary type for the representation of the const view of a row of a ThermoVector instance
class ThermoVectorConstRow
{
public:
    ThermoVectorConstRow(const ThermoVector& properties, unsigned irow);
    const double& val;
    const double& ddt;
    const double& ddp;
    operator ThermoScalar();
};

/// An auxiliary type for the representation of the view of a block of a ThermoVector instance
class ThermoVectorRows
{
public:
    ThermoVectorRows(ThermoVector& vector, unsigned irow, unsigned nrows);
    ThermoVectorRows(ThermoVector& vector, unsigned irow, unsigned icol, unsigned nrows, unsigned ncols);
    auto operator=(const ThermoVectorRows& block) -> ThermoVectorRows&;
    auto operator=(const ThermoVector& vector) -> ThermoVectorRows&;
    decltype(std::declval<Vector>().segment(0, 0)) val;
    decltype(std::declval<Vector>().segment(0, 0)) ddt;
    decltype(std::declval<Vector>().segment(0, 0)) ddp;
};

/// An auxiliary type for the representation of the const view of a block of a ThermoVector instance
class ThermoVectorConstRows
{
public:
    ThermoVectorConstRows(const ThermoVector& vector, unsigned irow, unsigned nrows);
    ThermoVectorConstRows(const ThermoVector& vector, unsigned irow, unsigned icol, unsigned nrows, unsigned ncols);
    decltype(std::declval<const Vector>().segment(0, 0)) val;
    decltype(std::declval<const Vector>().segment(0, 0)) ddt;
    decltype(std::declval<const Vector>().segment(0, 0)) ddp;
};

/// A type used to define the function signature for the calculation of many thermodynamic properties.
/// @see ThermoVector, ThermoScalar, ThermoScalarFunction
using ThermoVectorFunction = std::function<ThermoVector(double, double)>;

/// Compare two ThermoVector instances for equality
auto operator==(const ThermoVector& l, const ThermoVector& r) -> bool;

/// Unary addition operator for a ThermoVector instance
auto operator+(const ThermoVector& l) -> ThermoVector;

/// Unary subtraction operator for a ThermoVector instance
auto operator-(const ThermoVector& l) -> ThermoVector;

/// Add two ThermoVector instances
auto operator+(const ThermoVector& l, const ThermoVector& r) -> ThermoVector;

/// Subtract two ThermoVector instances
auto operator-(const ThermoVector& l, const ThermoVector& r) -> ThermoVector;

/// Left-multiply a ThermoVector instance by a scalar
auto operator*(double scalar, const ThermoVector& r) -> ThermoVector;

/// Right-multiply a ThermoVector instance by a scalar
auto operator*(const ThermoVector& l, double scalar) -> ThermoVector;

/// Left-multiply a ThermoVector instance by a ThermoScalar instance
auto operator*(const ThermoScalar& scalar, const ThermoVector& r) -> ThermoVector;

/// Right-multiply a ThermoVector instance by a ThermoScalar instance
auto operator*(const ThermoVector& l, const ThermoScalar& scalar) -> ThermoVector;

/// Left-divide a ThermoVector instance by a scalar
auto operator/(double scalar, const ThermoVector& r) -> ThermoVector;

/// Right-divide a ThermoVector instance by a scalar
auto operator/(const ThermoVector& l, double scalar) -> ThermoVector;

/// Left-divide a ThermoVector instance by a ThermoScalar
auto operator/(const ThermoScalar& scalar, const ThermoVector& r) -> ThermoVector;

/// Right-divide a ThermoVector instance by a ThermoScalar
auto operator/(const ThermoVector& l, const ThermoScalar& scalar) -> ThermoVector;

/// Divide a ThermoVector instance by another
auto operator/(const ThermoVector& l, const ThermoVector& r) -> ThermoVector;

/// Multiply two ThermoVector instances component-wise
auto operator%(const ThermoVector& l, const ThermoVector& r) -> ThermoVector;

/// Return the power of a ThermoVector instance
auto pow(const ThermoVector& l, double power) -> ThermoVector;

/// Return the natural exponential of a ThermoVector instance
auto exp(const ThermoVector& l) -> ThermoVector;

/// Return the natural log of a ThermoVector instance
auto log(const ThermoVector& l) -> ThermoVector;

/// Return the log10 of a ThermoVector instance
auto log10(const ThermoVector& l) -> ThermoVector;

/// Return the sum of the components of a ThermoVector instance
auto sum(const ThermoVector& l) -> ThermoScalar;

} // namespace Reaktoro
