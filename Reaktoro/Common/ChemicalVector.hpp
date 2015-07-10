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
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

// Forward declarations
class ChemicalScalar;
class ChemicalVectorBlock;
class ChemicalVectorBlockConst;
class ChemicalVectorRow;
class ChemicalVectorRowConst;
class ChemicalVectorRows;
class ChemicalVectorRowsConst;
class ChemicalVectorRowsCols;
class ChemicalVectorRowsColsConst;
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
    /// Return a ChemicalVector instance from a vector of amounts (in units of mol)
    static auto Composition(const Vector& n) -> ChemicalVector;

    /// Return a ChemicalVector instance from a list of amounts (in units of mol)
    static auto Composition(std::initializer_list<double> n) -> ChemicalVector;

    /// Construct a default ChemicalVector instance
    ChemicalVector();

    /// Construct a ChemicalVector instance with given number of species.
    ChemicalVector(unsigned nspecies);

    /// Construct a ChemicalVector instance with given dimensions.
    ChemicalVector(unsigned nrows, unsigned nspecies);

    /// Construct a ChemicalVector instance with given data members.
    /// @param val The vector value of the chemical property
    /// @param ddt The partial temperature derivatives of the vector chemical property
    /// @param ddp The partial pressure derivative of the vector chemical property
    /// @param ddn The partial molar derivatives of the vector chemical property
    ChemicalVector(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn);

    /// Construct a ChemicalVector instance with given rows.
    ChemicalVector(const ChemicalVectorRows& rows);

    /// Construct a ChemicalVector instance with given const rows.
    ChemicalVector(const ChemicalVectorRowsConst& rows);

    /// Construct a ChemicalVector instance with given block.
    ChemicalVector(const ChemicalVectorRowsCols& block);

    /// Construct a ChemicalVector instance with given const block.
    ChemicalVector(const ChemicalVectorRowsColsConst& block);

    /// Return a reference of a row of this ChemicalVector instance
    auto row(unsigned irow) -> ChemicalVectorRow;

    /// Return a reference of a row of this ChemicalVector instance
    auto row(unsigned irow, unsigned icol, unsigned ncols) -> ChemicalVectorRow;

    /// Return a const reference of a row of this ChemicalVector instance
    auto row(unsigned irow) const -> ChemicalVectorRowConst;

    /// Return a const reference of a row of this ChemicalVector instance
    auto row(unsigned irow, unsigned icol, unsigned ncols) const -> ChemicalVectorRowConst;

    /// Return a reference of a sequence of rows of this ChemicalVector instance
    auto rows(unsigned irow, unsigned nrows) -> ChemicalVectorBlock;

    /// Return a reference of a sequence of rows of this ChemicalVector instance
    auto rows(unsigned irow, unsigned icol, unsigned nrows, unsigned ncols) -> ChemicalVectorBlock;

    /// Return a const reference of a sequence of rows of this ChemicalVector instance
    auto rows(unsigned irow, unsigned nrows) const -> ChemicalVectorBlockConst;

    /// Return a const reference of a sequence of rows of this ChemicalVector instance
    auto rows(unsigned irow, unsigned icol, unsigned nrows, unsigned ncols) const -> ChemicalVectorBlockConst;

    /// Return a reference of some rows of this ChemicalVector instance
    auto rows(const Indices& irows) -> ChemicalVectorRows;

    /// Return a const reference of some rows of this ChemicalVector instance
    auto rows(const Indices& irows) const -> ChemicalVectorRowsConst;

    /// Return a reference of some rows and cols of this ChemicalVector instance
    auto rows(const Indices& irows, const Indices& icols) -> ChemicalVectorRowsCols;

    /// Return a const reference of some rows and cols of this ChemicalVector instance
    auto rows(const Indices& irows, const Indices& icols) const -> ChemicalVectorRowsColsConst;

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

    /// Return a reference of a row of this ChemicalVector instance
    auto operator[](unsigned irow) -> ChemicalVectorRow;

    /// Return a const reference of a row of this ChemicalVector instance
    auto operator[](unsigned irow) const -> ChemicalVectorRowConst;

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
    ChemicalVectorRow(ChemicalVector& vector, unsigned irow, unsigned icol, unsigned ncols);

    auto operator=(const ChemicalVectorRow& row) -> ChemicalVectorRow&;
    auto operator=(const ChemicalScalar& scalar) -> ChemicalVectorRow&;

    auto operator+=(const ChemicalVectorRow& other) -> ChemicalVectorRow&;
    auto operator+=(const ChemicalVectorRowConst& other) -> ChemicalVectorRow&;
    auto operator+=(const ChemicalScalar& other) -> ChemicalVectorRow&;
    auto operator+=(const ThermoScalar& other) -> ChemicalVectorRow&;

    auto operator-=(const ChemicalVectorRow& other) -> ChemicalVectorRow&;
    auto operator-=(const ChemicalVectorRowConst& other) -> ChemicalVectorRow&;
    auto operator-=(const ChemicalScalar& other) -> ChemicalVectorRow&;
    auto operator-=(const ThermoScalar& other) -> ChemicalVectorRow&;

    double& val;
    double& ddt;
    double& ddp;
    decltype(std::declval<Matrix>().row(0).segment(0, 0)) ddn;
};

/// An auxiliary type for the representation of the const view of a row of a ChemicalVector instance
class ChemicalVectorRowConst
{
public:
    ChemicalVectorRowConst(const ChemicalVector& vector, unsigned irow);
    ChemicalVectorRowConst(const ChemicalVector& vector, unsigned irow, unsigned icol, unsigned ncols);
    const double& val;
    const double& ddt;
    const double& ddp;
    decltype(std::declval<const Matrix>().row(0).segment(0, 0)) ddn;
};

/// An auxiliary type for the representation of the view of rows of a ChemicalVector instance
class ChemicalVectorRows
{
public:
    ChemicalVectorRows(ChemicalVector& vector, const Indices& irows);
    auto operator=(const ChemicalVector& vector) -> ChemicalVectorRows&;
    MatrixViewRows<Vector> val;
    MatrixViewRows<Vector> ddt;
    MatrixViewRows<Vector> ddp;
    MatrixViewRows<Matrix> ddn;
};

/// An auxiliary type for the representation of the const view of rows of a ChemicalVector instance
class ChemicalVectorRowsConst
{
public:
    ChemicalVectorRowsConst(const ChemicalVector& vector, const Indices& irows);
    MatrixViewRowsConst<Vector> val;
    MatrixViewRowsConst<Vector> ddt;
    MatrixViewRowsConst<Vector> ddp;
    MatrixViewRowsConst<Matrix> ddn;
};

/// An auxiliary type for the representation of the view of rows and cols of a ChemicalVector instance
class ChemicalVectorRowsCols
{
public:
    ChemicalVectorRowsCols(ChemicalVector& vector, const Indices& irows, const Indices& icols);
    auto operator=(const ChemicalVector& vector) -> ChemicalVectorRowsCols&;
    MatrixViewRows<Vector> val;
    MatrixViewRows<Vector> ddt;
    MatrixViewRows<Vector> ddp;
    MatrixViewRowsCols<Matrix> ddn;
};

/// An auxiliary type for the representation of the view of const rows and cols of a ChemicalVector instance
class ChemicalVectorRowsColsConst
{
public:
    ChemicalVectorRowsColsConst(const ChemicalVector& vector, const Indices& irows, const Indices& icols);
    MatrixViewRowsConst<Vector> val;
    MatrixViewRowsConst<Vector> ddt;
    MatrixViewRowsConst<Vector> ddp;
    MatrixViewRowsColsConst<Matrix> ddn;
};

/// An auxiliary type for the representation of the view of a block of a ChemicalVector instance
class ChemicalVectorBlock
{
public:
    ChemicalVectorBlock(ChemicalVector& vector, unsigned irow, unsigned nrows);
    ChemicalVectorBlock(ChemicalVector& vector, unsigned irow, unsigned icol, unsigned nrows, unsigned ncols);
    auto operator=(const ChemicalVectorBlock& block) -> ChemicalVectorBlock&;
    auto operator=(const ChemicalVector& vector) -> ChemicalVectorBlock&;
    decltype(std::declval<Vector>().segment(0, 0)) val;
    decltype(std::declval<Vector>().segment(0, 0)) ddt;
    decltype(std::declval<Vector>().segment(0, 0)) ddp;
    decltype(std::declval<Matrix>().block(0, 0, 0, 0)) ddn;
};

/// An auxiliary type for the representation of the const view of a block of a ChemicalVector instance
class ChemicalVectorBlockConst
{
public:
    ChemicalVectorBlockConst(const ChemicalVector& vector, unsigned irow, unsigned nrows);
    ChemicalVectorBlockConst(const ChemicalVector& vector, unsigned irow, unsigned icol, unsigned nrows, unsigned ncols);
    decltype(std::declval<const Vector>().segment(0, 0)) val;
    decltype(std::declval<const Vector>().segment(0, 0)) ddt;
    decltype(std::declval<const Vector>().segment(0, 0)) ddp;
    decltype(std::declval<const Matrix>().block(0, 0, 0, 0)) ddn;
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

/// Left-divide a ChemicalVector instance by a scalar
auto operator/(double scalar, const ChemicalVector& r) -> ChemicalVector;

/// Right-divide a ChemicalVector instance by a scalar
auto operator/(const ChemicalVector& l, double scalar) -> ChemicalVector;

/// Left-divide a ChemicalVector instance by a ChemicalScalar
auto operator/(const ChemicalScalar& scalar, const ChemicalVector& r) -> ChemicalVector;

/// Right-divide a ChemicalVector instance by a ChemicalScalar
auto operator/(const ChemicalVector& l, const ChemicalScalar& scalar) -> ChemicalVector;

/// Divide a ChemicalVector instance by another
auto operator/(const ChemicalVector& l, const ChemicalVector& r) -> ChemicalVector;

/// Multiply two ChemicalVector instances component-wise
auto operator%(const ChemicalVector& l, const ChemicalVector& r) -> ChemicalVector;

/// Return the power of a ChemicalVector instance
auto pow(const ChemicalVector& l, double power) -> ChemicalVector;

/// Return the natural exponential of a ChemicalVector instance
auto exp(const ChemicalVector& l) -> ChemicalVector;

/// Return the natural log of a ChemicalVector instance
auto log(const ChemicalVector& l) -> ChemicalVector;

/// Return the log10 of a ChemicalVector instance
auto log10(const ChemicalVector& l) -> ChemicalVector;

/// Return the sum of the components of a ChemicalVector instance
auto sum(const ChemicalVector& l) -> ChemicalScalar;

} // namespace Reaktoro
