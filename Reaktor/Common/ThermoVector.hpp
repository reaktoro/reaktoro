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
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
class  ThermoScalar;
struct ThermoVectorConstRow;
struct ThermoVectorRow;

/// A type that defines a vector thermodynamic quantity.
/// A ThermoVector instance not only holds the value of the
/// thermodynamic quantity, but also is partial temperature,
/// pressure and molar derivatives.
/// @see ThermoScalar
class ThermoVector
{
public:
    // Forward declaration
    struct Row;
    struct ConstRow;

	/// Construct a default ThermoVector instance
    ThermoVector();

    /// Construct a ThermoVector instance with given dimensions
    /// @param nrows The number of rows of the vector quantities
    /// @param nrows The number of columns of the matrix quantities
    ThermoVector(unsigned nrows, unsigned ncols);

	/// Construct a ThermoVector instance with given data members
    /// @param val The vector value of the thermodynamic quantity
    /// @param ddt The partial temperature derivatives of the vector thermodynamic quantity
    /// @param ddp The partial pressure derivative of the vector thermodynamic quantity
    /// @param ddn The partial molar derivatives of the vector thermodynamic quantity
    ThermoVector(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn);

    /// Get the vector value of the thermodynamic quantity
    auto val() const -> const Vector&;

    /// Get the partial temperature derivatives of the vector thermodynamic quantity
    auto ddt() const -> const Vector&;

    /// Get the partial pressure derivative of the vector thermodynamic quantity
    auto ddp() const -> const Vector&;

    /// Get the partial molar derivatives of the vector thermodynamic quantity
    auto ddn() const -> const Matrix&;

    /// Get a reference of a row of this ThermoVector instance
    auto row(unsigned irow) -> ThermoVectorRow;

    /// Get a const reference of a row of this ThermoVector instance
    auto row(unsigned irow) const -> ThermoVectorConstRow;

private:
    /// The vector value of the thermodynamic quantity
    Vector m_val;

    /// The partial temperature derivatives of the vector thermodynamic quantity
    Vector m_ddt;

    /// The partial pressure derivative of the vector thermodynamic quantity
    Vector m_ddp;

    /// The partial molar derivatives of the vector thermodynamic quantity
    Matrix m_ddn;
};

/// An auxiliary type for the representation of the view of a row of a ThermoVector instance
struct ThermoVectorRow
{
    ThermoVectorRow(const ThermoVector& vector, unsigned irow);
    auto operator=(const ThermoScalar& scalar) -> ThermoVectorRow&;
    VectorRow val;
    VectorRow ddt;
    VectorRow ddp;
    MatrixRow ddn;
};

/// An auxiliary type for the representation of the const view of a row of a ThermoVector instance
struct ThermoVectorConstRow
{
    ThermoVectorConstRow(const ThermoVector& vector, unsigned irow);
    const VectorRow val;
    const VectorRow ddt;
    const VectorRow ddp;
    const MatrixRow ddn;
};

/// Compares two ThermoVector instances for equality
auto operator==(const ThermoVector& l, const ThermoVector& r) -> bool;

} // namespace Reaktor
