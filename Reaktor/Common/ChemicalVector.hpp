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
class  ChemicalScalar;
struct ChemicalVectorConstRow;
struct ChemicalVectorRow;

/// A type that defines a vector thermodynamic quantity.
/// A ChemicalVector instance not only holds the value of the
/// thermodynamic quantity, but also is partial temperature,
/// pressure and molar derivatives.
/// @see ChemicalScalar
class ChemicalVector
{
public:
    // Forward declaration
    struct Row;
    struct ConstRow;

	/// Construct a default ChemicalVector instance
    ChemicalVector();

    /// Construct a ChemicalVector instance with given dimensions
    /// @param nrows The number of rows of the vector quantities
    /// @param nrows The number of columns of the matrix quantities
    ChemicalVector(unsigned nrows, unsigned ncols);

	/// Construct a ChemicalVector instance with given data members
    /// @param val The vector value of the thermodynamic quantity
    /// @param ddt The partial temperature derivatives of the vector thermodynamic quantity
    /// @param ddp The partial pressure derivative of the vector thermodynamic quantity
    /// @param ddn The partial molar derivatives of the vector thermodynamic quantity
    ChemicalVector(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn);

    /// Get the vector value of the thermodynamic quantity
    auto val() const -> const Vector&;

    /// Get the partial temperature derivatives of the vector thermodynamic quantity
    auto ddt() const -> const Vector&;

    /// Get the partial pressure derivative of the vector thermodynamic quantity
    auto ddp() const -> const Vector&;

    /// Get the partial molar derivatives of the vector thermodynamic quantity
    auto ddn() const -> const Matrix&;

    /// Get a reference of a row of this ChemicalVector instance
    auto row(unsigned irow) -> ChemicalVectorRow;

    /// Get a const reference of a row of this ChemicalVector instance
    auto row(unsigned irow) const -> ChemicalVectorConstRow;

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

/// An auxiliary type for the representation of the view of a row of a ChemicalVector instance
struct ChemicalVectorRow
{
    ChemicalVectorRow(const ChemicalVector& vector, unsigned irow);
    auto operator=(const ChemicalScalar& scalar) -> ChemicalVectorRow&;
    VectorRow val;
    VectorRow ddt;
    VectorRow ddp;
    MatrixRow ddn;
};

/// An auxiliary type for the representation of the const view of a row of a ChemicalVector instance
struct ChemicalVectorConstRow
{
    ChemicalVectorConstRow(const ChemicalVector& vector, unsigned irow);
    const VectorRow val;
    const VectorRow ddt;
    const VectorRow ddp;
    const MatrixRow ddn;
};

/// Compares two ChemicalVector instances for equality
auto operator==(const ChemicalVector& l, const ChemicalVector& r) -> bool;

} // namespace Reaktor
