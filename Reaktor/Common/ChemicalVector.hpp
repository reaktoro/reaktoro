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
class  ChemicalScalar;
struct ChemicalVectorConstRow;
struct ChemicalVectorRow;

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

    /// Get the vector value of the chemical property
    auto val() const -> const Vector&;

    /// Get the partial temperature derivatives of the vector chemical property
    auto ddt() const -> const Vector&;

    /// Get the partial pressure derivative of the vector chemical property
    auto ddp() const -> const Vector&;

    /// Get the partial molar derivatives of the vector chemical property
    auto ddn() const -> const Matrix&;

    /// Get a reference of a row of this ChemicalVector instance
    auto row(unsigned irow) -> ChemicalVectorRow;

    /// Get a const reference of a row of this ChemicalVector instance
    auto row(unsigned irow) const -> ChemicalVectorConstRow;

    friend class ChemicalVectorRow;
    friend class ChemicalVectorConstRow;

private:
    /// The vector value of the chemical property
    Vector m_val;

    /// The partial temperature derivatives of the vector chemical property
    Vector m_ddt;

    /// The partial pressure derivative of the vector chemical property
    Vector m_ddp;

    /// The partial molar derivatives of the vector chemical property
    Matrix m_ddn;
};

/// An auxiliary type for the representation of the view of a row of a ChemicalVector instance
struct ChemicalVectorRow
{
    ChemicalVectorRow(ChemicalVector& vector, unsigned irow);
    auto operator=(const ChemicalScalar& scalar) -> ChemicalVectorRow&;
    double& val;
    double& ddt;
    double& ddp;
    decltype(std::declval<Matrix>().row(0)) ddn;
};

/// An auxiliary type for the representation of the const view of a row of a ChemicalVector instance
struct ChemicalVectorConstRow
{
    ChemicalVectorConstRow(const ChemicalVector& vector, unsigned irow);
    const double& val;
    const double& ddt;
    const double& ddp;
    decltype(std::declval<const Matrix>().row(0)) ddn;
};

/// Compares two ChemicalVector instances for equality
auto operator==(const ChemicalVector& l, const ChemicalVector& r) -> bool;

/// A type used to define the function signature for the calculation of a vector of chemical properties.
/// @see ChemicalVector, ChemicalScalarFunction
typedef std::function<ChemicalVector(double, double, const Vector&)> ChemicalVectorFunction;

} // namespace Reaktor
