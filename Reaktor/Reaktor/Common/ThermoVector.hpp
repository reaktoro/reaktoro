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
class ThermoScalar;
class ThermoVectorRow;

/// A type that defines a vector thermodynamic quantity
///
/// A ThermoVector instance not only holds the value of the
/// thermodynamic quantity, but also is partial temperature,
/// pressure and molar derivatives.
///
/// @see ThermoScalar
class ThermoVector
{
public:
	/// Construct a default ThermoVector instance
	ThermoVector();

	/// Construct a ThermoVector instance with given data members
    ThermoVector(const Vector& val, const Vector& ddt, const Vector& ddp, const Matrix& ddn);

    /// Construct a ThermoVector instance with given number of rows and species
    /// @param nrows The number rows that the vector thermodynamic quantity has
    /// @param nspecies The number of species for which the partial molar derivatives are calculated
    ThermoVector(unsigned nrows, unsigned nspecies);

    /// Get a view of the data on a given row of this ThermoVector instance
    /// @param irow The index of the row
    auto row(const Index& irow) -> ThermoVectorRow;

    /// Get a const view of the data on a given row of this ThermoVector instance
    /// @param irow The index of the row
    auto row(const Index& irow) const -> const ThermoVectorRow;

    /// Return a zero ThermoVector instance
    /// @param nrows The number rows that the vector thermodynamic quantity has
    /// @param nspecies The number of species for which the partial molar derivatives are calculated
    static auto zero(unsigned nrows, unsigned nspecies) -> ThermoVector;

    /// The value (vector) of the thermodynamic quantity
    Vector val;

    /// The partial derivative of the thermodynamic quantity w.r.t. temperature (in units of K)
    Vector ddt;

    /// The partial derivative of the thermodynamic quantity w.r.t. pressure (in units of Pa)
    Vector ddp;

    /// The partial derivative of the thermodynamic quantity w.r.t. composition (in units of mol)
    Matrix ddn;
};

/// A type that defines a view of a row of a vector thermodynamic quantity
/// @see ThermoScalar, ThermoVector
class ThermoVectorRow
{
public:
    /// Construct a ThermoVectorRow instance
    /// @param vector The thermodynamic vector quantity and its partial derivatives
    /// @param irow The index of the row of the thermodynamic vector quantity
    ThermoVectorRow(const ThermoVector& vector, unsigned irow);

    /// Assign a ThermoScalar instance to this row of a ThermoVector instance
    /// @param scalar The thermodynamic scalar quantity
    auto operator=(const ThermoScalar& scalar) -> ThermoVectorRow&;

    /// The view of a row of a thermodynamic quantity
    VectorRow val;

    /// The view of a row of the partial derivative of the thermodynamic quantity w.r.t. temperature (in units of K)
    VectorRow ddt;

    /// The view of a row of the partial derivative of the thermodynamic quantity w.r.t. pressure (in units of Pa)
    VectorRow ddp;

    /// The view of a row of the partial derivative of the thermodynamic quantity w.r.t. composition (in units of mol)
    MatrixRow ddn;
};

} // namespace Reaktor
