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
#include <Reaktor/Common/Vector.hpp>
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

// Forward declarations
class ThermoScalar;
class ThermoVectorRow;

/// A type that defines a vector thermodynamic quantity with its partial temperature, pressure and molar derivatives
class ThermoVector
{
public:
    /// Construct a default ThermoVector instance
    ThermoVector();

    /// Set the vector value of the thermodynamic quantity
    auto val(const Vector& val) -> ThermoVector&;

    /// Set the partial derivative of the thermodynamic quantity with respect to temperature (temperature in units of K)
    auto ddt(const Vector& ddt) -> ThermoVector&;

    /// Set the partial derivative of the thermodynamic quantity with respect to pressure (pressure in units of Pa)
    auto ddp(const Vector& ddp) -> ThermoVector&;

    /// Set the partial derivative of the thermodynamic quantity with respect to composition (composition in units of mol)
    auto ddn(const Matrix& ddn) -> ThermoVector&;

    /// Get the vector value of the thermodynamic quantity
    auto val() const -> const Vector&;

    /// Get the partial derivative of the thermodynamic quantity with respect to temperature (temperature in units of K)
    auto ddt() const -> const Vector&;

    /// Get the partial derivative of the thermodynamic quantity with respect to pressure (pressure in units of Pa)
    auto ddp() const -> const Vector&;

    /// Get the partial derivative of the thermodynamic quantity with respect to composition (composition in units of mol)
    auto ddn() const -> const Matrix&;

    /// Get a view of a row of the thermodynamic vector quantity
    /// @param irow The index of the row
    auto row(const Index& irow) -> ThermoVectorRow;

    /// Get a view of a row of the thermodynamic vector quantity
    /// @param irow The index of the row
    auto row(const Index& irow) const -> const ThermoVectorRow;

private:
    /// The vector value of the thermodynamic quantity
    Vector m_val;

    /// The partial derivative of the thermodynamic quantity with respect to temperature (temperature in units of K)
    Vector m_ddt;

    /// The partial derivative of the thermodynamic quantity with respect to pressure (pressure in units of Pa)
    Vector m_ddp;

    /// The partial derivative of the thermodynamic quantity with respect to composition (composition in units of mol)
    Matrix m_ddn;
};

/// A type that defines a view of a row of vector thermodynamic quantity with its partial temperature, pressure and molar derivatives
class ThermoVectorRow
{
public:
    /// Construct a default ThermoVectorRow instance
    ThermoVectorRow();

    /// Construct a ThermoVectorRow instance
    /// @param vector The thermodynamic vector quantity and its partial derivatives
    /// @param irow The index of the row of the thermodynamic vector quantity
    ThermoVectorRow(const ThermoVector& vector, unsigned irow);

    /// Assign a ThermoScalar instance to this row of a ThermoVector instance
    /// @param scalar The thermodynamic scalar quantity
    auto operator=(const ThermoScalar& scalar) -> ThermoVectorRow&;

    /// Get the vector value of the thermodynamic quantity
    auto val() const -> const VectorRow&;

    /// Get the partial derivative of the thermodynamic quantity with respect to temperature (temperature in units of K)
    auto ddt() const -> const VectorRow&;

    /// Get the partial derivative of the thermodynamic quantity with respect to pressure (pressure in units of Pa)
    auto ddp() const -> const VectorRow&;

    /// Get the partial derivative of the thermodynamic quantity with respect to composition (composition in units of mol)
    auto ddn() const -> const MatrixRow&;

private:
    /// The vector value of the thermodynamic quantity
    VectorRow m_val;

    /// The partial derivative of the thermodynamic quantity with respect to temperature (temperature in units of K)
    VectorRow m_ddt;

    /// The partial derivative of the thermodynamic quantity with respect to pressure (pressure in units of Pa)
    VectorRow m_ddp;

    /// The partial derivative of the thermodynamic quantity with respect to composition (composition in units of mol)
    MatrixRow m_ddn;
};

/// Get the vector value of the thermodynamic quantity
inline auto val(const ThermoVector& vector) -> const Vector&
{
    return vector.val();
}

/// Get the partial derivative of the thermodynamic quantity with respect to temperature (temperature in units of K)
inline auto ddt(const ThermoVector& vector) -> const Vector&
{
    return vector.ddt();
}

/// Get the partial derivative of the thermodynamic quantity with respect to pressure (pressure in units of Pa)
inline auto ddp(const ThermoVector& vector) -> const Vector&
{
    return vector.ddp();
}

/// Get the partial derivative of the thermodynamic quantity with respect to composition (composition in units of mol)
inline auto ddn(const ThermoVector& vector) -> const Matrix&
{
    return vector.ddn();
}

} // namespace Reaktor
