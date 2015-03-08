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
#include <Reaktor/Common/Matrix.hpp>

namespace Reaktor {

// Forward declarations
struct ChemicalVectorRow;
struct ChemicalVectorConstRow;

/// A type that defines a scalar chemical property.
/// A chemical property means here any property that depends on
/// temperature, pressure and composition. A ChemicalScalar instance
/// not only holds the value of the chemical property, but also is partial
/// temperature, pressure and molar derivatives.
/// @see ChemicalVector
class ChemicalScalar
{
public:
    /// Construct a default ChemicalScalar instance
    ChemicalScalar();

    /// Construct a ChemicalScalar instance
    /// @param val The scalar value of the thermodynamic quantity
    /// @param ddt The partial temperature derivative of the thermodynamic quantity
    /// @param ddp The partial pressure derivative of the thermodynamic quantity
    /// @param ddn The partial molar derivatives of the thermodynamic quantity
    ChemicalScalar(double val, double ddt, double ddp, const Vector& ddn);

    /// Construct a ChemicalScalar instance from a ChemicalVectorRow instance
    /// @param row The row of a ChemicalVector instance
    ChemicalScalar(const ChemicalVectorRow& row);

    /// Construct a ChemicalScalar instance from a ChemicalVectorConstRow instance
    /// @param row The row of a const ChemicalVector instance
    ChemicalScalar(const ChemicalVectorConstRow& row);

    /// Get the scalar value of the thermodynamic quantity
    auto val() const -> double;

    /// Get the partial temperature derivative of the thermodynamic quantity
    auto ddt() const -> double;

    /// Get the partial pressure derivative of the thermodynamic quantity
    auto ddp() const -> double;

    /// Get the partial molar derivatives of the thermodynamic quantity
    auto ddn() const -> const Vector&;

    /// Assign a row of a ChemicalVector instance to this ChemicalScalar instance
    auto operator=(const ChemicalVectorRow& row) -> ChemicalScalar&;

    /// Assign a row of a ChemicalVector instance to this ChemicalScalar instance
    auto operator=(const ChemicalVectorConstRow& row) -> ChemicalScalar&;

private:
    /// The scalar value of the chemical property
    double m_val;

    /// The partial temperature derivative of the chemical property
    double m_ddt;

    /// The partial pressure derivative of the chemical property
    double m_ddp;

    /// The partial molar derivatives of the chemical property
    Vector m_ddn;
};

/// Compares two ChemicalScalar instances for equality
auto operator==(const ChemicalScalar& l, const ChemicalScalar& r) -> bool;

/// A type used to define the function signature for the calculation of a chemical property.
/// @see ChemicalScalar, ChemicalVectorFunction
using ChemicalScalarFunction = std::function<ChemicalScalar(double, double, const Vector&)>;

} // namespace Reaktor
