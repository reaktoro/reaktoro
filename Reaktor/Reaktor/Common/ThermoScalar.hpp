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
#include <Reaktor/Common/Vector.hpp>

namespace Reaktor {

// Forward declarations
struct ThermoVectorRow;
struct ThermoVectorConstRow;

/// A type that defines a scalar thermodynamic quantity.
/// A ThermoScalar instance not only holds the value of the
/// thermodynamic quantity, but also is partial temperature,
/// pressure and molar derivatives.
/// @see ThermoVector
class ThermoScalar
{
public:
	/// Construct a default ThermoScalar instance
    ThermoScalar();

    /// Construct a ThermoScalar instance
    /// @param val The scalar value of the thermodynamic quantity
	/// @param ddn The partial molar derivatives of the thermodynamic quantity
    ThermoScalar(double val, const Vector& ddn);
    
    /// Construct a ThermoScalar instance
    /// @param val The scalar value of the thermodynamic quantity
	/// @param ddt The partial temperature derivative of the thermodynamic quantity
	/// @param ddp The partial pressure derivative of the thermodynamic quantity 
    ThermoScalar(double val, double ddt, double ddp);

    /// Construct a ThermoScalar instance
    /// @param val The scalar value of the thermodynamic quantity
	/// @param ddt The partial temperature derivative of the thermodynamic quantity
	/// @param ddp The partial pressure derivative of the thermodynamic quantity 
	/// @param ddn The partial molar derivatives of the thermodynamic quantity
    ThermoScalar(double val, double ddt, double ddp, const Vector& ddn);

    /// Get the scalar value of the thermodynamic quantity
	auto val() const -> double;

	/// Get the partial temperature derivative of the thermodynamic quantity
	auto ddt() const -> double;

	/// Get the partial pressure derivative of the thermodynamic quantity 
	auto ddp() const -> double;

	/// Get the partial molar derivatives of the thermodynamic quantity
	auto ddn() const -> const Vector&;

	auto operator=(const ThermoVectorRow& row) -> ThermoScalar&;

	auto operator=(const ThermoVectorConstRow& row) -> ThermoScalar&;

private:
	/// The scalar value of the thermodynamic quantity
	double m_val;

	/// The partial temperature derivative of the thermodynamic quantity
	double m_ddt;

	/// The partial pressure derivative of the thermodynamic quantity 
	double m_ddp;

	/// The partial molar derivatives of the thermodynamic quantity
	Vector m_ddn;
};

} // namespace Reaktor
