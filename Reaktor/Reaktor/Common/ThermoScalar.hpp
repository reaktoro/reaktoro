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

namespace Reaktor {

// Forward declarations
class ThermoVectorRow;

/// A type that defines a scalar thermodynamic quantity
///
/// A ThermoScalar instance not only holds the value of the
/// thermodynamic quantity, but also is partial temperature,
/// pressure and molar derivatives.
///
/// @see ThermoVector
class ThermoScalar
{
public:
	/// Construct a default ThermoScalar instance
    ThermoScalar();

    /// Construct a ThermoScalar instance with given data members
    ThermoScalar(double val, double ddt, double ddp, const Vector& ddn);

	/// Construct a ThermoScalar instance with given number of species
	/// @param nspecies The number of species for which the partial molar derivatives are calculated
    explicit ThermoScalar(unsigned nspecies);

	/// Construct a ThermoScalar instance from a row of a ThermoVector instance
	/// @param res The ThermoVectorRow instance from which the ThermoScalar instance is built
	ThermoScalar(const ThermoVectorRow& row);

	/// Return a zero ThermoScalar instance
	/// @param nspecies The number of species for which the partial molar derivatives are calculated
	static auto zero(unsigned nspecies) -> ThermoScalar;

	/// The value (scalar) of the thermodynamic quantity
	double val;

	/// The partial derivative of the thermodynamic quantity w.r.t. temperature (in units of K)
	double ddt;

	/// The partial derivative of the thermodynamic quantity w.r.t. pressure (in units of Pa)
	double ddp;

	/// The partial derivative of the thermodynamic quantity w.r.t. composition (in units of mol)
	Vector ddn;
};

} // namespace Reaktor
