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
class ThermoVectorBlock;
class ThermoVectorConstBlock;
class ThermoVectorConstRow;
class ThermoVectorRow;

/// Describe the thermodynamic properties and their partial temperature and pressure
/// derivatives of a collection of species or reactions.
/// ingroup Common
class ThermoVector
{
public:
    /// Construct a ThermoVector instance with given dimension
    /// @param nrows The number of rows of the vector quantities
    ThermoVector(unsigned nrows);

    /// Construct a ThermoVector instance
    /// @param val The values of the thermodynamic properties
    /// @param ddt The partial temperature derivatives of the thermodynamic properties
    /// @param ddp The partial pressure derivatives of the thermodynamic properties
    ThermoVector(const Vector& val, const Vector& ddt, const Vector& ddp);

    /// Return a reference of a row of this ThermoVector instance
    auto row(unsigned irow) -> ThermoVectorRow;

    /// Return a const reference of a row of this ThermoVector instance
    auto row(unsigned irow) const -> ThermoVectorConstRow;

    /// Return a reference of a block of this ThermoVector instance
    auto block(unsigned irow, unsigned nrows) -> ThermoVectorBlock;

    /// Return a const reference of a block of this ThermoVector instance
    auto block(unsigned irow, unsigned nrows) const -> ThermoVectorConstBlock;

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
};

/// An auxiliary type for the representation of the const view of a row of a ThermoVector instance
class ThermoVectorConstRow
{
public:
    ThermoVectorConstRow(const ThermoVector& properties, unsigned irow);
    const double& val;
    const double& ddt;
    const double& ddp;
};

/// An auxiliary type for the representation of the view of a block of a ThermoVector instance
class ThermoVectorBlock
{
public:
    ThermoVectorBlock(ThermoVector& vector, unsigned irow, unsigned nrows);
    ThermoVectorBlock(ThermoVector& vector, unsigned irow, unsigned icol, unsigned nrows, unsigned ncols);
    auto operator=(const ThermoVectorBlock& block) -> ThermoVectorBlock&;
    auto operator=(const ThermoVector& vector) -> ThermoVectorBlock&;
    decltype(std::declval<Vector>().segment(0, 0)) val;
    decltype(std::declval<Vector>().segment(0, 0)) ddt;
    decltype(std::declval<Vector>().segment(0, 0)) ddp;
};

/// An auxiliary type for the representation of the const view of a block of a ThermoVector instance
class ThermoVectorConstBlock
{
public:
    ThermoVectorConstBlock(const ThermoVector& vector, unsigned irow, unsigned nrows);
    ThermoVectorConstBlock(const ThermoVector& vector, unsigned irow, unsigned icol, unsigned nrows, unsigned ncols);
    decltype(std::declval<const Vector>().segment(0, 0)) val;
    decltype(std::declval<const Vector>().segment(0, 0)) ddt;
    decltype(std::declval<const Vector>().segment(0, 0)) ddp;
};

/// Compare two ThermoVector instances for equality
auto operator==(const ThermoVector& l, const ThermoVector& r) -> bool;

/// A type used to define the function signature for the calculation of many thermodynamic properties.
/// @see ThermoVector, ThermoScalar, ThermoScalarFunction
using ThermoVectorFunction = std::function<ThermoVector(double, double)>;

} // namespace Reaktoro
