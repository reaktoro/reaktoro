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

// Eigen includes
#include <Reaktor/eigen/Core>

// Reaktor includes
#include <Reaktor/Common/Index.hpp>

namespace Reaktor {

/// Define an alias to the vector type of the Eigen library
using Vector = Eigen::VectorXd;

/// A type used to define a view of some rows of a Vector instance
class VectorViewRows
{
public:
    /// Construct a VectorView instance
    /// @param vec The vector for which this view is defined
    /// @param irows The indices of the rows of this vector view
    VectorViewRows(Vector& vec, const Indices& irows);

    /// Assign a Vector instance to this VectorView instance
    auto operator=(const Vector& other) -> VectorViewRows&;

    /// Convert this VectorView instance into a Vector instance
    operator Vector() const;

private:
    /// The vector for which this view is defined
    Vector& vec;

    /// The indices of the rows of this vector view
    const Indices& irows;
};

/// A type used to define a const view of some rows of a Vector instance
class VectorViewRowsConst
{
public:
    /// Construct a VectorViewConst instance
    /// @param vec The vector for which this view is defined
    /// @param irows The indices of the rows of this vector view
    VectorViewRowsConst(const Vector& vec, const Indices& irows);

    /// Convert this VectorViewConst instance into a Vector instance
    operator Vector() const;

private:
    /// The vector for which this view is defined
    const Vector& vec;

    /// The indices of the rows of this vector view
    const Indices& irows;
};

/// Return an expression of a zero vector
/// @param rows The number of rows
/// @return The expression of a zero vector
auto zeros(unsigned rows) -> decltype(Vector::Zero(rows));

/// Return an expression of a vector with entries equal to one
/// @param rows The number of rows
/// @return The expression of a vector with entries equal to one
auto ones(unsigned rows) -> decltype(Vector::Ones(rows));

/// Return a view of some rows of a vector
/// @param vec The vector for which the view is created
/// @param irows The indices of the rows of the vector
auto rows(Vector& vec, const Indices& irows) -> VectorViewRows;

/// Return a const view of some rows of a vector
/// @param vec The vector for which the view is created
/// @param irows The indices of the rows of the vector
auto rows(const Vector& vec, const Indices& irows) -> VectorViewRowsConst;

} // namespace Reaktor

#include "Vector.hxx"
