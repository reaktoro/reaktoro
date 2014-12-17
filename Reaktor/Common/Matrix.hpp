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

/// Define an alias to the matrix type of the Eigen library
using Matrix = Eigen::MatrixXd;

/// A type used to define a view of some rows of a matrix instance
class MatrixViewRows
{
public:
    /// Construct a MatrixView instance
    /// @param mat The matrix for which this view is defined
    /// @param irows The indices of the rows of this matrix view
    MatrixViewRows(Matrix& mat, const Indices& irows);

    /// Assign a matrix instance to this MatrixView instance
    auto operator=(const Matrix& other) -> MatrixViewRows&;

    /// Convert this MatrixViewRows instance into a matrix instance
    operator Matrix() const;

private:
    /// The matrix for which this view is defined
    Matrix& mat;

    /// The indices of the rows of this matrix view
    const Indices& irows;
};

/// A type used to define a const view of some rows of a matrix instance
class MatrixViewRowsConst
{
public:
    /// Construct a MatrixViewConst instance
    /// @param mat The matrix for which this view is defined
    /// @param irows The indices of the rows of this matrix view
    MatrixViewRowsConst(const Matrix& mat, const Indices& irows);

    /// Convert this MatrixViewConst instance into a matrix instance
    operator Matrix() const;

private:
    /// The matrix for which this view is defined
    const Matrix& mat;

    /// The indices of the rows of this matrix view
    const Indices& irows;
};

/// A type used to define a view of some columns of a matrix instance
class MatrixViewCols
{
public:
    /// Construct a MatrixViewCols instance
    /// @param mat The matrix for which this view is defined
    /// @param icols The indices of the columns of this matrix view
    MatrixViewCols(Matrix& mat, const Indices& icols);

    /// Assign a matrix instance to this MatrixViewCols instance
    auto operator=(const Matrix& other) -> MatrixViewCols&;

    /// Convert this MatrixViewCols instance into a matrix instance
    operator Matrix() const;

private:
    /// The matrix for which this view is defined
    Matrix& mat;

    /// The indices of the columns of this matrix view
    const Indices& icols;
};

/// A type used to define a const view of some columns of a matrix instance
class MatrixViewColsConst
{
public:
    /// Construct a MatrixViewColsConst instance
    /// @param mat The matrix for which this view is defined
    /// @param icols The indices of the columns of this matrix view
    MatrixViewColsConst(const Matrix& mat, const Indices& icols);

    /// Convert this MatrixViewColsConst instance into a matrix instance
    operator Matrix() const;

private:
    /// The matrix for which this view is defined
    const Matrix& mat;

    /// The indices of the columns of this matrix view
    const Indices& icols;
};

/// A type used to define a view of some rows and columns of a matrix instance
class MatrixViewRowsCols
{
public:
    /// Construct a MatrixViewRowsCols instance
    /// @param mat The matrix for which this view is defined
    /// @param irows The indices of the rows of this matrix view
    /// @param icols The indices of the columns of this matrix view
    MatrixViewRowsCols(Matrix& mat, const Indices& irows, const Indices& icols);

    /// Assign a matrix instance to this MatrixViewCols instance
    auto operator=(const Matrix& other) -> MatrixViewRowsCols&;

    /// Convert this MatrixViewRowsCols instance into a matrix instance
    operator Matrix() const;

private:
    /// The matrix for which this view is defined
    Matrix& mat;

    /// The indices of the rows of this matrix view
    const Indices& irows;

    /// The indices of the columns of this matrix view
    const Indices& icols;
};

/// A type used to define a const view of some rows and columns of a matrix instance
class MatrixViewRowsColsConst
{
public:
    /// Construct a MatrixViewRowsColsConst instance
    /// @param mat The matrix for which this view is defined
    /// @param irows The indices of the rows of this matrix view
    /// @param icols The indices of the columns of this matrix view
    MatrixViewRowsColsConst(const Matrix& mat, const Indices& irows, const Indices& icols);

    /// Convert this MatrixViewRowsColsConst instance into a matrix instance
    operator Matrix() const;

private:
    /// The matrix for which this view is defined
    const Matrix& mat;

    /// The indices of the rows of this matrix view
    const Indices& irows;

    /// The indices of the columns of this matrix view
    const Indices& icols;
};

/// Return an expression of an identity matrix
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of an identity matrix
auto identity(unsigned rows, unsigned cols) -> decltype(Matrix::Identity(rows, cols));

/// Return an expression of a zero matrix
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a zero matrix
auto zeros(unsigned rows, unsigned cols) -> decltype(Matrix::Zero(rows, cols));

/// Return an expression of a matrix with entries equal to one
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a matrix with entries equal to one
auto ones(unsigned rows, unsigned cols) -> decltype(Matrix::Ones(rows, cols));

/// Return a view of a sequence of rows of a matrix
/// @param start The row index of the start of the sequence
/// @param num The number of rows in the sequence
template<typename MatrixType>
auto rows(MatrixType&& mat, unsigned start, unsigned num) -> decltype(mat.middleRows(start, num));

/// Return a view of some rows of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
auto rows(Matrix& mat, const Indices& irows) -> MatrixViewRows;

/// Return a const view of some rows of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
auto rows(const Matrix& mat, const Indices& irows) -> MatrixViewRowsConst;

/// Return a view of a sequence of columns of a matrix
/// @param start The column index of the start of the sequence
/// @param num The number of columns in the sequence
template<typename MatrixType>
auto cols(MatrixType&& mat, unsigned start, unsigned num) -> decltype(mat.middleCols(start, num));

/// Return a view of some columns of a matrix
/// @param mat The matrix for which the view is created
/// @param icols The indices of the columns of the matrix
auto cols(Matrix& mat, const Indices& icols) -> MatrixViewCols;

/// Return a const view of some columns of a matrix
/// @param mat The matrix for which the view is created
/// @param icols The indices of the columns of the matrix
auto cols(const Matrix& mat, const Indices& icols) -> MatrixViewColsConst;

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
auto submatrix(Matrix& mat, const Indices& irows, const Indices& icols) -> MatrixViewRowsCols;

/// Return a const view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
auto submatrix(const Matrix& mat, const Indices& irows, const Indices& icols) -> MatrixViewRowsColsConst;

} // namespace Reaktor
