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

// Eigen includes
#include <Reaktoro/Eigen/Core>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>

namespace Reaktoro {

/// Define an alias to the vector type of the Eigen library
using Vector = Eigen::VectorXd;

/// Define an alias to the matrix type of the Eigen library
using Matrix = Eigen::MatrixXd;

template<typename Derived> class MatrixViewRows;
template<typename Derived> class MatrixViewRowsConst;

/// A type used to define a view of some rows of a matrix instance
template<typename Derived>
class MatrixViewRows
{
public:
    /// Construct a MatrixViewRows instance
    /// @param mat The matrix for which this view is defined
    /// @param irows The indices of the rows of this matrix view
    MatrixViewRows(Eigen::MatrixBase<Derived>& mat, const Indices& irows);

    /// Assign a MatrixViewRows instance to this MatrixViewRows instance
    auto operator=(const MatrixViewRows& other) -> MatrixViewRows&;

    /// Assign a MatrixViewRowsConst instance to this MatrixViewRows instance
    template<typename DerivedOther>
    auto operator=(const MatrixViewRowsConst<DerivedOther>& other) -> MatrixViewRows&;

    /// Assign a matrix instance to this MatrixViewRows instance
    template<typename DerivedOther>
    auto operator=(const Eigen::MatrixBase<DerivedOther>& other) -> MatrixViewRows&;

    /// Assign a MatrixViewRows instance to a scalar
    auto operator=(const typename Derived::Scalar& scalar) -> MatrixViewRows&;

    /// Transfer the data under this view to another matrix
    template<typename DerivedOther>
    auto to(Eigen::MatrixBase<DerivedOther>& other) const -> void;

    /// Convert this MatrixViewRows instance into a matrix instance
    operator Derived() const;

    /// The matrix for which this view is defined
    Eigen::MatrixBase<Derived>& mat;

    /// The indices of the rows of this matrix view
    const Indices& irows;
};

/// A type used to define a const view of some rows of a matrix instance
template<typename Derived>
class MatrixViewRowsConst
{
public:
    /// Construct a MatrixViewConst instance
    /// @param mat The matrix for which this view is defined
    /// @param irows The indices of the rows of this matrix view
    MatrixViewRowsConst(const Eigen::MatrixBase<Derived>& mat, const Indices& irows);

    /// Transfer the data under this view to another matrix
    template<typename DerivedOther>
    auto to(Eigen::MatrixBase<DerivedOther>& other) const -> void;

    /// Convert this MatrixViewConst instance into a matrix instance
    operator Derived() const;

    /// The matrix for which this view is defined
    const Eigen::MatrixBase<Derived>& mat;

    /// The indices of the rows of this matrix view
    const Indices& irows;
};

/// A type used to define a view of some columns of a matrix instance
template<typename Derived>
class MatrixViewCols
{
public:
    /// Construct a MatrixViewCols instance
    /// @param mat The matrix for which this view is defined
    /// @param icols The indices of the columns of this matrix view
    MatrixViewCols(Eigen::MatrixBase<Derived>& mat, const Indices& icols);

    /// Assign a matrix instance to this MatrixViewCols instance
    template<typename DerivedOther>
    auto operator=(const Eigen::MatrixBase<DerivedOther>& other) -> MatrixViewCols&;

    /// Transfer the data under this view to another matrix
    template<typename DerivedOther>
    auto to(Eigen::MatrixBase<DerivedOther>& other) const -> void;

    /// Convert this MatrixViewCols instance into a matrix instance
    operator Derived() const;

    /// The matrix for which this view is defined
    Eigen::MatrixBase<Derived>& mat;

    /// The indices of the columns of this matrix view
    const Indices& icols;
};

/// A type used to define a const view of some columns of a matrix instance
template<typename Derived>
class MatrixViewColsConst
{
public:
    /// Construct a MatrixViewColsConst instance
    /// @param mat The matrix for which this view is defined
    /// @param icols The indices of the columns of this matrix view
    MatrixViewColsConst(const Eigen::MatrixBase<Derived>& mat, const Indices& icols);

    /// Transfer the data under this view to another matrix
    template<typename DerivedOther>
    auto to(Eigen::MatrixBase<DerivedOther>& other) const -> void;

    /// Convert this MatrixViewColsConst instance into a matrix instance
    operator Derived() const;

    /// The matrix for which this view is defined
    const Eigen::MatrixBase<Derived>& mat;

    /// The indices of the columns of this matrix view
    const Indices& icols;
};

/// A type used to define a view of some rows and columns of a matrix instance
template<typename Derived>
class MatrixViewRowsCols
{
public:
    /// Construct a MatrixViewRowsCols instance
    /// @param mat The matrix for which this view is defined
    /// @param irows The indices of the rows of this matrix view
    /// @param icols The indices of the columns of this matrix view
    MatrixViewRowsCols(Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols);

    /// Assign a matrix instance to this MatrixViewCols instance
    template<typename DerivedOther>
    auto operator=(const Eigen::MatrixBase<DerivedOther>& other) -> MatrixViewRowsCols&;

    /// Transfer the data under this view to another matrix
    template<typename DerivedOther>
    auto to(Eigen::MatrixBase<DerivedOther>& other) const -> void;

    /// Convert this MatrixViewRowsCols instance into a matrix instance
    operator Derived() const;

    /// The matrix for which this view is defined
    Eigen::MatrixBase<Derived>& mat;

    /// The indices of the rows of this matrix view
    const Indices& irows;

    /// The indices of the columns of this matrix view
    const Indices& icols;
};

/// A type used to define a const view of some rows and columns of a matrix instance
template<typename Derived>
class MatrixViewRowsColsConst
{
public:
    /// Construct a MatrixViewRowsColsConst instance
    /// @param mat The matrix for which this view is defined
    /// @param irows The indices of the rows of this matrix view
    /// @param icols The indices of the columns of this matrix view
    MatrixViewRowsColsConst(const Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols);

    /// Transfer the data under this view to another matrix
    template<typename DerivedOther>
    auto to(Eigen::MatrixBase<DerivedOther>& other) const -> void;

    /// Convert this MatrixViewRowsColsConst instance into a matrix instance
    operator Derived() const;

    /// The matrix for which this view is defined
    const Eigen::MatrixBase<Derived>& mat;

    /// The indices of the rows of this matrix view
    const Indices& irows;

    /// The indices of the columns of this matrix view
    const Indices& icols;
};

/// Return an expression of a zero vector
/// @param rows The number of rows
/// @return The expression of a zero vector
auto zeros(unsigned rows) -> decltype(Vector::Zero(rows));

/// Return an expression of a vector with entries equal to one
/// @param rows The number of rows
/// @return The expression of a vector with entries equal to one
auto ones(unsigned rows) -> decltype(Vector::Ones(rows));

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

/// Return an expression of an identity matrix
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of an identity matrix
auto identity(unsigned rows, unsigned cols) -> decltype(Matrix::Identity(rows, cols));

/// Return a view of a sequence of rows of a matrix
/// @param start The row index of the start of the sequence
/// @param num The number of rows in the sequence
template<typename MatrixType>
auto rows(MatrixType&& mat, unsigned start, unsigned num) -> decltype(mat.middleRows(start, num));

/// Return a view of some rows of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
template<typename Derived>
auto rows(Eigen::MatrixBase<Derived>& mat, const Indices& irows) -> MatrixViewRows<Derived>;

/// Return a const view of some rows of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
template<typename Derived>
auto rows(const Eigen::MatrixBase<Derived>& mat, const Indices& irows) -> MatrixViewRowsConst<Derived>;

/// Return a view of a sequence of columns of a matrix
/// @param start The column index of the start of the sequence
/// @param num The number of columns in the sequence
template<typename MatrixType>
auto cols(MatrixType&& mat, unsigned start, unsigned num) -> decltype(mat.middleCols(start, num));

/// Return a view of some columns of a matrix
/// @param mat The matrix for which the view is created
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto cols(Eigen::MatrixBase<Derived>& mat, const Indices& icols) -> MatrixViewCols<Derived>;

/// Return a const view of some columns of a matrix
/// @param mat The matrix for which the view is created
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto cols(const Eigen::MatrixBase<Derived>& mat, const Indices& icols) -> MatrixViewColsConst<Derived>;

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename MatrixType>
auto block(MatrixType&& mat, unsigned irow, unsigned icol, unsigned nrows, unsigned ncols) -> decltype(mat.block(irow, icol, nrows, ncols));

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto submatrix(Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols) -> MatrixViewRowsCols<Derived>;

/// Return a const view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto submatrix(const Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols) -> MatrixViewRowsColsConst<Derived>;

/// Return the transpose of the matrix
template<typename Derived>
auto tr(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.transpose());

/// Return the component-wise inverse of the matrix
template<typename Derived>
auto inv(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseInverse());

/// Return a diagonal matrix representation of a vector
template<typename Derived>
auto diag(const Eigen::MatrixBase<Derived>& vec) -> decltype(vec.asDiagonal());

/// Return a vector representation of the diagonal of a matrix
template<typename MatrixType>
auto diagonal(MatrixType&& mat) -> decltype(mat.diagonal());

/// Return the Lp norm of a matrix
template<int p, typename Derived>
auto norm(const Eigen::MatrixBase<Derived>& mat) -> double;

/// Return the L2 norm of a matrix
template<typename Derived>
auto norm(const Eigen::MatrixBase<Derived>& mat) -> double;

/// Return the L-inf norm of a matrix
template<typename Derived>
auto norminf(const Eigen::MatrixBase<Derived>& mat) -> double;

/// Return the sum of the components of a matrix
template<typename Derived>
auto sum(const Eigen::DenseBase<Derived>& mat) -> typename Derived::Scalar;

/// Return the dot product of two matrices
template<typename DerivedLHS, typename DerivedRHS>
auto dot(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.dot(rhs));

/// Return the minimum component of a matrix
template<typename Derived>
auto min(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.minCoeff());

/// Return the component-wise minimum of two matrices
template<typename DerivedLHS, typename DerivedRHS>
auto min(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseMin(rhs));

/// Return the maximum component of a matrix
template<typename Derived>
auto max(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.maxCoeff());

/// Return the component-wise maximum of two matrices
template<typename DerivedLHS, typename DerivedRHS>
auto max(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseMax(rhs));

/// Return the component-wise absolute entries of a matrix
template<typename Derived>
auto abs(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseAbs());

/// Return the component-wise square root of a matrix
template<typename Derived>
auto sqrt(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseSqrt());

/// Return the component-wise exponential of a matrix
template<typename Derived>
auto pow(const Eigen::MatrixBase<Derived>& mat, double power) -> decltype(mat.array().pow(power));

/// Return the component-wise natural exponential of a matrix
template<typename Derived>
auto exp(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().exp());

/// Return the component-wise natural log of a matrix
template<typename Derived>
auto log(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().log());

/// Return the component-wise multiplication of two matrices
template<typename DerivedLHS, typename DerivedRHS>
auto operator%(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseProduct(rhs));

/// Return the component-wise division of two matrices
template<typename DerivedLHS, typename DerivedRHS>
auto operator/(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseQuotient(rhs));

/// Return the component-wise division of two matrices
template<typename Derived>
auto operator/(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype(scalar*mat.cwiseInverse());

/// Return the component-wise division of two matrices
template<typename Derived>
auto operator+(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype((scalar + mat.array()).matrix());

/// Return the component-wise division of two matrices
template<typename Derived>
auto operator+(const Eigen::MatrixBase<Derived>& mat, const typename Derived::Scalar& scalar) -> decltype((scalar + mat.array()).matrix());

/// Return the component-wise division of two matrices
template<typename Derived>
auto operator-(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype((scalar - mat.array()).matrix());

/// Return the component-wise division of two matrices
template<typename Derived>
auto operator-(const Eigen::MatrixBase<Derived>& mat, const typename Derived::Scalar& scalar) -> decltype((mat.array() - scalar).matrix());

} // namespace Reaktoro

#include "Matrix.hxx"
