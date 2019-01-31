// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

#define EIGEN_MATRIX_PLUGIN "Reaktoro/Math/EigenMatrixPlugin.hpp"

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/Core>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>

namespace Reaktoro {

using Vector           = Eigen::VectorXd;                             /// < Alias to Eigen type VectorXd.
using VectorRef        = Eigen::Ref<Eigen::VectorXd>;                 /// < Alias to Eigen type Ref<VectorXd>.
using VectorStridedRef = Eigen::Ref<Vector, 0, Eigen::InnerStride<>>; /// < Alias to Eigen type Ref<VectorXd>.
using VectorConstRef   = Eigen::Ref<const Eigen::VectorXd>;           /// < Alias to Eigen type Ref<const VectorXd>.
using VectorMap        = Eigen::Map<Eigen::VectorXd>;                 /// < Alias to Eigen type Map<VectorXd>.
using VectorConstMap   = Eigen::Map<const Eigen::VectorXd>;           /// < Alias to Eigen type Map<const VectorXd>.

using RowVector           = Eigen::RowVectorXd;                                      /// < Alias to Eigen type RowVectorXd.
using RowVectorRef        = Eigen::Ref<Eigen::RowVectorXd>;                          /// < Alias to Eigen type Ref<RowVectorXd>.
using RowVectorStridedRef = Eigen::Ref<Eigen::RowVectorXd, 0, Eigen::InnerStride<>>; /// < Alias to Eigen type Ref<RowVectorXd>.
using RowVectorConstRef   = Eigen::Ref<const Eigen::RowVectorXd>;                    /// < Alias to Eigen type Ref<const RowVectorXd>.
using RowVectorMap        = Eigen::Map<Eigen::RowVectorXd>;                          /// < Alias to Eigen type Map<VectorXd>.
using RowVectorConstMap   = Eigen::Map<const Eigen::RowVectorXd>;                    /// < Alias to Eigen type Map<const VectorXd>.

using Matrix         = Eigen::MatrixXd;                   ///< Alias to Eigen type MatrixXd.
using MatrixRef      = Eigen::Ref<Eigen::MatrixXd>;       ///< Alias to Eigen type Ref<MatrixXd>.
using MatrixConstRef = Eigen::Ref<const Eigen::MatrixXd>; ///< Alias to Eigen type Ref<const MatrixXd>.
using MatrixMap      = Eigen::Map<Eigen::MatrixXd>;       ///< Alias to Eigen type Map<MatrixXd>.
using MatrixConstMap = Eigen::Map<const Eigen::MatrixXd>; ///< Alias to Eigen type Map<const MatrixXd>.

using Vector = Eigen::VectorXd; /// Alias to Eigen type Eigen::VectorXd.
using VectorXd = Eigen::VectorXd; /// Alias to Eigen type Eigen::VectorXd.
using VectorXi = Eigen::VectorXi; /// Alias to Eigen type Eigen::VectorXi.

using VectorRef = Eigen::Ref<VectorXd>; ///< Alias to Eigen type Eigen::Ref<VectorXd>.
using VectorXdRef = Eigen::Ref<VectorXd>; ///< Alias to Eigen type Eigen::Ref<VectorXd>.
using VectorXiRef = Eigen::Ref<VectorXi>; ///< Alias to Eigen type Eigen::Ref<VectorXi>.

using VectorConstRef = Eigen::Ref<const VectorXd>; ///< Alias to Eigen type Eigen::Ref<const VectorXd>.
using VectorXdConstRef = Eigen::Ref<const VectorXd>; ///< Alias to Eigen type Eigen::Ref<const VectorXd>.
using VectorXiConstRef = Eigen::Ref<const VectorXi>; ///< Alias to Eigen type Eigen::Ref<const VectorXi>.

using VectorMap = Eigen::Map<VectorXd>; ///< Alias to Eigen type Eigen::Map<VectorXd>.
using VectorXdMap = Eigen::Map<VectorXd>; ///< Alias to Eigen type Eigen::Map<VectorXd>.
using VectorXiMap = Eigen::Map<VectorXi>; ///< Alias to Eigen type Eigen::Map<VectorXi>.

using VectorConstMap = Eigen::Map<const VectorXd>; ///< Alias to Eigen type Eigen::Map<const VectorXd>.
using VectorXdConstMap = Eigen::Map<const VectorXd>; ///< Alias to Eigen type Eigen::Map<const VectorXd>.
using VectorXiConstMap = Eigen::Map<const VectorXi>; ///< Alias to Eigen type Eigen::Map<const VectorXi>.

using Matrix = Eigen::MatrixXd; /// Alias to Eigen type Eigen::MatrixXd.
using MatrixXd = Eigen::MatrixXd; /// Alias to Eigen type Eigen::MatrixXd.
using MatrixXi = Eigen::MatrixXi; /// Alias to Eigen type Eigen::MatrixXi.

using MatrixRef = Eigen::Ref<MatrixXd>; ///< Alias to Eigen type Eigen::Ref<MatrixXd>.
using MatrixXdRef = Eigen::Ref<MatrixXd>; ///< Alias to Eigen type Eigen::Ref<MatrixXd>.
using MatrixXiRef = Eigen::Ref<MatrixXi>; ///< Alias to Eigen type Eigen::Ref<MatrixXi>.

using MatrixConstRef = Eigen::Ref<const MatrixXd>; ///< Alias to Eigen type Eigen::Ref<const MatrixXd>.
using MatrixXdConstRef = Eigen::Ref<const MatrixXd>; ///< Alias to Eigen type Eigen::Ref<const MatrixXd>.
using MatrixXiConstRef = Eigen::Ref<const MatrixXi>; ///< Alias to Eigen type Eigen::Ref<const MatrixXi>.

using MatrixMap = Eigen::Map<MatrixXd>; ///< Alias to Eigen type Eigen::Map<MatrixXd>.
using MatrixXdMap = Eigen::Map<MatrixXd>; ///< Alias to Eigen type Eigen::Map<MatrixXd>.
using MatrixXiMap = Eigen::Map<MatrixXi>; ///< Alias to Eigen type Eigen::Map<MatrixXi>.

using MatrixConstMap = Eigen::Map<const MatrixXd>; ///< Alias to Eigen type Eigen::Map<const MatrixXd>.
using MatrixXdConstMap = Eigen::Map<const MatrixXd>; ///< Alias to Eigen type Eigen::Map<const MatrixXd>.
using MatrixXiConstMap = Eigen::Map<const MatrixXi>; ///< Alias to Eigen type Eigen::Map<const MatrixXi>.

/// Define an alias to a permutation matrix type of the Eigen library
using PermutationMatrix = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>;

/// Return an expression of a zero vector
/// @param rows The number of rows
/// @return The expression of a zero vector
auto zeros(Index rows) -> decltype(Vector::Zero(rows));

/// Return an expression of a vector with entries equal to one
/// @param rows The number of rows
/// @return The expression of a vector with entries equal to one
auto ones(Index rows) -> decltype(Vector::Ones(rows));

/// Return an expression of a vector with random entries
/// @param rows The number of rows
/// @return The expression of a vector with random entries equal to one
auto random(Index rows) -> decltype(Vector::Random(rows));

/// Return a linearly spaced vector
/// @param rows The number of rows
/// @param start The start of the sequence
/// @param stop The stop of the sequence
/// @return The expression of a vector with linearly spaced entries
auto linspace(Index rows, double start, double stop) -> decltype(Vector::LinSpaced(rows, start, stop));

/// Return an expression of a unit vector
/// @param rows The number of rows
/// @param i The index at which the component is one
/// @return The expression of a unit vector
auto unit(Index rows, Index i) -> decltype(Vector::Unit(rows, i));

/// Return an expression of a zero matrix
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a zero matrix
auto zeros(Index rows, Index cols) -> decltype(Matrix::Zero(rows, cols));

/// Return an expression of a matrix with entries equal to one
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a matrix with entries equal to one
auto ones(Index rows, Index cols) -> decltype(Matrix::Ones(rows, cols));

/// Return an expression of a matrix with random entries
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a matrix with random entries
auto random(Index rows, Index cols) -> decltype(Matrix::Random(rows, cols));

/// Return an expression of an identity matrix
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of an identity matrix
auto identity(Index rows, Index cols) -> decltype(Matrix::Identity(rows, cols));

/// Return a view of a sequence of rows of a matrix
/// @param start The row index of the start of the sequence
/// @param num The number of rows in the sequence
template<typename Derived>
auto rows(Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleRows(start, num));

/// Return a view of a sequence of rows of a matrix
/// @param start The row index of the start of the sequence
/// @param num The number of rows in the sequence
template<typename Derived>
auto rows(const Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleRows(start, num));

/// Return a view of some rows of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
template<typename Derived, typename Indices>
auto rows(Eigen::MatrixBase<Derived>& mat, const Indices& irows) -> decltype(mat(irows, Eigen::all));

/// Return a const view of some rows of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
template<typename Derived, typename Indices>
auto rows(const Eigen::MatrixBase<Derived>& mat, const Indices& irows) -> decltype(mat(irows, Eigen::all));

/// Return a view of a sequence of columns of a matrix
/// @param start The column index of the start of the sequence
/// @param num The number of columns in the sequence
template<typename Derived>
auto cols(Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleCols(start, num));

/// Return a view of a sequence of columns of a matrix
/// @param start The column index of the start of the sequence
/// @param num The number of columns in the sequence
template<typename Derived>
auto cols(const Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleCols(start, num));

/// Return a view of some columns of a matrix
/// @param mat The matrix for which the view is created
/// @param icols The indices of the columns of the matrix
template<typename Derived, typename Indices>
auto cols(Eigen::MatrixBase<Derived>& mat, const Indices& icols) -> decltype(mat(Eigen::all, icols));

/// Return a const view of some columns of a matrix
/// @param mat The matrix for which the view is created
/// @param icols The indices of the columns of the matrix
template<typename Derived, typename Indices>
auto cols(const Eigen::MatrixBase<Derived>& mat, const Indices& icols) -> decltype(mat(Eigen::all, icols));

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto segment(Eigen::MatrixBase<Derived>& vec, Index irow, Index nrows) -> decltype(vec.segment(irow, nrows));

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto segment(const Eigen::MatrixBase<Derived>& vec, Index irow, Index nrows) -> decltype(vec.segment(irow, nrows));

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto block(Eigen::MatrixBase<Derived>& mat, Index irow, Index icol, Index nrows, Index ncols) -> decltype(mat.block(irow, icol, nrows, ncols));

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived>
auto block(const Eigen::MatrixBase<Derived>& mat, Index irow, Index icol, Index nrows, Index ncols) -> decltype(mat.block(irow, icol, nrows, ncols));

/// Return a view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived, typename Indices>
auto submatrix(Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols) -> decltype(mat(irows, icols));

/// Return a const view of some rows and columns of a matrix
/// @param mat The matrix for which the view is created
/// @param irows The indices of the rows of the matrix
/// @param icols The indices of the columns of the matrix
template<typename Derived, typename Indices>
auto submatrix(const Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols) -> decltype(mat(irows, icols));

/// Return a block mapped view of a matrix.
/// @param mat The matrix from which the mapped view is created.
/// @param row The index of the row at which the view starts.
/// @param col The index of the column at which the view starts.
/// @param nrows The number of rows of the mapped view.
/// @param ncols The number of columns of the mapped view.
template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
auto blockmap(Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& mat, Index row, Index col, Index nrows, Index ncols) -> Eigen::Map<Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>, Eigen::Unaligned, Eigen::Stride<Rows,Cols>>
{
    Eigen::Stride<Rows,Cols> stride(mat.outerStride(), mat.innerStride());
    return Eigen::Map<Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>, Eigen::Unaligned, Eigen::Stride<MaxRows,MaxCols>>(
        mat.block(row, col, nrows, ncols).data(), nrows, ncols, stride);
}

/// Return a const block mapped view of a matrix.
/// @param mat The matrix from which the mapped view is created.
/// @param row The index of the row at which the view starts.
/// @param col The index of the column at which the view starts.
/// @param nrows The number of rows of the mapped view.
/// @param ncols The number of columns of the mapped view.
template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
auto blockmap(const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& mat, Index row, Index col, Index nrows, Index ncols) -> Eigen::Map<const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>, Eigen::Unaligned, Eigen::Stride<Rows,Cols>>
{
    Eigen::Stride<Rows,Cols> stride(mat.outerStride(), mat.innerStride());
    return Eigen::Map<const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>, Eigen::Unaligned, Eigen::Stride<MaxRows,MaxCols>>(
        mat.block(row, col, nrows, ncols).data(), nrows, ncols, stride);
}

/// Return a mapped view of a sequence of rows of a matrix.
/// @param mat The matrix from which the mapped view is created.
/// @param row The index of the row at which the view starts.
/// @param nrows The number of rows of the mapped view.
template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
auto rowsmap(Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& mat, Index row, Index nrows) -> Eigen::Map<Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>, Eigen::Unaligned, Eigen::Stride<Rows,Cols>>
{
    return blockmap(mat, row, 0, nrows, mat.cols());
}

/// Return a const mapped view of a sequence of rows of a matrix.
/// @param mat The matrix from which the mapped view is created.
/// @param row The index of the row at which the view starts.
/// @param nrows The number of rows of the mapped view.
template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
auto rowsmap(const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& mat, Index row, Index nrows) -> Eigen::Map<const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>, Eigen::Unaligned, Eigen::Stride<Rows,Cols>>
{
    return blockmap(mat, row, 0, nrows, mat.cols());
}

/// Return a mapped view of a sequence of columns of a matrix.
/// @param mat The matrix from which the mapped view is created.
/// @param row The index of the row at which the view starts.
/// @param col The index of the column at which the view starts.
/// @param nrows The number of rows of the mapped view.
/// @param ncols The number of columns of the mapped view.
template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
auto colsmap(Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& mat, Index col, Index ncols) -> Eigen::Map<Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>, Eigen::Unaligned, Eigen::Stride<Rows,Cols>>
{
    return blockmap(mat, 0, col, mat.rows(), ncols);
}

/// Return a const mapped view of a sequence of columns of a matrix.
/// @param mat The matrix from which the mapped view is created.
/// @param row The index of the row at which the view starts.
/// @param col The index of the column at which the view starts.
/// @param nrows The number of rows of the mapped view.
/// @param ncols The number of columns of the mapped view.
template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
auto colsmap(const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>& mat, Index col, Index ncols) -> Eigen::Map<const Eigen::Matrix<Scalar,Rows,Cols,Options,MaxRows,MaxCols>, Eigen::Unaligned, Eigen::Stride<Rows,Cols>>
{
    return blockmap(mat, 0, col, mat.rows(), ncols);
}

/// Return the transpose of the matrix
template<typename Derived>
auto tr(Eigen::MatrixBase<Derived>& mat) -> decltype(mat.transpose());

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
template<typename Derived>
auto diagonal(Eigen::MatrixBase<Derived>& mat) -> decltype(mat.diagonal());

/// Return a vector representation of the diagonal of a matrix
template<typename Derived>
auto diagonal(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.diagonal());

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
auto exp(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().exp().matrix());

/// Return the component-wise natural log of a matrix
template<typename Derived>
auto log(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().log().matrix());

/// Return the component-wise log10 of a matrix
template<typename Derived>
auto log10(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().log10().matrix());

} // namespace Reaktoro

#include "Matrix.hxx"
