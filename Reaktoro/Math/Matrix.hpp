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

// Eigen includes
#include <Reaktoro/deps/eigen3/Eigen/Core>

// Reaktoro includes
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/Real.hpp>

namespace Reaktoro {

///< Convenient alias to Eigen type.
template<typename T>
using VectorX = Eigen::Matrix<T, -1, 1, 0, -1, 1>;

///< Convenient alias to Eigen type.
template<typename T>
using MatrixX = Eigen::Matrix<T, -1, -1, 0, -1, -1>;

///< Convenient alias to Eigen type.
template<typename T>
using ArrayX = Eigen::Array<T, -1, 1,  0,  -1, 1>;

///< Convenient alias to Eigen type.
template<typename T>
using ArrayXX = Eigen::Array<T, -1, -1,  0,  -1, -1>;

///< Convenient alias to Eigen type.
template<typename T>
using RowVectorX = Eigen::Matrix<T, 1, -1, 0, 1, -1>;

//---------------------------------------------------------------------------------------------------------------------
// == VECTOR TYPE ALIASES ==
//---------------------------------------------------------------------------------------------------------------------

using VectorXr                = VectorX<real>;                                       ///< Convenient alias to Eigen type.
using VectorXrRef             = Eigen::Ref<VectorXr>;                                ///< Convenient alias to Eigen type.
using VectorXrConstRef        = Eigen::Ref<const VectorXr>;                          ///< Convenient alias to Eigen type.
using VectorXrMap             = Eigen::Map<VectorXr>;                                ///< Convenient alias to Eigen type.
using VectorXrConstMap        = Eigen::Map<const VectorXr>;                          ///< Convenient alias to Eigen type.
using VectorXrStridedRef      = Eigen::Ref<VectorXr, 0, Eigen::InnerStride<>>;       ///< Convenient alias to Eigen type.
using VectorXrStridedConstRef = Eigen::Ref<const VectorXr, 0, Eigen::InnerStride<>>; ///< Convenient alias to Eigen type.

using VectorXl                = VectorX<long>;                                       ///< Convenient alias to Eigen type.
using VectorXlRef             = Eigen::Ref<VectorXl>;                                ///< Convenient alias to Eigen type.
using VectorXlConstRef        = Eigen::Ref<const VectorXl>;                          ///< Convenient alias to Eigen type.
using VectorXlMap             = Eigen::Map<VectorXl>;                                ///< Convenient alias to Eigen type.
using VectorXlConstMap        = Eigen::Map<const VectorXl>;                          ///< Convenient alias to Eigen type.
using VectorXlStridedRef      = Eigen::Ref<VectorXl, 0, Eigen::InnerStride<>>;       ///< Convenient alias to Eigen type.
using VectorXlStridedConstRef = Eigen::Ref<const VectorXl, 0, Eigen::InnerStride<>>; ///< Convenient alias to Eigen type.

using VectorXd                = Eigen::VectorXd;                                     ///< Convenient alias to Eigen type.
using VectorXdRef             = Eigen::Ref<VectorXd>;                                ///< Convenient alias to Eigen type.
using VectorXdConstRef        = Eigen::Ref<const VectorXd>;                          ///< Convenient alias to Eigen type.
using VectorXdMap             = Eigen::Map<VectorXd>;                                ///< Convenient alias to Eigen type.
using VectorXdConstMap        = Eigen::Map<const VectorXd>;                          ///< Convenient alias to Eigen type.
using VectorXdStridedRef      = Eigen::Ref<VectorXd, 0, Eigen::InnerStride<>>;       ///< Convenient alias to Eigen type.
using VectorXdStridedConstRef = Eigen::Ref<const VectorXd, 0, Eigen::InnerStride<>>; ///< Convenient alias to Eigen type.


//---------------------------------------------------------------------------------------------------------------------
// == ARRAY TYPE ALIASES ==
//---------------------------------------------------------------------------------------------------------------------

using ArrayXr                = ArrayX<real>;                                       ///< Convenient alias to Eigen type.
using ArrayXrRef             = Eigen::Ref<ArrayXr>;                                ///< Convenient alias to Eigen type.
using ArrayXrConstRef        = Eigen::Ref<const ArrayXr>;                          ///< Convenient alias to Eigen type.
using ArrayXrMap             = Eigen::Map<ArrayXr>;                                ///< Convenient alias to Eigen type.
using ArrayXrConstMap        = Eigen::Map<const ArrayXr>;                          ///< Convenient alias to Eigen type.
using ArrayXrStridedRef      = Eigen::Ref<ArrayXr, 0, Eigen::InnerStride<>>;       ///< Convenient alias to Eigen type.
using ArrayXrStridedConstRef = Eigen::Ref<const ArrayXr, 0, Eigen::InnerStride<>>; ///< Convenient alias to Eigen type.

using ArrayXl                = ArrayX<long>;                                       ///< Convenient alias to Eigen type.
using ArrayXlRef             = Eigen::Ref<ArrayXl>;                                ///< Convenient alias to Eigen type.
using ArrayXlConstRef        = Eigen::Ref<const ArrayXl>;                          ///< Convenient alias to Eigen type.
using ArrayXlMap             = Eigen::Map<ArrayXl>;                                ///< Convenient alias to Eigen type.
using ArrayXlConstMap        = Eigen::Map<const ArrayXl>;                          ///< Convenient alias to Eigen type.
using ArrayXlStridedRef      = Eigen::Ref<ArrayXl, 0, Eigen::InnerStride<>>;       ///< Convenient alias to Eigen type.
using ArrayXlStridedConstRef = Eigen::Ref<const ArrayXl, 0, Eigen::InnerStride<>>; ///< Convenient alias to Eigen type.

using ArrayXd                = Eigen::ArrayXd;                                     ///< Convenient alias to Eigen type.
using ArrayXdRef             = Eigen::Ref<ArrayXd>;                                ///< Convenient alias to Eigen type.
using ArrayXdConstRef        = Eigen::Ref<const ArrayXd>;                          ///< Convenient alias to Eigen type.
using ArrayXdMap             = Eigen::Map<ArrayXd>;                                ///< Convenient alias to Eigen type.
using ArrayXdConstMap        = Eigen::Map<const ArrayXd>;                          ///< Convenient alias to Eigen type.
using ArrayXdStridedRef      = Eigen::Ref<ArrayXd, 0, Eigen::InnerStride<>>;       ///< Convenient alias to Eigen type.
using ArrayXdStridedConstRef = Eigen::Ref<const ArrayXd, 0, Eigen::InnerStride<>>; ///< Convenient alias to Eigen type.

//---------------------------------------------------------------------------------------------------------------------
// == MATRIX TYPE ALIASES ==
//---------------------------------------------------------------------------------------------------------------------

using MatrixXr                = MatrixX<real>;              ///< Convenient alias to Eigen type.
using MatrixXrRef             = Eigen::Ref<MatrixXr>;       ///< Convenient alias to Eigen type.
using MatrixXrConstRef        = Eigen::Ref<const MatrixXr>; ///< Convenient alias to Eigen type.
using MatrixXrMap             = Eigen::Map<MatrixXr>;       ///< Convenient alias to Eigen type.
using MatrixXrConstMap        = Eigen::Map<const MatrixXr>; ///< Convenient alias to Eigen type.

using MatrixXd                = Eigen::MatrixXd;            ///< Convenient alias to Eigen type.
using MatrixXdRef             = Eigen::Ref<MatrixXd>;       ///< Convenient alias to Eigen type.
using MatrixXdConstRef        = Eigen::Ref<const MatrixXd>; ///< Convenient alias to Eigen type.
using MatrixXdMap             = Eigen::Map<MatrixXd>;       ///< Convenient alias to Eigen type.
using MatrixXdConstMap        = Eigen::Map<const MatrixXd>; ///< Convenient alias to Eigen type.

//---------------------------------------------------------------------------------------------------------------------
// == ROW VECTOR TYPE ALIASES ==
//---------------------------------------------------------------------------------------------------------------------

using RowVectorXr                = RowVectorX<real>;              ///< Convenient alias to Eigen type.
using RowVectorXrRef             = Eigen::Ref<RowVectorXr>;       ///< Convenient alias to Eigen type.
using RowVectorXrConstRef        = Eigen::Ref<const RowVectorXr>; ///< Convenient alias to Eigen type.
using RowVectorXrMap             = Eigen::Map<RowVectorXr>;       ///< Convenient alias to Eigen type.
using RowVectorXrConstMap        = Eigen::Map<const RowVectorXr>; ///< Convenient alias to Eigen type.

using RowVectorXl                = RowVectorX<long>;              ///< Convenient alias to Eigen type.
using RowVectorXlRef             = Eigen::Ref<RowVectorXl>;       ///< Convenient alias to Eigen type.
using RowVectorXlConstRef        = Eigen::Ref<const RowVectorXl>; ///< Convenient alias to Eigen type.
using RowVectorXlMap             = Eigen::Map<RowVectorXl>;       ///< Convenient alias to Eigen type.
using RowVectorXlConstMap        = Eigen::Map<const RowVectorXl>; ///< Convenient alias to Eigen type.

using RowVectorXd                = Eigen::RowVectorXd;            ///< Convenient alias to Eigen type.
using RowVectorXdRef             = Eigen::Ref<RowVectorXd>;       ///< Convenient alias to Eigen type.
using RowVectorXdConstRef        = Eigen::Ref<const RowVectorXd>; ///< Convenient alias to Eigen type.
using RowVectorXdMap             = Eigen::Map<RowVectorXd>;       ///< Convenient alias to Eigen type.
using RowVectorXdConstMap        = Eigen::Map<const RowVectorXd>; ///< Convenient alias to Eigen type.

//---------------------------------------------------------------------------------------------------------------------
// == OTHER TYPE ALIASES ==
//---------------------------------------------------------------------------------------------------------------------

/// Define an alias to a permutation matrix type of the Eigen library
using PermutationMatrix = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>;

//---------------------------------------------------------------------------------------------------------------------
// == FUNCTION ALIASES ==
//---------------------------------------------------------------------------------------------------------------------

/// Return an expression of a zero vector
/// @param rows The number of rows
/// @return The expression of a zero vector
auto zeros(Index rows) -> decltype(VectorXd::Zero(rows));

/// Return an expression of a vector with entries equal to one
/// @param rows The number of rows
/// @return The expression of a vector with entries equal to one
auto ones(Index rows) -> decltype(VectorXd::Ones(rows));

/// Return an expression of a vector with random entries
/// @param rows The number of rows
/// @return The expression of a vector with random entries equal to one
auto random(Index rows) -> decltype(VectorXd::Random(rows));

/// Return a linearly spaced vector
/// @param rows The number of rows
/// @param start The start of the sequence
/// @param stop The stop of the sequence
/// @return The expression of a vector with linearly spaced entries
auto linspace(Index rows, double start, double stop) -> decltype(VectorXd::LinSpaced(rows, start, stop));

/// Return an expression of a unit vector
/// @param rows The number of rows
/// @param i The index at which the component is one
/// @return The expression of a unit vector
auto unit(Index rows, Index i) -> decltype(VectorXd::Unit(rows, i));

/// Return an expression of a zero matrix
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a zero matrix
auto zeros(Index rows, Index cols) -> decltype(MatrixXd::Zero(rows, cols));

/// Return an expression of a matrix with entries equal to one
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a matrix with entries equal to one
auto ones(Index rows, Index cols) -> decltype(MatrixXd::Ones(rows, cols));

/// Return an expression of a matrix with random entries
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of a matrix with random entries
auto random(Index rows, Index cols) -> decltype(MatrixXd::Random(rows, cols));

/// Return an expression of an identity matrix
/// @param rows The number of rows
/// @param cols The number of columns
/// @return The expression of an identity matrix
auto identity(Index rows, Index cols) -> decltype(MatrixXd::Identity(rows, cols));

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
