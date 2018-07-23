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

namespace Reaktoro {

inline auto zeros(Index rows) -> decltype(Vector::Zero(rows))
{
    return Vector::Zero(rows);
}

inline auto ones(Index rows) -> decltype(Vector::Ones(rows))
{
    return Vector::Ones(rows);
}

inline auto constants(Index rows, double val) -> decltype(Vector::Constant(rows, val))
{
    return Vector::Constant(rows, val);
}

inline auto random(Index rows) -> decltype(Vector::Random(rows))
{
    return Vector::Random(rows);
}

inline auto linspace(Index rows, double start, double stop) -> decltype(Vector::LinSpaced(rows, start, stop))
{
    return Vector::LinSpaced(rows, start, stop);
}

inline auto unit(Index rows, Index i) -> decltype(Vector::Unit(rows, i))
{
    return Vector::Unit(rows, i);
}

inline auto zeros(Index rows, Index cols) -> decltype(Matrix::Zero(rows, cols))
{
    return Matrix::Zero(rows, cols);
}

inline auto ones(Index rows, Index cols) -> decltype(Matrix::Ones(rows, cols))
{
    return Matrix::Ones(rows, cols);
}

inline auto random(Index rows, Index cols) -> decltype(Matrix::Random(rows, cols))
{
    return Matrix::Random(rows, cols);
}

inline auto identity(Index rows, Index cols) -> decltype(Matrix::Identity(rows, cols))
{
    return Matrix::Identity(rows, cols);
}

template<typename Derived>
auto rows(Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleRows(start, num))
{
    return mat.middleRows(start, num);
}

template<typename Derived>
auto rows(const Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleRows(start, num))
{
    return mat.middleRows(start, num);
}

template<typename Derived, typename Indices>
auto rows(Eigen::MatrixBase<Derived>& mat, const Indices& irows) -> Eigen::MatrixRowsView<Derived, Indices>
{
    return {mat, irows};
}

template<typename Derived, typename Indices>
auto rows(const Eigen::MatrixBase<Derived>& mat, const Indices& irows) -> Eigen::MatrixRowsViewConst<Derived, Indices>
{
    return {mat, irows};
}

template<typename Derived>
auto cols(Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleCols(start, num))
{
    return mat.middleCols(start, num);
}

template<typename Derived>
auto cols(const Eigen::MatrixBase<Derived>& mat, Index start, Index num) -> decltype(mat.middleCols(start, num))
{
    return mat.middleCols(start, num);
}

template<typename Derived, typename Indices>
auto cols(Eigen::MatrixBase<Derived>& mat, const Indices& icols) -> Eigen::MatrixColsView<Derived, Indices>
{
    return {mat, icols};
}

template<typename Derived, typename Indices>
auto cols(const Eigen::MatrixBase<Derived>& mat, const Indices& icols) -> Eigen::MatrixColsViewConst<Derived, Indices>
{
    return {mat, icols};
}

template<typename Derived>
auto segment(Eigen::MatrixBase<Derived>& vec, Index irow, Index nrows) -> decltype(vec.segment(irow, nrows))
{
    return vec.segment(irow, nrows);
}

template<typename Derived>
auto segment(const Eigen::MatrixBase<Derived>& vec, Index irow, Index nrows) -> decltype(vec.segment(irow, nrows))
{
    return vec.segment(irow, nrows);
}

template<typename Derived>
auto block(Eigen::MatrixBase<Derived>& mat, Index irow, Index icol, Index nrows, Index ncols) -> decltype(mat.block(irow, icol, nrows, ncols))
{
    return mat.block(irow, icol, nrows, ncols);
}

template<typename Derived>
auto block(const Eigen::MatrixBase<Derived>& mat, Index irow, Index icol, Index nrows, Index ncols) -> decltype(mat.block(irow, icol, nrows, ncols))
{
    return mat.block(irow, icol, nrows, ncols);
}

template<typename Derived, typename Indices>
auto submatrix(Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols) -> Eigen::MatrixSubView<Derived, Indices>
{
    return {mat, irows, icols};
}

template<typename Derived, typename Indices>
auto submatrix(const Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols) -> Eigen::MatrixSubViewConst<Derived, Indices>
{
    return {mat, irows, icols};
}

template<typename Derived>
auto tr(Eigen::MatrixBase<Derived>& mat) -> decltype(mat.transpose())
{
    return mat.transpose();
}

template<typename Derived>
auto tr(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.transpose())
{
    return mat.transpose();
}

template<typename Derived>
auto inv(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseInverse())
{
    return mat.cwiseInverse();
}

template<typename Derived>
auto diag(const Eigen::MatrixBase<Derived>& vec) -> decltype(vec.asDiagonal())
{
    return vec.asDiagonal();
}

template<typename Derived>
auto diagonal(Eigen::MatrixBase<Derived>& mat) -> decltype(mat.diagonal())
{
    return mat.diagonal();
}

template<typename Derived>
auto diagonal(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.diagonal())
{
    return mat.diagonal();
}

template<int p, typename Derived>
auto norm(const Eigen::MatrixBase<Derived>& mat) -> double
{
    return mat.template lpNorm<p>();
}

template<typename Derived>
auto norm(const Eigen::MatrixBase<Derived>& mat) -> double
{
    return mat.norm();
}

template<typename Derived>
auto norminf(const Eigen::MatrixBase<Derived>& mat) -> double
{
    return mat.template lpNorm<Eigen::Infinity>();
}

template<typename Derived>
auto sum(const Eigen::DenseBase<Derived>& mat) -> typename Derived::Scalar
{
    return mat.sum();
}

template<typename DerivedLHS, typename DerivedRHS>
auto dot(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.dot(rhs))
{
    return lhs.dot(rhs);
}

template<typename Derived>
auto min(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.minCoeff())
{
    return mat.minCoeff();
}

template<typename DerivedLHS, typename DerivedRHS>
auto min(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseMin(rhs))
{
    return lhs.cwiseMin(rhs);
}

template<typename Derived>
auto max(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.maxCoeff())
{
    return mat.maxCoeff();
}

template<typename DerivedLHS, typename DerivedRHS>
auto max(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseMax(rhs))
{
    return lhs.cwiseMax(rhs);
}

template<typename DerivedLHS, typename DerivedRHS>
auto operator%(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseProduct(rhs))
{
    return lhs.cwiseProduct(rhs);
}

template<typename DerivedLHS, typename DerivedRHS>
auto operator/(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseQuotient(rhs))
{
    return lhs.cwiseQuotient(rhs);
}

template<typename Derived>
auto operator/(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype(scalar*mat.cwiseInverse())
{
    return scalar*mat.cwiseInverse();
}

template<typename Derived>
auto operator+(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype((scalar + mat.array()).matrix())
{
    return (scalar + mat.array()).matrix();
}

template<typename Derived>
auto operator+(const Eigen::MatrixBase<Derived>& mat, const typename Derived::Scalar& scalar) -> decltype((scalar + mat.array()).matrix())
{
    return (scalar + mat.array()).matrix();
}

template<typename Derived>
auto operator-(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype((scalar - mat.array()).matrix())
{
    return (scalar - mat.array()).matrix();
}

template<typename Derived>
auto operator-(const Eigen::MatrixBase<Derived>& mat, const typename Derived::Scalar& scalar) -> decltype((mat.array() - scalar).matrix())
{
    return (mat.array() - scalar).matrix();
}

template<typename Derived>
auto abs(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseAbs())
{
    return mat.cwiseAbs();
}

template<typename Derived>
auto sqrt(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseSqrt())
{
    return mat.cwiseSqrt();
}

template<typename Derived>
auto pow(const Eigen::MatrixBase<Derived>& mat, double power) -> decltype(mat.array().pow(power).matrix())
{
    return mat.array().pow(power).matrix();
}

template<typename Derived>
auto exp(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().exp().matrix())
{
    return mat.array().exp().matrix();
}

template<typename Derived>
auto log(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().log().matrix())
{
    return mat.array().log().matrix();
}

template<typename Derived>
auto log10(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().log10().matrix())
{
    return mat.array().log10().matrix();
}

} // namespace Reaktoro
