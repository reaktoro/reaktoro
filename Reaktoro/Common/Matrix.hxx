// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
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

namespace Reaktoro {

template<typename Derived>
MatrixViewRows<Derived>::MatrixViewRows(Eigen::MatrixBase<Derived>& mat, const Indices& irows)
: mat(mat), irows(irows) {}

template<typename Derived>
MatrixViewRowsConst<Derived>::MatrixViewRowsConst(const Eigen::MatrixBase<Derived>& mat, const Indices& irows)
: mat(mat), irows(irows) {}

template<typename Derived>
MatrixViewCols<Derived>::MatrixViewCols(Eigen::MatrixBase<Derived>& mat, const Indices& icols)
: mat(mat), icols(icols) {}

template<typename Derived>
MatrixViewColsConst<Derived>::MatrixViewColsConst(const Eigen::MatrixBase<Derived>& mat, const Indices& icols)
: mat(mat), icols(icols) {}

template<typename Derived>
MatrixViewRowsCols<Derived>::MatrixViewRowsCols(Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols)
: mat(mat), irows(irows), icols(icols) {}

template<typename Derived>
MatrixViewRowsColsConst<Derived>::MatrixViewRowsColsConst(const Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols)
: mat(mat), irows(irows), icols(icols) {}

template<typename Derived>
auto MatrixViewRows<Derived>::operator=(const MatrixViewRows& other) -> MatrixViewRows&
{
    for(Index i = 0; i < other.irows.size(); ++i)
        mat.row(irows[i]) = other.mat.row(other.irows[i]);
    return *this;
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewRows<Derived>::operator=(const MatrixViewRowsConst<DerivedOther>& other) -> MatrixViewRows&
{
    for(Index i = 0; i < other.irows.size(); ++i)
        mat.row(irows[i]) = other.mat.row(other.irows[i]);
    return *this;
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewRows<Derived>::operator=(const Eigen::MatrixBase<DerivedOther>& other) -> MatrixViewRows&
{
    for(int i = 0; i < other.rows(); ++i)
        mat.row(irows[i]) = other.row(i);
    return *this;
}

template<typename Derived>
auto MatrixViewRows<Derived>::operator=(const typename Derived::Scalar& scalar) -> MatrixViewRows&
{
    for(Index i : irows)
        mat.row(i).fill(scalar);
    return *this;
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewRows<Derived>::to(Eigen::MatrixBase<DerivedOther>& other) const -> void
{
    other.derived().resize(irows.size(), mat.cols());
    for(Index i = 0; i < irows.size(); ++i)
        other.row(i) = mat.row(irows[i]);
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewRowsConst<Derived>::to(Eigen::MatrixBase<DerivedOther>& other) const -> void
{
    other.derived().resize(irows.size(), mat.cols());
    for(Index i = 0; i < irows.size(); ++i)
        other.row(i) = mat.row(irows[i]);
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewCols<Derived>::to(Eigen::MatrixBase<DerivedOther>& other) const -> void
{
    other.derived().resize(mat.rows(), icols.size());
    for(Index i = 0; i < icols.size(); ++i)
        other.col(i) = mat.col(icols[i]);
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewColsConst<Derived>::to(Eigen::MatrixBase<DerivedOther>& other) const -> void
{
    other.derived().resize(mat.rows(), icols.size());
    for(Index i = 0; i < icols.size(); ++i)
        other.col(i) = mat.col(icols[i]);
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewRowsCols<Derived>::to(Eigen::MatrixBase<DerivedOther>& other) const -> void
{
    other.derived().resize(irows.size(), icols.size());
    for(Index i = 0; i < irows.size(); ++i)
        for(Index j = 0; j < icols.size(); ++j)
            other(i, j) = mat(irows[i], icols[j]);
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewRowsColsConst<Derived>::to(Eigen::MatrixBase<DerivedOther>& other) const -> void
{
    other.derived().resize(irows.size(), icols.size());
    for(Index i = 0; i < irows.size(); ++i)
        for(Index j = 0; j < icols.size(); ++j)
            other(i, j) = mat(irows[i], icols[j]);
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewCols<Derived>::operator=(const Eigen::MatrixBase<DerivedOther>& other) -> MatrixViewCols&
{
    for(int i = 0; i < other.cols(); ++i)
        mat.col(icols[i]) = other.col(i);
    return *this;
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewRowsCols<Derived>::operator=(const Eigen::MatrixBase<DerivedOther>& other) -> MatrixViewRowsCols&
{
    for(int i = 0; i < other.rows(); ++i)
        for(int j = 0; j < other.cols(); ++j)
            mat(irows[i], icols[j]) = other(i, j);
    return *this;
}

template<typename Derived>
MatrixViewRows<Derived>::operator Derived() const
{
    Derived res(irows.size(), mat.cols());
    for(int i = 0; i < res.rows(); ++i)
        res.row(i) = mat.row(irows[i]);
    return res;
}

template<typename Derived>
MatrixViewRowsConst<Derived>::operator Derived() const
{
    Derived res(irows.size(), mat.cols());
    for(int i = 0; i < res.rows(); ++i)
        res.row(i) = mat.row(irows[i]);
    return res;
}

template<typename Derived>
MatrixViewCols<Derived>::operator Derived() const
{
    Derived res(mat.rows(), icols.size());
    for(int i = 0; i < res.cols(); ++i)
        res.col(i) = mat.col(icols[i]);
    return res;
}

template<typename Derived>
MatrixViewColsConst<Derived>::operator Derived() const
{
    Derived res(mat.rows(), icols.size());
    for(int i = 0; i < res.cols(); ++i)
        res.col(i) = mat.col(icols[i]);
    return res;
}

template<typename Derived>
MatrixViewRowsCols<Derived>::operator Derived() const
{
    Derived res(irows.size(), icols.size());
    for(int i = 0; i < res.rows(); ++i)
        for(int j = 0; j < res.cols(); ++j)
            res(i, j) = mat(irows[i], icols[j]);
    return res;
}

template<typename Derived>
MatrixViewRowsColsConst<Derived>::operator Derived() const
{
    Derived res(irows.size(), icols.size());
    for(int i = 0; i < res.rows(); ++i)
        for(int j = 0; j < res.cols(); ++j)
            res(i, j) = mat(irows[i], icols[j]);
    return res;
}

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

inline auto identity(Index rows, Index cols) -> decltype(Matrix::Identity(rows, cols))
{
    return Matrix::Identity(rows, cols);
}

template<typename MatrixType>
inline auto rows(MatrixType&& mat, Index start, Index num) -> decltype(mat.middleRows(start, num))
{
    return mat.middleRows(start, num);
}

template<typename Derived>
inline auto rows(Eigen::MatrixBase<Derived>& mat, const Indices& irows) -> MatrixViewRows<Derived>
{
    return MatrixViewRows<Derived>(mat, irows);
}

template<typename Derived>
inline auto rows(const Eigen::MatrixBase<Derived>& mat, const Indices& irows) -> MatrixViewRowsConst<Derived>
{
    return MatrixViewRowsConst<Derived>(mat, irows);
}

template<typename MatrixType>
inline auto cols(MatrixType&& mat, Index start, Index num) -> decltype(mat.middleCols(start, num))
{
    return mat.middleCols(start, num);
}

template<typename Derived>
inline auto cols(Eigen::MatrixBase<Derived>& mat, const Indices& icols) -> MatrixViewCols<Derived>
{
    return MatrixViewCols<Derived>(mat, icols);
}

template<typename Derived>
inline auto cols(const Eigen::MatrixBase<Derived>& mat, const Indices& icols) -> MatrixViewColsConst<Derived>
{
    return MatrixViewColsConst<Derived>(mat, icols);
}

template<typename VectorType>
inline auto segment(VectorType&& vec, Index irow, Index nrows) -> decltype(vec.segment(irow, nrows))
{
    return vec.segment(irow, nrows);
}

template<typename MatrixType>
inline auto block(MatrixType&& mat, Index irow, Index icol, Index nrows, Index ncols) -> decltype(mat.block(irow, icol, nrows, ncols))
{
    return mat.block(irow, icol, nrows, ncols);
}

template<typename Derived>
inline auto submatrix(Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols) -> MatrixViewRowsCols<Derived>
{
    return MatrixViewRowsCols<Derived>(mat, irows, icols);
}

template<typename Derived>
inline auto submatrix(const Eigen::MatrixBase<Derived>& mat, const Indices& irows, const Indices& icols) -> MatrixViewRowsColsConst<Derived>
{
    return MatrixViewRowsColsConst<Derived>(mat, irows, icols);
}

template<typename Derived>
auto rowascol(Eigen::MatrixBase<Derived>& mat, Index irow) -> Eigen::Map<Eigen::Matrix<typename Derived::Scalar, -1, 1>, 0, Eigen::InnerStride<>>
{
    return Eigen::Map<Eigen::Matrix<typename Derived::Scalar, -1, 1>, 0, Eigen::InnerStride<>>(mat.row(irow).data(), mat.cols(), Eigen::InnerStride<>(mat.rows()));
}

template<typename Derived>
auto rowascol(const Eigen::MatrixBase<Derived>& mat, Index irow) -> Eigen::Map<const Eigen::Matrix<typename Derived::Scalar, -1, 1>, 0, Eigen::InnerStride<>>
{
    return Eigen::Map<const Eigen::Matrix<typename Derived::Scalar, -1, 1>, 0, Eigen::InnerStride<>>(mat.row(irow).data(), mat.cols(), Eigen::InnerStride<>(mat.rows()));
}

template<typename Derived>
inline auto tr(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.transpose())
{
    return mat.transpose();
}

template<typename Derived>
inline auto inv(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseInverse())
{
    return mat.cwiseInverse();
}

template<typename Derived>
inline auto diag(const Eigen::MatrixBase<Derived>& vec) -> decltype(vec.asDiagonal())
{
    return vec.asDiagonal();
}

template<typename MatrixType>
inline auto diagonal(MatrixType&& mat) -> decltype(mat.diagonal())
{
    return mat.diagonal();
}

template<int p, typename Derived>
inline auto norm(const Eigen::MatrixBase<Derived>& mat) -> double
{
    return mat.template lpNorm<p>();
}

template<typename Derived>
inline auto norm(const Eigen::MatrixBase<Derived>& mat) -> double
{
    return mat.norm();
}

template<typename Derived>
inline auto norminf(const Eigen::MatrixBase<Derived>& mat) -> double
{
    return mat.template lpNorm<Eigen::Infinity>();
}

template<typename Derived>
inline auto sum(const Eigen::DenseBase<Derived>& mat) -> typename Derived::Scalar
{
    return mat.sum();
}

template<typename DerivedLHS, typename DerivedRHS>
inline auto dot(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.dot(rhs))
{
    return lhs.dot(rhs);
}

template<typename Derived>
inline auto min(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.minCoeff())
{
    return mat.minCoeff();
}

template<typename DerivedLHS, typename DerivedRHS>
inline auto min(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseMin(rhs))
{
    return lhs.cwiseMin(rhs);
}

template<typename Derived>
inline auto max(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.maxCoeff())
{
    return mat.maxCoeff();
}

template<typename DerivedLHS, typename DerivedRHS>
inline auto max(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseMax(rhs))
{
    return lhs.cwiseMax(rhs);
}

template<typename Derived>
inline auto abs(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseAbs())
{
    return mat.cwiseAbs();
}

template<typename Derived>
inline auto sqrt(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.cwiseSqrt())
{
    return mat.cwiseSqrt();
}

template<typename Derived>
inline auto pow(const Eigen::MatrixBase<Derived>& mat, double power) -> decltype(mat.array().pow(power).matrix())
{
    return mat.array().pow(power).matrix();
}

template<typename Derived>
inline auto exp(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().exp().matrix())
{
    return mat.array().exp().matrix();
}

template<typename Derived>
inline auto log(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().log().matrix())
{
    return mat.array().log().matrix();
}

template<typename Derived>
inline auto log10(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().log10().matrix())
{
    return mat.array().log10().matrix();
}

template<typename DerivedLHS, typename DerivedRHS>
inline auto operator%(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseProduct(rhs))
{
    return lhs.cwiseProduct(rhs);
}

template<typename DerivedLHS, typename DerivedRHS>
inline auto operator/(const Eigen::MatrixBase<DerivedLHS>& lhs, const Eigen::MatrixBase<DerivedRHS>& rhs) -> decltype(lhs.cwiseQuotient(rhs))
{
    return lhs.cwiseQuotient(rhs);
}

template<typename Derived>
inline auto operator/(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype(scalar*mat.cwiseInverse())
{
    return scalar*mat.cwiseInverse();
}

template<typename Derived>
inline auto operator+(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype((scalar + mat.array()).matrix())
{
    return (scalar + mat.array()).matrix();
}

template<typename Derived>
inline auto operator+(const Eigen::MatrixBase<Derived>& mat, const typename Derived::Scalar& scalar) -> decltype((scalar + mat.array()).matrix())
{
    return (scalar + mat.array()).matrix();
}

template<typename Derived>
inline auto operator-(const typename Derived::Scalar& scalar, const Eigen::MatrixBase<Derived>& mat) -> decltype((scalar - mat.array()).matrix())
{
    return (scalar - mat.array()).matrix();
}

template<typename Derived>
inline auto operator-(const Eigen::MatrixBase<Derived>& mat, const typename Derived::Scalar& scalar) -> decltype((mat.array() - scalar).matrix())
{
    return (mat.array() - scalar).matrix();
}

} // namespace Reaktoro
