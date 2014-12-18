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

namespace Reaktor {

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
template<typename DerivedOther>
auto MatrixViewRows<Derived>::operator=(const Eigen::MatrixBase<DerivedOther>& other) -> MatrixViewRows&
{
    for(unsigned i = 0; i < other.rows(); ++i)
        mat.row(irows[i]) = other.row(i);
    return *this;
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewCols<Derived>::operator=(const Eigen::MatrixBase<DerivedOther>& other) -> MatrixViewCols&
{
    for(unsigned i = 0; i < other.cols(); ++i)
        mat.col(icols[i]) = other.col(i);
    return *this;
}

template<typename Derived>
template<typename DerivedOther>
auto MatrixViewRowsCols<Derived>::operator=(const Eigen::MatrixBase<DerivedOther>& other) -> MatrixViewRowsCols&
{
    for(unsigned i = 0; i < other.rows(); ++i)
        for(unsigned j = 0; j < other.cols(); ++j)
            mat(irows[i], icols[j]) = other(i, j);
    return *this;
}

template<typename Derived>
MatrixViewRows<Derived>::operator Derived() const
{
    Derived res(irows.size(), mat.cols());
    for(unsigned i = 0; i < res.rows(); ++i)
        res.row(i) = mat.row(irows[i]);
    return res;
}

template<typename Derived>
MatrixViewRowsConst<Derived>::operator Derived() const
{
    Derived res(irows.size(), mat.cols());
    for(unsigned i = 0; i < res.rows(); ++i)
        res.row(i) = mat.row(irows[i]);
    return res;
}

template<typename Derived>
MatrixViewCols<Derived>::operator Derived() const
{
    Derived res(mat.rows(), icols.size());
    for(unsigned i = 0; i < res.cols(); ++i)
        res.col(i) = mat.col(icols[i]);
    return res;
}

template<typename Derived>
MatrixViewColsConst<Derived>::operator Derived() const
{
    Derived res(mat.rows(), icols.size());
    for(unsigned i = 0; i < res.cols(); ++i)
        res.col(i) = mat.col(icols[i]);
    return res;
}

template<typename Derived>
MatrixViewRowsCols<Derived>::operator Derived() const
{
    Derived res(irows.size(), icols.size());
    for(unsigned i = 0; i < res.rows(); ++i)
        for(unsigned j = 0; j < res.cols(); ++j)
            res(i, j) = mat(irows[i], icols[j]);
    return res;
}

template<typename Derived>
MatrixViewRowsColsConst<Derived>::operator Derived() const
{
    Derived res(irows.size(), icols.size());
    for(unsigned i = 0; i < res.rows(); ++i)
        for(unsigned j = 0; j < res.cols(); ++j)
            res(i, j) = mat(irows[i], icols[j]);
    return res;
}

inline auto zeros(unsigned rows) -> decltype(Vector::Zero(rows))
{
    return Vector::Zero(rows);
}

inline auto ones(unsigned rows) -> decltype(Vector::Ones(rows))
{
    return Vector::Ones(rows);
}

inline auto zeros(unsigned rows, unsigned cols) -> decltype(Matrix::Zero(rows, cols))
{
    return Matrix::Zero(rows, cols);
}

inline auto ones(unsigned rows, unsigned cols) -> decltype(Matrix::Ones(rows, cols))
{
    return Matrix::Ones(rows, cols);
}

inline auto identity(unsigned rows, unsigned cols) -> decltype(Matrix::Identity(rows, cols))
{
    return Matrix::Identity(rows, cols);
}

template<typename MatrixType>
inline auto rows(MatrixType&& mat, unsigned start, unsigned num) -> decltype(mat.middleRows(start, num))
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
inline auto cols(MatrixType&& mat, unsigned start, unsigned num) -> decltype(mat.middleCols(start, num))
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

template<typename MatrixType>
inline auto block(MatrixType&& mat, unsigned irow, unsigned icol, unsigned nrows, unsigned ncols) -> decltype(mat.block(irow, icol, nrows, ncols))
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
inline auto norminf(const Eigen::DenseBase<Derived>& mat) -> double
{
    return mat.template lpNorm<Eigen::Infinity>();
}

template<typename Derived>
inline auto sum(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.sum())
{
    return mat.sum();
}

template<typename Derived>
inline auto sum(const Eigen::ArrayBase<Derived>& arr) -> decltype(arr.sum())
{
    return arr.sum();
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
inline auto sqrt(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().sqrt())
{
    return mat.array().sqrt();
}

template<typename Derived>
inline auto log(const Eigen::MatrixBase<Derived>& mat) -> decltype(mat.array().log())
{
    return mat.array().log();
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

} // namespace Reaktor
