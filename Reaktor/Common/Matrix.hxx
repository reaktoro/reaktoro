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

MatrixViewRows::MatrixViewRows(Matrix& mat, const Indices& irows)
: mat(mat), irows(irows) {}

MatrixViewRowsConst::MatrixViewRowsConst(const Matrix& mat, const Indices& irows)
: mat(mat), irows(irows) {}

MatrixViewCols::MatrixViewCols(Matrix& mat, const Indices& icols)
: mat(mat), icols(icols) {}

MatrixViewColsConst::MatrixViewColsConst(const Matrix& mat, const Indices& icols)
: mat(mat), icols(icols) {}

MatrixViewRowsCols::MatrixViewRowsCols(Matrix& mat, const Indices& irows, const Indices& icols)
: mat(mat), irows(irows), icols(icols) {}

MatrixViewRowsColsConst::MatrixViewRowsColsConst(const Matrix& mat, const Indices& irows, const Indices& icols)
: mat(mat), irows(irows), icols(icols) {}

auto MatrixViewRows::operator=(const Matrix& other) -> MatrixViewRows&
{
    for(unsigned i = 0; i < other.rows(); ++i)
        mat.row(irows[i]) = other.row(i);
    return *this;
}

auto MatrixViewCols::operator=(const Matrix& other) -> MatrixViewCols&
{
    for(unsigned i = 0; i < other.cols(); ++i)
        mat.col(icols[i]) = other.col(i);
    return *this;
}

auto MatrixViewCols::operator=(const Matrix& other) -> MatrixViewRowsCols&
{
    for(unsigned i = 0; i < other.rows(); ++i)
        for(unsigned j = 0; j < other.cols(); ++j)
            mat(irows[i], icols[j]) = other(i, j);
    return *this;
}

MatrixViewRows::operator Matrix() const
{
    Matrix res(irows.size(), mat.cols());
    for(unsigned i = 0; i < res.rows(); ++i)
        res.row(i) = mat.row(irows[i]);
    return res;
}

MatrixViewRowsConst::operator Matrix() const
{
    Matrix res(irows.size(), mat.cols());
    for(unsigned i = 0; i < res.rows(); ++i)
        res.row(i) = mat.row(irows[i]);
    return res;
}

MatrixViewCols::operator Matrix() const
{
    Matrix res(mat.rows(), icols.size());
    for(unsigned i = 0; i < res.cols(); ++i)
        res.col(i) = mat.col(icols[i]);
    return res;
}

MatrixViewColsConst::operator Matrix() const
{
    Matrix res(mat.rows(), icols.size());
    for(unsigned i = 0; i < res.cols(); ++i)
        res.col(i) = mat.col(icols[i]);
    return res;
}

MatrixViewRowsCols::operator Matrix() const
{
    Matrix res(irows.size(), icols.size());
    for(unsigned i = 0; i < res.rows(); ++i)
        for(unsigned j = 0; j < res.cols(); ++j)
            res(i, j) = mat(irows[i], icols[j]);
    return res;
}

MatrixViewRowsColsConst::operator Matrix() const
{
    Matrix res(irows.size(), icols.size());
    for(unsigned i = 0; i < res.rows(); ++i)
        for(unsigned j = 0; j < res.cols(); ++j)
            res(i, j) = mat(irows[i], icols[j]);
    return res;
}

inline auto identity(unsigned rows, unsigned cols) -> decltype(Matrix::Identity(rows, cols))
{
    return Matrix::Identity(rows, cols);
}

inline auto zeros(unsigned rows, unsigned cols) -> decltype(Matrix::Zero(rows, cols))
{
    return Matrix::Zero(rows, cols);
}

inline auto ones(unsigned rows, unsigned cols) -> decltype(Matrix::Ones(rows, cols))
{
    return Matrix::Ones(rows, cols);
}

template<typename MatrixType>
inline auto rows(MatrixType&& mat, unsigned start, unsigned num) -> decltype(mat.middleRows(start, num))
{
    return mat.middleRows(start, num);
}

inline auto rows(Matrix& mat, const Indices& irows) -> MatrixViewRows
{
    return MatrixViewRows(mat, irows);
}

inline auto rows(const Matrix& mat, const Indices& irows) -> MatrixViewRowsConst
{
    return MatrixViewRowsConst(mat, irows);
}

template<typename MatrixType>
inline auto cols(MatrixType&& mat, unsigned start, unsigned num) -> decltype(mat.middleCols(start, num))
{
    return mat.middleCols(start, num);
}

inline auto cols(Matrix& mat, const Indices& icols) -> MatrixViewCols
{
    return MatrixViewCols(mat, icols);
}

inline auto cols(const Matrix& mat, const Indices& icols) -> MatrixViewColsConst
{
    return MatrixViewColsConst(mat, icols);
}

inline auto submatrix(Matrix& mat, const Indices& irows, const Indices& icols) -> MatrixViewRowsCols
{
    return MatrixViewRowsCols(mat, irows, icols);
}

inline auto submatrix(const Matrix& mat, const Indices& irows, const Indices& icols) -> MatrixViewRowsColsConst
{
    return MatrixViewRowsColsConst(mat, irows, icols);
}


} // namespace Reaktor
