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

VectorViewRows::VectorViewRows(Vector& vec, const Indices& irows)
: vec(vec), irows(irows) {}

VectorViewRowsConst::VectorViewRowsConst(const Vector& vec, const Indices& irows)
: vec(vec), irows(irows) {}

auto VectorViewRows::operator=(const Vector& other) -> VectorViewRows&
{
    for(unsigned i = 0; i < other.rows(); ++i)
        vec[irows[i]] = other[i];
    return *this;
}

VectorViewRows::operator Vector() const
{
    Vector res(irows.size());
    for(unsigned i = 0; i < res.size(); ++i)
        res[i] = vec[irows[i]];
    return res;
}

VectorViewRowsConst::operator Vector() const
{
    Vector res(irows.size());
    for(unsigned i = 0; i < res.size(); ++i)
        res[i] = vec[irows[i]];
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

inline auto rows(Vector& vec, const Indices& irows) -> VectorViewRows
{
    return VectorViewRows(vec, irows);
}

inline auto rows(const Vector& vec, const Indices& irows) -> VectorViewRowsConst
{
    return VectorViewRowsConst(vec, irows);
}

} // namespace Reaktor
