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

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

// Forward declaration
template<typename V, typename T, typename P>
class ThermoVectorBase;

/// A type that defines a vector of thermodynamic properties.
/// A thermo property means here any property that depends on
/// temperature and pressure. A ThermoVector instance
/// not only holds the values of the thermo properties, but also their
/// partial temperature and pressure derivatives.
/// @see ThermoScalar, ChemicalScalar, ChemicalVector
using ThermoVector = ThermoVectorBase<Vector,Vector,Vector>;

/// A type that defines a vector of thermodynamic properties.
using ThermoVectorRef = ThermoVectorBase<VectorRef,VectorRef,VectorRef>;

/// A type that defines a vector of thermodynamic properties.
using ThermoVectorConstRef = ThermoVectorBase<VectorConstRef,VectorConstRef,VectorConstRef>;

/// A template base class to represent a vector of thermodynamic scalars and their partial derivatives.
/// @see ThermoScalar, ThermoVector, ChemicalScalar, ChemicalVector
template<typename V, typename T, typename P>
class ThermoVectorBase
{
public:
    /// The vector of values of the thermodynamic properties.
    V val;

    /// The vector of partial temperature derivatives of the thermodynamic properties.
    T ddT;

    /// The vector of partial pressure derivatives of the thermodynamic properties.
    P ddP;

    /// Construct a default ThermoVectorBase instance.
    ThermoVectorBase()
    {}

    /// Construct a ThermoVectorBase instance with given number of rows.
    /// @param nrows The number of rows in the thermo vector
    explicit ThermoVectorBase(Index nrows)
    : ThermoVectorBase(zeros(nrows), zeros(nrows), zeros(nrows)) {}

    /// Construct a ThermoVectorBase instance with given number of rows and value.
    /// @param nrows The number of rows in the thermo vector
    /// @param val The constant value
    ThermoVectorBase(Index nrows, double val)
    : ThermoVectorBase(constants(nrows, val), zeros(nrows), zeros(nrows)) {}

    /// Construct a ChemicalVectorBase instance with given values and derivatives.
    /// @param val The vector of values of the thermo scalars
    /// @param ddT The vector of partial temperature derivatives of the thermo scalars
    /// @param ddP The vector of partial pressure derivatives of the thermo scalars
    ThermoVectorBase(const V& val, const T& ddT, const P& ddP)
    : val(val), ddT(ddT), ddP(ddP) {}

    /// Construct a ChemicalVectorBase instance from another.
    template<typename VR, typename TR, typename PR>
    ThermoVectorBase(ThermoVectorBase<VR,TR,PR>& other)
    : val(other.val), ddT(other.ddT), ddP(other.ddP)
    {}

    /// Construct a ChemicalVectorBase instance from another.
    template<typename VR, typename TR, typename PR>
    ThermoVectorBase(const ThermoVectorBase<VR,TR,PR>& other)
    : val(other.val), ddT(other.ddT), ddP(other.ddP)
    {}

    /// Return the number of rows in this ChemicalVectorBase instance.
    auto size() const -> Index
    {
        return val.size();
    }

    /// Resize this ChemicalVectorBase instance with new number of rows.
    /// @param nrows The new number of rows
    auto resize(Index nrows) -> void
    {
        val = zeros(nrows);
        ddT = zeros(nrows);
        ddP = zeros(nrows);
    }

    /// Assign a ThermoScalarBase instance to this.
    template<typename VR>
    auto fill(const ThermoScalarBase<VR>& other) -> void
    {
        val.fill(other.val);
        ddT.fill(other.ddT);
        ddP.fill(other.ddP);
    }

    /// Assign a scalarsto this.
    auto fill(double value) -> void
	{
    	val.fill(value);
    	ddT.fill(0.0);
    	ddP.fill(0.0);
	}

    /// Assign a ThermoVectorBase instance to this ThermoVectorBase instance.
    template<typename VR, typename TR, typename PR>
    auto operator=(const ThermoVectorBase<VR,TR,PR>& other) -> ThermoVectorBase&
    {
        val = other.val;
        ddT = other.ddT;
        ddP = other.ddP;
        return *this;
    }

    /// Assign a ThermoScalarBase instance to this ThermoVectorBase instance.
    template<typename VR>
    auto operator=(const ThermoScalarBase<VR>& other) -> ThermoVectorBase&
    {
        val.fill(other.val);
        ddT.fill(other.ddT);
        ddP.fill(other.ddP);
        return *this;
    }

    /// Assign a scalar to this ThermoVectorBase instance.
    auto operator=(double other) -> ThermoVectorBase&
    {
        val.fill(other);
        ddT.fill(0.0);
        ddP.fill(0.0);
        return *this;
    }

    /// Assign-addition of a ThermoVectorBase instance.
    template<typename VR, typename TR, typename PR>
    auto operator+=(const ThermoVectorBase<VR,TR,PR>& other) -> ThermoVectorBase&
    {
        val += other.val;
        ddT += other.ddT;
        ddP += other.ddP;
        return *this;
    }

    /// Assign-addition of a ThermoScalarBase instance.
    template<typename VR>
    auto operator+=(const ThermoScalarBase<VR>& other) -> ThermoVectorBase&
    {
        val.array() += other.val;
        ddT.array() += other.ddT;
        ddP.array() += other.ddP;
        return *this;
    }

    /// Assign-addition of a scalar.
    auto operator+=(double scalar) -> ThermoVectorBase&
    {
        val.array() += scalar;
        return *this;
    }

    /// Assign-subtraction of a ThermoVectorBase instance.
    template<typename VR, typename TR, typename PR>
    auto operator-=(const ThermoVectorBase<VR,TR,PR>& other) -> ThermoVectorBase&
    {
        val -= other.val;
        ddT -= other.ddT;
        ddP -= other.ddP;
        return *this;
    }

    /// Assign-subtraction of a ThermoVectorBase instance.
    template<typename VR>
    auto operator-=(const ThermoScalarBase<VR>& other) -> ThermoVectorBase&
    {
        val.array() -= other.val;
        ddT.array() -= other.ddT;
        ddP.array() -= other.ddP;
        return *this;
    }

    /// Assign-subtraction of a scalar.
    auto operator-=(double other) -> ThermoVectorBase&
    {
        val.array() -= other;
        return *this;
    }

    /// Assign-multiplication of a ThermoVectorBase instance.
    template<typename VR, typename TR, typename PR>
    auto operator*=(const ThermoVectorBase<VR,TR,PR>& other) -> ThermoVectorBase&
    {
        ddT = diag(ddT) * other.val + diag(val) * other.ddT;
        ddP = diag(ddP) * other.val + diag(val) * other.ddP;
        val = diag(val) * other.val;
        return *this;
    }

    /// Assign-multiplication of a ThermoVectorBase instance.
    auto operator*=(double scalar) -> ThermoVectorBase&
    {
        val *= scalar;
        ddT *= scalar;
        ddP *= scalar;
        return *this;
    }

    /// Assign-division of a ThermoVectorBase instance.
    auto operator/=(double scalar) -> ThermoVectorBase&
    {
        *this *= 1.0/scalar;
        return *this;
    }

    /// Return a ChemicalScalarBase with reference to the thermo scalar in a given row.
    auto operator[](Index irow) -> ThermoScalarBase<double&>
    {
        return {val[irow], ddT[irow], ddP[irow]};
    }

    /// Return a ChemicalScalarBase with const reference to the thermo scalar in a given row.
    auto operator[](Index irow) const -> ThermoScalarBase<double>
    {
        return {val[irow], ddT[irow], ddP[irow]};
    }

    /// Return a view of an interval of the ThermoVectorBase instance.
    /// @param irow The index of the row starting the view.
    /// @param nrows The number of rows in the view.
    auto view(Index irow, Index nrows) -> ThermoVectorRef
    {
        return {rowsmap(val, irow, nrows), rowsmap(ddT, irow, nrows), rowsmap(ddP, irow, nrows)};
    }

    /// Return a view of an interval of the ThermoVectorBase instance.
    /// @param irow The index of the row starting the view.
    /// @param nrows The number of rows in the view.
    auto view(Index irow, Index nrows) const -> ThermoVectorConstRef
    {
        return {rowsmap(val, irow, nrows), rowsmap(ddT, irow, nrows), rowsmap(ddP, irow, nrows)};
    }

    /// Explicitly converts this ThermoVector instance into a Vector.
    explicit operator Vector() const
    {
        return val;
    }
};

/// Return a ThermoVectorBase expression representing zeros with same dimension of given vector.
template<typename V, typename T, typename P>
auto zeros(const ThermoVectorBase<V,T,P>& v) -> ThermoVectorBase<decltype(zeros(0)), decltype(zeros(0)), decltype(zeros(0))>
{
    const Index n = v.size();
    return {zeros(n), zeros(n), zeros(n)};
}

/// Return a ThermoVectorBase expression representing ones with same dimension of given vector.
template<typename V, typename T, typename P>
auto ones(const ThermoVectorBase<V,T,P>& v) -> ThermoVectorBase<decltype(ones(0)), decltype(zeros(0)), decltype(zeros(0))>
{
    const Index n = v.size();
    return {ones(n), zeros(n), zeros(n)};
}

template<typename V, typename T, typename P>
auto operator<<(std::ostream& out, const ThermoVectorBase<V,T,P>& a) -> std::ostream&
{
    out << a.val;
    return out;
}

template<typename V, typename T, typename P>
auto operator+(const ThermoVectorBase<V,T,P>& l) -> ThermoVectorBase<V,T,P>
{
    return l;
}

template<typename V, typename T, typename P>
auto operator-(const ThermoVectorBase<V,T,P>& l) -> ThermoVectorBase<decltype(-l.val), decltype(-l.ddT), decltype(-l.ddP)>
{
    return {-l.val, -l.ddT, -l.ddP};
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator+(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ThermoVectorBase<decltype(l.val + r.val), decltype(l.ddT + r.ddT), decltype(l.ddP + r.ddP)>
{
    return {l.val + r.val, l.ddT + r.ddT, l.ddP + r.ddP};
}

template<typename V, typename T, typename P>
auto operator+(const ThermoVectorBase<V,T,P>& l, VectorConstRef r) -> ThermoVectorBase<decltype(l.val + r),T,P>
{
    return {l.val + r, l.ddT, l.ddP};
}

template<typename V, typename T, typename P>
auto operator+(VectorConstRef l, const ThermoVectorBase<V,T,P>& r) -> decltype(r + l)
{
    return r + l;
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator-(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ThermoVectorBase<decltype(l.val - r.val), decltype(l.ddT - r.ddT), decltype(l.ddP - r.ddP)>
{
    return {l.val - r.val, l.ddT - r.ddT, l.ddP - r.ddP};
}

template<typename V, typename T, typename P>
auto operator-(const ThermoVectorBase<V,T,P>& l, VectorConstRef r) -> ThermoVectorBase<decltype(l.val - r),T,P>
{
    return {l.val - r, l.ddT, l.ddP};
}

template<typename V, typename T, typename P>
auto operator-(VectorConstRef l, const ThermoVectorBase<V,T,P>& r) -> decltype(-(r - l))
{
    return -(r - l);
}

template<typename V, typename T, typename P>
auto operator*(double l, const ThermoVectorBase<V,T,P>& r) -> ThermoVectorBase<decltype(l * r.val), decltype(l * r.ddT), decltype(l * r.ddP)>
{
    return {l * r.val, l * r.ddT, l * r.ddP};
}

template<typename V, typename T, typename P>
auto operator*(const ThermoVectorBase<V,T,P>& l, double r) -> decltype(r * l)
{
    return r * l;
}

template<typename VL, typename V, typename T, typename P>
auto operator*(const ThermoScalarBase<VL>& l, const ThermoVectorBase<V,T,P>& r) -> ThermoVectorBase<decltype(l.val * r.val), decltype(l.val * r.ddT + l.ddT * r.val), decltype(l.val * r.ddP + l.ddP * r.val)>
{
    return {l.val * r.val, l.val * r.ddT + l.ddT * r.val, l.val * r.ddP + l.ddP * r.val};
}

template<typename V, typename T, typename P, typename VR>
auto operator*(const ThermoVectorBase<V,T,P>& l, const ThermoScalarBase<VR>& r) -> decltype(r * l)
{
    return r * l;
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator%(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ThermoVectorBase<decltype(diag(l.val) * r.val), decltype(diag(l.val) * r.ddT + diag(r.val) * l.ddT), decltype(diag(l.val) * r.ddP + diag(r.val) * l.ddP)>
{
    return {diag(l.val) * r.val,
            diag(l.val) * r.ddT + diag(r.val) * l.ddT,
            diag(l.val) * r.ddP + diag(r.val) * l.ddP};
}

template<typename V, typename T, typename P>
auto operator%(VectorConstRef l, const ThermoVectorBase<V,T,P>& r) -> ThermoVectorBase<decltype(diag(l) * r.val), decltype(diag(l) * r.ddT), decltype(diag(l) * r.ddP)>
{
    return {diag(l) * r.val,
            diag(l) * r.ddT,
            diag(l) * r.ddP};
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator%(const ThermoVectorBase<VL,TL,PL>& l, VectorConstRef r) -> decltype(r % l)
{
    return r % l;
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator/(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ThermoVector
{
    const Vector tmp = 1.0/(r.val % r.val);
    return {l.val/r.val,
            diag(tmp) * (diag(r.val) * l.ddT - diag(l.val) * r.ddT),
            diag(tmp) * (diag(r.val) * l.ddP - diag(l.val) * r.ddP)};
}

template<typename V, typename T, typename P, typename VR>
auto operator/(const ThermoVectorBase<V,T,P>& l, const ThermoScalarBase<VR>& r) -> ThermoVectorBase<decltype(l.val/r.val), decltype(double() * (l.ddT * r.val - l.val * r.ddT)), decltype(double() * (l.ddP * r.val - l.val * r.ddP))>
{
    const double tmp = 1.0/(r.val * r.val);
    return {l.val/r.val, tmp * (l.ddT * r.val - l.val * r.ddT), tmp * (l.ddP * r.val - l.val * r.ddP)};
}

template<typename V, typename T, typename P>
auto operator/(const ThermoVectorBase<V,T,P>& l, double r) -> decltype((1.0/r) * l)
{
    return (1.0/r) * l;
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator<(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> bool
{
    return l.val < r.val;
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator<=(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> bool
{
    return l.val <= r.val;
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator>(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> bool
{
    return l.val > r.val;
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator>=(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> bool
{
    return l.val >= r.val;
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator==(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> bool
{
    return l.val == r.val;
}

template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator!=(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> bool
{
    return l.val == r.val;
}

/// Return a reference of a row of this ThermoVectorBase instance.
template<typename V, typename T, typename P>
auto row(ThermoVectorBase<V,T,P>& vec, Index irow) -> ThermoScalarBase<double&>
{
    return {vec.val[irow], vec.ddT[irow], vec.ddP[irow]};
}

/// Return a const reference of a row of this ThermoVectorBase instance.
template<typename V, typename T, typename P>
auto row(const ThermoVectorBase<V,T,P>& vec, Index irow) -> ThermoScalarBase<const double&>
{
    return {vec.val[irow], vec.ddT[irow], vec.ddP[irow]};
}

/// Return a reference of a sequence of rows of this ThermoVectorBase instance.
template<typename V, typename T, typename P>
auto rows(ThermoVectorBase<V,T,P>& vec, Index irow, Index nrows) -> ThermoVectorBase<decltype(vec.val.segment(irow, nrows)), decltype(vec.ddT.segment(irow, nrows)), decltype(vec.ddP.segment(irow, nrows))>
{
    return {vec.val.segment(irow, nrows), vec.ddT.segment(irow, nrows), vec.ddP.segment(irow, nrows)};
}

/// Return a const reference of a sequence of rows of this ThermoVectorBase instance.
template<typename V, typename T, typename P>
auto rows(const ThermoVectorBase<V,T,P>& vec, Index irow, Index nrows) -> ThermoVectorBase<decltype(vec.val.segment(irow, nrows)), decltype(vec.ddT.segment(irow, nrows)), decltype(vec.ddP.segment(irow, nrows))>
{
    return {vec.val.segment(irow, nrows), vec.ddT.segment(irow, nrows), vec.ddP.segment(irow, nrows)};
}

/// Return a reference of some rows of this ThermoVectorBase instance.
template<typename V, typename T, typename P>
auto rows(ThermoVectorBase<V,T,P>& vec, const Indices& irows) -> ThermoVectorBase<decltype(Reaktoro::rows(vec.val, irows)), decltype(Reaktoro::rows(vec.ddT, irows)), decltype(Reaktoro::rows(vec.ddP, irows))>
{
    return {Reaktoro::rows(vec.val, irows), Reaktoro::rows(vec.ddT, irows), Reaktoro::rows(vec.ddP, irows)};
}

/// Return a const reference of some rows of this ThermoVectorBase instance.
template<typename V, typename T, typename P>
auto rows(const ThermoVectorBase<V,T,P>& vec, const Indices& irows) -> ThermoVectorBase<decltype(Reaktoro::rows(vec.val, irows)), decltype(Reaktoro::rows(vec.ddT, irows)), decltype(Reaktoro::rows(vec.ddP, irows))>
{
    return {Reaktoro::rows(vec.val, irows), Reaktoro::rows(vec.ddT, irows), Reaktoro::rows(vec.ddP, irows)};
}

template<typename V, typename T, typename P>
auto abs(const ThermoVectorBase<V,T,P>& l) -> ThermoVector
{
    const Vector tmp1 = abs(l.val);
    const Vector tmp2 = l.val/tmp1;
    return {tmp1, diag(tmp2) * l.ddT, diag(tmp2) * l.ddP};
}

template<typename V, typename T, typename P>
auto sqrt(const ThermoVectorBase<V,T,P>& l) -> ThermoVector
{
    const Vector tmp1 = sqrt(l.val);
    const Vector tmp2 = 0.5 * tmp1/l.val;
    return {tmp1, diag(tmp2) * l.ddT, diag(tmp2) * l.ddP};
}

template<typename V, typename T, typename P>
auto pow(const ThermoVectorBase<V,T,P>& l, double power) -> ThermoVector
{
    const Vector tmp1 = pow(l.val, power);
    const Vector tmp2 = power * tmp1/l.val;
    return {tmp1, diag(tmp2) * l.ddT, diag(tmp2) * l.ddP};
}

template<typename V, typename T, typename P>
auto exp(const ThermoVectorBase<V,T,P>& l) -> ThermoVector
{
    const Vector tmp = exp(l.val);
    return {tmp, diag(tmp) * l.ddT, diag(tmp) * l.ddP};
}

template<typename V, typename T, typename P>
auto log(const ThermoVectorBase<V,T,P>& l) -> ThermoVector
{
    const Vector tmp1 = log(l.val);
    const Vector tmp2 = 1.0/l.val;
    return {tmp1, diag(tmp2) * l.ddT, diag(tmp2) * l.ddP};
}

template<typename V, typename T, typename P>
auto log10(const ThermoVectorBase<V,T,P>& l) -> ThermoVector
{
    const double log10e = 0.4342944819032518166679324;
    const Vector tmp1 = log10e*log(l.val);
    const Vector tmp2 = log10e/l.val;
    return {tmp1, diag(tmp2) * l.ddT, diag(tmp2) * l.ddP, diag(tmp2) * l.ddn};
}

} // namespace Reaktoro
