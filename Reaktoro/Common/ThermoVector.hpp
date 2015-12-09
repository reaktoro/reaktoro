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

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>

namespace Reaktoro {

/// A template base class to represent a vector of thermodynamic scalars and their partial derivatives.
/// @see ThermoScalar, ThermoVector, ChemicalScalar, ChemicalVector
template<typename V, typename T, typename P>
class ThermoVectorBase
{
public:
    /// The vector of values of the thermodynamic properties.
    V val;

    /// The vector of partial temperature derivatives of the thermodynamic properties.
    T ddt;

    /// The vector of partial pressure derivatives of the thermodynamic properties.
    P ddp;

    /// Return a ThermoVector with zeros and zero derivatives.
    /// @param nrows The number of rows in the thermo vector
    static auto Zero(Index nrows) -> ThermoVectorBase
    {
        return ThermoVectorBase<V,T,P>(nrows, 0.0);
    }

    /// Return a ThermoVector with ones and zero derivatives.
    /// @param nrows The number of rows in the thermo vector
    static auto One(Index nrows) -> ThermoVectorBase
    {
        return ThermoVectorBase<V,T,P>(nrows, 1.0);
    }

    /// Return a ThermoVector with a given constant and zero derivatives.
    /// @param nrows The number of rows in the thermo vector
    static auto Constant(Index nrows, double val) -> ThermoVectorBase
    {
        return ThermoVectorBase<V,T,P>(nrows, val);
    }

    /// Construct a default ThermoVectorBase instance.
    ThermoVectorBase()
    {}

    /// Construct a ThermoVectorBase instance with given number of rows.
    /// @param nrows The number of rows in the thermo vector
    ThermoVectorBase(Index nrows)
    : ThermoVectorBase(zeros(nrows), zeros(nrows), zeros(nrows)) {}

    /// Construct a ThermoVectorBase instance with given number of rows and value.
    /// @param nrows The number of rows in the thermo vector
    /// @param val The constant value
    ThermoVectorBase(Index nrows, double val)
    : ThermoVectorBase(constants(nrows, val), zeros(nrows), zeros(nrows)) {}

    /// Construct a ChemicalVectorBase instance with given values and derivatives.
    /// @param val The vector of values of the thermo scalars
    /// @param ddt The vector of partial temperature derivatives of the thermo scalars
    /// @param ddp The vector of partial pressure derivatives of the thermo scalars
    ThermoVectorBase(const V& val, const T& ddt, const P& ddp)
    : val(val), ddt(ddt), ddp(ddp) {}

    /// Construct a ChemicalVectorBase instance from another.
    template<typename VR, typename TR, typename PR>
    ThermoVectorBase(const ThermoVectorBase<VR,TR,PR>& other)
    : val(other.val), ddt(other.ddt), ddp(other.ddp)
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
        val.resize(nrows);
        ddt.resize(nrows);
        ddp.resize(nrows);
    }

    /// Assign a ThermoVectorBase instance to this ThermoVectorBase instance.
    template<typename VR, typename TR, typename PR>
    auto operator=(const ThermoVectorBase<VR,TR,PR>& other) -> ThermoVectorBase&
    {
        val = other.val;
        ddt = other.ddt;
        ddp = other.ddp;
        return *this;
    }

    /// Assign a ThermoScalarBase instance to this ThermoVectorBase instance.
    template<typename VR>
    auto operator=(const ThermoScalarBase<VR>& other) -> ThermoVectorBase&
    {
        val.fill(other.val);
        ddt.fill(other.ddt);
        ddp.fill(other.ddp);
        return *this;
    }

    /// Assign a scalar to this ThermoVectorBase instance.
    auto operator=(double other) -> ThermoVectorBase&
    {
        val.fill(other);
        ddt.fill(0.0);
        ddp.fill(0.0);
        return *this;
    }

    /// Assign-addition of a ThermoVectorBase instance.
    template<typename VR, typename TR, typename PR>
    auto operator+=(const ThermoVectorBase<VR,TR,PR>& other) -> ThermoVectorBase&
    {
        val += other.val;
        ddt += other.ddt;
        ddp += other.ddp;
        return *this;
    }

    /// Assign-addition of a ThermoScalarBase instance.
    template<typename VR>
    auto operator+=(const ThermoScalarBase<VR>& other) -> ThermoVectorBase&
    {
        val.array() += other.val;
        ddt.array() += other.ddt;
        ddp.array() += other.ddp;
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
        ddt -= other.ddt;
        ddp -= other.ddp;
        return *this;
    }

    /// Assign-subtraction of a ThermoScalar instance.
    template<typename VR>
    auto operator-=(const ThermoScalarBase<VR>& other) -> ThermoVectorBase&
    {
        val.array() -= other.val;
        ddt.array() -= other.ddt;
        ddp.array() -= other.ddp;
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
        ddt = diag(ddt) * other.val + diag(val) * other.ddt;
        ddp = diag(ddp) * other.val + diag(val) * other.ddp;
        val = diag(val) * other.val;
        return *this;
    }

    /// Assign-multiplication of a ThermoVectorBase instance.
    auto operator*=(double scalar) -> ThermoVectorBase&
    {
        val *= scalar;
        ddt *= scalar;
        ddp *= scalar;
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
        return {val[irow], ddt[irow], ddp[irow]};
    }

    /// Return a ChemicalScalarBase with const reference to the thermo scalar in a given row.
    auto operator[](Index irow) const -> ThermoScalarBase<const double&>
    {
        return {val[irow], ddt[irow], ddp[irow]};
    }

    /// Return a reference of a row of this ThermoVectorBase instance.
    auto row(Index irow) -> ThermoScalarBase<double&>
    {
        return {val[irow], ddt[irow], ddp[irow]};
    }

    /// Return a const reference of a row of this ThermoVectorBase instance.
    auto row(Index irow) const -> ThermoScalarBase<const double&>
    {
        return {val[irow], ddt[irow], ddp[irow]};
    }

    /// Return a reference of a sequence of rows of this ThermoVectorBase instance.
    auto rows(Index irow, Index nrows) -> ThermoVectorBase<decltype(val.segment(irow, nrows)), decltype(ddt.segment(irow, nrows)), decltype(ddp.segment(irow, nrows))>
    {
        return {val.segment(irow, nrows), ddt.segment(irow, nrows), ddp.segment(irow, nrows)};
    }

    /// Return a const reference of a sequence of rows of this ThermoVectorBase instance.
    auto rows(Index irow, Index nrows) const -> ThermoVectorBase<decltype(val.segment(irow, nrows)), decltype(ddt.segment(irow, nrows)), decltype(ddp.segment(irow, nrows))>
    {
        return {val.segment(irow, nrows), ddt.segment(irow, nrows), ddp.segment(irow, nrows)};
    }

    /// Return a reference of some rows of this ThermoVectorBase instance.
    auto rows(const Indices& irows) -> ThermoVectorBase<decltype(Reaktoro::rows(val, irows)), decltype(Reaktoro::rows(ddt, irows)), decltype(Reaktoro::rows(ddp, irows))>
    {
        return {Reaktoro::rows(val, irows), Reaktoro::rows(ddt, irows), Reaktoro::rows(ddp, irows)};
    }

    /// Return a const reference of some rows of this ThermoVectorBase instance.
    auto rows(const Indices& irows) const -> ThermoVectorBase<decltype(Reaktoro::rows(val, irows)), decltype(Reaktoro::rows(ddt, irows)), decltype(Reaktoro::rows(ddp, irows))>
    {
        return {Reaktoro::rows(val, irows), Reaktoro::rows(ddt, irows), Reaktoro::rows(ddp, irows)};
    }
};

/// A type that defines a vector thermo property.
/// A thermo property means here any property that depends on
/// temperature and pressure. A ThermoVector instance
/// not only holds the values of the thermo properties, but also their
/// partial temperature and pressure derivatives.
/// @see ThermoScalar, ChemicalScalar, ChemicalVector
using ThermoVector = ThermoVectorBase<Vector,Vector,Vector>;

/// Unary addition operator for a ThermoScalar instance
template<typename V, typename T, typename P>
auto operator+(const ThermoVectorBase<V,T,P>& l) -> ThermoVectorBase<V,T,P>
{
    return l;
}

/// Unary subtraction operator for a ThermoScalar instance
template<typename V, typename T, typename P>
auto operator-(const ThermoVectorBase<V,T,P>& l) -> ThermoVectorBase<decltype(-l.val), decltype(-l.ddt), decltype(-l.ddp)>
{
    return {-l.val, -l.ddt, -l.ddp};
}

/// Add two ThermoScalar instances
template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator+(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ThermoVectorBase<decltype(l.val + r.val), decltype(l.ddt + r.ddt), decltype(l.ddp + r.ddp)>
{
    return {l.val + r.val, l.ddt + r.ddt, l.ddp + r.ddp};
}

/// Add a ThermoScalar instance and a scalar
template<typename V, typename T, typename P>
auto operator+(const ThermoVectorBase<V,T,P>& l, const Vector& r) -> ThermoVectorBase<decltype(l.val + r),T,P>
{
    return {l.val + r, l.ddt, l.ddp};
}

/// Add a scalar and a ThermoScalar instance
template<typename V, typename T, typename P>
auto operator+(const Vector& l, const ThermoVectorBase<V,T,P>& r) -> decltype(r + l)
{
    return r + l;
}

/// Subtract two ThermoScalar instances
template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator-(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ThermoVectorBase<decltype(l.val - r.val), decltype(l.ddt - r.ddt), decltype(l.ddp - r.ddp)>
{
    return {l.val - r.val, l.ddt - r.ddt, l.ddp - r.ddp};
}

/// Subtract a ThermoScalar instance and a scalar
template<typename V, typename T, typename P>
auto operator-(const ThermoVectorBase<V,T,P>& l, const Vector& r) -> ThermoVectorBase<decltype(l.val - r),T,P>
{
    return {l.val - r, l.ddt, l.ddp};
}

/// Subtract a scalar and a ThermoScalar instance
template<typename V, typename T, typename P>
auto operator-(const Vector& l, const ThermoVectorBase<V,T,P>& r) -> decltype(-(r - l))
{
    return -(r - l);
}

/// Left-multiply a ThermoScalar instance by a scalar
template<typename V, typename T, typename P>
auto operator*(double l, const ThermoVectorBase<V,T,P>& r) -> ThermoVectorBase<decltype(l * r.val), decltype(l * r.ddt), decltype(l * r.ddp)>
{
    return {l * r.val, l * r.ddt, l * r.ddp};
}

/// Right-multiply a ThermoScalar instance by a scalar
template<typename V, typename T, typename P>
auto operator*(const ThermoVectorBase<V,T,P>& l, double r) -> decltype(r * l)
{
    return r * l;
}

/// Left-multiply a ThermoScalar instance by a ThermoScalar
template<typename VL, typename V, typename T, typename P>
auto operator*(const ThermoScalarBase<VL>& l, const ThermoVectorBase<V,T,P>& r) -> ThermoVectorBase<decltype(l.val * r.val), decltype(l.val * r.ddt + l.ddt * r.val), decltype(l.val * r.ddp + l.ddp * r.val)>
{
    return {l.val * r.val, l.val * r.ddt + l.ddt * r.val, l.val * r.ddp + l.ddp * r.val};
}

/// Right-multiply a ThermoScalar instance by a ThermoScalar
template<typename V, typename T, typename P, typename VR>
auto operator*(const ThermoVectorBase<V,T,P>& l, const ThermoScalarBase<VR>& r) -> decltype(r * l)
{
    return r * l;
}

/// Multiply two ThermoVector instances
template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator%(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ThermoVectorBase<decltype(diag(l.val) * r.val), decltype(diag(l.val) * r.ddt + diag(r.val) * l.ddt), decltype(diag(l.val) * r.ddp + diag(r.val) * l.ddp)>
{
    return {diag(l.val) * r.val,
            diag(l.val) * r.ddt + diag(r.val) * l.ddt,
            diag(l.val) * r.ddp + diag(r.val) * l.ddp};
}

/// Left-multiply a ThermoVector instance by a Vector instance
template<typename V, typename T, typename P>
auto operator%(const Vector& l, const ThermoVectorBase<V,T,P>& r) -> ThermoVectorBase<decltype(diag(l) * r.val), decltype(diag(l) * r.ddt), decltype(diag(l) * r.ddp)>
{
    return {diag(l) * r.val,
            diag(l) * r.ddt,
            diag(l) * r.ddp};
}

/// Right-multiply a ThermoVector instance by a Vector instance
template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator%(const ThermoVectorBase<VL,TL,PL>& l, const Vector& r) -> decltype(r % l)
{
    return r % l;
}

/// Divide a ThermoVector instance by another
template<typename VL, typename TL, typename PL, typename VR, typename TR, typename PR>
auto operator/(const ThermoVectorBase<VL,TL,PL>& l, const ThermoVectorBase<VR,TR,PR>& r) -> ThermoVector
{
    const Vector tmp = 1.0/(r.val % r.val);
    return {l.val/r.val,
            diag(tmp) * (diag(r.val) * l.ddt - diag(l.val) * r.ddt),
            diag(tmp) * (diag(r.val) * l.ddp - diag(l.val) * r.ddp)};
}

/// Right-divide a ThermoVector instance by a ThermoScalar
template<typename V, typename T, typename P, typename VR>
auto operator/(const ThermoVectorBase<V,T,P>& l, const ThermoScalarBase<VR>& r) -> ThermoVectorBase<decltype(l.val/r.val), decltype(double() * (l.ddt * r.val - l.val * r.ddt)), decltype(double() * (l.ddp * r.val - l.val * r.ddp))>
{
    const double tmp = 1.0/(r.val * r.val);
    return {l.val/r.val, tmp * (l.ddt * r.val - l.val * r.ddt), tmp * (l.ddp * r.val - l.val * r.ddp)};
}

/// Right-divide a ThermoVector instance by a scalar
template<typename V, typename T, typename P>
auto operator/(const ThermoVectorBase<V,T,P>& l, double r) -> decltype((1.0/r) * l)
{
    return (1.0/r) * l;
}

/// Return the square root of a ThermoScalar instance
template<typename V, typename T, typename P>
auto sqrt(const ThermoVectorBase<V,T,P>& l) -> ThermoVector
{
    const Vector tmp1 = sqrt(l.val);
    const Vector tmp2 = 0.5 * tmp1/l.val;
    return {tmp1, diag(tmp2) * l.ddt, diag(tmp2) * l.ddp};
}

/// Return the power of a ThermoScalar instance
template<typename V, typename T, typename P>
auto pow(const ThermoVectorBase<V,T,P>& l, double power) -> ThermoVector
{
    const Vector tmp1 = pow(l.val, power);
    const Vector tmp2 = power * tmp1/l.val;
    return {tmp1, diag(tmp2) * l.ddt, diag(tmp2) * l.ddp};
}

/// Return the natural exponential of a ThermoScalar instance
template<typename V, typename T, typename P>
auto exp(const ThermoVectorBase<V,T,P>& l) -> ThermoVector
{
    const Vector tmp = exp(l.val);
    return {tmp, diag(tmp) * l.ddt, diag(tmp) * l.ddp};
}

/// Return the natural log of a ThermoScalar instance
template<typename V, typename T, typename P>
auto log(const ThermoVectorBase<V,T,P>& l) -> ThermoVector
{
    const Vector tmp1 = log(l.val);
    const Vector tmp2 = 1.0/l.val;
    return {tmp1, diag(tmp2) * l.ddt, diag(tmp2) * l.ddp};
}

/// Return the log10 of a ThermoScalar instance
template<typename V, typename T, typename P>
auto log10(const ThermoVectorBase<V,T,P>& l) -> ThermoVector
{
    const double log10e = 0.4342944819032518;
    const Vector tmp1 = log10e*log(l.val);
    const Vector tmp2 = log10e/l.val;
    return {tmp1, diag(tmp2) * l.ddt, diag(tmp2) * l.ddp, diag(tmp2) * l.ddn};
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

/// Output a ThermoScalar instance
inline auto operator<<(std::ostream& out, const ThermoVector& vector) -> std::ostream&
{
    out << vector.val;
    return out;
}

} // namespace Reaktoro
