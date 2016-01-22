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

// C++ includes
#include <cmath>
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Common/Matrix.hpp>

namespace Reaktoro {

/// A template base class to represent a thermodynamic scalar and its partial derivatives.
/// A *thermodynamic property* is a quantity that depends on temperature and pressure.
/// @see ThermoScalar, ChemicalScalar, ThermoVector
template<typename V>
class ThermoScalarBase
{
public:
    /// The value of the thermodynamic property.
    V val;

    /// The partial temperature derivative of the thermodynamic property.
    V ddt;

    /// The partial pressure derivative of the thermodynamic property.
    V ddp;

    /// Construct a default ThermoScalar instance
    ThermoScalarBase()
    : ThermoScalarBase(0.0) {}

    /// Construct a custom ThermoScalarBase instance with given value only.
    /// @param val The value of the thermodynamic property
    explicit ThermoScalarBase(double val)
    : ThermoScalarBase(val, 0.0, 0.0) {}

    /// Construct a custom ThermoScalarBase instance with given value and derivatives.
    /// @param val The value of the thermodynamic property
    /// @param ddt The partial temperature derivative of the thermodynamic property
    /// @param ddp The partial pressure derivative of the thermodynamic property
    ThermoScalarBase(const V& val, const V& ddt, const V& ddp)
    : val(val), ddt(ddt), ddp(ddp) {}

    /// Construct a copy of a ThermoScalar instance.
    template<typename VR>
    ThermoScalarBase(const ThermoScalarBase<VR>& other)
    : val(other.val), ddt(other.ddt), ddp(other.ddp) {}

    /// Assign another ThermoScalarBase instance to this ThermoScalarBase instance.
    template<typename VR>
    auto operator=(const ThermoScalarBase<VR>& other) -> ThermoScalarBase&
    {
        val = other.val;
        ddt = other.ddt;
        ddp = other.ddp;
        return *this;
    }

    /// Assign a scalar to this ThermoScalarBase instance.
    auto operator=(double other) -> ThermoScalarBase&
    {
        val = other;
        ddt = 0.0;
        ddp = 0.0;
        return *this;
    }

    /// Assign-addition of a ThermoScalar instance
    template<typename VR>
    auto operator+=(const ThermoScalarBase<VR>& other) -> ThermoScalarBase&
    {
        val += other.val;
        ddt += other.ddt;
        ddp += other.ddp;
        return *this;
    }

    /// Assign-subtraction of a ThermoScalar instance
    template<typename VR>
    auto operator-=(const ThermoScalarBase<VR>& other) -> ThermoScalarBase&
    {
        val -= other.val;
        ddt -= other.ddt;
        ddp -= other.ddp;
        return *this;
    }

    /// Assign-multiplication of a ThermoScalar instance
    template<typename VR>
    auto operator*=(const ThermoScalarBase<VR>& other) -> ThermoScalarBase&
    {
        ddt  = ddt * other.val + val * other.ddt;
        ddp  = ddp * other.val + val * other.ddp;
        val *= other.val;
        return *this;
    }

    /// Assign-division of a ThermoScalar instance
    template<typename VR>
    auto operator/=(const ThermoScalarBase<VR>& other) -> ThermoScalarBase&
    {
        const double tmp1 = 1.0/other.val;
        const double tmp2 = tmp1 * tmp1;
        ddt  = (ddt * other.val - val * other.ddt) * tmp2;
        ddp  = (ddp * other.val - val * other.ddp) * tmp2;
        val *= tmp1;
        return *this;
    }

    /// Assign-addition of a scalar
    auto operator+=(double other) -> ThermoScalarBase&
    {
        val += other;
        return *this;
    }

    /// Assign-subtraction of a scalar
    auto operator-=(double other) -> ThermoScalarBase&
    {
        val -= other;
        return *this;
    }

    /// Assign-multiplication of a ThermoScalar instance
    auto operator*=(double other) -> ThermoScalarBase&
    {
        val *= other;
        ddt *= other;
        ddp *= other;
        return *this;
    }

    /// Assign-division of a ThermoScalar instance
    auto operator/=(double other) -> ThermoScalarBase&
    {
        *this *= 1.0/other;
        return *this;
    }

    /// Explicitly converts this ThermoScalar instance into a double.
    explicit operator double()
    {
        return val;
    }
};

/// A type that defines a scalar thermo property.
/// A thermo property means here any property that depends on
/// temperature and pressure. A ThermoScalar instance not only holds
/// the value of the thermo property, but also is partial
/// temperature and pressure derivatives.
/// @see ChemicalVector
using ThermoScalar = ThermoScalarBase<double>;

/// A type that describes temperature in units of K
class Temperature : public ThermoScalar
{
public:
    /// Construct a default Temperature instance
    Temperature() : Temperature(0.0) {}

    /// Construct a Temperature instance with given value
    Temperature(double val) : ThermoScalarBase(val, 1.0, 0.0) {}

    /// Convert this Temperature instance into a double
//    operator double() { return val; }
};

/// A type that describes pressure in units of Pa
class Pressure : public ThermoScalar
{
public:
    /// Construct a default Pressure instance
    Pressure() : Pressure(0.0) {}

    /// Construct a Pressure instance with given value
    Pressure(double val) : ThermoScalarBase(val, 0.0, 1.0) {}

    /// Convert this Pressure instance into a double
//    operator double() { return val; }
};

/// Unary addition operator for a ThermoScalar instance
template<typename V>
inline auto operator+(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    return l;
}

/// Add two ThermoScalar instances
template<typename VL, typename VR>
inline auto operator+(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> ThermoScalarBase<double>
{
    return {l.val + r.val, l.ddt + r.ddt, l.ddp + r.ddp};
}

template<typename V>
inline auto operator+(double l, const ThermoScalarBase<V>& r) -> ThermoScalarBase<double>
{
    return {l + r.val, r.ddt, r.ddp};
}

template<typename V>
inline auto operator+(const ThermoScalarBase<V>& l, double r) -> ThermoScalarBase<double>
{
    return r + l;
}

/// Unary subtraction operator for a ThermoScalar instance
template<typename V>
inline auto operator-(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    return {-l.val, -l.ddt, -l.ddp};
}

/// Subtract two ThermoScalar instances
template<typename VL, typename VR>
inline auto operator-(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> ThermoScalarBase<double>
{
    return {l.val - r.val, l.ddt - r.ddt, l.ddp - r.ddp};
}

/// Right-subtract a ThermoScalar instance by a scalar
template<typename V>
inline auto operator-(const ThermoScalarBase<V>& l, double r) -> ThermoScalarBase<double>
{
    return {l.val - r, l.ddt, l.ddp};
}

/// Left-subtract a ThermoScalar instance by a scalar
template<typename V>
inline auto operator-(double l, const ThermoScalarBase<V>& r) -> ThermoScalarBase<double>
{
    return {l - r.val, -r.ddt, -r.ddp};
}

/// Multiply two ThermoScalar instances
template<typename VL, typename VR>
inline auto operator*(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> ThermoScalarBase<double>
{
    return {l.val * r.val, l.val * r.ddt + l.ddt * r.val, l.val * r.ddp + l.ddp * r.val};
}

/// Left-multiply a ThermoScalar instance by a scalar
template<typename V>
inline auto operator*(double l, const ThermoScalarBase<V>& r) -> ThermoScalarBase<double>
{
    return {l * r.val, l * r.ddt, l * r.ddp};
}

/// Right-multiply a ThermoScalar instance by a scalar
template<typename V>
inline auto operator*(const ThermoScalarBase<V>& l, double r) -> ThermoScalarBase<double>
{
    return r * l;
}

/// Divide a ThermoScalar instance by another
template<typename VL, typename VR>
inline auto operator/(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> ThermoScalarBase<double>
{
    const double tmp1 = 1.0/r.val;
    const double tmp2 = tmp1 * tmp1;
    return {tmp1 * l.val, tmp2 * (l.ddt * r.val - l.val * r.ddt), tmp2 * (l.ddp * r.val - l.val * r.ddp)};
}

/// Left-divide a ThermoScalar instance by a scalar
template<typename V>
inline auto operator/(double l, const ThermoScalarBase<V>& r) -> ThermoScalarBase<double>
{
    const double tmp1 = 1.0/r.val;
    const double tmp2 = -l*tmp1*tmp1;
    return {tmp1 * l, tmp2 * r.ddt, tmp2 * r.ddp};
}

/// Right-divide a ThermoScalar instance by a scalar
template<typename V>
inline auto operator/(const ThermoScalarBase<V>& l, double r) -> ThermoScalarBase<double>
{
    return (1.0/r) * l;
}

/// Return the square root of a ThermoScalar instance
template<typename V>
inline auto sqrt(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    const double tmp1 = std::sqrt(l.val);
    const double tmp2 = 0.5 * tmp1/l.val;
    return {tmp1, tmp2 * l.ddt, tmp2 * l.ddp};
}

/// Return the power of a ThermoScalar instance
template<typename V>
inline auto pow(const ThermoScalarBase<V>& l, double power) -> ThermoScalarBase<double>
{
    const double tmp1 = std::pow(l.val, power);
    const double tmp2 = power * tmp1/l.val;
    return {tmp1, tmp2 * l.ddt, tmp2 * l.ddp};
}

/// Return the power of a ThermoScalar instance
template<typename VL, typename VR>
inline auto pow(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& power) -> ThermoScalarBase<double>
{
    const double logl = std::log(l.val);
    const double powl = std::pow(l.val, power.val);
    const double tmp = power.val/l.val;
    return {powl, powl * (logl * power.ddt + tmp * l.ddt), powl * (logl * power.ddp + tmp * l.ddp)};
}

/// Return the natural exponential of a ThermoScalar instance
template<typename V>
inline auto exp(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    const double tmp = std::exp(l.val);
    return {tmp, tmp * l.ddt, tmp * l.ddp};
}

/// Return the natural log of a ThermoScalar instance
template<typename V>
inline auto log(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    const double tmp1 = std::log(l.val);
    const double tmp2 = 1.0/l.val;
    return {tmp1, tmp2 * l.ddt, tmp2 * l.ddp};
}

/// Return the log10 of a ThermoScalar instance
template<typename V>
inline auto log10(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    const double ln10 = 2.302585092994046;
    return log(l)/ln10;
}

/// Return true if a ThermoScalar instance is less than another
template<typename VL, typename VR>
inline auto operator<(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val < r.val;
}

/// Return true if a ThermoScalar instance is less or equal than another
template<typename VL, typename VR>
inline auto operator<=(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val <= r.val;
}

/// Return true if a ThermoScalar instance is greater than another
template<typename VL, typename VR>
inline auto operator>(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val > r.val;
}

/// Return true if a ThermoScalar instance is greater or equal than another
template<typename VL, typename VR>
inline auto operator>=(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val >= r.val;
}

/// Return true if a ThermoScalar instance is equal to another
template<typename VL, typename VR>
inline auto operator==(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val == r.val;
}

/// Return true if a ThermoScalar instance is not equal to another
template<typename VL, typename VR>
inline auto operator!=(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val != r.val;
}

/// Return true if a scalar is less than a ThermoScalar instance
template<typename V>
inline auto operator<(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l < r.val;
}

/// Return true if a ThermoScalar instance is less than a scalar
template<typename V>
inline auto operator<(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val < r;
}

/// Return true if a scalar is less or equal than a ThermoScalar instance
template<typename V>
inline auto operator<=(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l <= r.val;
}

/// Return true if a ThermoScalar instance is less or equal than a scalar
template<typename V>
inline auto operator<=(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val <= r;
}

/// Return true if a scalar is greater than a ThermoScalar instance
template<typename V>
inline auto operator>(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l > r.val;
}

/// Return true if a ThermoScalar is greater than a scalar
template<typename V>
inline auto operator>(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val > r;
}

/// Return true if a scalar is greater or equal than a ThermoScalar instance
template<typename V>
inline auto operator>=(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l >= r.val;
}

/// Return true if a ThermoScalar instance is greater or equal than a scalar
template<typename V>
inline auto operator>=(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val >= r;
}

/// Return true if a scalar is equal to a ThermoScalar instance
template<typename V>
inline auto operator==(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l == r.val;
}

/// Return true if a ThermoScalar instance is equal to a scalar
template<typename V>
inline auto operator==(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val == r;
}

/// Return true if a scalar is not equal to a ThermoScalar instance
template<typename V>
inline auto operator!=(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l != r.val;
}

/// Return true if a ThermoScalar instance is not equal to a scalar
template<typename V>
inline auto operator!=(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val != r;
}

/// Output a ThermoScalar instance
template<typename V>
inline auto operator<<(std::ostream& out, const ThermoScalarBase<V>& scalar) -> std::ostream&
{
    out << scalar.val;
    return out;
}

} // namespace Reaktoro
