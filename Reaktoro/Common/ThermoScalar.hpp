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

// C++ includes
#include <cmath>
#include <iostream>

// Reaktoro includes
#include <Reaktoro/Math/Matrix.hpp>

namespace Reaktoro {

// Forward declaration
template<typename V>
class ThermoScalarBase;

/// A type that defines a scalar thermo property.
/// A thermo property means here any property that depends on
/// temperature and pressure. A ThermoScalar instance not only holds
/// the value of the thermo property, but also is partial
/// temperature and pressure derivatives.
/// @see ChemicalVector
using ThermoScalar = ThermoScalarBase<double>;

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
    V ddT;

    /// The partial pressure derivative of the thermodynamic property.
    V ddP;

    /// Construct a default ThermoScalar instance
    ThermoScalarBase()
    : ThermoScalarBase(0.0) {}

    /// Construct a custom ThermoScalarBase instance with given value only.
    /// @param val The value of the thermodynamic property
    explicit ThermoScalarBase(double val)
    : ThermoScalarBase(val, 0.0, 0.0) {}

    /// Construct a custom ThermoScalarBase instance with given value and derivatives.
    /// @param val The value of the thermodynamic property
    /// @param ddT The partial temperature derivative of the thermodynamic property
    /// @param ddP The partial pressure derivative of the thermodynamic property
    ThermoScalarBase(const V& val, const V& ddT, const V& ddP)
    : val(val), ddT(ddT), ddP(ddP) {}

    /// Construct a copy of a ThermoScalar instance.
    template<typename VR>
    ThermoScalarBase(const ThermoScalarBase<VR>& other)
    : val(other.val), ddT(other.ddT), ddP(other.ddP) {}

    /// Assign another ThermoScalarBase instance to this ThermoScalarBase instance.
    template<typename VR>
    auto operator=(const ThermoScalarBase<VR>& other) -> ThermoScalarBase&
    {
        val = other.val;
        ddT = other.ddT;
        ddP = other.ddP;
        return *this;
    }

    /// Assign a scalar to this ThermoScalarBase instance.
    auto operator=(double other) -> ThermoScalarBase&
    {
        val = other;
        ddT = 0.0;
        ddP = 0.0;
        return *this;
    }

    /// Assign-addition of a ThermoScalar instance
    template<typename VR>
    auto operator+=(const ThermoScalarBase<VR>& other) -> ThermoScalarBase&
    {
        val += other.val;
        ddT += other.ddT;
        ddP += other.ddP;
        return *this;
    }

    /// Assign-subtraction of a ThermoScalar instance
    template<typename VR>
    auto operator-=(const ThermoScalarBase<VR>& other) -> ThermoScalarBase&
    {
        val -= other.val;
        ddT -= other.ddT;
        ddP -= other.ddP;
        return *this;
    }

    /// Assign-multiplication of a ThermoScalar instance
    template<typename VR>
    auto operator*=(const ThermoScalarBase<VR>& other) -> ThermoScalarBase&
    {
        ddT  = ddT * other.val + val * other.ddT;
        ddP  = ddP * other.val + val * other.ddP;
        val *= other.val;
        return *this;
    }

    /// Assign-division of a ThermoScalar instance
    template<typename VR>
    auto operator/=(const ThermoScalarBase<VR>& other) -> ThermoScalarBase&
    {
        const double tmp1 = 1.0/other.val;
        const double tmp2 = tmp1 * tmp1;
        ddT  = (ddT * other.val - val * other.ddT) * tmp2;
        ddP  = (ddP * other.val - val * other.ddP) * tmp2;
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
        ddT *= other;
        ddP *= other;
        return *this;
    }

    /// Assign-division of a ThermoScalar instance
    auto operator/=(double other) -> ThermoScalarBase&
    {
        *this *= 1.0/other;
        return *this;
    }

    /// Explicitly converts this ThermoScalar instance into a double.
    explicit operator double() const
    {
        return val;
    }
};

/// A type that describes temperature in units of K
class Temperature : public ThermoScalar
{
public:
    /// Construct a default Temperature instance
    Temperature() : Temperature(0.0) {}

    /// Construct a Temperature instance with given value
    Temperature(double val) : ThermoScalarBase(val, 1.0, 0.0) {}

    /// Converts this Temperature instance into a double.
    explicit operator double() const
    {
        return val;
    }
};

/// A type that describes pressure in units of Pa
class Pressure : public ThermoScalar
{
public:
    /// Construct a default Pressure instance
    Pressure() : Pressure(0.0) {}

    /// Construct a Pressure instance with given value
    Pressure(double val) : ThermoScalarBase(val, 0.0, 1.0) {}

    /// Converts this Pressure instance into a double.
    explicit operator double() const
    {
        return val;
    }
};

template<typename V>
auto operator+(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    return l;
}

template<typename VL, typename VR>
auto operator+(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> ThermoScalarBase<double>
{
    return {l.val + r.val, l.ddT + r.ddT, l.ddP + r.ddP};
}

template<typename V>
auto operator+(double l, const ThermoScalarBase<V>& r) -> ThermoScalarBase<double>
{
    return {l + r.val, r.ddT, r.ddP};
}

template<typename V>
auto operator+(const ThermoScalarBase<V>& l, double r) -> ThermoScalarBase<double>
{
    return r + l;
}

template<typename V>
auto operator-(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    return {-l.val, -l.ddT, -l.ddP};
}

template<typename VL, typename VR>
auto operator-(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> ThermoScalarBase<double>
{
    return {l.val - r.val, l.ddT - r.ddT, l.ddP - r.ddP};
}

template<typename V>
auto operator-(const ThermoScalarBase<V>& l, double r) -> ThermoScalarBase<double>
{
    return {l.val - r, l.ddT, l.ddP};
}

template<typename V>
auto operator-(double l, const ThermoScalarBase<V>& r) -> ThermoScalarBase<double>
{
    return {l - r.val, -r.ddT, -r.ddP};
}

template<typename VL, typename VR>
auto operator*(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> ThermoScalarBase<double>
{
    return {l.val * r.val, l.val * r.ddT + l.ddT * r.val, l.val * r.ddP + l.ddP * r.val};
}

template<typename V>
auto operator*(double l, const ThermoScalarBase<V>& r) -> ThermoScalarBase<double>
{
    return {l * r.val, l * r.ddT, l * r.ddP};
}

template<typename V>
auto operator*(const ThermoScalarBase<V>& l, double r) -> ThermoScalarBase<double>
{
    return r * l;
}

template<typename VL, typename VR>
auto operator/(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> ThermoScalarBase<double>
{
    const double tmp1 = 1.0/r.val;
    const double tmp2 = tmp1 * tmp1;
    return {tmp1 * l.val, tmp2 * (l.ddT * r.val - l.val * r.ddT), tmp2 * (l.ddP * r.val - l.val * r.ddP)};
}

template<typename V>
auto operator/(double l, const ThermoScalarBase<V>& r) -> ThermoScalarBase<double>
{
    const double tmp1 = 1.0/r.val;
    const double tmp2 = -l*tmp1*tmp1;
    return {tmp1 * l, tmp2 * r.ddT, tmp2 * r.ddP};
}

template<typename V>
auto operator/(const ThermoScalarBase<V>& l, double r) -> ThermoScalarBase<double>
{
    return (1.0/r) * l;
}

template<typename VL, typename VR>
auto operator<(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val < r.val;
}

template<typename VL, typename VR>
auto operator<=(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val <= r.val;
}

template<typename VL, typename VR>
auto operator>(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val > r.val;
}

template<typename VL, typename VR>
auto operator>=(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val >= r.val;
}

template<typename VL, typename VR>
auto operator==(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val == r.val;
}

template<typename VL, typename VR>
auto operator!=(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& r) -> bool
{
    return l.val != r.val;
}

template<typename V>
auto operator<(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l < r.val;
}

template<typename V>
auto operator<(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val < r;
}

template<typename V>
auto operator<=(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l <= r.val;
}

template<typename V>
auto operator<=(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val <= r;
}

template<typename V>
auto operator>(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l > r.val;
}

template<typename V>
auto operator>(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val > r;
}

template<typename V>
auto operator>=(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l >= r.val;
}

template<typename V>
auto operator>=(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val >= r;
}

template<typename V>
auto operator==(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l == r.val;
}

template<typename V>
auto operator==(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val == r;
}

template<typename V>
auto operator!=(double l, const ThermoScalarBase<V>& r) -> bool
{
    return l != r.val;
}

template<typename V>
auto operator!=(const ThermoScalarBase<V>& l, double r) -> bool
{
    return l.val != r;
}

template<typename V>
auto operator<<(std::ostream& out, const ThermoScalarBase<V>& scalar) -> std::ostream&
{
    out << scalar.val;
    return out;
}

template<typename V>
auto abs(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    if(l.val == 0.0) return {};
    const double tmp1 = std::abs(l.val);
    const double tmp2 = l.val/tmp1;
    return {tmp1, tmp2 * l.ddT, tmp2 * l.ddP};
}

template<typename V>
auto sqrt(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    if(l.val == 0.0) return {};
    const double tmp1 = std::sqrt(l.val);
    const double tmp2 = 0.5 * tmp1/l.val;
    return {tmp1, tmp2 * l.ddT, tmp2 * l.ddP};
}

template<typename V>
auto pow(const ThermoScalarBase<V>& l, double power) -> ThermoScalarBase<double>
{
    if(l.val == 0.0) return {};
    const double tmp1 = std::pow(l.val, power - 1);
    const double tmp2 = power * tmp1;
    return {tmp1 * l.val, tmp2 * l.ddT, tmp2 * l.ddP};
}

template<typename VL, typename VR>
auto pow(const ThermoScalarBase<VL>& l, const ThermoScalarBase<VR>& power) -> ThermoScalarBase<double>
{
    if(l.val == 0.0) return {};
    const double logl = std::log(l.val);
    const double powl = std::pow(l.val, power.val);
    const double tmp = power.val/l.val;
    return {powl, powl * (logl * power.ddT + tmp * l.ddT), powl * (logl * power.ddP + tmp * l.ddP)};
}

template<typename V>
auto exp(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    const double tmp = std::exp(l.val);
    return {tmp, tmp * l.ddT, tmp * l.ddP};
}

template<typename V>
auto log(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    const double tmp1 = std::log(l.val);
    const double tmp2 = 1.0/l.val;
    return {tmp1, tmp2 * l.ddT, tmp2 * l.ddP};
}

template<typename V>
auto log10(const ThermoScalarBase<V>& l) -> ThermoScalarBase<double>
{
    const double ln10 = 2.302585092994046;
    return log(l)/ln10;
}

} // namespace Reaktoro
