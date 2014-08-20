/*
 * Constant.hpp
 *
 *  Created on: 28 Oct 2013
 *      Author: allan
 */

#pragma once

// Units++ includes
#include "ConversionUtils.hpp"
#include "StringUtils.hpp"
#include "Unit.hpp"

namespace units {

template<typename Unit>
class Constant;

using Mass                = Constant<g>;
using Length              = Constant<m>;
using Time                = Constant<s>;
using Current             = Constant<A>;
using Temperature         = Constant<K>;
using Amount              = Constant<mol>;
using LuminousIntensity   = Constant<cd>;
using Angle               = Constant<radian>;
using SolidAngle          = Constant<sr>;
using Area                = Constant<m2>;
using Volume              = Constant<m3>;
using Velocity            = Constant<mps>;
using Acceleration        = Constant<mpss>;
using Frequency           = Constant<Hz>;
using FlowRate            = Constant<gps>;
using Force               = Constant<N>;
using Pressure            = Constant<Pa>;
using Energy              = Constant<J>;
using Power               = Constant<W>;
using Charge              = Constant<C>;
using Voltage             = Constant<V>;
using Resistance          = Constant<ohm>;
using Conductance         = Constant<S>;
using Capacitance         = Constant<F>;
using MagneticFlux        = Constant<Wb>;
using MagneticFluxDensity = Constant<T>;
using Inductance          = Constant<H>;
using Molality            = Constant<molal>;

template<typename Unit>
class Constant
{
public:
    Constant()
    : val(0)
    {}

    template<typename OtherUnit>
    Constant(const Constant<OtherUnit>& other)
    : val(other.in(Unit()))
    {}

    Constant(double val, const std::string& unit)
    : val(convert(val, unit, symbol(Unit())))
    {}

    explicit Constant(double val)
    : val(val)
    {}

    double value() const
    {
        return val;
    }

    double operator()() const
    {
        return val;
    }

    inline auto operator+() -> Constant&
    {
        return *this;
    }

    inline auto operator-() -> Constant&
    {
        val *= -1;
        return *this;
    }

    template<typename OtherUnit>
    inline auto in(OtherUnit) const -> double
    {
        return convert(Unit(), OtherUnit(), val);
    }

    template<typename OtherUnit>
    inline auto operator+=(const Constant<OtherUnit>& c) -> Constant&
    {
        val += c.in(Unit());
        return *this;
    }

    template<typename OtherUnit>
    inline auto operator-=(const Constant<OtherUnit>& c) -> Constant&
    {
        val -= c.in(Unit());
        return *this;
    }

    inline auto operator*=(double scalar) -> Constant&
    {
        val *= scalar;
        return *this;
    }

    inline auto operator/=(double scalar) -> Constant&
    {
        val /= scalar;
        return *this;
    }

private:
    double val;
};

template<typename U>
std::ostream& operator<<(std::ostream& out, const Constant<U>& c)
{
    out << c.value() << " " << symbol(U()); return out;
}

template<typename U>
inline constexpr auto operator*(double val, Unit<U> unit) -> Constant<U>
{
    return Constant<U>(val);
}

template<typename U1, typename U2>
inline constexpr auto operator*(const Constant<U1>& c, Unit<U2> unit) -> Constant<Mult<U1, U2>>
{
    return Constant<Mult<U1, U2>>(c.value());
}

template<typename U1, typename U2>
inline constexpr auto operator/(const Constant<U1>& c, Unit<U2> unit) -> Constant<Mult<U1, Inv<U2>>>
{
    return Constant<Mult<U1, Inv<U2>>>(c.value());
}

template<typename U1, typename U2>
inline constexpr auto operator+(const Constant<U1>& c1, const Constant<U2>& c2) -> Constant<U1>
{
    return Constant<U1>(c1.value() + c2.in(U1()));
}

template<typename U1, typename U2>
inline constexpr auto operator-(const Constant<U1>& c1, const Constant<U2>& c2) -> Constant<U1>
{
    return Constant<U1>(c1.value() - c2.in(U1()));
}

template<typename U1, typename U2>
inline constexpr auto operator*(const Constant<U1>& c1, const Constant<U2>& c2) -> Constant<Mult<U1, U2>>
{
    return Constant<Mult<U1, U2>>(c1.value() * c2.value());
}

template<typename U1, typename U2>
inline constexpr auto operator/(const Constant<U1>& c1, const Constant<U2>& c2) -> Constant<Mult<U1, Inv<U2>>>
{
    return Constant<Mult<U1, Inv<U2>>>(c1.value() / c2.value());
}

template<typename U>
inline constexpr auto operator*(double val, const Constant<U>& c) -> Constant<U>
{
    return Constant<U>(val * c.value());
}

template<typename U>
inline constexpr auto operator*(const Constant<U>& c, double val) -> Constant<U>
{
    return Constant<U>(val * c.value());
}

template<typename U>
inline constexpr auto operator/(const Constant<U>& c, double val) -> Constant<U>
{
    return Constant<U>(c.value() / val);
}

template<typename U1, typename U2>
inline constexpr auto operator==(const Constant<U1>& c1, const Constant<U2>& c2) -> bool
{
    return c1.value() == c2.in(U1());
}

template<typename U1, typename U2>
inline constexpr auto operator!=(const Constant<U1>& c1, const Constant<U2>& c2) -> bool
{
    return c1.value() != c2.in(U1());
}

template<typename U1, typename U2>
inline constexpr auto operator<(const Constant<U1>& c1, const Constant<U2>& c2) -> bool
{
    return c1.value() < c2.in(U1());
}

template<typename U1, typename U2>
inline constexpr auto operator<=(const Constant<U1>& c1, const Constant<U2>& c2) -> bool
{
    return c1.value() <= c2.in(U1());
}

template<typename U1, typename U2>
inline constexpr auto operator>(const Constant<U1>& c1, const Constant<U2>& c2) -> bool
{
    return c1.value() > c2.in(U1());
}

template<typename U1, typename U2>
inline constexpr auto operator>=(const Constant<U1>& c1, const Constant<U2>& c2) -> bool
{
    return c1.value() >= c2.in(U1());
}

} /* namespace units */
