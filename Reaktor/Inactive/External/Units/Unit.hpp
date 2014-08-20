/*
 * Unit.hpp
 *
 *  Created on: 28 Oct 2013
 *      Author: allan
 */

#pragma once

// C++ includes
#include <cstdint>
#include <string>

namespace units {

using Integer = std::intmax_t;

template<typename U>
struct Unit
{};

template<typename U, Integer num, Integer den = 1>
struct Add : Unit<Add<U, num, den>>
{};

template<typename U, Integer num, Integer den = 1>
struct Scale : Unit<Scale<U, num, den>>
{};

template<typename U>
struct Inv : Unit<Inv<U>>
{};

template<typename U1, typename U2>
struct Mult : Unit<Mult<U1, U2>>
{};

template<Integer num, Integer den = 1>
struct Ratio {};

template<typename U, Integer num, Integer den>
auto operator+(const Unit<U>& u, const Ratio<num, den>& r) -> Add<U, num, den>
{ return {}; }

template<typename U, Integer num, Integer den>
auto operator+(const Ratio<num, den>& r, const Unit<U>& u) -> Add<U, num, den>
{ return {}; }

template<typename U, Integer num, Integer den>
auto operator-(const Unit<U>& u, const Ratio<num, den>& r) -> Add<U, -num, den>
{ return {}; }

template<typename A, Integer num, Integer den>
auto operator*(const Ratio<num, den>& r, const Unit<A>& a) -> Scale<A, num, den>
{ return {}; }

template<typename A, Integer num, Integer den>
auto operator*(const Unit<A>& a, const Ratio<num, den>& r) -> Scale<A, num, den>
{ return {}; }

template<typename A, Integer num, Integer den>
auto operator/(const Ratio<num, den>& r, const Unit<A>& a) -> Scale<Inv<A>, num, den>
{ return {}; }

template<typename A, Integer num, Integer den>
auto operator/(const Unit<A>& a, const Ratio<num, den>& r) -> Scale<A, den, num>
{ return {}; }

template<typename A, typename B>
auto operator*(const Unit<A>& a, const Unit<B>& b) -> Mult<A, B>
{ return {}; }

template<typename A, typename B>
auto operator/(const Unit<A>& a, const Unit<B>& b) -> Mult<A, Inv<B>>
{ return {}; }

// Convenience Types
template<typename U> using atto  = Scale<U, 1, 1000000000000000000>; // atto  = 10^-18
template<typename U> using femto = Scale<U, 1, 1000000000000000>;    // femto = 10^-15
template<typename U> using pico  = Scale<U, 1, 1000000000000>;       // pico  = 10^-12
template<typename U> using nano  = Scale<U, 1, 1000000000>;          // nano  = 10^-9
template<typename U> using micro = Scale<U, 1, 1000000>;             // micro = 10^-6
template<typename U> using milli = Scale<U, 1, 1000>;                // milli = 10^-3
template<typename U> using centi = Scale<U, 1, 100>;                 // centi = 10^-2
template<typename U> using deci  = Scale<U, 1, 10>;                  // deci  = 10^-1
template<typename U> using deca  = Scale<U, 10, 1>;                  // deca  = 10^1
template<typename U> using hecto = Scale<U, 100, 1>;                 // hecto = 10^2
template<typename U> using kilo  = Scale<U, 1000, 1>;                // kilo  = 10^3
template<typename U> using mega  = Scale<U, 1000000, 1>;             // mega  = 10^6
template<typename U> using giga  = Scale<U, 1000000000, 1>;          // giga  = 10^9
template<typename U> using tera  = Scale<U, 1000000000000, 1>;       // tera  = 10^12
template<typename U> using peta  = Scale<U, 1000000000000000, 1>;    // peta  = 10^15
template<typename U> using exa   = Scale<U, 1000000000000000000, 1>; // exa   = 10^18
using one = Ratio<1>;

/**
 * Definition of the base units
 */
namespace base {

struct Gram      : Unit<Gram>      {}; // Mass
struct Meter     : Unit<Meter>     {}; // Length
struct Second    : Unit<Second>    {}; // Time
struct Ampere    : Unit<Ampere>    {}; // Current
struct Kelvin    : Unit<Kelvin>    {}; // Temperature
struct Mol       : Unit<Mol>       {}; // Amount of Substance
struct Candela   : Unit<Candela>   {}; // Luminous Intensity
struct Radian    : Unit<Radian>    {}; // Angle
struct Steradian : Unit<Steradian> {}; // Solid Angle

} /* namespace base */

#define DefineUnit(name, expression) \
    using name = expression; \
    inline std::string symbol(name) { return #name; }

#define DefineAlias(alias, name) \
    using alias = name;

#define add(unit, num, den) \
    Add<unit, num, den>

#define scale(unit, num, den) \
    Scale<unit, num, den>

// Mass Units
DefineUnit(g            , base::Gram)
DefineUnit(kg           , kilo<g>)
DefineUnit(mg           , milli<g>)
DefineUnit(ug           , micro<g>)
DefineUnit(ng           , nano<g>)
DefineUnit(tonne        , kilo<kg>)
DefineUnit(lb           , scale(g, 45359237, 100000))
DefineUnit(amu          , scale(atto<kg>, 1660538921, 1))
DefineUnit(ton          , scale(lb, 2000, 1))
DefineAlias(gram        , g)
DefineAlias(grams       , g)
DefineAlias(kilogram    , g)
DefineAlias(kilograms   , g)
DefineAlias(lbm         , lb)

// Length Units
DefineUnit(m            , base::Meter)
DefineUnit(km           , kilo<m>)
DefineUnit(cm           , centi<m>)
DefineUnit(mm           , milli<m>)
DefineUnit(um           , micro<m>)
DefineUnit(nm           , nano<m>)
DefineUnit(Angstrom     , deci<nm>)
DefineUnit(in           , scale(mm, 254, 10))
DefineUnit(mil          , milli<in>)
DefineUnit(ft           , scale(in, 12, 1))
DefineUnit(yd           , scale(ft, 3, 1))
DefineUnit(mi           , scale(ft, 5280, 1))
DefineUnit(nauticalmile , scale(m, 1852, 1))
DefineUnit(lightyear    , scale(m, 9460730472580800, 1))
DefineUnit(au           , scale(m, 149597870700, 1))
DefineAlias(feet        , ft)
DefineAlias(foot        , ft)
DefineAlias(inch        , in)
DefineAlias(inches      , in)
DefineAlias(meter       , m)
DefineAlias(meters      , m)
DefineAlias(mile        , mi)
DefineAlias(miles       , mi)

// Time Units
DefineUnit(s            , base::Second)
DefineUnit(ms           , milli<s>)
DefineUnit(us           , micro<s>)
DefineUnit(ns           , nano<s>)
DefineUnit(minute       , scale(s, 60, 1))
DefineUnit(hour         , scale(s, 3600, 1))
DefineUnit(day          , scale(hour, 24, 1))
DefineUnit(week         , scale(day, 7, 1))
DefineUnit(month        , scale(day, 3043685, 100000))
DefineUnit(year         , scale(day, 3652422, 10000))
DefineAlias(second      , s)
DefineAlias(seconds     , s)
DefineAlias(minutes     , minute)
DefineAlias(hours       , hour)
DefineAlias(days        , day)
DefineAlias(weeks       , week)
DefineAlias(months      , month)
DefineAlias(years       , year)

// Current Units
DefineUnit(A            , base::Ampere)
DefineUnit(MA           , mega<A>)
DefineUnit(kA           , kilo<A>)
DefineUnit(mA           , milli<A>)
DefineUnit(uA           , micro<A>)
DefineUnit(nA           , nano<A>)
DefineUnit(pA           , pico<A>)
DefineAlias(ampere      , A)
DefineAlias(amperes     , A)

// Temperature Units
DefineUnit(K            , base::Kelvin)
DefineUnit(degC         , add(K, -27315, 100))
DefineUnit(degF         , add(scale(degC, 9, 5), 32, 1))
DefineUnit(degR         , add(degF, 45967, 100))
DefineAlias(kelvin      , K)
DefineAlias(celsius     , degC)
DefineAlias(fahrenheit  , degF)
DefineAlias(rankine     , degR)

// Amount of Substance Units
DefineUnit(mol          , base::Mol)
DefineUnit(mmol         , milli<mol>)
DefineAlias(mole        , mol)
DefineAlias(moles       , mol)

// Luminous Intensity Units
DefineUnit(cd           , base::Candela)
DefineAlias(candela     , cd)

// Angle Units
DefineUnit(radian       , base::Radian)
DefineUnit(mradian      , milli<radian>)
DefineUnit(uradian      , micro<radian>)
DefineUnit(deg          , scale(radian, 1745329251994330, 100000000000000000))
DefineAlias(degree      , deg)
DefineAlias(degrees     , deg)

// Solid Angle Units
DefineUnit(sr           , base::Steradian)
DefineAlias(steradian   , sr)

// Area Units
DefineUnit(m2           , decltype(m()*m()))
DefineUnit(mm2          , decltype(mm()*mm()))
DefineUnit(cm2          , decltype(cm()*cm()))
DefineUnit(in2          , decltype(in()*in()))
DefineUnit(ft2          , decltype(ft()*ft()))
DefineUnit(hectare      , scale(m2, 10000, 1))
DefineUnit(acre         , scale(m2, 40468726, 10000))
DefineUnit(barn         , deci<nano<atto<m2>>>)
DefineAlias(hectares    , hectare)
DefineAlias(acres       , acre)
DefineAlias(barns       , barn)

// Volume Units
DefineUnit(m3           , decltype(m()*m()*m()))
DefineUnit(cm3          , decltype(cm()*cm()*cm()))
DefineUnit(mm3          , decltype(mm()*mm()*mm()))
DefineUnit(in3          , decltype(in()*in()*in()))
DefineUnit(ft3          , decltype(ft()*ft()*ft()))
DefineUnit(l            , milli<m3>)
DefineUnit(ml           , milli<l>)
DefineUnit(ul           , micro<l>)
DefineUnit(gal          , scale(l, 454609, 100000))
DefineUnit(quart        , scale(gal, 1, 4))
DefineAlias(cc          , cm3)
DefineAlias(gallon      , gal)
DefineAlias(gallons     , gal)
DefineAlias(liter       , l)
DefineAlias(liters      , l)

// Velocity Units
DefineUnit(mps          , decltype(m()/s()))
DefineUnit(cmps         , decltype(cm()/s()))
DefineUnit(mph          , decltype(mile()/hour()))
DefineUnit(knot         , decltype(nauticalmile()/hour()))
DefineUnit(fps          , decltype(foot()/s()))

// Acceleration Units
DefineUnit(mpss         , decltype(mps()/s()))
DefineUnit(fpss         , decltype(fps()/s()))
DefineUnit(gravity      , scale(mpss, 980665, 100000))

// Frequency Units
DefineUnit(Hz           , decltype(one()/second()))
DefineUnit(kHz          , kilo<Hz>)
DefineUnit(MHz          , mega<Hz>)
DefineUnit(GHz          , giga<Hz>)
DefineUnit(rpm          , decltype(one()/minute()))
DefineAlias(hertz       , Hz)

// Flow Rate Units
DefineUnit(gps          , decltype(gal()/s()))
DefineUnit(gpm          , decltype(gal()/minute()))
DefineUnit(gph          , decltype(gal()/hour()))
DefineUnit(cfs          , decltype(ft3()/s()))
DefineUnit(cfm          , decltype(ft3()/minute()))
DefineUnit(lps          , decltype(liter()/s()))
DefineUnit(lpm          , decltype(liter()/minute()))

// Force Units
DefineUnit(N            , decltype(kg()*m()/s()/s()))
DefineUnit(lbf          , scale(N, 44482216, 10000000))
DefineUnit(dyne         , scale(N, 1, 100000))

DefineAlias(dynes       , dyne)
DefineAlias(newton      , N)
DefineAlias(newtons     , N)

// Pressure Units
DefineUnit(Pa           , decltype(N()/m2()))
DefineUnit(kPa          , kilo<Pa>)
DefineUnit(MPa          , mega<Pa>)
DefineUnit(GPa          , giga<Pa>)
DefineUnit(atm          , scale(Pa, 101325, 1))
DefineUnit(mmHg         , scale(Pa, 13332239, 100000))
DefineUnit(inHg         , scale(Pa, 33863886, 10000))
DefineUnit(psi          , decltype(lbf()/in2()))
DefineUnit(kpsi         , kilo<psi>)
DefineUnit(Mpsi         , mega<psi>)
DefineUnit(psf          , decltype(lbf()/ft2()))
DefineUnit(bar          , scale(Pa, 100000, 1))
DefineUnit(torr         , scale(Pa, 133322368, 1000000))
DefineUnit(inH2O        , scale(Pa, 24908891, 100000))
DefineUnit(ftH2O        , scale(inH2O, 12, 1))
DefineAlias(pascal      , Pa)
DefineAlias(pascals     , Pa)

// Energy Units
DefineUnit(J            , decltype(N()*m()))
DefineUnit(erg          , decltype(dyne()*cm()))
DefineUnit(cal          , scale(J, 4184, 1000))
DefineUnit(kcal         , kilo<cal>)
DefineUnit(eV           , scale(deci<atto<J>>, 16021765, 10000000))
DefineUnit(keV          , kilo<eV>)
DefineUnit(MeV          , mega<eV>)
DefineUnit(GeV          , giga<eV>)
DefineUnit(btu          , scale(J, 10550559, 10000))
DefineAlias(calorie     , cal)
DefineAlias(calories    , cal)
DefineAlias(joule       , J)
DefineAlias(joules      , J)

// Power Units
DefineUnit(W            , decltype(J()/s()))
DefineAlias(watts       , W)

// Charge Units
DefineUnit(C            , decltype(A()*s()))
DefineAlias(coulomb     , C)
DefineAlias(coulombs    , C)

// Voltage Units
DefineUnit(V            , decltype(J()/C()))
DefineAlias(volt        , V)
DefineAlias(volts       , V)

// Resistance Units
DefineUnit(ohm          , decltype(V()/A()))
DefineAlias(ohms        , ohm)

// Conductance Units
DefineUnit(S            , decltype(one()/ohm()))
DefineAlias(siemens     , S)

// Capacitance Units
DefineUnit(F            , decltype(C()/V()))
DefineUnit(pF           , pico<F>)
DefineUnit(nF           , nano<F>)
DefineUnit(uF           , micro<F>)
DefineUnit(mF           , milli<F>)
DefineAlias(farad       , F)
DefineAlias(farads      , F)

// Magnetic Flux Units
DefineUnit(Wb           , decltype(V()*s()))
DefineAlias(weber       , Wb)
DefineAlias(webers      , Wb)

// Magnetic Flux Density Units
DefineUnit(T            , decltype(Wb()/m()/m()))
DefineUnit(gauss        , deci<milli<T>>)
DefineAlias(tesla       , T)
DefineAlias(teslas      , T)

// Inductance Units
DefineUnit(H            , decltype(Wb()/A()))
DefineUnit(mH           , milli<H>)
DefineUnit(uH           , micro<H>)
DefineAlias(henry       , H)

// Molality Units
DefineUnit(molal        , decltype(mol()/kg()))
DefineUnit(mmolal       , milli<molal>)
DefineUnit(umolal       , micro<molal>)

} /* namespace units */

#define unit(x) units::x()

#undef DefineUnit
#undef DefineAlias
#undef scale
#undef add
