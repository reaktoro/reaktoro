/*
 * ConversionUtils.hpp
 *
 *  Created on: 28 Oct 2013
 *      Author: allan
 */

#pragma once

// Units++ includes
#include "Unit.hpp"

namespace units {

template<typename U1, typename U2>
constexpr double convert(Unit<U1> from, Unit<U2> to, double value);

namespace internal {

template<typename...>
struct List {};

template<typename>
struct Size {};

template<typename, typename>
struct Join {};

template<typename, typename>
struct Concat {};

template<>
struct Size<List<>> { constexpr static int value = 0; };

template<typename Arg, typename... Args>
struct Size<List<Arg, Args...>> { constexpr static int value = 1 + Size<List<Args...>>::value; };

template<typename Entry, typename... Args>
struct Join<Entry, List<Args...>> { using type = List<Entry, Args...>; };

template<typename... Args1, typename... Args2>
struct Concat<List<Args1...>, List<Args2...>> { using type = List<Args1..., Args2...>; };

template<typename, typename>
struct Contains {};

template<typename Entry>
struct Contains<Entry, List<>> { constexpr static bool value = false; };

template<typename Entry, typename... Args>
struct Contains<Entry, List<Entry, Args...>> { constexpr static bool value = true; };

template<typename Entry, typename Arg, typename... Args>
struct Contains<Entry, List<Arg, Args...>> : Contains<Entry, List<Args...>> {};

template<typename, typename>
struct Contained {};

template<typename... Args1>
struct Contained<List<Args1...>, List<>> { constexpr static bool value = false; };

template<typename... Args2>
struct Contained<List<>, List<Args2...>> { constexpr static bool value = true; };

template<typename Arg1, typename... Args1, typename... Args2>
struct Contained<List<Arg1, Args1...>, List<Args2...>>
{
    constexpr static bool value =
        Contains<Arg1, List<Args2...>>::value and
            Contained<List<Args1...>, List<Args2...>>::value;
};

template<typename, typename>
struct Equal {};

template<typename... Args1, typename... Args2>
struct Equal<List<Args1...>, List<Args2...>>
{
    constexpr static bool res1 = Contained<List<Args1...>, List<Args2...>>::value;
    constexpr static bool res2 = Size<List<Args1...>>::value == Size<List<Args2...>>::value;
    constexpr static bool value = res1 and res2;
};

template<typename>
struct Unique {};

template<>
struct Unique<List<>> { using type = List<>; };

template<typename Arg, typename... Args>
struct Unique<List<Arg, Args...>>
{
    using res1 = typename Unique<List<Args...>>::type;
    using res2 = typename Join<Arg, res1>::type;
    using type = typename std::conditional<Contains<Arg, List<Args...>>::value, res1, res2>::type;
};

template<template <typename> class Function, typename List>
struct Map
{};

template<template <typename> class Function, typename Arg>
struct Map<Function, List<Arg>>
{
    using type = List<typename Function<Arg>::type>;
};

template<template <typename> class Function, typename Arg, typename... Args>
struct Map<Function, List<Arg, Args...>>
{
    using type = typename Concat<List<typename Function<Arg>::type>,
        typename Map<Function, List<Args...>>::type>::type;
};

template<template <typename> class Condition, typename ListType>
struct Remove
{};

template<template <typename> class Condition>
struct Remove<Condition, List<>>
{
    using type = List<>;
};

template<template <typename> class Condition, typename Arg, typename... Args>
struct Remove<Condition, List<Arg, Args...>>
{
    using cond1 = typename Remove<Condition, List<Args...>>::type;
    using cond2 = typename Concat<List<Arg>, typename Remove<Condition, List<Args...>>::type>::type;
    using type = typename std::conditional<Condition<Arg>::value, cond1, cond2>::type;
};

template<typename U>
struct BaseUnits
{
    template<typename W>
    struct ListUnits { using type = List<W>; };

    template<typename W>
    struct ListUnits<Unit<W>> { using type = typename ListUnits<W>::type; };

    template<typename W, Integer num, Integer den>
    struct ListUnits<Add<W, num, den>> { using type = typename ListUnits<W>::type; };

    template<typename W, Integer num, Integer den>
    struct ListUnits<Scale<W, num, den>> { using type = typename ListUnits<W>::type; };

    template<typename W>
    struct ListUnits<Inv<W>> { using type = typename ListUnits<W>::type; };

    template<typename W1, typename W2>
    struct ListUnits<Mult<W1, W2>>
    {
        using type = typename Concat<typename ListUnits<W1>::type, typename ListUnits<W2>::type>::type;
    };

    using type = typename Unique<typename ListUnits<U>::type>::type;
};

template<typename Base, typename U>
struct BaseDimension { constexpr static int value = 0; };

template<typename Base>
struct BaseDimension<Base, Base> { constexpr static int value = 1; };

template<typename Base, typename U, Integer num, Integer den>
struct BaseDimension<Base, Add<U, num, den>> : BaseDimension<Base, U> {};

template<typename Base, typename U, Integer num, Integer den>
struct BaseDimension<Base, Scale<U, num, den>> : BaseDimension<Base, U> {};

template<typename Base, typename U>
struct BaseDimension<Base, Inv<U>>
{
    constexpr static int value = -BaseDimension<Base, U>::value;
};

template<typename Base, typename U1, typename U2>
struct BaseDimension<Base, Mult<U1, U2>>
{
    constexpr static int value = BaseDimension<Base, U1>::value + BaseDimension<Base, U2>::value;
};

template<typename Unit, int dim>
struct Tuple {};

template<typename U>
struct Dimension
{
    template<typename Base>
    struct TupleBaseDimension { using type = Tuple<Base, BaseDimension<Base, U>::value>; };

    template<typename>
    struct HasZeroDimension { constexpr static bool value = false; };

    template<typename T, int dim>
    struct HasZeroDimension<Tuple<T, dim>> { constexpr static bool value = dim == 0; };

    // Map a list of base units into a list of pairs (unit, dimension)
    using tuples = typename Map<TupleBaseDimension, typename BaseUnits<U>::type>::type;

    // Remove the tuples with zero-dimension units
    using type = typename Remove<HasZeroDimension, tuples>::type;
};

template<typename U>
struct Factor
{
    constexpr static double value = 1.0;
};

template<typename U>
struct Factor<Unit<U>>
{
    constexpr static double value = Factor<U>::value;
};

template<typename U, Integer num, Integer den>
struct Factor<Scale<U, num, den>>
{
    constexpr static double value = Factor<U>::value * double(num)/den;
};

template<typename U>
struct Factor<Inv<U>>
{
    constexpr static double value = 1.0/Factor<U>::value;
};

template<typename U1, typename U2>
struct Factor<Mult<U1, U2>>
{
    constexpr static double value = Factor<U1>::value * Factor<U2>::value;
};

template<typename U>
constexpr double factor(const Unit<U>& u)
{
    return Factor<U>::value;
}

template<typename U1, typename U2>
constexpr double auxconvert(Unit<U1> from, Unit<U2> to, double value)
{
    return factor(U1())/factor(U2()) * value;
}

template<typename U1, typename U2, Integer num, Integer den>
constexpr double auxconvert(Add<U1, num, den> from, Unit<U2> to, double value)
{
    return auxconvert(U1(), U2(), value - double(num)/den);
}

template<typename U1, typename U2, Integer num, Integer den>
constexpr double auxconvert(Unit<U1> from, Add<U2, num, den> to, double value)
{
    return auxconvert(U1(), U2(), value) + double(num)/den;
}

} /* namespace internal*/

template<typename U1, typename U2>
struct Convertible
{
    using dim1 = typename internal::Dimension<U1>::type;
    using dim2 = typename internal::Dimension<U2>::type;
    constexpr static bool value = internal::Equal<dim1, dim2>::value;
};

template<typename U1, typename U2>
constexpr double convert(Unit<U1> from, Unit<U2> to, double value)
{
    static_assert(Convertible<U1, U2>::value,
        "*** Error *** conversion of units must have matched dimensions");
    return internal::auxconvert(U1(), U2(), value);
}

} /* namespace units */
