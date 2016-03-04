/*********************************************************************************
 * This file is part of CUTE.
 *
 * CUTE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CUTE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CUTE.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2013 Peter Sommerlad
 *
 *********************************************************************************/

#ifndef CUTE_DETERMINE_TRAITS_H_
#define CUTE_DETERMINE_TRAITS_H_
#include "cute_determine_version.h"
#if defined(USE_STD11)
#include <type_traits>
#elif defined(USE_TR1)
#include <tr1/type_traits>
#else
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/make_signed.hpp>
#endif
#if defined(USE_STD11)
	namespace impl_place_for_traits = std;
#elif defined(USE_TR1)
	namespace impl_place_for_traits = std::tr1;
#else
	namespace impl_place_for_traits = boost;
#endif

#endif /* CUTE_DETERMINE_TRAITS_H_ */
