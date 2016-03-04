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
 * Copyright 2007-2011 Peter Sommerlad, Emanuel Graf
 *
 *********************************************************************************/

#ifndef CUTE_EQUALS_H_
#define CUTE_EQUALS_H_
#include "cute_base.h"
#include "cute_diff_values.h"
#include "cute_determine_traits.h"
#include <cmath>
#include <limits>
#include <algorithm>


namespace cute {
	namespace cute_do_equals {
		// provide some template meta programming tricks to select "correct" comparison for floating point and integer types
		template <typename ExpectedValue, typename ActualValue, typename DeltaValue>
		bool do_equals_floating_with_delta(ExpectedValue const &expected
				,ActualValue const &actual
				,DeltaValue const &delta) {
			return std::abs(delta)  >= std::abs(expected-actual);
		}
		template <typename ExpectedValue, typename ActualValue, bool select_non_floating_point_type>
		bool do_equals_floating(ExpectedValue const &expected
					,ActualValue const &actual,const impl_place_for_traits::integral_constant<bool, select_non_floating_point_type>&){
			return expected==actual; // normal case for most types uses operator==!
		}
		template <typename ExpectedValue, typename ActualValue>
		bool do_equals_floating(ExpectedValue const &expected
					,ActualValue const &actual,const impl_place_for_traits::true_type&){
			const ExpectedValue automatic_delta_masking_last_significant_digit=(10*std::numeric_limits<ExpectedValue>::epsilon())*expected;
			return do_equals_floating_with_delta(expected,actual,automatic_delta_masking_last_significant_digit);
		}
		// TMP-overload dispatch for floating points 2 bool params --> 4 overloads
		template <typename ExpectedValue, typename ActualValue, bool select_non_integral_type>
		bool do_equals(ExpectedValue const &expected
					,ActualValue const &actual
					,const impl_place_for_traits::integral_constant<bool, select_non_integral_type>&exp_is_integral
					,const impl_place_for_traits::integral_constant<bool, select_non_integral_type>&act_is_integral){
			return do_equals_floating(expected,actual,impl_place_for_traits::is_floating_point<ExpectedValue>());
		}
		template <typename ExpectedValue, typename ActualValue, bool select_non_integral_type>
		bool do_equals(ExpectedValue const &expected
					,ActualValue const &actual
					,const impl_place_for_traits::false_type&,const impl_place_for_traits::false_type&){
			return do_equals_floating(expected,actual,impl_place_for_traits::is_floating_point<ActualValue>());
		}
		template <typename ExpectedValue, typename ActualValue, bool select_non_integral_type>
		bool do_equals(ExpectedValue const &expected
					,ActualValue const &actual
					,const impl_place_for_traits::integral_constant<bool, select_non_integral_type>&exp_is_integral
					,const impl_place_for_traits::true_type&){
			return do_equals_floating(expected,actual,impl_place_for_traits::is_floating_point<ExpectedValue>());
		}
		template <typename ExpectedValue, typename ActualValue, bool select_non_integral_type>
		bool do_equals(ExpectedValue const &expected
					,ActualValue const &actual
					,const impl_place_for_traits::true_type&,const impl_place_for_traits::integral_constant<bool, select_non_integral_type>&act_is_integral){
			return do_equals_floating(expected,actual,impl_place_for_traits::is_floating_point<ActualValue>());
		}
		// can I get rid of the following complexity by doing a do_equals_integral
		// parameterized by is_signed<ExpectedValue>==is_signed<ActualValue> or nofBits<A> < nofBits<B>


		// this is an optimization for avoiding if and sign-extend overhead if both int types are the same as below
		template <typename IntType>
		bool do_equals(IntType const &expected
					,IntType const &actual
					,const impl_place_for_traits::true_type&,const impl_place_for_traits::true_type&){
			return expected==actual;
		}
		// bool cannot be made signed, therefore we need the following three overloads, also to avoid ambiguity
		template <typename IntType>
		bool do_equals(bool const &expected
					,IntType const &actual
					,const impl_place_for_traits::true_type&,const impl_place_for_traits::true_type&){
			return expected== !(!actual); // warning from VS
		}
		template <typename IntType>
		bool do_equals(IntType const &expected
					,bool const &actual
					,const impl_place_for_traits::true_type&,const impl_place_for_traits::true_type&){
			return !(!expected)==actual; // warning from VS
		}
		// do not forget the inline on a non-template overload!
		// this overload is needed to actually avoid ambiguity for comparing bool==bool as a best match
		inline bool do_equals(bool const &expected
				      ,bool const &actual
				      , const impl_place_for_traits::true_type&,const impl_place_for_traits::true_type&){
			return expected==actual;
		}
		// overload for char const *, my test case failed because VC++ doesn't use string constant folding like g++/clang
		// a feature where we should do string comparison
		inline bool do_equals(char const *const &expected
				      ,char const *const &actual
				      , const impl_place_for_traits::false_type&,const impl_place_for_traits::false_type&){
			return std::string(expected) == actual;
		}
		template <typename IntegralType>
		size_t nof_bits(IntegralType const &){
			return std::numeric_limits<IntegralType>::digits;
		}
#if defined(USE_STD11)||defined(USE_TR1)
		template <typename ExpectedValue, typename ActualValue>
		bool do_equals_integral(ExpectedValue const &expected
				,ActualValue const &actual
				,const impl_place_for_traits::true_type&,const impl_place_for_traits::true_type&){
			if (nof_bits(expected) < nof_bits(actual))
						return static_cast<ActualValue>(expected) == actual;
			else
						return expected == static_cast<ExpectedValue>(actual);
			return false;
		}
		template <typename ExpectedValue, typename ActualValue>
		bool do_equals_integral(ExpectedValue const &expected
				,ActualValue const &actual
				,const impl_place_for_traits::false_type&,const impl_place_for_traits::true_type&){
//TODO complicated case, one signed one unsigned type. since it is about equality casting to the longer should work?
			if (sizeof(ExpectedValue) >	sizeof(ActualValue))
				return expected==static_cast<ExpectedValue>(actual);
			else
				return static_cast<ActualValue>(expected) == actual;
		}
		template <typename ExpectedValue, typename ActualValue>
		bool do_equals_integral(ExpectedValue const &expected
				,ActualValue const &actual
				,const impl_place_for_traits::true_type&,const impl_place_for_traits::false_type&){
//TODO
			if (sizeof(ExpectedValue) < sizeof(ActualValue))
				return static_cast<ActualValue>(expected)==	actual;
			else
				return expected == static_cast<ExpectedValue>(actual);
		}
		template <typename ExpectedValue, typename ActualValue>
		bool do_equals_integral(ExpectedValue const &expected
				,ActualValue const &actual
				,const impl_place_for_traits::false_type&,const impl_place_for_traits::false_type&){
			if (nof_bits(expected) < nof_bits(actual))
						return static_cast<ActualValue>(expected) == actual;
			else
						return expected == static_cast<ExpectedValue>(actual);
			return false;
		}
#endif
		// will not work if either type is bool, therefore the overloads above.
		template <typename ExpectedValue, typename ActualValue>
		bool do_equals(ExpectedValue const &expected
					,ActualValue const &actual
					,const impl_place_for_traits::true_type&,const impl_place_for_traits::true_type&){
#if defined(USE_STD11) || defined(USE_TR1)
			return do_equals_integral(expected,actual,
					impl_place_for_traits::is_signed<ExpectedValue>()
					,impl_place_for_traits::is_signed<ActualValue>());
#else
//TODO: replace the following code with a dispatcher on signed/unsigned
			typedef typename impl_place_for_traits::make_signed<ExpectedValue>::type ex_s;
			typedef typename impl_place_for_traits::make_signed<ActualValue>::type ac_s;
				// need to sign extend with the longer type, should work...
			    // might be done through more template meta prog tricks....but...
				if (nof_bits(expected) < nof_bits(actual))
					return static_cast<ac_s>(expected) == static_cast<ac_s>(actual);
				else
					return static_cast<ex_s>(expected) == static_cast<ex_s>(actual);
#endif
		}
	} // namespace equals_impl
	template <typename ExpectedValue, typename ActualValue>
	void assert_equal(ExpectedValue const &expected
				,ActualValue const &actual
				,std::string const &msg
				,char const *file
				,int line) {
		// unfortunately there is no is_integral_but_not_bool_or_enum
		typedef typename impl_place_for_traits::is_integral<ExpectedValue> exp_integral;
		typedef typename impl_place_for_traits::is_integral<ActualValue> act_integral;
		if (cute_do_equals::do_equals(expected,actual,exp_integral(),act_integral()))
			return;
		throw test_failure(msg + diff_values(expected,actual),file,line);
	}

	template <typename ExpectedValue, typename ActualValue, typename DeltaValue>
	void assert_equal_delta(ExpectedValue const &expected
				,ActualValue const &actual
				,DeltaValue const &delta
				,std::string const &msg
				,char const *file
				,int line) {
		if (cute_do_equals::do_equals_floating_with_delta(expected,actual,delta)) return;
		throw test_failure(msg + diff_values(expected,actual),file,line);
	}

}

#define ASSERT_EQUALM(msg,expected,actual) cute::assert_equal((expected),(actual),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg),__FILE__,__LINE__)
#define ASSERT_EQUAL(expected,actual) ASSERT_EQUALM(#expected " == " #actual, (expected),(actual))
#define ASSERT_EQUAL_DELTAM(msg,expected,actual,delta) cute::assert_equal_delta((expected),(actual),(delta),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg),__FILE__,__LINE__)
#define ASSERT_EQUAL_DELTA(expected,actual,delta) ASSERT_EQUAL_DELTAM(#expected " == " #actual " with error " #delta  ,(expected),(actual),(delta))
#endif /*CUTE_EQUALS_H_*/
