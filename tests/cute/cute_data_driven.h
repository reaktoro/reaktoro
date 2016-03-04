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

#ifndef CUTE_DATA_DRIVEN_H_
#define CUTE_DATA_DRIVEN_H_

#include "cute_base.h"
#include "cute_equals.h"
#include "cute_relops.h"

#define DDTM(msg) (cute::test_failure(cute::cute_to_string::backslashQuoteTabNewline(msg),__FILE__,__LINE__))
#define DDT() DDTM("")

#define ASSERT_DDTM(cond,msg,failure) do{ \
	if (!(cond)) \
		throw cute::test_failure(CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg)+(failure).reason,\
				(failure).filename.c_str(),(failure).lineno);\
} while (false)
#define ASSERT_DDT(cond,failure) ASSERT_DDTM(cond,#cond,(failure))

#define ASSERT_EQUAL_DDTM(msg,expected,actual,failure) cute::assert_equal((expected),(actual),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg)\
		+(failure).reason,(failure).filename.c_str(),(failure).lineno)
#define ASSERT_EQUAL_DDT(expected,actual,failure) ASSERT_EQUAL_DDTM(#expected " == " #actual, (expected),(actual),(failure))
#define ASSERT_EQUAL_DELTA_DDTM(msg,expected,actual,delta,failure) cute::assert_equal_delta((expected),(actual),(delta),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg)\
		+(failure).reason,(failure).filename.c_str(),(failure).lineno)
#define ASSERT_EQUAL_DELTA_DDT(expected,actual,delta,failure) ASSERT_EQUAL_DELTA_DDTM(#expected " == " #actual " with error " #delta  ,(expected),(actual),(delta),(failure))

#define ASSERT_LESS_DDTM(msg,left,right,failure) cute::assert_relop<std::less>((left),(right),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg)\
		+(failure).reason,(failure).filename.c_str(),(failure).lineno)
#define ASSERT_LESS_DDT(left,right,failure) ASSERT_LESS_DDTM(#left " < " #right, (left),(right),(failure))

#define ASSERT_LESS_EQUAL_DDTM(msg,left,right,failure) cute::assert_relop<std::less_equal>((left),(right),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg)\
		+(failure).reason,(failure).filename.c_str(),(failure).lineno)
#define ASSERT_LESS_EQUAL_DDT(left,right,failure) ASSERT_LESS_EQUAL_DDTM(#left " <= " #right, (left),(right),(failure))

#define ASSERT_GREATER_DDTM(msg,left,right,failure) cute::assert_relop<std::greater>((left),(right),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg)\
		+(failure).reason,(failure).filename.c_str(),(failure).lineno)
#define ASSERT_GREATER_DDT(left,right,failure) ASSERT_GREATER_DDTM(#left " > " #right, (left),(right),(failure))

#define ASSERT_GREATER_EQUAL_DDTM(msg,left,right,failure) cute::assert_relop<std::greater_equal>((left),(right),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg)\
		+(failure).reason,(failure).filename.c_str(),(failure).lineno)
#define ASSERT_GREATER_EQUAL_DDT(left,right,failure) ASSERT_GREATER_EQUAL_DDTM(#left " >= " #right, (left),(right),(failure))

#define ASSERT_NOT_EQUAL_TO_DDTM(msg,left,right,failure) cute::assert_relop<std::not_equal_to>((left),(right),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg)\
		+(failure).reason,(failure).filename.c_str(),(failure).lineno)
#define ASSERT_NOT_EQUAL_TO_DDT(left,right,failure) ASSERT_NOT_EQUAL_TO_DDTM(#left " != " #right, (left),(right),(failure))



#endif /* CUTE_DATA_DRIVEN_H_ */
