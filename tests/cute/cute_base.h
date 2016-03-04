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
 * Copyright 2007-2013 Peter Sommerlad
 *
 *********************************************************************************/

#ifndef CUTE_BASE_H_
#define CUTE_BASE_H_
#include <string>
#include "cute_to_string.h"
#include "cute_determine_version.h"
namespace cute{
	struct test_failure {
		std::string reason;
		std::string filename;
		int lineno;

		test_failure(std::string const &r,char const *f, int line)
		:reason(r),filename(f),lineno(line)
		{ 	}
		char const * what() const { return reason.c_str(); }
	};
}
#if defined(USE_STD11) && !defined(_MSC_VER)
#define CUTE_FUNCNAME_PREFIX std::string(__func__)+": "
#else
#if defined( _MSC_VER ) || defined(__GNUG__)
//#if defined(USE_STD11)
//#define CUTE_FUNCNAME_PREFIX
// workaround doesn't work namespace { char const __FUNCTION__ []="lambda";}
// MSC can not use lambdas outside of function bodies for tests.
//#endif
// use -D CUTE_FUNCNAME_PREFIX if you want to use non-local lambdas with test macros
#if !defined(CUTE_FUNCNAME_PREFIX)
#define CUTE_FUNCNAME_PREFIX std::string(__FUNCTION__)+": "
#endif
#else // could provide #defines for other compiler-specific extensions... need to experiment, i.e., MS uses __FUNCTION__
#define CUTE_FUNCNAME_PREFIX std::string("")
#endif
#endif
#define ASSERTM(msg,cond) do { if (!(cond)) throw cute::test_failure(CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg),__FILE__,__LINE__);} while(false)
#define ASSERT(cond) ASSERTM(#cond,cond)
#define FAIL() ASSERTM("FAIL()",false)
#define FAILM(msg) ASSERTM(msg,false)
#endif /*CUTE_BASE_H_*/
