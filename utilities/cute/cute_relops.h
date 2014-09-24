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
#ifndef CUTE_RELOPS_H_
#define CUTE_RELOPS_H_
#include "cute_base.h"
#include "cute_diff_values.h"
#include <functional>
namespace cute {

namespace cute_relops_detail{
	template <typename LeftValue, typename RightValue>
	std::string compare_values(LeftValue const &left
						, RightValue const & right){
		return cute::diff_values(left,right,"left","right");
	}
}
template <template<typename> class RELOP,typename TL, typename TR>
void assert_relop(TL const &left
		, TR const &right
		,std::string const &msg
		,char const *file
		,int line) {
	if (RELOP<TL>()(left,right)) return;
	throw test_failure(msg + cute_relops_detail::compare_values(left,right),file,line);
}

}
#define ASSERT_LESSM(msg,left,right) cute::assert_relop<std::less>((left),(right),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg),__FILE__,__LINE__)
#define ASSERT_LESS(left,right) ASSERT_LESSM(#left " < " #right, (left),(right))
#define ASSERT_LESS_EQUALM(msg,left,right) cute::assert_relop<std::less_equal>((left),(right),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg),__FILE__,__LINE__)
#define ASSERT_LESS_EQUAL(left,right) ASSERT_LESS_EQUALM(#left " <= " #right, (left),(right))
#define ASSERT_GREATERM(msg,left,right) cute::assert_relop<std::greater>((left),(right),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg),__FILE__,__LINE__)
#define ASSERT_GREATER(left,right) ASSERT_GREATERM(#left " > " #right, (left),(right))
#define ASSERT_GREATER_EQUALM(msg,left,right) cute::assert_relop<std::greater_equal>((left),(right),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg),__FILE__,__LINE__)
#define ASSERT_GREATER_EQUAL(left,right) ASSERT_GREATER_EQUALM(#left " >= " #right, (left),(right))
#define ASSERT_NOT_EQUAL_TOM(msg,left,right) cute::assert_relop<std::not_equal_to>((left),(right),\
		CUTE_FUNCNAME_PREFIX+cute::cute_to_string::backslashQuoteTabNewline(msg),__FILE__,__LINE__)
#define ASSERT_NOT_EQUAL_TO(left,right) ASSERT_NOT_EQUAL_TOM(#left " != " #right, (left),(right))


#endif /* CUTE_RELOPS_H_ */
