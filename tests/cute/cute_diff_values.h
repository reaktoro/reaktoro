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

#ifndef CUTE_DIFF_VALUES_H_
#define CUTE_DIFF_VALUES_H_
#include "cute_to_string.h"
namespace cute {
	// you could provide your own overload for diff_values for your app-specific types
	// be sure to use tabs as given below, then the CUTE eclipse plug-in will parse correctly
	template <typename ExpectedValue, typename ActualValue>
	std::string diff_values(ExpectedValue const &expected
						, ActualValue const & actual
						, char const *left="expected"
						, char const *right="but was"){
		// construct a simple message...to be parsed by IDE support
		std::string res;
		res += ' ';
		res += left;
		res+=":\t" + cute_to_string::backslashQuoteTabNewline(cute_to_string::to_string(expected))+'\t';
		res += right;
		res +=":\t"+cute_to_string::backslashQuoteTabNewline(cute_to_string::to_string(actual))+'\t';
		return res;
	}
}

#endif /* CUTE_DIFF_VALUES_H_ */
