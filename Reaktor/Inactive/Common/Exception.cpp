/*
 * Reaktor is a C++ library for computational reaction modelling.
 *
 * Copyright (C) 2014 Allan Leal
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "Exception.hpp"

// C++ includes
#include <algorithm>
#include <sstream>

namespace Reaktor {
namespace internal {

/**
 * Creates the location string from the file name and line number
 *
 * The result of this function on the file `/home/user/git/Reaktor/Reaktor/Reaktor/Core/ChemicalSystem.cpp`
 * will be `Reaktor/Core/ChemicalSystem.cpp`.
 *
 * @param file The full path to the file
 * @param line The line number
 */
std::string location(const std::string& file, int line)
{
    std::string str = "Reaktor/";
    auto pos = std::find_end(file.begin(), file.end(), str.begin(), str.end()) - file.begin();
	std::stringstream ss;
	ss << file.substr(pos, file.size() - pos) << ":" << line;
    return ss.str();
}

std::string message(const Exception& exception, const std::string& file, int line)
{
    std::string error = exception.error.str();
    std::string reason = exception.reason.str();
    std::string location = internal::location(file, line);
    unsigned length = std::max(error.size(), std::max(reason.size(), location.size())) + 25;
    std::string bar(length, '*');
    std::stringstream message;
    message << std::endl;
    message << bar << std::endl;
    message << "*** Error: " << error << std::endl;
    message << "*** Reason: " << reason << std::endl;
    message << "*** Location:  This error was encountered in " << location << "." << std::endl;
    message << bar << std::endl;
    message << std::endl;
    return message.str();
}

} /* namespace internal */
} // namespace Reaktor
