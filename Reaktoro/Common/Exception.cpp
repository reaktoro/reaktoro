// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "Exception.hpp"

// C++ includes
#include <algorithm>
#include <sstream>

namespace Reaktoro {
namespace internal {

/// Creates the location string from the file name and line number.
/// The result of this function on the file `/home/user/git/Reaktoro/Reaktoro/Reaktoro/Core/Species.cpp`
/// will be `Reaktoro/Core/Species.cpp`.
/// @param file The full path to the file
/// @param line The line number
std::string location(const std::string& file, int line)
{
    std::string str = "Reaktoro/";
    auto pos = std::find_end(file.begin(), file.end(), str.begin(), str.end()) - file.begin();
	std::stringstream ss;
	ss << file.substr(pos, file.size() - pos) << ":" << line;
    return ss.str();
}

std::string message(const Exception& exception, const std::string& file, int line)
{
    std::string error = exception.error.str();
    std::string reason = exception.reason.str();
    std::string loc = location(file, line);
    unsigned length = std::max(error.size(), std::max(reason.size(), loc.size())) + 25;
    std::string bar(length, '*');
    std::stringstream message;
    message << std::endl;
    message << bar << std::endl;
    message << "*** Error: " << error << std::endl;
    message << "*** Reason: " << reason << std::endl;
    message << "*** Location:  This error was encountered in " << loc << "." << std::endl;
    message << bar << std::endl;
    message << std::endl;
    return message.str();
}

} // namespace internal
} // namespace Reaktoro
