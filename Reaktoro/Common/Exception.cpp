// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

std::string message(const Exception& exception, const std::string& file, int line)
{
    std::string error = exception.error.str();
    std::string reason = exception.reason.str();
    std::stringstream message;
    message << std::endl;
    message << "*******************************************************************************" << std::endl;
    message << "* Error: " << error << std::endl;
    message << "* Reason: " << reason << std::endl;
    message << "* Location: " << file << ":" << line << std::endl;
    message << "*******************************************************************************" << std::endl;
    message << std::endl;
    return message.str();
}

} // namespace internal
} // namespace Reaktoro
