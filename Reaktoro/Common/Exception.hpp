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

#pragma once

// C++ includes
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Reaktoro includes
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {

/// Provides a convenient way to initialized an exception with helpful error messages.
/// **Usage**
/// Below we demonstrate a convenient way to raise an exception using an
/// Exception instance and the @ref raise macro.
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// using namespace Reaktoro;
/// Exception exception;
/// exception.error << "Cannot calculate the actitivy of species " << species.name() << ".";
/// exception.reason << "The species " << species.name() << " does not exist in the chemical system.";
/// RaiseError(exception);
/// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// @see raise
/// @ingroup Common
struct Exception
{
    /// The error message to be displayed when the exception is raised
    std::stringstream error;

    /// The reason message to be displayed when the exception is raised
    std::stringstream reason;
};

namespace internal {

/// Create a message error from an Exception instance
/// @param exception The exception instance to be raised
/// @param file The name of the file where the exception was raised
/// @param line The number of the line where the exception was raised
/// @return The message error of the exception
/// @see Exception
std::string message(const Exception& exception, const std::string& file, int line);

} // namespace internal

/// Define a macro to raise a runtime exception from an Exception instance.
/// @see Exception
/// @ingroup Common
#define RaiseError(exception) \
    throw std::runtime_error(Reaktoro::internal::message(exception, __FILE__, __LINE__));

/// Define a macro to raise a runtime exception from a error string and a reason string.
/// @see Exception
/// @ingroup Common
#define RuntimeError(errorstr, reasonstr) \
    { \
        Reaktoro::Exception exception; \
        exception.error << errorstr; \
        exception.reason << reasonstr; \
        RaiseError(exception); \
    }

/// Define a macro to raise a runtime exception from a error string and a reason string.
/// @see Exception
/// @ingroup Common
#define Assert(condition, errorstr, reasonstr) \
    { \
        if(!(condition)) { \
            Reaktoro::Exception exception; \
            exception.error << errorstr; \
            exception.reason << reasonstr; \
            RaiseError(exception); \
        } \
    }

/// Issue a warning message if condition is true.
template<typename... Args>
auto warning(bool condition, Args... items) -> void
{
    if(condition)
        std::cerr << "\033[1;33m***WARNING***\n" << str(items...) << "\n\033[0m";
}

/// Raise a runtime error if condition is true.
template<typename... Args>
auto error(bool condition, Args... items) -> void
{
    if(condition)
        throw std::runtime_error(str("\033[1;31m***ERROR***\n", str(items...), "\n\033[0m"));
}

} // namespace Reaktoro
