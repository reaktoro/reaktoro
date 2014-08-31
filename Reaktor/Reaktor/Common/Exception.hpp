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

#pragma once

// C++ includes
#include <sstream>
#include <stdexcept>
#include <string>

namespace Reaktor {

/**
 * Provides a convenient way to initialized an exception with helpful error messages
 *
 * **Usage**
 *
 * Below we demonstrate a convenient way to raise an exception using an
 * Exception instance and the @ref raise macro.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * using namespace Reaktor;
 * Exception exception;
 * exception.error << "Cannot calculate the actitivy of species " << species.name() << ".";
 * exception.reason << "The species " << species.name() << " does not exist in the chemical system.";
 * raise(exception);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * @see raise
 * @ingroup Common
 */
struct Exception
{
    /// The error message to be displayed when the exception is raised
    std::stringstream error;

    /// The reason message to be displayed when the exception is raised
    std::stringstream reason;
};

namespace internal {

/**
 * Creates a message error from an Exception instance
 * @param exception The exception instance to be raised
 * @param file The name of the file where the exception was raised
 * @param line The number of the line where the exception was raised
 * @return The message error of the exception
 * @see Exception
 */
std::string message(const Exception& exception, const std::string& file, int line);

} /* namespace internal */

/**
 * Defines a macro to raise a runtime exception from an Exception instance
 *
 * @see Exception
 * @ingroup Common
 */
#define raise(exception) \
    throw std::runtime_error(internal::message(exception, __FILE__, __LINE__));

/**
 * Defines a macro to raise a runtime exception from a error string and a reason string
 *
 * @see Exception
 * @ingroup Common
 */
#define error(errorstr, reasonstr) \
    { \
        Exception exception; \
        exception.error << errorstr; \
        exception.reason << reasonstr; \
        raise(exception); \
    }

} // namespace Reaktor
