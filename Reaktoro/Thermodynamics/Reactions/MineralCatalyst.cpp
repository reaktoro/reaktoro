// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "MineralCatalyst.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {
namespace internal {

inline auto invalidCatalystError(const std::string& catalyst) -> void
{
    Exception exception;
    exception.error << "Cannot set the mineral catalyst with given catalyst string: " << catalyst << ".";
    exception.reason << "The provided catalyst string is not formated correctly.";
    RaiseError(exception);
}

inline auto checkCatalystQuantity(const std::string& quantity) -> void
{
    if(quantity != "a" and quantity != "activity" and quantity != "p" and quantity != "pressure")
    {
        Exception exception;
        exception.error << "Cannot set the mineral catalyst with given catalyst quantity: " << quantity << ".";
        exception.reason << "The provided catalyst quantity is not supported.";
        RaiseError(exception);
    }
}

} /* namespace internal */

MineralCatalyst::MineralCatalyst()
{}

MineralCatalyst::MineralCatalyst(const std::string& species, const std::string& quantity, double power)
: species(species), quantity(quantity), power(power)
{
    internal::checkCatalystQuantity(quantity);
}

MineralCatalyst::MineralCatalyst(const std::string& catalyst)
{
    // Split the string in two parts: before and after the '=' sign
    std::vector<std::string> words = split(catalyst, "= ");

    // Check if the split of the string result in only two new strings
    if(words.size() != 2) internal::invalidCatalystError(catalyst);

    // Set the power of the catalyst with the second word
    power = tofloat(words[1]);

    // Split the first word in two parts separated by the characters '[' and ']'
    words = split(words[0], "[]");

    // Set the species and quantity members
    species = words[1];
    quantity = words[0];

    // Check if the provided catalyser quantity is valid
    internal::checkCatalystQuantity(quantity);
}

} // namespace Reaktoro
