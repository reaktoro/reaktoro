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

#include "MineralCatalyst.hpp"

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/StringUtils.hpp>

namespace Reaktoro {
namespace internal {

inline auto invalidCatalystError(std::string catalyst) -> void
{
    Exception exception;
    exception.error << "Cannot set the mineral catalyst with given catalyst string: " << catalyst << ".";
    exception.reason << "The provided catalyst string is not formated correctly.";
    RaiseError(exception);
}

inline auto checkCatalystQuantity(std::string quantity) -> void
{
    if(quantity != "a" && quantity != "activity" && quantity != "p" && quantity != "pressure")
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

MineralCatalyst::MineralCatalyst(std::string species, std::string quantity, double power)
: species(species), quantity(quantity), power(power)
{
    internal::checkCatalystQuantity(quantity);
}

MineralCatalyst::MineralCatalyst(std::string catalyst)
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

auto operator<(const MineralCatalyst& lhs, const MineralCatalyst& rhs) -> bool
{
    return lhs.species < rhs.species;
}

auto operator==(const MineralCatalyst& lhs, const MineralCatalyst& rhs) -> bool
{
    return lhs.species == rhs.species;
}

} // namespace Reaktoro
