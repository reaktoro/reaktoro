// Reaktor is a C++ library for computational reaction modelling.
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

// Cute includes
#include <cute/cute.h>

// Reaktor includes
#include <Reaktor/Core/Species.hpp>
using namespace Reaktor;

auto testSpecies() -> void
{
    Species species;
    species.setName("ABC-");
    species.setCharge(-1.0);
    species.setElements({{"A", 1}, {"B", 1}, {"C", 1}});
    species.setFormula("ABC");
    species.setMolarMass(100.0);

    ASSERT_EQUAL("ABC-", species.name());
    ASSERT_EQUAL(-1.0, species.charge());
    ASSERT_EQUAL(3, species.elements().size());
    ASSERT_EQUAL(1, species.elements().at("A"));
    ASSERT_EQUAL(1, species.elements().at("B"));
    ASSERT_EQUAL(1, species.elements().at("C"));
    ASSERT_EQUAL("ABC", species.formula());
    ASSERT_EQUAL(100.0, species.molarMass());
}

int main(int argc, char **argv)
{
    testSpecies();
}
