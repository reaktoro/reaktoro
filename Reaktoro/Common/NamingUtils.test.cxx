// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Common/NamingUtils.hpp>
using namespace Reaktoro;

namespace Check {

/// Checks if method Reaktoro::splitSpeciesNameSuffix produces correct results.
auto splitSpeciesNameSuffix(std::string substance, std::string name, std::string suffix)
{
    const auto [name_, suffix_] = Reaktoro::splitSpeciesNameSuffix(substance);
    CHECK( name == name_ );
    CHECK( suffix == suffix_ );
};

} // namespace Check

TEST_CASE("Testing NamingUtils module", "[NamingUtils]")
{
    Check::splitSpeciesNameSuffix("H2O(aq)", "H2O", "aq");
    Check::splitSpeciesNameSuffix("CO2(g)", "CO2", "g");
    Check::splitSpeciesNameSuffix("CaCO3(s, calcite)", "CaCO3", "s, calcite");
    Check::splitSpeciesNameSuffix("Methane(aq)", "Methane", "aq");
    Check::splitSpeciesNameSuffix("H2O", "H2O", "");
    Check::splitSpeciesNameSuffix("Fe[3+]", "Fe[3+]", "");
    Check::splitSpeciesNameSuffix("Fe[3+](aq)", "Fe[3+]", "aq");
    Check::splitSpeciesNameSuffix("Ca++(aq)", "Ca++", "aq");
    Check::splitSpeciesNameSuffix("Ca++", "Ca++", "");
}
