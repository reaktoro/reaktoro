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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/Utils.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing CoreUtils", "[CoreUtils]")
{
    ChemicalSystem system = test::createChemicalSystem();

    const auto mmH2O = system.species().get("H2O(aq)").molarMass(); // kg/mol
    const auto iH2O = system.species().index("H2O(aq)");

    real actual = {};
    real expected = {};

    //======================================================================
    // Testing computeSpeciesAmount
    //======================================================================
    actual = detail::computeSpeciesAmount(system, iH2O, 1.0, "mol");
    expected = 1.0;
    CHECK( actual == expected );

    actual = detail::computeSpeciesAmount(system, iH2O, 1.0, "mmol");
    expected = 0.001;
    CHECK( actual == expected );

    actual = detail::computeSpeciesAmount(system, iH2O, 1.0, "kg");
    expected = 1.0 / mmH2O;
    CHECK( actual == Approx(expected) );

    actual = detail::computeSpeciesAmount(system, iH2O, 1.0, "g");
    expected = 0.001 / mmH2O;
    CHECK( actual == Approx(expected) );
}
