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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Equilibrium/EquilibriumRestrictions.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumRestrictions", "[EquilibriumRestrictions]")
{
    ChemicalSystem system = test::createChemicalSystem();

    const auto idx = [&](auto speciesname) { return system.species().index(speciesname); };

    EquilibriumRestrictions restrictions(system);

    restrictions.cannotReact("SiO2(s)");

    restrictions.cannotIncrease("CaCO3(s)");
    restrictions.cannotIncrease("CO2(g)");

    restrictions.cannotIncreaseAbove("CO(g)", 2.0, "g");
    restrictions.cannotIncreaseAbove("HCl(aq)", 1.0, "umol");

    restrictions.cannotDecrease("MgCO3(s)");
    restrictions.cannotDecrease("CH4(g)");

    restrictions.cannotDecreaseBelow("NaCl(s)", 0.1, "mmol");
    restrictions.cannotDecreaseBelow("HCl(aq)", 0.1, "umol");

    CHECK( restrictions.speciesCannotIncrease().size() == 3 );
    CHECK( restrictions.speciesCannotIncrease().count(idx("SiO2(s)")) );
    CHECK( restrictions.speciesCannotIncrease().count(idx("CaCO3(s)")) );
    CHECK( restrictions.speciesCannotIncrease().count(idx("CO2(g)")) );

    CHECK( restrictions.speciesCannotIncreaseAbove().size() == 2 );
    CHECK( restrictions.speciesCannotIncreaseAbove().at(idx("CO(g)")) == Approx(2.0/28.0104) ); // molar mass of CO is 28.0104 g/mol
    CHECK( restrictions.speciesCannotIncreaseAbove().at(idx("HCl(aq)")) == Approx(1.0e-6) );

    CHECK( restrictions.speciesCannotDecrease().size() == 3 );
    CHECK( restrictions.speciesCannotDecrease().count(idx("SiO2(s)")) );
    CHECK( restrictions.speciesCannotDecrease().count(idx("MgCO3(s)")) );
    CHECK( restrictions.speciesCannotDecrease().count(idx("CH4(g)")) );

    CHECK( restrictions.speciesCannotDecreaseBelow().size() == 2 );
    CHECK( restrictions.speciesCannotDecreaseBelow().at(idx("NaCl(s)")) == Approx(0.1e-3) );
    CHECK( restrictions.speciesCannotDecreaseBelow().at(idx("HCl(aq)")) == Approx(0.1e-6) );

    restrictions.cannotIncreaseAbove("HCl(aq)", 100.0, "mol");
    restrictions.cannotDecreaseBelow("HCl(aq)", 1.0e-6, "mol");

    CHECK( restrictions.speciesCannotIncreaseAbove().at(idx("HCl(aq)")) == Approx(100.0) ); // check the previous change in bounds of HCl(aq) is correct
    CHECK( restrictions.speciesCannotDecreaseBelow().at(idx("HCl(aq)")) == Approx(1.0e-6) ); // check the previous change in bounds of HCl(aq) is correct

    restrictions.canReactFreely("SiO2(s)");

    restrictions.canIncreaseFreely("CaCO3(s)");
    restrictions.canIncreaseFreely("CO2(g)");
    restrictions.canIncreaseFreely("CO(g)");
    restrictions.canIncreaseFreely("HCl(aq)");

    restrictions.canDecreaseFreely("MgCO3(s)");
    restrictions.canDecreaseFreely("CH4(g)");
    restrictions.canDecreaseFreely("NaCl(s)");
    restrictions.canDecreaseFreely("HCl(aq)");

    CHECK_NOTHROW( restrictions.canIncreaseFreely("H2O(aq)") ); // no bounds were set for H2O(aq) before, so ensure this call does not raise any error
    CHECK_NOTHROW( restrictions.canDecreaseFreely("H2O(aq)") ); // no bounds were set for H2O(aq) before, so ensure this call does not raise any error

    CHECK( restrictions.speciesCannotIncrease().size() == 0 );
    CHECK( restrictions.speciesCannotIncreaseAbove().size() == 0 );
    CHECK( restrictions.speciesCannotDecrease().size() == 0 );
    CHECK( restrictions.speciesCannotDecreaseBelow().size() == 0 );
}
