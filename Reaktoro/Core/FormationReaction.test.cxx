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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/FormationReaction.hpp>
#include <Reaktoro/Core/Species.hpp>
#include <Reaktoro/Singletons/PeriodicTable.hpp>
using namespace Reaktoro;


TEST_CASE("Testing FormationReaction class", "[FormationReaction]")
{
    // FORMATION REACTIONS CONSIDERED IN THE TESTS BELOW
    //    A + 2B = C   ---   (0th level of recursion when computing standard thermo props)
    //    B + 3C = D   ---   (1st level of recursion when computing standard thermo props)
    //    C - 2D = E   ---   (2nd level of recursion when computing standard thermo props)

    const auto R = universalGasConstant;

    const auto lgK_C = 1.234;
    const auto lgK_D = 2.345;
    const auto lgK_E = 3.456;

    const auto dH0_C = 123.4;
    const auto dH0_D = 234.5;
    const auto dH0_E = 345.6;

    const auto A = Species()
        .withName("A")
        .withStandardGibbsEnergy(0.0);

    const auto B = Species()
        .withName("B")
        .withStandardGibbsEnergy(0.0);

    const auto C = Species()
        .withName("C")
        .withFormationReaction(
            FormationReaction()
                .withProduct("C")
                .withReactants({{A, 1}, {B, 2}})
                .withEquilibriumConstant(lgK_C)
                .withEnthalpyChange(dH0_C));

    const auto D = Species()
        .withName("D")
        .withFormationReaction(
            FormationReaction()
                .withProduct("D")
                .withReactants({{B, 1}, {C, 3}})
                .withEquilibriumConstant(lgK_D)
                .withEnthalpyChange(dH0_D));

    const auto E = Species()
        .withName("E")
        .withFormationReaction(
            FormationReaction()
                .withProduct("E")
                .withReactants({{C, 1}, {D, -2}})
                .withEquilibriumConstant(lgK_E)
                .withEnthalpyChange(dH0_E));

    REQUIRE(  A.reaction().product() == ""         );
    REQUIRE(  A.reaction().reactants().size() == 0 );
    REQUIRE( !A.reaction().equilibriumConstantFn() );
    REQUIRE( !A.reaction().enthalpyChangeFn()      );
    REQUIRE( !A.reaction().standardGibbsEnergyFn() );
    REQUIRE( !A.reaction().standardEnthalpyFn()    );

    REQUIRE(  B.reaction().product() == ""         );
    REQUIRE(  B.reaction().reactants().size() == 0 );
    REQUIRE( !B.reaction().equilibriumConstantFn() );
    REQUIRE( !B.reaction().enthalpyChangeFn()      );
    REQUIRE( !B.reaction().standardGibbsEnergyFn() );
    REQUIRE( !B.reaction().standardEnthalpyFn()    );

    REQUIRE( C.reaction().product() == "C"                      );
    REQUIRE( C.reaction().reactants().size() == 2               );
    REQUIRE( C.reaction().reactants().at(0).first.name() == "A" );
    REQUIRE( C.reaction().reactants().at(1).first.name() == "B" );
    REQUIRE( C.reaction().reactants().at(0).second == 1.0       );
    REQUIRE( C.reaction().reactants().at(1).second == 2.0       );
    REQUIRE( C.reaction().equilibriumConstantFn()               );
    REQUIRE( C.reaction().enthalpyChangeFn()                    );
    REQUIRE( C.reaction().standardGibbsEnergyFn()               );
    REQUIRE( C.reaction().standardEnthalpyFn()                  );
    REQUIRE( C.reaction().stoichiometry("A") == 1.0             );
    REQUIRE( C.reaction().stoichiometry("B") == 2.0             );

    REQUIRE( D.reaction().product() == "D"                      );
    REQUIRE( D.reaction().reactants().size() == 2               );
    REQUIRE( D.reaction().reactants().at(0).first.name() == "B" );
    REQUIRE( D.reaction().reactants().at(1).first.name() == "C" );
    REQUIRE( D.reaction().reactants().at(0).second == 1.0       );
    REQUIRE( D.reaction().reactants().at(1).second == 3.0       );
    REQUIRE( D.reaction().equilibriumConstantFn()               );
    REQUIRE( D.reaction().enthalpyChangeFn()                    );
    REQUIRE( D.reaction().standardGibbsEnergyFn()               );
    REQUIRE( D.reaction().standardEnthalpyFn()                  );
    REQUIRE( D.reaction().stoichiometry("B") == 1.0             );
    REQUIRE( D.reaction().stoichiometry("C") == 3.0             );

    REQUIRE( E.reaction().product() == "E"                      );
    REQUIRE( E.reaction().reactants().size() == 2               );
    REQUIRE( E.reaction().reactants().at(0).first.name() == "C" );
    REQUIRE( E.reaction().reactants().at(1).first.name() == "D" );
    REQUIRE( E.reaction().reactants().at(0).second ==  1.0      );
    REQUIRE( E.reaction().reactants().at(1).second == -2.0      );
    REQUIRE( E.reaction().equilibriumConstantFn()               );
    REQUIRE( E.reaction().enthalpyChangeFn()                    );
    REQUIRE( E.reaction().standardGibbsEnergyFn()               );
    REQUIRE( E.reaction().standardEnthalpyFn()                  );
    REQUIRE( E.reaction().stoichiometry("C") ==  1.0            );
    REQUIRE( E.reaction().stoichiometry("D") == -2.0            );

    const auto T = 300.0;
    const auto P = 1.0e5;

    const auto G0_A = 0.0;
    const auto G0_B = 0.0;
    const auto G0_C = G0_A + 2*G0_B - R*T*ln10*lgK_C;
    const auto G0_D = G0_B + 3*G0_C - R*T*ln10*lgK_D;
    const auto G0_E = G0_C - 2*G0_D - R*T*ln10*lgK_E;

    const auto H0_A = 0.0;
    const auto H0_B = 0.0;
    const auto H0_C = H0_A + 2*H0_B + dH0_C;
    const auto H0_D = H0_B + 3*H0_C + dH0_D;
    const auto H0_E = H0_C - 2*H0_D + dH0_E;

    // Convenient functions that evaluate a reaction thermo prop at given temperature and pressure
    auto G0 = [](auto reaction, auto T, auto P) { return reaction.standardGibbsEnergyFn()(T, P); };
    auto H0 = [](auto reaction, auto T, auto P) { return reaction.standardEnthalpyFn()(T, P); };

    REQUIRE( G0(C.reaction(), T, P)  == Approx(G0_C) );
    REQUIRE( G0(D.reaction(), T, P)  == Approx(G0_D) );
    REQUIRE( G0(E.reaction(), T, P)  == Approx(G0_E) );

    REQUIRE( H0(C.reaction(), T, P)  == Approx(H0_C) );
    REQUIRE( H0(D.reaction(), T, P)  == Approx(H0_D) );
    REQUIRE( H0(E.reaction(), T, P)  == Approx(H0_E) );
}
