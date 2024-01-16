// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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

    const auto dH0_C =   0.0;
    const auto dH0_D = 234.5;
    const auto dH0_E = 345.6;

    const auto dCp0_C =   0.0;
    const auto dCp0_D = 345.6;
    const auto dCp0_E = 457.7;

    const auto V0_C = 16.324;
    const auto V0_D = 17.435;
    const auto V0_E = 18.546;

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
                .withReactants({{A, 1}, {B, 2}})
                .withProductStandardVolume(V0_C)
                .withEquilibriumConstant(lgK_C)
            );

    const auto D = Species()
        .withName("D")
        .withFormationReaction(
            FormationReaction()
                .withReactants({{B, 1}, {C, 3}})
                .withProductStandardVolume(V0_D)
                .withReactionStandardThermoModel(
                    [=](ReactionStandardThermoProps& res, ReactionStandardThermoModelArgs args) {
                        const auto& [T, P, dV0] = args;
                        res.dG0  = -R*T*ln10*lgK_D;
                        res.dH0  = dH0_D;
                        res.dCp0 = dCp0_D;
                        return res;
                    })
            );

    const auto E = Species()
        .withName("E")
        .withFormationReaction(
            FormationReaction()
                .withReactants({{C, 1}, {D, -2}})
                .withProductStandardVolume(V0_E)
                .withReactionStandardThermoModel(
                    [=](ReactionStandardThermoProps& res, ReactionStandardThermoModelArgs args) {
                        const auto& [T, P, dV0] = args;
                        res.dG0  = -R*T*ln10*lgK_E;
                        res.dH0  = dH0_E;
                        res.dCp0 = dCp0_E;
                        return res;
                    })
            );

    REQUIRE(  A.reaction().reactants().size() == 0 );
    REQUIRE( !A.reaction().reactionThermoModel() );
    REQUIRE(  A.reaction().initialized() == false );

    REQUIRE_THROWS( A.reaction().createStandardThermoModel() );

    REQUIRE(  B.reaction().reactants().size() == 0 );
    REQUIRE( !B.reaction().reactionThermoModel()   );
    REQUIRE(  B.reaction().initialized() == false );

    REQUIRE_THROWS( B.reaction().createStandardThermoModel() );

    REQUIRE( C.reaction().reactants().size() == 2               );
    REQUIRE( C.reaction().reactants().at(0).first.name() == "A" );
    REQUIRE( C.reaction().reactants().at(1).first.name() == "B" );
    REQUIRE( C.reaction().reactants().at(0).second == 1.0       );
    REQUIRE( C.reaction().reactants().at(1).second == 2.0       );
    REQUIRE( C.reaction().reactionThermoModel()                 );
    REQUIRE( C.reaction().createStandardThermoModel()           );
    REQUIRE( C.reaction().stoichiometry("A") == 1.0             );
    REQUIRE( C.reaction().stoichiometry("B") == 2.0             );
    REQUIRE( C.reaction().initialized() == true                 );

    REQUIRE( D.reaction().reactants().size() == 2               );
    REQUIRE( D.reaction().reactants().at(0).first.name() == "B" );
    REQUIRE( D.reaction().reactants().at(1).first.name() == "C" );
    REQUIRE( D.reaction().reactants().at(0).second == 1.0       );
    REQUIRE( D.reaction().reactants().at(1).second == 3.0       );
    REQUIRE( D.reaction().reactionThermoModel()                 );
    REQUIRE( D.reaction().createStandardThermoModel()           );
    REQUIRE( D.reaction().stoichiometry("B") == 1.0             );
    REQUIRE( D.reaction().stoichiometry("C") == 3.0             );
    REQUIRE( D.reaction().initialized() == true                 );

    REQUIRE( E.reaction().reactants().size() == 2               );
    REQUIRE( E.reaction().reactants().at(0).first.name() == "C" );
    REQUIRE( E.reaction().reactants().at(1).first.name() == "D" );
    REQUIRE( E.reaction().reactants().at(0).second ==  1.0      );
    REQUIRE( E.reaction().reactants().at(1).second == -2.0      );
    REQUIRE( E.reaction().reactionThermoModel()                 );
    REQUIRE( E.reaction().createStandardThermoModel()           );
    REQUIRE( E.reaction().stoichiometry("C") ==  1.0            );
    REQUIRE( E.reaction().stoichiometry("D") == -2.0            );
    REQUIRE( E.reaction().initialized() == true                 );

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

    const auto Cp0_A = 0.0;
    const auto Cp0_B = 0.0;
    const auto Cp0_C = Cp0_A + 2*Cp0_B + dCp0_C;
    const auto Cp0_D = Cp0_B + 3*Cp0_C + dCp0_D;
    const auto Cp0_E = Cp0_C - 2*Cp0_D + dCp0_E;

    // Convenient functions that evaluate a reaction thermo prop at given temperature and pressure
    auto G0  = [](auto reaction, auto T, auto P) { return reaction.createStandardThermoModel()(T, P).G0;  };
    auto H0  = [](auto reaction, auto T, auto P) { return reaction.createStandardThermoModel()(T, P).H0;  };
    auto Cp0 = [](auto reaction, auto T, auto P) { return reaction.createStandardThermoModel()(T, P).Cp0; };

    REQUIRE( G0(C.reaction(), T, P)  == Approx(G0_C) );
    REQUIRE( G0(D.reaction(), T, P)  == Approx(G0_D) );
    REQUIRE( G0(E.reaction(), T, P)  == Approx(G0_E) );

    REQUIRE( H0(C.reaction(), T, P)  == Approx(H0_C) );
    REQUIRE( H0(D.reaction(), T, P)  == Approx(H0_D) );
    REQUIRE( H0(E.reaction(), T, P)  == Approx(H0_E) );

    REQUIRE( Cp0(C.reaction(), T, P)  == Approx(Cp0_C) );
    REQUIRE( Cp0(D.reaction(), T, P)  == Approx(Cp0_D) );
    REQUIRE( Cp0(E.reaction(), T, P)  == Approx(Cp0_E) );
}
