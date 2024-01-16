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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Core/Reaction.hpp>
using namespace Reaktoro;

namespace test {

/// Return a mock Database object for test reasons (imported).
extern auto createDatabase() -> Database;

} // namespace test

TEST_CASE("Testing Reaction class", "[Reaction]")
{
    Reaction reaction;

    reaction = reaction.withName("Dolomite");
    reaction = reaction.withEquation("CaCO3(s) = Ca++ + CO3--");
    reaction = reaction.withRateModel([](ChemicalProps const& props) -> ReactionRate { return 1.0; });

    REQUIRE( reaction.name() == "Dolomite" );
    REQUIRE( reaction.equation().size() == 3 );
    REQUIRE( reaction.equation().coefficient("CaCO3(s)") == -1 );
    REQUIRE( reaction.equation().coefficient("Ca++") == 1 );
    REQUIRE( reaction.equation().coefficient("CO3--") == 1 );
    REQUIRE( reaction.rateModel().initialized() );

    //-------------------------------------------------------------------------
    // TESTING METHOD: Reaction::props
    //-------------------------------------------------------------------------
    Database db = test::createDatabase();

    const auto T  = 10.0;
    const auto P  = 20.0;
    const auto RT = universalGasConstant * T;

    ReactionThermoProps rprops;

    auto G0  = [&](auto name) { return db.species().get(name).props(T, P).G0;  };
    auto H0  = [&](auto name) { return db.species().get(name).props(T, P).H0;  };
    auto V0  = [&](auto name) { return db.species().get(name).props(T, P).V0;  };
    auto VT0 = [&](auto name) { return db.species().get(name).props(T, P).VT0; };
    auto VP0 = [&](auto name) { return db.species().get(name).props(T, P).VP0; };
    auto Cp0 = [&](auto name) { return db.species().get(name).props(T, P).Cp0; };
    auto Cv0 = [&](auto name) { return db.species().get(name).props(T, P).Cv0; };
    auto U0  = [&](auto name) { return db.species().get(name).props(T, P).U0;  };
    auto S0  = [&](auto name) { return db.species().get(name).props(T, P).S0;  };
    auto A0  = [&](auto name) { return db.species().get(name).props(T, P).A0;  };
    auto lgK = [&](auto dG0) { return -dG0/(RT*ln10); };

    reaction = db.reaction("H2O(aq) = H+(aq) + OH-(aq)");
    rprops = reaction.props(T, P);

    REQUIRE( rprops.dG0  == Approx(  G0("H+(aq)") +  G0("OH-(aq)") -  G0("H2O(aq)") ) );
    REQUIRE( rprops.dH0  == Approx(  H0("H+(aq)") +  H0("OH-(aq)") -  H0("H2O(aq)") ) );
    REQUIRE( rprops.dV0  == Approx(  V0("H+(aq)") +  V0("OH-(aq)") -  V0("H2O(aq)") ) );
    REQUIRE( rprops.dVT0 == Approx( VT0("H+(aq)") + VT0("OH-(aq)") - VT0("H2O(aq)") ) );
    REQUIRE( rprops.dVP0 == Approx( VP0("H+(aq)") + VP0("OH-(aq)") - VP0("H2O(aq)") ) );
    REQUIRE( rprops.dCp0 == Approx( Cp0("H+(aq)") + Cp0("OH-(aq)") - Cp0("H2O(aq)") ) );
    REQUIRE( rprops.dCv0 == Approx( Cv0("H+(aq)") + Cv0("OH-(aq)") - Cv0("H2O(aq)") ) );
    REQUIRE( rprops.dU0  == Approx(  U0("H+(aq)") +  U0("OH-(aq)") -  U0("H2O(aq)") ) );
    REQUIRE( rprops.dS0  == Approx(  S0("H+(aq)") +  S0("OH-(aq)") -  S0("H2O(aq)") ) );
    REQUIRE( rprops.dA0  == Approx(  A0("H+(aq)") +  A0("OH-(aq)") -  A0("H2O(aq)") ) );
    REQUIRE( rprops.lgK  == Approx( lgK(rprops.dG0) ) );

    reaction = db.reaction("CaMg(CO3)2(s) + 2*H+(aq) = Ca++(aq) + Mg++(aq) + 2*HCO3-(aq)");
    rprops = reaction.props(T, P);

    REQUIRE( rprops.dG0  == Approx(  G0("Ca++(aq)") +  G0("Mg++(aq)") + 2* G0("HCO3-(aq)") -  G0("CaMg(CO3)2(s)") - 2* G0("H+(aq)") ) );
    REQUIRE( rprops.dH0  == Approx(  H0("Ca++(aq)") +  H0("Mg++(aq)") + 2* H0("HCO3-(aq)") -  H0("CaMg(CO3)2(s)") - 2* H0("H+(aq)") ) );
    REQUIRE( rprops.dV0  == Approx(  V0("Ca++(aq)") +  V0("Mg++(aq)") + 2* V0("HCO3-(aq)") -  V0("CaMg(CO3)2(s)") - 2* V0("H+(aq)") ) );
    REQUIRE( rprops.dVT0 == Approx( VT0("Ca++(aq)") + VT0("Mg++(aq)") + 2*VT0("HCO3-(aq)") - VT0("CaMg(CO3)2(s)") - 2*VT0("H+(aq)") ) );
    REQUIRE( rprops.dVP0 == Approx( VP0("Ca++(aq)") + VP0("Mg++(aq)") + 2*VP0("HCO3-(aq)") - VP0("CaMg(CO3)2(s)") - 2*VP0("H+(aq)") ) );
    REQUIRE( rprops.dCp0 == Approx( Cp0("Ca++(aq)") + Cp0("Mg++(aq)") + 2*Cp0("HCO3-(aq)") - Cp0("CaMg(CO3)2(s)") - 2*Cp0("H+(aq)") ) );
    REQUIRE( rprops.dCv0 == Approx( Cv0("Ca++(aq)") + Cv0("Mg++(aq)") + 2*Cv0("HCO3-(aq)") - Cv0("CaMg(CO3)2(s)") - 2*Cv0("H+(aq)") ) );
    REQUIRE( rprops.dU0  == Approx(  U0("Ca++(aq)") +  U0("Mg++(aq)") + 2* U0("HCO3-(aq)") -  U0("CaMg(CO3)2(s)") - 2* U0("H+(aq)") ) );
    REQUIRE( rprops.dS0  == Approx(  S0("Ca++(aq)") +  S0("Mg++(aq)") + 2* S0("HCO3-(aq)") -  S0("CaMg(CO3)2(s)") - 2* S0("H+(aq)") ) );
    REQUIRE( rprops.dA0  == Approx(  A0("Ca++(aq)") +  A0("Mg++(aq)") + 2* A0("HCO3-(aq)") -  A0("CaMg(CO3)2(s)") - 2* A0("H+(aq)") ) );
    REQUIRE( rprops.lgK  == Approx( lgK(rprops.dG0) ) );
}
