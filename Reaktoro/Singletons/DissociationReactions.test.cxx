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
#include <Reaktoro/Singletons/DissociationReactions.hpp>
using namespace Reaktoro;

TEST_CASE("Testing DissociationReactions", "[DissociationReactions]")
{
    REQUIRE( DissociationReactions::size() == DissociationReactions::reactions().size() );

    REQUIRE( DissociationReactions::get("NaCl").has_value() );
    REQUIRE( DissociationReactions::get("CaCO3").has_value() );
    REQUIRE( DissociationReactions::get("Ca(CO3)").has_value() );

    REQUIRE( DissociationReactions::coefficient("NaCl", "Na+")   == 1.0 );
    REQUIRE( DissociationReactions::coefficient("NaCl", "Na[+]") == 1.0 );

    REQUIRE( DissociationReactions::coefficient("NaCl", "Cl-")   == 1.0 );
    REQUIRE( DissociationReactions::coefficient("NaCl", "Cl[-]") == 1.0 );

    REQUIRE( DissociationReactions::coefficient("NaCl", "H+")    == 0.0 );
    REQUIRE( DissociationReactions::coefficient("NaCl", "H[+]")  == 0.0 );

    REQUIRE( DissociationReactions::coefficient("CaCl2", "Ca++")   == 1.0 );
    REQUIRE( DissociationReactions::coefficient("CaCl2", "Ca+2")   == 1.0 );
    REQUIRE( DissociationReactions::coefficient("CaCl2", "Ca[2+]") == 1.0 );

    REQUIRE( DissociationReactions::coefficient("CaCl2", "Cl-") == 2.0 );
    REQUIRE( DissociationReactions::coefficient("CaCl2", "Cl[-]") == 2.0 );

    REQUIRE( DissociationReactions::coefficient("CaCl2", "H+") == 0.0 );
    REQUIRE( DissociationReactions::coefficient("CaCl2", "H[+]") == 0.0 );

    DissociationReactions::append({ "NaKCl2", {{1, "Na+"}, {1, "K+"}, {2, "Cl-"}} });

    REQUIRE( DissociationReactions::get("NaKClCl").has_value() );

    DissociationReactions::append({ "NaCl", {} }); // Remove ions in NaCl dissociation reactions

    REQUIRE( DissociationReactions::get("NaCl").has_value() );
    REQUIRE( DissociationReactions::get("NaCl").value().ions.size() == 0 );
}
