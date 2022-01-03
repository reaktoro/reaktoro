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
#include <Reaktoro/Extensions/Nasa/NasaDatabase.hpp>
using namespace Reaktoro;

TEST_CASE("Testing NasaDatabase module", "[NasaDatabase]")
{
    NasaDatabase db("cea");

    CHECK( db.species().size() == 2085 );

    CHECK( db.species().findWithName("ALF+") );
    CHECK( db.species().findWithName("B3H9") );
    CHECK( db.species().findWithName("C3S2") );
    CHECK( db.species().findWithName("Zr+") );
    CHECK( db.species().findWithName("Ti+") );
    CHECK( db.species().findWithName("CO2") );
    CHECK( db.species().findWithName("PbBr2(cr)") );
    CHECK( db.species().findWithName("Jet-A(g)") );

    const auto Tr = 298.15; // in K
    const auto Pr = 1.0e5;  // in Pa

    CHECK(db.species().getWithName("CO2").props(Tr, Pr).H0 == Approx( -393510.000) );
    CHECK(db.species().getWithName("Zr+").props(Tr, Pr).H0 == Approx( 1246246.292) );
    CHECK(db.species().getWithName("Ti+").props(Tr, Pr).H0 == Approx( 1137624.029) );
}
