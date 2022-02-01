// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/Types.hpp>
#include <Reaktoro/Common/StringUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing StringUtils", "[StringUtils]")
{
    Strings words;

    //-------------------------------------------------------------------------
    // TESTING METHOD: stringfy
    //-------------------------------------------------------------------------
    REQUIRE( stringfy(".", "a", "b") == "a.b");
    REQUIRE( stringfy(".", 1, 2, "abc") == "1.2.abc");

    words = {"alpha", "beta", "gamma"};

    //-------------------------------------------------------------------------
    // TESTING METHOD: str
    //-------------------------------------------------------------------------
    REQUIRE( str(words) == "alpha, beta, gamma");
    REQUIRE( str(1, 2, " abc") == "12 abc");

    //-------------------------------------------------------------------------
    // TESTING METHOD: replace
    //-------------------------------------------------------------------------
    REQUIRE( replace("CARBON DIOXIDE", " ", "-")   == "CARBON-DIOXIDE" );
    REQUIRE( replace("CARBON DIOXIDE", "O", "X")   == "CARBXN DIXXIDE" );
    REQUIRE( replace("CARBON DIOXIDE", "", "XYZ")  == "CARBON DIOXIDE" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: lowercase
    //-------------------------------------------------------------------------
    REQUIRE( lowercase("ABC 123 DEF")  == "abc 123 def" );
    REQUIRE( lowercase("AbC 123 dEf")  == "abc 123 def" );
    REQUIRE( lowercase("abc 123 def")  == "abc 123 def" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: uppercase
    //-------------------------------------------------------------------------
    REQUIRE( uppercase("ABC 123 DEF")  == "ABC 123 DEF" );
    REQUIRE( uppercase("AbC 123 dEf")  == "ABC 123 DEF" );
    REQUIRE( uppercase("abc 123 def")  == "ABC 123 DEF" );

    //-------------------------------------------------------------------------
    // TESTING METHOD: trimleft
    //-------------------------------------------------------------------------
    REQUIRE( trimleft("ABC"    )  == "ABC"   );
    REQUIRE( trimleft("ABC "   )  == "ABC "  );
    REQUIRE( trimleft("ABC  "  )  == "ABC  " );
    REQUIRE( trimleft(" ABC"   )  == "ABC"   );
    REQUIRE( trimleft("  ABC"  )  == "ABC"   );
    REQUIRE( trimleft(" ABC "  )  == "ABC "  );
    REQUIRE( trimleft("  ABC  ")  == "ABC  " );
    REQUIRE( trimright("       ")  == ""     );
    REQUIRE( trimright(""       )  == ""     );

    //-------------------------------------------------------------------------
    // TESTING METHOD: trimright
    //-------------------------------------------------------------------------
    REQUIRE( trimright("ABC"    )  == "ABC"   );
    REQUIRE( trimright("ABC "   )  == "ABC"   );
    REQUIRE( trimright("ABC  "  )  == "ABC"   );
    REQUIRE( trimright(" ABC"   )  == " ABC"  );
    REQUIRE( trimright("  ABC"  )  == "  ABC" );
    REQUIRE( trimright(" ABC "  )  == " ABC"  );
    REQUIRE( trimright("  ABC  ")  == "  ABC" );
    REQUIRE( trimright("       ")  == ""      );
    REQUIRE( trimright(""       )  == ""      );

    //-------------------------------------------------------------------------
    // TESTING METHOD: trim
    //-------------------------------------------------------------------------
    REQUIRE( trimright("ABC"    )  == "ABC"   );
    REQUIRE( trimright("ABC "   )  == "ABC"   );
    REQUIRE( trimright("ABC  "  )  == "ABC"   );
    REQUIRE( trimright(" ABC"   )  == " ABC"  );
    REQUIRE( trimright("  ABC"  )  == "  ABC" );
    REQUIRE( trimright(" ABC "  )  == " ABC"  );
    REQUIRE( trimright("  ABC  ")  == "  ABC" );
    REQUIRE( trimright("       ")  == ""      );
    REQUIRE( trimright(""       )  == ""      );

    //-------------------------------------------------------------------------
    // TESTING METHOD: split
    //-------------------------------------------------------------------------
    REQUIRE( split("ABC;DEF GHI", ""  )  == Strings{"ABC;DEF GHI"} );
    REQUIRE( split("ABC;DEF GHI", " " )  == Strings{"ABC;DEF", "GHI"} );
    REQUIRE( split("ABC;DEF GHI", ";" )  == Strings{"ABC", "DEF GHI"} );
    REQUIRE( split("ABC;DEF GHI", " ;")  == Strings{"ABC", "DEF", "GHI"} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: split with trim transform
    //-------------------------------------------------------------------------
    REQUIRE( split(  "ABC;DEF.GHI",   ";.", trim)  == Strings{"ABC", "DEF", "GHI"} );
    REQUIRE( split(  "ABC;DEF.GHI ",  ";.", trim)  == Strings{"ABC", "DEF", "GHI"} );
    REQUIRE( split( " ABC;DEF.GHI",   ";.", trim)  == Strings{"ABC", "DEF", "GHI"} );
    REQUIRE( split("  ABC;DEF.GHI  ", ";.", trim)  == Strings{"ABC", "DEF", "GHI"} );

    //-------------------------------------------------------------------------
    // TESTING METHOD: join
    //-------------------------------------------------------------------------
    REQUIRE( join({"ABC", "DEF", "GHI"})       == "ABC DEF GHI" );
    REQUIRE( join({"ABC", "DEF", "GHI"}, ".")  == "ABC.DEF.GHI" );
    REQUIRE( join({"ABC"})                     == "ABC"         );
    REQUIRE( join({"ABC"}, ".")                == "ABC"         );
    REQUIRE( join({})                          == ""            );

    //-------------------------------------------------------------------------
    // TESTING METHOD: tofloat
    //-------------------------------------------------------------------------
    REQUIRE( tofloat("123") == 123.0 );
    REQUIRE( tofloat("123.") == 123.0 );
    REQUIRE( tofloat("123.4") == 123.4 );
    REQUIRE( tofloat("123.45") == 123.45 );

    REQUIRE( tofloat("0") == 0.0 );
    REQUIRE( tofloat("0.") == 0.0 );
    REQUIRE( tofloat("0.123") == 0.123 );

    REQUIRE( tofloat("1e5") == 1.0e+5 );
    REQUIRE( tofloat("1.0e+5") == 1.0e+5 );
    REQUIRE( tofloat("1.0e-5") == 1.0e-5 );

    REQUIRE( tofloat("1E5") == 1.0e+5 );
    REQUIRE( tofloat("1.0E+5") == 1.0e+5 );
    REQUIRE( tofloat("1.0E-5") == 1.0e-5 );

    //-------------------------------------------------------------------------
    // TESTING METHOD: makeunique
    //-------------------------------------------------------------------------
    words = {"alpha", "beta", "gamma", "beta", "gamma", "beta", "gamma", "alpha"};
    words = makeunique(words, "!");

    REQUIRE( words[0] == "alpha"   );
    REQUIRE( words[1] == "beta"    );
    REQUIRE( words[2] == "gamma"   );
    REQUIRE( words[3] == "beta!"   );
    REQUIRE( words[4] == "gamma!"  );
    REQUIRE( words[5] == "beta!!"  );
    REQUIRE( words[6] == "gamma!!" );
    REQUIRE( words[7] == "alpha!"  );
}
