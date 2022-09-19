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
#include <Reaktoro/Common/Table.hpp>
using namespace Reaktoro;

TEST_CASE("Testing basic usage of Table", "[Table]")
{
    Table table;

    table.column("Floats").appendFloat(0.0);
    table.column("Floats").appendFloat(1.0);
    table.column("Floats").appendFloat(2.0);
    table.column("Floats").appendFloat(3.0);
    table.column("Floats").appendFloat(4.0);

    table.column("Integers").appendInteger(2);
    table.column("Integers").appendInteger(3);
    table.column("Integers").appendInteger(4);
    table.column("Integers").appendInteger(5);

    table.column("Strings").appendString("Hello");
    table.column("Strings").appendString("World");
    table.column("Strings").appendString("!");

    table.column("Booleans").appendBoolean(true);
    table.column("Booleans").appendBoolean(false);

    CHECK( table.column("Floats").rows()   == 5 );
    CHECK( table.column("Integers").rows() == 4 );
    CHECK( table.column("Strings").rows()  == 3 );
    CHECK( table.column("Booleans").rows() == 2 );

    CHECK( table.column("Floats").dataType()   == TableColumn::DataType::Float );
    CHECK( table.column("Integers").dataType() == TableColumn::DataType::Integer );
    CHECK( table.column("Strings").dataType()  == TableColumn::DataType::String );
    CHECK( table.column("Booleans").dataType() == TableColumn::DataType::Boolean );

    CHECK( table.column("Floats").floats()     == Deque<double>{0.0, 1.0, 2.0, 3.0, 4.0} );
    CHECK( table.column("Integers").integers() == Deque<long>{2, 3, 4, 5} );
    CHECK( table.column("Strings").strings()   == Deque<String>{"Hello", "World", "!"} );
    CHECK( table.column("Booleans").booleans() == Deque<bool>{true, false} );

    CHECK_THROWS( table.column("Floats").integers() );
    CHECK_THROWS( table.column("Floats").strings() );
    CHECK_THROWS( table.column("Floats").booleans() );

    CHECK_THROWS( table.column("Integers").floats() );
    CHECK_THROWS( table.column("Integers").strings() );
    CHECK_THROWS( table.column("Integers").booleans() );

    CHECK_THROWS( table.column("Strings").floats() );
    CHECK_THROWS( table.column("Strings").integers() );
    CHECK_THROWS( table.column("Strings").booleans() );

    CHECK_THROWS( table.column("Booleans").floats() );
    CHECK_THROWS( table.column("Booleans").integers() );
    CHECK_THROWS( table.column("Booleans").strings() );

    CHECK_NOTHROW( table.column("Floats").append(5) ); // integer 5 can be inserted with append in a column of floats!
    CHECK_NOTHROW( table.column("Floats").append(6.0) ); // float 5.0 can be inserted with append in a column of floats - no problem!
    CHECK_THROWS( table.column("Floats").appendInteger(5) ); // integer 5 cannot be inserted with appendInteger in a column of floats!
    CHECK_THROWS( table.column("Floats").append("NOT OK") );
    CHECK_THROWS( table.column("Floats").append(true) );

    CHECK_NOTHROW( table.column("Integers").append(6) ); // an integer can be appended to a column of integers - no problem!
    CHECK_THROWS( table.column("Integers").append(1.0) ); // a float cannot be appended to a column of integers!
    CHECK_THROWS( table.column("Integers").append("NOT OK") );
    CHECK_THROWS( table.column("Integers").append(false) );

    CHECK_NOTHROW( table.column("Strings").append("Star") ); // a string can be appended to a column of strings - no problem!
    CHECK_THROWS( table.column("Strings").append(1.0) );
    CHECK_THROWS( table.column("Strings").append(3) );
    CHECK_THROWS( table.column("Strings").append(false) );

    CHECK_NOTHROW( table.column("Booleans").append(true) ); // a boolean can be appended to a column of booleans - no problem!
    CHECK_THROWS( table.column("Booleans").append(0.0) );
    CHECK_THROWS( table.column("Booleans").append(8) );
    CHECK_THROWS( table.column("Booleans").append("False") );
}

TEST_CASE("Testing convenient usage of Table", "[Table]")
{
    Table table;

    table.column("Floats") << 0.0;
    table.column("Floats") << 1.0;
    table.column("Floats") << 2.0;
    table.column("Floats") << 3.0;
    table.column("Floats") << 4.0;

    table.column("OtherFloats") << 2 << 3 << 4 << 5; // note operator<< will append integers as floats - check docs for reason!
    table.column("Strings") << "Hello" << String("World") << "!";
    table.column("Booleans") << true << false;

    CHECK( table.column("Floats").rows()      == 5 );
    CHECK( table.column("OtherFloats").rows() == 4 );
    CHECK( table.column("Strings").rows()     == 3 );
    CHECK( table.column("Booleans").rows()    == 2 );

    CHECK( table.column("Floats").dataType()      == TableColumn::DataType::Float );
    CHECK( table.column("OtherFloats").dataType() == TableColumn::DataType::Float );
    CHECK( table.column("Strings").dataType()     == TableColumn::DataType::String );
    CHECK( table.column("Booleans").dataType()    == TableColumn::DataType::Boolean );

    CHECK( table.column("Floats").floats()      == Deque<double>{0.0, 1.0, 2.0, 3.0, 4.0} );
    CHECK( table.column("OtherFloats").floats() == Deque<double>{2.0, 3.0, 4.0, 5.0} );
    CHECK( table.column("Strings").strings()    == Deque<String>{"Hello", "World", "!"} );
    CHECK( table.column("Booleans").booleans()  == Deque<bool>{true, false} );

    CHECK( table.dump() == "" );
}
