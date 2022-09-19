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

TEST_CASE("Testing class Table", "[Table]")
{
    SECTION("Testing basic usage of Table and TableColumn")
    {
        Table table;

        //----------------------------------------------------------------------------------------------------
        // Checking methods TableColumn::appendFloat|appendInteger|appendString|appendBoolean
        //----------------------------------------------------------------------------------------------------

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

        //----------------------------------------------------------------------------------------------------
        // Checking methods Table::rows and Table::cols
        //----------------------------------------------------------------------------------------------------

        CHECK(table.rows() == 5); // the number of rows in column Floats
        CHECK(table.cols() == 4); // the number of columns registered above

        //----------------------------------------------------------------------------------------------------
        // Checking method TableColumn::rows
        //----------------------------------------------------------------------------------------------------

        CHECK( table.column("Floats").rows()   == 5 );
        CHECK( table.column("Integers").rows() == 4 );
        CHECK( table.column("Strings").rows()  == 3 );
        CHECK( table.column("Booleans").rows() == 2 );

        //----------------------------------------------------------------------------------------------------
        // Checking methods TableColumn::dataType
        //----------------------------------------------------------------------------------------------------

        CHECK( table.column("Floats").dataType()   == TableColumn::DataType::Float );
        CHECK( table.column("Integers").dataType() == TableColumn::DataType::Integer );
        CHECK( table.column("Strings").dataType()  == TableColumn::DataType::String );
        CHECK( table.column("Booleans").dataType() == TableColumn::DataType::Boolean );

        //----------------------------------------------------------------------------------------------------
        // Checking conversion methods TableColumn::floats|integers|strings|booleans
        //----------------------------------------------------------------------------------------------------

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

        //----------------------------------------------------------------------------------------------------
        // Checking method TableColumn::append
        //----------------------------------------------------------------------------------------------------

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

        //----------------------------------------------------------------------------------------------------
        // Checking operator Table::operator[]
        //----------------------------------------------------------------------------------------------------

        // Check it on column of floats
        CHECK( table["Floats"][0] == 0.0 );
        CHECK( table["Floats"][1] == 1.0 );
        CHECK( table["Floats"][2] == 2.0 );
        CHECK( table["Floats"][3] == 3.0 );
        CHECK( table["Floats"][4] == 4.0 );

        // Check it on columns of non-floats (Table::operator[] only works with column of floats!)
        CHECK_THROWS( table["Integers"] );
        CHECK_THROWS( table["Strings"] );
        CHECK_THROWS( table["Booleans"] );

        //----------------------------------------------------------------------------------------------------
        // Checking method Table::dump
        //----------------------------------------------------------------------------------------------------
        CHECK( table.dump() ==
            " Floats | Integers | Strings | Booleans\n"
            "0.00000 |        2 |   Hello |        1\n"
            "1.00000 |        3 |   World |        0\n"
            "2.00000 |        4 |       ! |        1\n"
            "3.00000 |        5 |    Star |         \n"
            "4.00000 |        6 |         |         \n"
            "5.00000 |          |         |         \n"
            "6.00000 |          |         |         ");

        table.column("Column name with space").append(1.0);
        table.column("Column name with space").append(2.0);

        CHECK( table.dump() ==
            " Floats | Integers | Strings | Booleans | Column name with space\n"
            "0.00000 |        2 |   Hello |        1 |                1.00000\n"
            "1.00000 |        3 |   World |        0 |                2.00000\n"
            "2.00000 |        4 |       ! |        1 |                       \n"
            "3.00000 |        5 |    Star |          |                       \n"
            "4.00000 |        6 |         |          |                       \n"
            "5.00000 |          |         |          |                       \n"
            "6.00000 |          |         |          |                       ");

        Table::OutputOptions opts;
        opts.delimiter = " ";
        opts.precision = 4;
        opts.scientific = true;

        CHECK( table.dump(opts) ==
            "    Floats Integers Strings Booleans \"Column name with space\"\n"
            "0.0000e+00        2   Hello        1               1.0000e+00\n"
            "1.0000e+00        3   World        0               2.0000e+00\n"
            "2.0000e+00        4       !        1                         \n"
            "3.0000e+00        5    Star                                  \n"
            "4.0000e+00        6                                          \n"
            "5.0000e+00                                                   \n"
            "6.0000e+00                                                   ");
    }

    SECTION("Testing a more convenient usage of Table using operator<<")
    {
        Table table;

        //----------------------------------------------------------------------------------------------------
        // Checking operator TableColumn::operator<<
        //----------------------------------------------------------------------------------------------------

        table.column("Floats") << 10.0 << 20 << 30.0 << 40 << 50.0 << 60 << 70.0; // integers 20, 40, 60 will be cast to floats
        table.column("Integers") << 1 << 2 << 3 << 4 << 5;
        table.column("Strings") << "Hello" << String("World") << "!";
        table.column("Booleans") << true << false;

        CHECK(table.rows() == 7);
        CHECK(table.cols() == 4);

        CHECK( table.column("Floats").rows()   == 7 );
        CHECK( table.column("Integers").rows() == 5 );
        CHECK( table.column("Strings").rows()  == 3 );
        CHECK( table.column("Booleans").rows() == 2 );

        CHECK( table.column("Floats").dataType()   == TableColumn::DataType::Float );
        CHECK( table.column("Integers").dataType() == TableColumn::DataType::Integer );
        CHECK( table.column("Strings").dataType()  == TableColumn::DataType::String );
        CHECK( table.column("Booleans").dataType() == TableColumn::DataType::Boolean );

        CHECK( table.column("Floats").floats()     == Deque<double>{10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0} );
        CHECK( table.column("Integers").integers() == Deque<long>{1, 2, 3, 4, 5} );
        CHECK( table.column("Strings").strings()   == Deque<String>{"Hello", "World", "!"} );
        CHECK( table.column("Booleans").booleans() == Deque<bool>{true, false} );

        //----------------------------------------------------------------------------------------------------
        // Checking method TableColumn::dump
        //----------------------------------------------------------------------------------------------------

        CHECK( table.dump() ==
            " Floats | Integers | Strings | Booleans\n"
            "10.0000 |        1 |   Hello |        1\n"
            "20.0000 |        2 |   World |        0\n"
            "30.0000 |        3 |       ! |         \n"
            "40.0000 |        4 |         |         \n"
            "50.0000 |        5 |         |         \n"
            "60.0000 |          |         |         \n"
            "70.0000 |          |         |         ");
    }
}
