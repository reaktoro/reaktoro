# Reaktoro is a unified framework for modeling chemically reactive systems.
#
# Copyright © 2014-2022 Allan Leal
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library. If not, see <http://www.gnu.org/licenses/>.


from reaktoro import *
import pytest


def testTable():

    # =======================================================================================================
    # Testing basic usage of Table and TableColumn
    # =======================================================================================================
    table = Table()

    #----------------------------------------------------------------------------------------------------
    # Checking methods TableColumn.appendFloat|appendInteger|appendString|appendBoolean
    #----------------------------------------------------------------------------------------------------

    table.column("Floats").appendFloat(0.0)
    table.column("Floats").appendFloat(1.0)
    table.column("Floats").appendFloat(2.0)
    table.column("Floats").appendFloat(3.0)
    table.column("Floats").appendFloat(4.0)

    table.column("Integers").appendInteger(2)
    table.column("Integers").appendInteger(3)
    table.column("Integers").appendInteger(4)
    table.column("Integers").appendInteger(5)

    table.column("Strings").appendString("Hello")
    table.column("Strings").appendString("World")
    table.column("Strings").appendString("!")

    table.column("Booleans").appendBoolean(True)
    table.column("Booleans").appendBoolean(False)

    #----------------------------------------------------------------------------------------------------
    # Checking methods Table.rows and Table.cols
    #----------------------------------------------------------------------------------------------------

    assert table.rows() == 5  # the number of rows in column Floats
    assert table.cols() == 4  # the number of columns registered above

    #----------------------------------------------------------------------------------------------------
    # Checking method TableColumn.rows
    #----------------------------------------------------------------------------------------------------

    assert table.column("Floats").rows()   == 5
    assert table.column("Integers").rows() == 4
    assert table.column("Strings").rows()  == 3
    assert table.column("Booleans").rows() == 2

    #----------------------------------------------------------------------------------------------------
    # Checking methods TableColumn.dataType
    #----------------------------------------------------------------------------------------------------

    assert table.column("Floats").dataType()   == TableColumn.DataType.Float
    assert table.column("Integers").dataType() == TableColumn.DataType.Integer
    assert table.column("Strings").dataType()  == TableColumn.DataType.String
    assert table.column("Booleans").dataType() == TableColumn.DataType.Boolean

    #----------------------------------------------------------------------------------------------------
    # Checking conversion methods TableColumn.floats|integers|strings|booleans
    #----------------------------------------------------------------------------------------------------

    assert table.column("Floats").floats()     == [0.0, 1.0, 2.0, 3.0, 4.0]
    assert table.column("Integers").integers() == [2, 3, 4, 5]
    assert table.column("Strings").strings()   == ["Hello", "World", "!"]
    assert table.column("Booleans").booleans() == [True, False]

    with pytest.raises(Exception): table.column("Floats").integers()
    with pytest.raises(Exception): table.column("Floats").strings()
    with pytest.raises(Exception): table.column("Floats").booleans()

    with pytest.raises(Exception): table.column("Integers").floats()
    with pytest.raises(Exception): table.column("Integers").strings()
    with pytest.raises(Exception): table.column("Integers").booleans()

    with pytest.raises(Exception): table.column("Strings").floats()
    with pytest.raises(Exception): table.column("Strings").integers()
    with pytest.raises(Exception): table.column("Strings").booleans()

    with pytest.raises(Exception): table.column("Booleans").floats()
    with pytest.raises(Exception): table.column("Booleans").integers()
    with pytest.raises(Exception): table.column("Booleans").strings()

    #----------------------------------------------------------------------------------------------------
    # Checking method TableColumn.append
    #----------------------------------------------------------------------------------------------------

    try: table.column("Floats").append(5)  # integer 5 can be inserted with append in a column of floats!
    except Exception: pytest.fail("An exception was not expected here.")

    try: table.column("Floats").append(6.0)  # float 5.0 can be inserted with append in a column of floats - no problem!
    except Exception: pytest.fail("An exception was not expected here.")

    with pytest.raises(Exception): table.column("Floats").appendInteger(5)  # integer 5 cannot be inserted with appendInteger in a column of floats!
    with pytest.raises(Exception): table.column("Floats").append("NOT OK")
    with pytest.raises(Exception): table.column("Floats").append(True)

    try: table.column("Integers").append(6)  # an integer can be appended to a column of integers - no problem!
    except Exception: pytest.fail("An exception was not expected here.")

    with pytest.raises(Exception): table.column("Integers").append(1.0)  # a float cannot be appended to a column of integers!
    with pytest.raises(Exception): table.column("Integers").append("NOT OK")
    with pytest.raises(Exception): table.column("Integers").append(False)

    try: table.column("Strings").append("Star")  # a string can be appended to a column of strings - no problem!
    except Exception: pytest.fail("An exception was not expected here.")

    with pytest.raises(Exception): table.column("Strings").append(1.0)
    with pytest.raises(Exception): table.column("Strings").append(3)
    with pytest.raises(Exception): table.column("Strings").append(False)

    try: table.column("Booleans").append(True)  # a boolean can be appended to a column of booleans - no problem!
    except Exception: pytest.fail("An exception was not expected here.")

    with pytest.raises(Exception): table.column("Booleans").append(0.0)
    with pytest.raises(Exception): table.column("Booleans").append(8)
    with pytest.raises(Exception): table.column("Booleans").append("False")

    #----------------------------------------------------------------------------------------------------
    # Checking operator Table.operator[]
    #----------------------------------------------------------------------------------------------------

    # Check it on column of floats
    assert table["Floats"][0] == 0.0
    assert table["Floats"][1] == 1.0
    assert table["Floats"][2] == 2.0
    assert table["Floats"][3] == 3.0
    assert table["Floats"][4] == 4.0

    # Check it on columns of non-floats (Table.operator[] only works with column of floats!)
    with pytest.raises(Exception): table["Integers"]
    with pytest.raises(Exception): table["Strings"]
    with pytest.raises(Exception): table["Booleans"]

    #----------------------------------------------------------------------------------------------------
    # Checking method Table.dump
    #----------------------------------------------------------------------------------------------------

    assert table.dump() == (
        " Floats | Integers | Strings | Booleans\n"
        "0.00000 |        2 |   Hello |        1\n"
        "1.00000 |        3 |   World |        0\n"
        "2.00000 |        4 |       ! |        1\n"
        "3.00000 |        5 |    Star |         \n"
        "4.00000 |        6 |         |         \n"
        "5.00000 |          |         |         ")

    table.column("Column name with space").append(1.0)
    table.column("Column name with space").append(2.0)

    assert table.dump() == (
        " Floats | Integers | Strings | Booleans | Column name with space\n"
        "0.00000 |        2 |   Hello |        1 |                1.00000\n"
        "1.00000 |        3 |   World |        0 |                2.00000\n"
        "2.00000 |        4 |       ! |        1 |                       \n"
        "3.00000 |        5 |    Star |          |                       \n"
        "4.00000 |        6 |         |          |                       \n"
        "5.00000 |          |         |          |                       ")

    opts = Table.OutputOptions()
    opts.delimiter = " | "
    opts.precision = 4
    opts.scientific = True

    assert table.dump(opts) == (
        "    Floats | Integers | Strings | Booleans | Column name with space\n"
        "0.0000e+00 |        2 |   Hello |        1 |             1.0000e+00\n"
        "1.0000e+00 |        3 |   World |        0 |             2.0000e+00\n"
        "2.0000e+00 |        4 |       ! |        1 |                       \n"
        "3.0000e+00 |        5 |    Star |          |                       \n"
        "4.0000e+00 |        6 |         |          |                       \n"
        "5.0000e+00 |          |         |          |                       ")

    # =======================================================================================================
    # Testing a more convenient usage of Table using operator<<
    # =======================================================================================================
    table = Table()

    #----------------------------------------------------------------------------------------------------
    # Checking operator TableColumn.operator<<
    #----------------------------------------------------------------------------------------------------

    table.column("Floats") << 10.0 << 20 << 30.0 << 40 << 50.0 << 60 << 70.0  # integers 20, 40, 60 will be cast to floats
    table.column("Integers") << 1 << 2 << 3 << 4 << 5
    table.column("Strings") << "Hello" << "World" << "!"
    table.column("Booleans") << True << False

    assert table.rows() == 7
    assert table.cols() == 4

    assert table.column("Floats").rows()   == 7
    assert table.column("Integers").rows() == 5
    assert table.column("Strings").rows()  == 3
    assert table.column("Booleans").rows() == 2

    assert table.column("Floats").dataType()   == TableColumn.DataType.Float
    assert table.column("Integers").dataType() == TableColumn.DataType.Integer
    assert table.column("Strings").dataType()  == TableColumn.DataType.String
    assert table.column("Booleans").dataType() == TableColumn.DataType.Boolean

    assert table.column("Floats").floats()     == [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
    assert table.column("Integers").integers() == [1, 2, 3, 4, 5]
    assert table.column("Strings").strings()   == ["Hello", "World", "!"]
    assert table.column("Booleans").booleans() == [True, False]

    #----------------------------------------------------------------------------------------------------
    # Checking method TableColumn.dump
    #----------------------------------------------------------------------------------------------------

    assert table.dump() == (
        " Floats | Integers | Strings | Booleans\n"
        "10.0000 |        1 |   Hello |        1\n"
        "20.0000 |        2 |   World |        0\n"
        "30.0000 |        3 |       ! |         \n"
        "40.0000 |        4 |         |         \n"
        "50.0000 |        5 |         |         \n"
        "60.0000 |          |         |         ")
