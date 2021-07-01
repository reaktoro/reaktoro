// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2021 Allan Leal
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
#include <Reaktoro/Common/Memoization.hpp>
using namespace Reaktoro;

TEST_CASE("Testing Memoization - memoizeLast with return", "[Memoization]")
{
    int counter = 0; // a counter for how many times f1 below has been fully evaluated

    auto f1 = [&](double x, int y, real z)
    {
        ++counter;
        return x * y * z;
    };

    auto f2 = memoizeLast(f1); // f2 is the memoized version of f1

    CHECK( f1(2.0, 4, 3.0) == f2(2.0, 4, 3.0) );

    CHECK( counter == 2 ); // incremented once in f1, incremented once in f2

    f2(2.0, 4, 3.0); // evaluate f2 with same arguments - counter should remain the same

    CHECK( counter == 2 ); // no increment in last f2 call

    Memoization::disable(); // disable memoization

    CHECK( Memoization::isEnabled() == false );
    CHECK( Memoization::isDisabled() == true );

    f2(2.0, 4, 3.0); // memoization is disabled, so counter will be incremented, even though same arguments used here again

    CHECK( counter == 3 ); // one increment in last f2 call

    Memoization::enable(); // enable memoization

    CHECK( Memoization::isEnabled() == true );
    CHECK( Memoization::isDisabled() == false );

    f2(2.0, 4, 3.0); // memoization is enabled and using same arguments here - cached result should be returned without increment in counter

    CHECK( counter == 3 ); // no increment in last f2 call

    CHECK( f1(3.0, 7, 1.0) == f2(3.0, 7, 1.0) );

    CHECK( counter == 5 ); // two increments above, in f1 and f2, because of different arguments
}

/// The result of the function tested in the next test case.
struct DummyResult
{
    real r;
    real s;
};

TEST_CASE("Testing Memoization - memoizeLast without return", "[Memoization]")
{
    int counter = 0; // a counter for how many times f1 below has been fully evaluated

    auto f1 = [&](DummyResult& res, double x, int y, real z) -> void
    {
        ++counter;
        res.r = x + y + z;
        res.s = x * y * z;
    };

    auto f2 = memoizeLastUsingRef(f1); // f2 is the memoized version of f1

    DummyResult res1;
    DummyResult res2;

    f1(res1, 2.0, 4, 3.0);
    f2(res2, 2.0, 4, 3.0);

    CHECK( res1.r == res2.r );
    CHECK( res1.s == res2.s );

    CHECK( counter == 2 ); // incremented once in f1, incremented once in f2

    f2(res2, 2.0, 4, 3.0); // evaluate f2 with same arguments - counter should remain the same

    CHECK( counter == 2 ); // no increment in last f2 call

    Memoization::disable(); // disable memoization

    CHECK( Memoization::isEnabled() == false );
    CHECK( Memoization::isDisabled() == true );

    f2(res2, 2.0, 4, 3.0); // memoization is disabled, so counter will be incremented, even though same arguments used here again

    CHECK( counter == 3 ); // one increment in last f2 call

    Memoization::enable(); // enable memoization

    CHECK( Memoization::isEnabled() == true );
    CHECK( Memoization::isDisabled() == false );

    f2(res2, 2.0, 4, 3.0); // memoization is enabled and using same arguments here - cached result should be returned without increment in counter

    CHECK( counter == 3 ); // no increment in last f2 call

    f1(res1, 3.0, 7, 1.0);
    f2(res2, 3.0, 7, 1.0);

    CHECK( res1.r == res2.r );
    CHECK( res1.s == res2.s );

    CHECK( counter == 5 ); // two increments above, in f1 and f2, because of different arguments
}
