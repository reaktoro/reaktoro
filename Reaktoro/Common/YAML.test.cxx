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
#include <Reaktoro/Common/YAML.hpp>
using namespace Reaktoro;

const std::string element = R"(
symbol: H
name: Hydrogen
atomicNumber: 1
atomicWeight: 0.001007940
electronegativity: 2.20
tags: [ group1 ]
)";

TEST_CASE("Testing YAML", "[YAML]")
{
    yaml y = yaml::parse(element);


    CHECK( y["symbol"].as<std::string>() == "H" );
    CHECK( y["name"].as<std::string>() == "Hydrogen" );
    CHECK( y["atomicNumber"].as<int>() == 1 );
    CHECK( y["atomicWeight"].as<double>() == 0.001007940 );
    CHECK( y["electronegativity"].as<double>() == 2.20 );
    CHECK( y["tags"].size() == 1 );
    CHECK( y["tags"][0].as<std::string>() == "group1" );
}

TEST_CASE("Testing YAML for NaN", "[YAML]")
{
    yaml y = yaml::parse("[.nan, .nan, .nan]");

    auto vals = y.as<Vec<double>>();

    CHECK(vals[0] != vals[0]);
    CHECK(vals[1] != vals[1]);
    CHECK(vals[2] != vals[2]);

    y = vals;

    vals = y.as<Vec<double>>();

    CHECK(vals[0] != vals[0]);
    CHECK(vals[1] != vals[1]);
    CHECK(vals[2] != vals[2]);
}
