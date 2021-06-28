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

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    auto species = Species()
            .withName("AB2C3+2(aq)")
            .withFormula("AB2C3+2")
            .withSubstance("AB2C3+2")
            .withElements({{A, 1}, {B, 2}, {C, 3}})
            .withCharge(2.0)
            .withAggregateState(AggregateState::Aqueous)
            .withTags({"tag1", "tag2", "tag3"})
            .withAttachedData(String{"SomeData"})
    ;
    std::out << "Species name: " <<  species.name() << std::endl;
    std::out << "Species name: " <<  species.name() << std::endl;

/*
    SECTION("Testing attributes of the chemical species")
    {
        CHECK( species.name() == "AB2C3+2(aq)" );
        CHECK( species.formula() == "AB2C3+2" );
        CHECK( species.substance() == "AB2C3+2" );
        CHECK( species.elements().size() == 3 );
        CHECK( species.elements().coefficient("A") == 1 );
        CHECK( species.elements().coefficient("B") == 2 );
        CHECK( species.elements().coefficient("C") == 3 );
        CHECK( species.charge() == 2.0 );
        CHECK( species.tags().size() == 3.0 );
        CHECK( species.tags().at(0) == "tag1" );
        CHECK( species.tags().at(1) == "tag2" );
        CHECK( species.tags().at(2) == "tag3" );
        CHECK( species.attachedData().has_value() );
        CHECK( species.attachedData().type() == typeid(String) );
        CHECK( std::any_cast<String>(species.attachedData()) == "SomeData" );
    }
*/

}