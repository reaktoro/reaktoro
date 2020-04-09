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
#include <Reaktoro/Core/AggregateState.hpp>
using namespace Reaktoro;

TEST_CASE("Testing AggregateState module", "[AggregateState]")
{
    CHECK( parseAggregateState("g")   == AggregateState::Gas              );
    CHECK( parseAggregateState("l")   == AggregateState::Liquid           );
    CHECK( parseAggregateState("s")   == AggregateState::Solid            );
    CHECK( parseAggregateState("pl")  == AggregateState::Plasma           );
    CHECK( parseAggregateState("cd")  == AggregateState::CondensedPhase   );
    CHECK( parseAggregateState("fl")  == AggregateState::Fluid            );
    CHECK( parseAggregateState("lc")  == AggregateState::LiquidCrystal    );
    CHECK( parseAggregateState("cr")  == AggregateState::CrystallineSolid );
    CHECK( parseAggregateState("am")  == AggregateState::AmorphousSolid   );
    CHECK( parseAggregateState("vit") == AggregateState::Vitreous         );
    CHECK( parseAggregateState("ads") == AggregateState::Adsorbed         );
    CHECK( parseAggregateState("mon") == AggregateState::Monomeric        );
    CHECK( parseAggregateState("pol") == AggregateState::Polymeric        );
    CHECK( parseAggregateState("ss")  == AggregateState::SolidSolution    );
    CHECK( parseAggregateState("ex")  == AggregateState::IonExchange      );
    CHECK( parseAggregateState("aq")  == AggregateState::Aqueous          );
    CHECK( parseAggregateState("xy")  == AggregateState::Undefined        );

    CHECK( identifyAggregateState("XYZ(g)")       == AggregateState::Gas              );
    CHECK( identifyAggregateState("XYZ(l)")       == AggregateState::Liquid           );
    CHECK( identifyAggregateState("XYZ(s)")       == AggregateState::Solid            );
    CHECK( identifyAggregateState("XYZ(s, xyz)")  == AggregateState::Solid            );
    CHECK( identifyAggregateState("XYZ(pl)")      == AggregateState::Plasma           );
    CHECK( identifyAggregateState("XYZ(cd)")      == AggregateState::CondensedPhase   );
    CHECK( identifyAggregateState("XYZ(fl)")      == AggregateState::Fluid            );
    CHECK( identifyAggregateState("XYZ(lc)")      == AggregateState::LiquidCrystal    );
    CHECK( identifyAggregateState("XYZ(cr)")      == AggregateState::CrystallineSolid );
    CHECK( identifyAggregateState("XYZ(am)")      == AggregateState::AmorphousSolid   );
    CHECK( identifyAggregateState("XYZ(vit)")     == AggregateState::Vitreous         );
    CHECK( identifyAggregateState("XYZ(ads)")     == AggregateState::Adsorbed         );
    CHECK( identifyAggregateState("XYZ(mon)")     == AggregateState::Monomeric        );
    CHECK( identifyAggregateState("XYZ(pol)")     == AggregateState::Polymeric        );
    CHECK( identifyAggregateState("XYZ(ss)")      == AggregateState::SolidSolution    );
    CHECK( identifyAggregateState("XYZ(ex)")      == AggregateState::IonExchange      );
    CHECK( identifyAggregateState("XYZ(aq)")      == AggregateState::Aqueous          );
    CHECK( identifyAggregateState("XYZ(aq, uvw)") == AggregateState::Aqueous          );

    CHECK( identifyAggregateState("XYZ-")         == AggregateState::Aqueous          );
    CHECK( identifyAggregateState("XYZ--")        == AggregateState::Aqueous          );
    CHECK( identifyAggregateState("XYZ---")       == AggregateState::Aqueous          );
    CHECK( identifyAggregateState("XYZ-2")        == AggregateState::Aqueous          );
    CHECK( identifyAggregateState("XYZ-3")        == AggregateState::Aqueous          );

    CHECK( identifyAggregateState("XYZ+")         == AggregateState::Aqueous          );
    CHECK( identifyAggregateState("XYZ++")        == AggregateState::Aqueous          );
    CHECK( identifyAggregateState("XYZ+++")       == AggregateState::Aqueous          );
    CHECK( identifyAggregateState("XYZ+2")        == AggregateState::Aqueous          );
    CHECK( identifyAggregateState("XYZ+3")        == AggregateState::Aqueous          );
    CHECK( identifyAggregateState("XYZ[2-]")      == AggregateState::Aqueous          );

    CHECK( identifyAggregateState("XYZ-(pl)")     == AggregateState::Plasma           );
    CHECK( identifyAggregateState("XYZ--(pl)")    == AggregateState::Plasma           );
    CHECK( identifyAggregateState("XYZ---(pl)")   == AggregateState::Plasma           );
    CHECK( identifyAggregateState("XYZ-2(pl)")    == AggregateState::Plasma           );
    CHECK( identifyAggregateState("XYZ-3(pl)")    == AggregateState::Plasma           );

    CHECK( identifyAggregateState("XYZ+(pl)")     == AggregateState::Plasma           );
    CHECK( identifyAggregateState("XYZ++(pl)")    == AggregateState::Plasma           );
    CHECK( identifyAggregateState("XYZ+++(pl)")   == AggregateState::Plasma           );
    CHECK( identifyAggregateState("XYZ[3+](pl)")  == AggregateState::Plasma           );
    CHECK( identifyAggregateState("XYZ+2(pl)")    == AggregateState::Plasma           );
    CHECK( identifyAggregateState("XYZ+3(pl)")    == AggregateState::Plasma           );

    CHECK( identifyAggregateState("XYZ")          == AggregateState::Undefined        );
    CHECK( identifyAggregateState("XYZ(xy)")      == AggregateState::Undefined        );
}
