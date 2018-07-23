// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2018 Allan Leal
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

// Reaktoro includes
#include <Reaktoro/Interpreter/Utils.hpp>
using namespace Reaktoro;

// cute includes
#include <cute/cute.h>
#include <cute/cute_runner.h>
#include <cute/ide_listener.h>

auto test_preprocess() -> void
{
    std::string input = R"xyz(
KineticPath StateFinal:
   InitialCondition: StateIC
   KineticSpecies: Calcite Dolomite
   Duration: 24 hours
   Plot:
       Name: MolalityCa
       x: t
       y: molality element=Ca units=mmolal
       xlabel: t [hour]
       ylabel: Concentration [mmolal]
       ytitles: Ca
       Key: right center
   Plot:
       Name: Calcite
       x: t
       y: mass species=Calcite units=g
       xlabel: t [hour]
       ylabel: Concentration [mmolal]
       ytitles: Calcite
       Key: right center
)xyz";

    std::string output = R"xyz(
- KineticPath StateFinal:
   - InitialCondition: StateIC
   - KineticSpecies: Calcite Dolomite
   - Duration: 24 hours
   - Plot:
       - Name: MolalityCa
       - x: t
       - y: molality element=Ca units=mmolal
       - xlabel: t [hour]
       - ylabel: Concentration [mmolal]
       - ytitles: Ca
       - Key: right center
   - Plot:
       - Name: Calcite
       - x: t
       - y: mass species=Calcite units=g
       - xlabel: t [hour]
       - ylabel: Concentration [mmolal]
       - ytitles: Calcite
       - Key: right center
)xyz";

    ASSERT_EQUAL(output, preprocess(input));
}
int main()
{
    cute::suite s;

    s += CUTE(test_preprocess);

    cute::ide_listener<> lis;
    cute::makeRunner(lis)(s, "TestParseUtils");
}
