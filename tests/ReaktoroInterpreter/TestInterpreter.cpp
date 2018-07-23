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
#include <Reaktoro/Interpreter/Interpreter.hpp>
using namespace Reaktoro;

// cute includes
#include <cute/cute.h>
#include <cute/cute_runner.h>
#include <cute/ide_listener.h>

auto test_Equilibrium1() -> void
{
    const std::string s = R"xyz(
Equilibrium: 
   Temperature: 30 celsius
   Pressure: 10 bar
   Recipe: 1 kg H2O; 1 mmol NaCl
)xyz";

    Interpreter interpreter;
    interpreter.execute(s);
}

int main()
{
    cute::suite s;

    s += CUTE(test_Equilibrium1);

    cute::ide_listener<> lis;
    cute::makeRunner(lis)(s, "TestParseUtils");
}
