// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2015 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

// Cute includes
#include <cute/cute.h>
#include <cute/cute_runner.h>
#include <cute/ide_listener.h>

// Reaktoro includes
#include <tests/Activity/TestActivity.hpp>
#include <tests/Common/TestCommon.hpp>
#include <tests/Core/TestCore.hpp>
#include <tests/Math/TestMath.hpp>
#include <tests/Optimization/TestOptimization.hpp>
using namespace Reaktoro;

int main(int argc, char **argv)
{
    cute::suite s;

    s += testSuiteCommon();
    s += testSuiteCore();
    s += testSuiteMath();
    s += testSuiteActivity();
    s += testSuiteOptimization();

    cute::ide_listener<> lis;
    cute::makeRunner(lis)(s, "Reaktoro tests");
}
