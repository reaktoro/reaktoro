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

#include "TestActivity.hpp"

// Reaktoro includes
#include <tests/Activity/TestActivityUtils.hpp>
#include <tests/Activity/TestAqueousActivityModel.hpp>
#include <tests/Activity/TestGaseousActivity.hpp>
#include <tests/Activity/TestMineralActivity.hpp>

namespace Reaktoro {

auto testSuiteActivity() -> cute::suite
{
    cute::suite s;

    s += testSuiteActivityUtils();
    s += testSuiteAqueousActivity();
    s += testSuiteGaseousActivity();
    s += testSuiteMineralActivity();

    return s;
}

} // namespace Reaktoro
