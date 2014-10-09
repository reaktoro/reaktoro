// Reaktor is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
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

#include "TestAlgorithmUtils.hpp"

// Reaktor includes
#include <Reaktor/Reaktor.hpp>

namespace Reaktor {
namespace {

auto test_acceptable() -> void
{

}

auto test_extend() -> void
{

}

auto test_largestStep() -> void
{

}

auto test_fractionToTheBoundary() -> void
{

}

auto test_lessThan() -> void
{

}

auto test_greaterThan() -> void
{

}

} // namespace

auto testSuiteAlgorithmUtils() -> cute::suite
{
    cute::suite s;

    s += CUTE(test_acceptable);
    s += CUTE(test_extend);
    s += CUTE(test_largestStep);
    s += CUTE(test_fractionToTheBoundary);
    s += CUTE(test_lessThan);
    s += CUTE(test_greaterThan);

    return s;
}

} // namespace Reaktor
