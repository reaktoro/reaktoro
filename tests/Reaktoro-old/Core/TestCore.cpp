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

#include "TestCore.hpp"

// Reaktoro includes
#include <tests/Core/TestSpecies.hpp>
#include <tests/Core/TestChemicalSystem.hpp>
#include <tests/Core/TestPhase.hpp>
#include <tests/Core/TestPartition.hpp>

namespace Reaktoro {

auto testSuiteCore() -> cute::suite
{
    cute::suite s;

    s += testSuiteSpecies();
    s += testSuitePhase();
    s += testSuiteChemicalSystem();
    s += testSuitePartition();

    return s;
}

} // namespace Reaktoro
