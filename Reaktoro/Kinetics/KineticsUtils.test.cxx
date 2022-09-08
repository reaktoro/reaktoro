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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
#include <Reaktoro/Kinetics/KineticsUtils.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing KineticsUtils", "[KineticsUtils]")
{
    SECTION("Testing method createEquilibriumSpecsForKinetics")
    {
        ChemicalSystem system = test::createChemicalSystem();

        EquilibriumSpecs specs(system);

        const auto Nr = system.reactions().size();

        WHEN("temperature and pressure are given")
        {
            specs.temperature();
            specs.pressure();

            specs = detail::createEquilibriumSpecsForKinetics(specs);

            CHECK( specs.numControlVariablesP() == Nr ); // the extent of reaction change variables Δξ (one for each reaction)

            auto Np = specs.numControlVariablesP();
            auto rconstraints = specs.assembleReactivityConstraints();

            CHECK( rconstraints.ids.size() == Nr );
            CHECK( rconstraints.Kn == system.stoichiometricMatrix().transpose() );
            CHECK( rconstraints.Kp == -identity(Nr, Nr) );
        }

        WHEN("volume and internal energy are given")
        {
            specs.volume();
            specs.internalEnergy();

            specs = detail::createEquilibriumSpecsForKinetics(specs);

            CHECK( specs.numControlVariablesP() == Nr + 2 ); // the extent of reaction change variables Δξ (one for each reaction), temperature and pressure

            auto Np = specs.numControlVariablesP();
            auto rconstraints = specs.assembleReactivityConstraints();

            CHECK( rconstraints.ids.size() == Nr );
            CHECK( rconstraints.Kn == system.stoichiometricMatrix().transpose() );
            CHECK( rconstraints.Kp.leftCols(2).isZero() );
            CHECK( rconstraints.Kp.rightCols(Nr) == -identity(Nr, Nr) );
        }
    }
}
