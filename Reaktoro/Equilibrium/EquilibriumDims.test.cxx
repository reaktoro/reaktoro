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
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ReactionEquation.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumDims", "[EquilibriumDims]")
{
    ChemicalSystem system = test::createChemicalSystem();

    EquilibriumConstraints constraints(system);

    constraints.control().temperature();
    constraints.control().pressure();

    constraints.until().internalEnergy(1.0, "kJ");

    constraints.preserve().volume();

    constraints.fix().pH(4.5);
    constraints.fix().fugacity("O2(g)", 1.0, "bar");

    constraints.prevent().fromIncreasing("NaCl(aq)");
    constraints.prevent().fromIncreasing("H2(g)");
    constraints.prevent().fromDecreasing("SiO2(s)");
    constraints.prevent().fromDecreasing("O2(g)");
    constraints.prevent().fromReacting("CaCO3(s)");
    constraints.prevent().fromReacting("CO2(g)");
    constraints.prevent().fromReacting("O2(aq) + H2(aq) = H2O(aq)");
    constraints.prevent().fromReacting("CO2(g) = CO2(aq)");

    auto dims = EquilibriumDims(constraints);

    const auto Ne  = system.elements().size() + 1;
    const auto Nn  = system.species().size();
    const auto Npe = 1; // 3x calls to until()
    const auto Npp = 1; // 4x calls to preserve()
    const auto Np  = Npe + Npp;
    const auto Nq  = 2; // 2x calls to fix()
    const auto Nir = 2; // 2x calls to prevent().fromReacting(reaction)
    const auto Ncv = 4; // 2x calls to control() and 2x calls to fix()
    const auto Nx  = Nn + Np + Nq;
    const auto Nc  = Ne + Nir;

    REQUIRE( dims.Ne  == Ne  );
    REQUIRE( dims.Nn  == Nn  );
    REQUIRE( dims.Npe == Npe );
    REQUIRE( dims.Npp == Npp );
    REQUIRE( dims.Np  == Np  );
    REQUIRE( dims.Nq  == Nq  );
    REQUIRE( dims.Nir == Nir );
    REQUIRE( dims.Ncv == Ncv );
    REQUIRE( dims.Nx  == Nx  );
    REQUIRE( dims.Nc  == Nc  );
}
