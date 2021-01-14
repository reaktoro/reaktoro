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
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumDims", "[EquilibriumDims]")
{
    ChemicalSystem system = test::createChemicalSystem();

    // EquilibriumConditions conditions(system);

    // conditions.internalEnergy(1.0, "kJ");
    // conditions.volume(1.0, "m3");

    // conditions.constantTemperature();
    // conditions.constantPressure();
    // conditions.constantEntropy();

    // conditions.pH(4.5);
    // conditions.fugacity("O2(g)", 1.0, "bar");

    // conditions.titrate("CO2");
    // conditions.titrate("CH4");
    // conditions.titrate("H2S");
    // conditions.titrate("H2");
    // conditions.titrate("CO");

    // conditions.cannotIncrease("NaCl(aq)");
    // conditions.cannotIncrease("H2(g)");
    // conditions.cannotDecrease("SiO2(s)");
    // conditions.cannotDecrease("O2(g)");
    // conditions.cannotReact("CaCO3(s)");
    // conditions.cannotReact("CO2(g)");
    // conditions.cannotReact("O2(aq) + H2(aq) = H2O(aq)");
    // conditions.cannotReact("CO2(g) = CO2(aq)");

    // auto dims = EquilibriumDims(conditions);

    // const auto Ne  = system.elements().size() + 1;
    // const auto Nn  = system.species().size();
    // const auto Npe = 2; // internalEnergy and volume - both equation constraints
    // const auto Npp = 3; // constantTemperature, constantPressure, constantEntropy - constant property constraints
    // const auto Np  = Npe + Npp;
    // const auto Nq  = 2; // pH and fugacity - both chemical potential constraints
    // const auto Nir = 2; // 2x calls to cannotReact(reaction)
    // const auto Ncv = 7; // 2x from the chemical potential constraints, H+ and O2, and 5x from titrants CO2, CH4, H2S, H2, CO
    // const auto Nx  = Nn + Np + Nq;
    // const auto Nc  = Ne + Nir;

    // REQUIRE( dims.Ne  == Ne  );
    // REQUIRE( dims.Nn  == Nn  );
    // REQUIRE( dims.Npe == Npe );
    // REQUIRE( dims.Npp == Npp );
    // REQUIRE( dims.Np  == Np  );
    // REQUIRE( dims.Nq  == Nq  );
    // REQUIRE( dims.Nir == Nir );
    // REQUIRE( dims.Ncv == Ncv );
    // REQUIRE( dims.Nx  == Nx  );
    // REQUIRE( dims.Nc  == Nc  );
}
