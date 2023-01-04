// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2021 Allan Leal
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

// C++ includes
#include <iomanip>

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Extensions/Nasa.hpp>
#include <Reaktoro/Equilibrium.hpp>
#include <Reaktoro/Kinetics.hpp>
using namespace Reaktoro;

TEST_CASE("Testing KineticsSolver", "[KineticsSolver]")
{
    NasaDatabase db("nasa-cea");

    auto ratefn = [](ChemicalProps const& props)
    {
        const auto k0 = 0.01;
        const auto nc = props.speciesAmount("C(gr)");
        return k0 * nc;
    };

    ChemicalSystem system(db,
        CondensedPhase("C(gr)"),
        GaseousPhase("O2 CO2"),
        GeneralReaction("C(gr) + O2 = CO2").setRateModel(ratefn)
    );

    ChemicalState state(system);
    state.set("C(gr)", 1.0, "mol");
    state.set("O2", 1.0, "mol");

    ChemicalProps props(state);

    SECTION("When temperature and pressure are constant")
    {
        EquilibriumSpecs specs = EquilibriumSpecs::TP(system);

        KineticsSolver solver(specs);

        KineticsOptions options;
        options.optima.output.active = true;
        solver.setOptions(options);

        const auto dt = 1.0;

        auto res = solver.solve(state, dt);

        REQUIRE( res.succeeded() );

        CHECK( state.speciesAmount("C(gr)") == Approx(0.990099) );
    }

    SECTION("When enthalpy and pressure are constant")
    {
        EquilibriumSpecs specs = EquilibriumSpecs::HP(system);

        KineticsSolver solver(specs);

        KineticsOptions options;
        options.optima.output.active = true;
        solver.setOptions(options);

        const auto dt = 1.0;

        EquilibriumConditions conditions(specs);
        conditions.pressure(props.pressure());
        conditions.enthalpy(props.enthalpy());

        auto res = solver.solve(state, dt, conditions);

        REQUIRE( res.succeeded() );

        CHECK( state.speciesAmount("C(gr)") == Approx(0.990099) );
    }
}
