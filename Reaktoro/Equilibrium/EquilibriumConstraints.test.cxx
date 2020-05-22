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
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumConstraints", "[EquilibriumConstraints]")
{
    ChemicalSystem system = test::createChemicalSystem();

    EquilibriumConstraints constraints(system);

    SECTION("Testing method EquilibriumConstraints::control")
    {
        // Reference to the controls data
        const auto& controls = constraints.data().controls;

        REQUIRE( controls.T == false );
        constraints.control().temperature();
        REQUIRE( controls.T == true );

        REQUIRE( controls.P == false );
        constraints.control().pressure();
        REQUIRE( controls.P == true );

        REQUIRE( constraints.numControlVariables() == 2 );

        REQUIRE_THROWS( constraints.control().titrationOf("HCl") );               // not supported yet!
        REQUIRE_THROWS( constraints.control().titrationOfEither("HCl", "NaOH") ); // not supported yet!

        // REQUIRE( controls.titrants.size() == 0 );
        // constraints.control().titrationOf("HCl");
        // REQUIRE( controls.titrants.size() == 1 );
        // REQUIRE( controls.titrants[0].str() == "HCl" );

        // constraints.control().titrationOf("NaOH");
        // REQUIRE( controls.titrants.size() == 2 );
        // REQUIRE( controls.titrants[0].str() == "HCl" );
        // REQUIRE( controls.titrants[1].str() == "NaOH" );
    }

    SECTION("Testing method EquilibriumConstraints::until")
    {
        // Reference to the functional equilibrium constraints data
        const auto& econstraints = constraints.data().econstraints;

        ChemicalProps props(system);

        const auto t = VectorXr{{ 1.0, 2.0 }};

        const auto titrants = Map<String, real>{{"HCl", 2.1}, {"NaOH", 3.1}};

        EquilibriumEquationArgs args{props, t, titrants};

        constraints.until().internalEnergy(1.0, "kJ");

        REQUIRE( econstraints.size() == 1 );
        REQUIRE( econstraints[0](args) == -1.0e+3 );

        constraints.until().volume(1.0, "mm3");

        REQUIRE( econstraints.size() == 2 );
        REQUIRE( econstraints[0](args) == -1.0e+3 );
        REQUIRE( econstraints[1](args) == -1.0e-9 );

        constraints.until().enthalpy(1.0, "J");

        REQUIRE( econstraints.size() == 3 );
        REQUIRE( econstraints[0](args) == -1.0e+3 );
        REQUIRE( econstraints[1](args) == -1.0e-9 );
        REQUIRE( econstraints[2](args) == -1.0 );

        REQUIRE( constraints.numEquationConstraints() == 3 );
    }

    SECTION("Testing method EquilibriumConstraints::preserve")
    {
        // Reference to the chemical property preservation constraints data
        const auto& pconstraints = constraints.data().pconstraints;

        ChemicalProps props(system);

        constraints.preserve().internalEnergy();

        REQUIRE( pconstraints.size() == 1 );
        REQUIRE( pconstraints[0](props) == 0.0 ); // TODO: Create a mock ChemicalProps object with non-zero properties to ensure the correct method is called.

        constraints.preserve().volume();

        REQUIRE( pconstraints.size() == 2 );
        REQUIRE( pconstraints[0](props) == 0.0 );
        REQUIRE( pconstraints[1](props) == 0.0 );

        constraints.preserve().enthalpy();

        REQUIRE( pconstraints.size() == 3 );
        REQUIRE( pconstraints[0](props) == 0.0 );
        REQUIRE( pconstraints[1](props) == 0.0 );
        REQUIRE( pconstraints[2](props) == 0.0 );

        REQUIRE( constraints.numPropertyPreservationConstraints() == 3 );
    }

    SECTION("Testing method EquilibriumConstraints::fix")
    {
        // Reference to the chemical potential constraints data
        const auto& uconstraints = constraints.data().uconstraints;

        ChemicalProps props(system);

        const auto R = universalGasConstant;

        const auto T = props.temperature();
        const auto P = props.pressure();

        REQUIRE( uconstraints.size() == 0 );

        constraints.fix().pH(4.5);

        REQUIRE( uconstraints.size() == 1 );
        REQUIRE( uconstraints[0].formula.equivalent("H+") );
        REQUIRE( uconstraints[0].fn );

        constraints.fix().fugacity("O2(g)", 1.0, "bar");

        REQUIRE( uconstraints.size() == 2 );
        REQUIRE( uconstraints[1].formula.equivalent("O2") );
        REQUIRE( uconstraints[1].fn );

        REQUIRE( constraints.numChemicalPotentialConstraints() == 2 );
    }

    SECTION("Testing method EquilibriumConstraints::prevent")
    {
        // Reference to the reactivity constraints data
        const auto& get = constraints.data().restrictions;

        // Return the index of a species with given name
        const auto idx = [&](auto name) { return system.species().index(name); };

        // Return the pairs of species index and stoichiometries in a reaction equation
        const auto equation = [&](auto reaction) -> Pairs<Index, double>
        {
            auto pairs = parseReactionEquation(reaction);
            return vectorize(pairs, RKT_LAMBDA(x, std::make_pair(idx(x.first), x.second)));
        };

        constraints.prevent().fromIncreasing("NaCl(aq)");
        constraints.prevent().fromIncreasing("H2(g)");

        constraints.prevent().fromDecreasing("SiO2(s)");
        constraints.prevent().fromDecreasing("O2(g)");

        constraints.prevent().fromReacting("CaCO3(s)");
        constraints.prevent().fromReacting("CO2(g)");

        constraints.prevent().fromReacting("O2(aq) + H2(aq) = H2O(aq)");
        constraints.prevent().fromReacting("CO2(g) = CO2(aq)");

        REQUIRE( get.species_cannot_increase.size() == 2         );
        REQUIRE( get.species_cannot_increase.count(idx("NaCl(aq)"))  );
        REQUIRE( get.species_cannot_increase.count(idx("H2(g)")) );

        REQUIRE( get.species_cannot_decrease.size() == 2           );
        REQUIRE( get.species_cannot_decrease.count(idx("SiO2(s)")) );
        REQUIRE( get.species_cannot_decrease.count(idx("O2(g)"))   );

        REQUIRE( get.species_cannot_react.size() == 2            );
        REQUIRE( get.species_cannot_react.count(idx("CaCO3(s)")) );
        REQUIRE( get.species_cannot_react.count(idx("CO2(g)"))   );

        REQUIRE( get.reactions_cannot_react.size() == 2                                 );
        REQUIRE( get.reactions_cannot_react[0] == equation("O2(aq) + H2(aq) = H2O(aq)") );
        REQUIRE( get.reactions_cannot_react[1] == equation("CO2(g) = CO2(aq)")          );

        REQUIRE( constraints.numInertReactions() == 2 );
    }
}
