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
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumConditions", "[EquilibriumConditions]")
{
    ChemicalSystem system = test::createChemicalSystem();

    // EquilibriumConditions conditions(system);

    // // The mock temperature, pressure, and species amounts used for the tests below.
    // const auto T = 1.0;
    // const auto P = 11.0;
    // const auto n = ArrayXr::Ones(system.species().size()).eval();

    // // The mock chemical properties of the system used for the tests below
    // ChemicalProps props(system);
    // props.update(T, P, n);

    // SECTION("Testing method EquilibriumConditions::temperature")
    // {
    //     const auto& details = conditions.details();

    //     REQUIRE( details.unknownT == true );
    //     REQUIRE( details.constantT == false );

    //     conditions.temperature(T);

    //     REQUIRE( details.unknownT == false );
    //     REQUIRE( details.constantT == false );

    //     conditions.constantTemperature();

    //     REQUIRE( details.unknownT == false );
    //     REQUIRE( details.constantT == true );
    // }

    // SECTION("Testing method EquilibriumConditions::pressure")
    // {
    //     const auto& details = conditions.details();

    //     REQUIRE( details.unknownP == true );
    //     REQUIRE( details.constantP == false );

    //     conditions.pressure(P);

    //     REQUIRE( details.unknownP == false );
    //     REQUIRE( details.constantP == false );

    //     conditions.constantPressure();

    //     REQUIRE( details.unknownP == false );
    //     REQUIRE( details.constantP == true );
    // }

    // SECTION("Testing method EquilibriumConditions::titrate")
    // {
    //     const auto& details = conditions.details();

    //     conditions.titrate("CO2");

    //     REQUIRE( details.titrants.back().str() == "CO2" );
    // }

    // SECTION("Testing method EquilibriumConditions::until")
    // {
    //     // Reference to the details data
    //     const auto& details = conditions.data().details;

    //     // Reference to the functional equilibrium conditions data
    //     const auto& econstraints = conditions.data().econstraints;

    //     const auto q = VectorXr{{ 1.0, 2.0 }};

    //     EquilibriumEquationArgs args{props, q, details.titrants};

    //     conditions.until().internalEnergy(1.0, "kJ");

    //     REQUIRE( econstraints.size() == 1 );
    //     REQUIRE( econstraints[0].fn(args) == props.internalEnergy() - 1.0e+3 );

    //     conditions.until().volume(1.0, "mm3");

    //     REQUIRE( econstraints.size() == 2 );
    //     REQUIRE( econstraints[0].fn(args) == props.internalEnergy() - 1.0e+3 );
    //     REQUIRE( econstraints[1].fn(args) == props.volume() - 1.0e-9 );

    //     conditions.until().enthalpy(1.0, "J");

    //     REQUIRE( econstraints.size() == 3 );
    //     REQUIRE( econstraints[0].fn(args) == props.internalEnergy() - 1.0e+3 );
    //     REQUIRE( econstraints[1].fn(args) == props.volume() - 1.0e-9 );
    //     REQUIRE( econstraints[2].fn(args) == props.enthalpy() - 1.0 );

    //     //---------------------------------------------------------------------
    //     // Impose an enthalpy constraint again with updated value and ensure
    //     // that the existing constraint is updated instead of creating a new one.
    //     //---------------------------------------------------------------------
    //     conditions.until().enthalpy(300.0, "J");

    //     // Ensure the number of conditions remain unchanged!
    //     REQUIRE( econstraints.size() == 3 );

    //     // Ensure the existing enthalpy constraint gets updated!
    //     REQUIRE( econstraints[2].fn(args) == props.enthalpy() - 300.0 );
    // }

    // SECTION("Testing method EquilibriumConditions::preserve")
    // {
    //     // Reference to the chemical property preservation conditions data
    //     const auto& pconstraints = conditions.data().pconstraints;

    //     conditions.preserve().internalEnergy();

    //     REQUIRE( pconstraints.size() == 1 );
    //     REQUIRE( pconstraints[0].fn(props) == props.internalEnergy() );

    //     conditions.preserve().volume();

    //     REQUIRE( pconstraints.size() == 2 );
    //     REQUIRE( pconstraints[0].fn(props) == props.internalEnergy() );
    //     REQUIRE( pconstraints[1].fn(props) == props.volume() );

    //     conditions.preserve().enthalpy();

    //     REQUIRE( pconstraints.size() == 3 );
    //     REQUIRE( pconstraints[0].fn(props) == props.internalEnergy() );
    //     REQUIRE( pconstraints[1].fn(props) == props.volume() );
    //     REQUIRE( pconstraints[2].fn(props) == props.enthalpy() );

    //     //---------------------------------------------------------------------
    //     // Impose the preservation constraint on volume again and ensure
    //     // that the existing constraint is updated instead of creating a new one.
    //     //---------------------------------------------------------------------
    //     conditions.preserve().volume();

    //     // Ensure the number of conditions remain unchanged!
    //     REQUIRE( pconstraints.size() == 3 );
    //     REQUIRE( pconstraints[0].fn(props) == props.internalEnergy() );
    //     REQUIRE( pconstraints[1].fn(props) == props.volume() );
    //     REQUIRE( pconstraints[2].fn(props) == props.enthalpy() );
    // }

    // SECTION("Testing method EquilibriumConditions::fix")
    // {
    //     // Reference to the chemical potential conditions data
    //     const auto& uconstraints = conditions.data().uconstraints;

    //     // The Species objects for H+(aq) and O2(g)
    //     const auto hplus = system.database().species().getWithName("H+(aq)");
    //     const auto o2g = system.database().species().getWithName("O2(g)");

    //     // The standard chemical potentials of H+(aq) and O2(g) at (T, P)
    //     auto u0hplus = hplus.props(T, P).G0;
    //     auto u0o2g = o2g.props(T, P).G0;

    //     // The universal gas constant (in J/(mol*K))
    //     const auto R = universalGasConstant;

    //     conditions.fix().pH(4.5);

    //     REQUIRE( uconstraints.size() == 1 );
    //     REQUIRE( uconstraints.back().formula.equivalent("H+") );
    //     REQUIRE( uconstraints.back().fn(T, P) == Approx(u0hplus + R*T*log(pow(10.0, -4.5))) );

    //     conditions.fix().fugacity("O2(g)", 10.0, "bar");

    //     REQUIRE( uconstraints.size() == 2 );
    //     REQUIRE( uconstraints.back().formula.equivalent("O2") );
    //     REQUIRE( uconstraints.back().fn(T, P) == Approx(u0o2g + R*T*log(10.0)) );

    //     //---------------------------------------------------------------------
    //     // Impose pH with an updated value. Ensure the existing chemical
    //     // potential constraint for H+ is used instead of creating a new one.
    //     //---------------------------------------------------------------------

    //     conditions.fix().pH(8.0);

    //     // Ensure the number of conditions remain unchanged!
    //     REQUIRE( uconstraints.size() == 2 );

    //     // Ensure the existing chemical potential constraint for H+ has been updated!
    //     REQUIRE( uconstraints[0].fn(T, P) == Approx(u0hplus + R*T*log(pow(10.0, -8.0))) );

    //     //---------------------------------------------------------------------
    //     // Impose activity of H+ and ensure the existing chemical potential
    //     // constraint for H+ is reused instead of creating a new one.
    //     //---------------------------------------------------------------------
    //     conditions.fix().activity("H+(aq)", 1e-6);

    //     // Ensure the number of conditions remain unchanged!
    //     REQUIRE( uconstraints.size() == 2 );

    //     // Ensure the existing chemical potential constraint for H+ has been updated!
    //     REQUIRE( uconstraints[0].fn(T, P) == Approx(u0hplus + R*T*log(1e-6)) );
    // }

    // SECTION("Testing method EquilibriumConditions::prevent")
    // {
    //     // Reference to the reactivity conditions data
    //     const auto& get = conditions.data().restrictions;

    //     // Return the index of a species with given name
    //     const auto idx = [&](auto name) { return system.species().index(name); };

    //     // Return the pairs of species index and stoichiometries in a reaction equation
    //     const auto equation = [&](auto reaction) -> Pairs<Index, double>
    //     {
    //         auto pairs = parseReactionEquation(reaction);
    //         return vectorize(pairs, RKT_LAMBDA(x, std::make_pair(idx(x.first), x.second)));
    //     };

    //     conditions.prevent().fromIncreasing("NaCl(aq)");
    //     conditions.prevent().fromIncreasing("H2(g)");

    //     conditions.prevent().fromDecreasing("SiO2(s)");
    //     conditions.prevent().fromDecreasing("O2(g)");

    //     conditions.prevent().fromReacting("CaCO3(s)");
    //     conditions.prevent().fromReacting("CO2(g)");

    //     conditions.prevent().fromReacting("O2(aq) + H2(aq) = H2O(aq)");
    //     conditions.prevent().fromReacting("CO2(g) = CO2(aq)");

    //     REQUIRE( get.species_cannot_increase.size() == 2         );
    //     REQUIRE( get.species_cannot_increase.count(idx("NaCl(aq)"))  );
    //     REQUIRE( get.species_cannot_increase.count(idx("H2(g)")) );

    //     REQUIRE( get.species_cannot_decrease.size() == 2           );
    //     REQUIRE( get.species_cannot_decrease.count(idx("SiO2(s)")) );
    //     REQUIRE( get.species_cannot_decrease.count(idx("O2(g)"))   );

    //     REQUIRE( get.species_cannot_react.size() == 2            );
    //     REQUIRE( get.species_cannot_react.count(idx("CaCO3(s)")) );
    //     REQUIRE( get.species_cannot_react.count(idx("CO2(g)"))   );

    //     REQUIRE( get.reactions_cannot_react.size() == 2                                 );
    //     REQUIRE( get.reactions_cannot_react[0] == equation("O2(aq) + H2(aq) = H2O(aq)") );
    //     REQUIRE( get.reactions_cannot_react[1] == equation("CO2(g) = CO2(aq)")          );
    // }

    // SECTION("Testing when the EquilibriumConditions object is locked.")
    // {
    //     conditions.control().temperature();

    //     conditions.until().volume(1.0, "m3");

    //     conditions.fix().pH(5.0);
    //     conditions.fix().fugacity("H2(g)", 5.0, "bar");

    //     conditions.prevent().fromReacting("H2O(aq) = H+(aq) + OH-(aq)");

    //     conditions.lock();

    //     REQUIRE_NOTHROW( conditions.until().volume(10.0, "m3")           );
    //     REQUIRE_NOTHROW( conditions.fix().pH(2.0)                        );
    //     REQUIRE_NOTHROW( conditions.fix().fugacity("H2(g)", 10.0, "bar") );
    //     REQUIRE_NOTHROW( conditions.prevent().fromReacting("CaCO3(s)")   );
    //     REQUIRE_NOTHROW( conditions.prevent().fromIncreasing("SiO2(s)")  );
    //     REQUIRE_NOTHROW( conditions.prevent().fromDecreasing("NaCl(s)")  );

    //     REQUIRE_THROWS( conditions.control()                                  );
    //     REQUIRE_THROWS( conditions.preserve()                                 );
    //     REQUIRE_THROWS( conditions.until().internalEnergy(1.0, "J")           );
    //     REQUIRE_THROWS( conditions.until().enthalpy(1.0, "J")                 );
    //     REQUIRE_THROWS( conditions.prevent().fromReacting("CO2(g) = CO2(aq)") );
    // }
}
