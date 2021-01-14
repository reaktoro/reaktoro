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
#include <Reaktoro/Equilibrium/EquilibriumSpecs.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumSpecs", "[EquilibriumSpecs]")
{
    // The mock chemical system for the tests below.
    ChemicalSystem system = test::createChemicalSystem();

    // // The mock temperature, pressure, and species amounts used for the tests below.
    // const auto T = 1.0;
    // const auto P = 11.0;
    // const auto n = ArrayXr::Ones(system.species().size()).eval();

    // // The mock chemical properties of the system used for the tests below.
    // ChemicalProps props(system);
    // props.update(T, P, n);

    // // The EquilibriumSpecs object tested below.
    // EquilibriumSpecs specs(system);

    // REQUIRE( specs.details().unknownT == true );
    // specs.temperature();
    // REQUIRE( specs.details().unknownT == false );

    // REQUIRE( specs.details().unknownP == true );
    // specs.pressure();
    // REQUIRE( specs.details().unknownP == false );

    // specs.volume();
    // specs.internalEnergy();
    // specs.enthalpy();
    // specs.gibbsEnergy();
    // specs.helmholtzEnergy();
    // specs.entropy();

    // REQUIRE( specs.details().econstraints.size() == 6 );
    // REQUIRE( specs.details().econstraints[0].id == "volume" );
    // REQUIRE( specs.details().econstraints[1].id == "internalEnergy" );
    // REQUIRE( specs.details().econstraints[2].id == "enthalpy" );
    // REQUIRE( specs.details().econstraints[3].id == "gibbsEnergy" );
    // REQUIRE( specs.details().econstraints[4].id == "helmholtzEnergy" );
    // REQUIRE( specs.details().econstraints[5].id == "entropy" );

    // Params params;
    // params.set("V", 1.0);
    // params.set("U", 2.0);
    // params.set("H", 3.0);
    // params.set("G", 4.0);
    // params.set("A", 5.0);
    // params.set("S", 6.0);

    // REQUIRE( specs.details().econstraints[0].fn(props, params) == props.volume() - 1.0 );
    // REQUIRE( specs.details().econstraints[1].fn(props, params) == props.internalEnergy() - 2.0 );
    // REQUIRE( specs.details().econstraints[2].fn(props, params) == props.enthalpy() - 3.0 );
    // REQUIRE( specs.details().econstraints[3].fn(props, params) == props.gibbsEnergy() - 4.0 );
    // REQUIRE( specs.details().econstraints[4].fn(props, params) == props.helmholtzEnergy() - 5.0 );
    // REQUIRE( specs.details().econstraints[5].fn(props, params) == props.entropy() - 6.0 );


    // SECTION("Testing method EquilibriumSpecs::temperature")
    // {
    //     const auto& details = specs.details();

    //     REQUIRE( details.unknownT == true );
    //     REQUIRE( details.constantT == false );

    //     specs.temperature(T);

    //     REQUIRE( details.unknownT == false );
    //     REQUIRE( details.constantT == false );

    //     specs.constantTemperature();

    //     REQUIRE( details.unknownT == false );
    //     REQUIRE( details.constantT == true );
    // }

    // SECTION("Testing method EquilibriumSpecs::pressure")
    // {
    //     const auto& details = specs.details();

    //     REQUIRE( details.unknownP == true );
    //     REQUIRE( details.constantP == false );

    //     specs.pressure(P);

    //     REQUIRE( details.unknownP == false );
    //     REQUIRE( details.constantP == false );

    //     specs.constantPressure();

    //     REQUIRE( details.unknownP == false );
    //     REQUIRE( details.constantP == true );
    // }

    // SECTION("Testing method EquilibriumSpecs::titrate")
    // {
    //     const auto& details = specs.details();

    //     specs.titrate("CO2");

    //     REQUIRE( details.titrants.back().str() == "CO2" );
    // }

    // SECTION("Testing method EquilibriumSpecs::until")
    // {
    //     // Reference to the details data
    //     const auto& details = specs.data().details;

    //     // Reference to the functional equilibrium specs data
    //     const auto& econstraints = specs.data().econstraints;

    //     const auto q = VectorXr{{ 1.0, 2.0 }};

    //     EquilibriumEquationArgs args{props, q, details.titrants};

    //     specs.until().internalEnergy(1.0, "kJ");

    //     REQUIRE( econstraints.size() == 1 );
    //     REQUIRE( econstraints[0].fn(args) == props.internalEnergy() - 1.0e+3 );

    //     specs.until().volume(1.0, "mm3");

    //     REQUIRE( econstraints.size() == 2 );
    //     REQUIRE( econstraints[0].fn(args) == props.internalEnergy() - 1.0e+3 );
    //     REQUIRE( econstraints[1].fn(args) == props.volume() - 1.0e-9 );

    //     specs.until().enthalpy(1.0, "J");

    //     REQUIRE( econstraints.size() == 3 );
    //     REQUIRE( econstraints[0].fn(args) == props.internalEnergy() - 1.0e+3 );
    //     REQUIRE( econstraints[1].fn(args) == props.volume() - 1.0e-9 );
    //     REQUIRE( econstraints[2].fn(args) == props.enthalpy() - 1.0 );

    //     //---------------------------------------------------------------------
    //     // Impose an enthalpy constraint again with updated value and ensure
    //     // that the existing constraint is updated instead of creating a new one.
    //     //---------------------------------------------------------------------
    //     specs.until().enthalpy(300.0, "J");

    //     // Ensure the number of specs remain unchanged!
    //     REQUIRE( econstraints.size() == 3 );

    //     // Ensure the existing enthalpy constraint gets updated!
    //     REQUIRE( econstraints[2].fn(args) == props.enthalpy() - 300.0 );
    // }

    // SECTION("Testing method EquilibriumSpecs::preserve")
    // {
    //     // Reference to the chemical property preservation specs data
    //     const auto& pconstraints = specs.data().pconstraints;

    //     specs.preserve().internalEnergy();

    //     REQUIRE( pconstraints.size() == 1 );
    //     REQUIRE( pconstraints[0].fn(props) == props.internalEnergy() );

    //     specs.preserve().volume();

    //     REQUIRE( pconstraints.size() == 2 );
    //     REQUIRE( pconstraints[0].fn(props) == props.internalEnergy() );
    //     REQUIRE( pconstraints[1].fn(props) == props.volume() );

    //     specs.preserve().enthalpy();

    //     REQUIRE( pconstraints.size() == 3 );
    //     REQUIRE( pconstraints[0].fn(props) == props.internalEnergy() );
    //     REQUIRE( pconstraints[1].fn(props) == props.volume() );
    //     REQUIRE( pconstraints[2].fn(props) == props.enthalpy() );

    //     //---------------------------------------------------------------------
    //     // Impose the preservation constraint on volume again and ensure
    //     // that the existing constraint is updated instead of creating a new one.
    //     //---------------------------------------------------------------------
    //     specs.preserve().volume();

    //     // Ensure the number of specs remain unchanged!
    //     REQUIRE( pconstraints.size() == 3 );
    //     REQUIRE( pconstraints[0].fn(props) == props.internalEnergy() );
    //     REQUIRE( pconstraints[1].fn(props) == props.volume() );
    //     REQUIRE( pconstraints[2].fn(props) == props.enthalpy() );
    // }

    // SECTION("Testing method EquilibriumSpecs::fix")
    // {
    //     // Reference to the chemical potential specs data
    //     const auto& uconstraints = specs.data().uconstraints;

    //     // The Species objects for H+(aq) and O2(g)
    //     const auto hplus = system.database().species().getWithName("H+(aq)");
    //     const auto o2g = system.database().species().getWithName("O2(g)");

    //     // The standard chemical potentials of H+(aq) and O2(g) at (T, P)
    //     auto u0hplus = hplus.props(T, P).G0;
    //     auto u0o2g = o2g.props(T, P).G0;

    //     // The universal gas constant (in J/(mol*K))
    //     const auto R = universalGasConstant;

    //     specs.fix().pH(4.5);

    //     REQUIRE( uconstraints.size() == 1 );
    //     REQUIRE( uconstraints.back().formula.equivalent("H+") );
    //     REQUIRE( uconstraints.back().fn(T, P) == Approx(u0hplus + R*T*log(pow(10.0, -4.5))) );

    //     specs.fix().fugacity("O2(g)", 10.0, "bar");

    //     REQUIRE( uconstraints.size() == 2 );
    //     REQUIRE( uconstraints.back().formula.equivalent("O2") );
    //     REQUIRE( uconstraints.back().fn(T, P) == Approx(u0o2g + R*T*log(10.0)) );

    //     //---------------------------------------------------------------------
    //     // Impose pH with an updated value. Ensure the existing chemical
    //     // potential constraint for H+ is used instead of creating a new one.
    //     //---------------------------------------------------------------------

    //     specs.fix().pH(8.0);

    //     // Ensure the number of specs remain unchanged!
    //     REQUIRE( uconstraints.size() == 2 );

    //     // Ensure the existing chemical potential constraint for H+ has been updated!
    //     REQUIRE( uconstraints[0].fn(T, P) == Approx(u0hplus + R*T*log(pow(10.0, -8.0))) );

    //     //---------------------------------------------------------------------
    //     // Impose activity of H+ and ensure the existing chemical potential
    //     // constraint for H+ is reused instead of creating a new one.
    //     //---------------------------------------------------------------------
    //     specs.fix().activity("H+(aq)", 1e-6);

    //     // Ensure the number of specs remain unchanged!
    //     REQUIRE( uconstraints.size() == 2 );

    //     // Ensure the existing chemical potential constraint for H+ has been updated!
    //     REQUIRE( uconstraints[0].fn(T, P) == Approx(u0hplus + R*T*log(1e-6)) );
    // }

    // SECTION("Testing method EquilibriumSpecs::prevent")
    // {
    //     // Reference to the reactivity specs data
    //     const auto& get = specs.data().restrictions;

    //     // Return the index of a species with given name
    //     const auto idx = [&](auto name) { return system.species().index(name); };

    //     // Return the pairs of species index and stoichiometries in a reaction equation
    //     const auto equation = [&](auto reaction) -> Pairs<Index, double>
    //     {
    //         auto pairs = parseReactionEquation(reaction);
    //         return vectorize(pairs, RKT_LAMBDA(x, std::make_pair(idx(x.first), x.second)));
    //     };

    //     specs.prevent().fromIncreasing("NaCl(aq)");
    //     specs.prevent().fromIncreasing("H2(g)");

    //     specs.prevent().fromDecreasing("SiO2(s)");
    //     specs.prevent().fromDecreasing("O2(g)");

    //     specs.prevent().fromReacting("CaCO3(s)");
    //     specs.prevent().fromReacting("CO2(g)");

    //     specs.prevent().fromReacting("O2(aq) + H2(aq) = H2O(aq)");
    //     specs.prevent().fromReacting("CO2(g) = CO2(aq)");

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

    // SECTION("Testing when the EquilibriumSpecs object is locked.")
    // {
    //     specs.control().temperature();

    //     specs.until().volume(1.0, "m3");

    //     specs.fix().pH(5.0);
    //     specs.fix().fugacity("H2(g)", 5.0, "bar");

    //     specs.prevent().fromReacting("H2O(aq) = H+(aq) + OH-(aq)");

    //     specs.lock();

    //     REQUIRE_NOTHROW( specs.until().volume(10.0, "m3")           );
    //     REQUIRE_NOTHROW( specs.fix().pH(2.0)                        );
    //     REQUIRE_NOTHROW( specs.fix().fugacity("H2(g)", 10.0, "bar") );
    //     REQUIRE_NOTHROW( specs.prevent().fromReacting("CaCO3(s)")   );
    //     REQUIRE_NOTHROW( specs.prevent().fromIncreasing("SiO2(s)")  );
    //     REQUIRE_NOTHROW( specs.prevent().fromDecreasing("NaCl(s)")  );

    //     REQUIRE_THROWS( specs.control()                                  );
    //     REQUIRE_THROWS( specs.preserve()                                 );
    //     REQUIRE_THROWS( specs.until().internalEnergy(1.0, "J")           );
    //     REQUIRE_THROWS( specs.until().enthalpy(1.0, "J")                 );
    //     REQUIRE_THROWS( specs.prevent().fromReacting("CO2(g) = CO2(aq)") );
    // }
}
