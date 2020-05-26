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
#include <Reaktoro/Equilibrium/EquilibriumConstraints.hpp>
#include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumProblem", "[EquilibriumProblem]")
{
    ChemicalSystem system = test::createChemicalSystem();

    EquilibriumConstraints constraints(system);

    constraints.control().temperature();
    constraints.control().pressure();

    REQUIRE_THROWS( EquilibriumProblem(constraints) ); // two functional constraints are needed for the two control variables

    constraints.until().volume(1.0, "m3");
    constraints.preserve().internalEnergy();

    constraints.fix().pH(3.0);
    constraints.fix().fugacity("O2(g)", 1.0, "bar");

    constraints.prevent().fromReacting("H2O(aq) = H+(aq) + OH-(aq)");
    constraints.prevent().fromReacting("CO2(g) = CO2(aq)");

    EquilibriumProblem problem(constraints);

    const auto N = system.species().size();
    const auto E = system.elements().size();

    REQUIRE( problem.numComponents() == E + 1 + 2 ); // elements, charge, inert reactions
    REQUIRE( problem.numVariables() == N + 4 ); // species, temperature, pressure, control H+ (from fixed pH), control O2 (from fixed fugacity)

    const auto Nx = problem.numVariables();

    //-------------------------------------------------------------------------
    // Testing EquilibriumProblem::conservationMatrix
    //-------------------------------------------------------------------------
    const auto C = problem.conservationMatrix();

    REQUIRE( C.rows() == problem.numComponents() );
    REQUIRE( C.cols() == problem.numVariables() );

    const auto Cx = MatrixXd // the expected conservation matrix!
    {{
        { 2,  1,  1,  2,  0,  0,  0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2,  2,  4,  0,  0,  0,  0,  0,  0,  1,  0 },
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  1,  0,  0,  0,  0,  0 },
        { 1,  0,  1,  0,  2,  0,  0,  0,  0,  1,  0,  0,  2,  3,  3,  0,  0,  2,  2,  2,  0,  1,  0,  1,  0,  3,  2,  0,  0,  0,  2 },
        { 0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0 },
        { 0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 },
        { 0,  1, -1,  0,  0,  1, -1,  0,  0,  0,  2,  2,  0, -1, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 },
        {-1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
        { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
    }};

    {
        INFO("C(assembled) = (transpose)\n" << C.transpose());
        INFO("C(expected) = (transpose)\n" << Cx.transpose());
        REQUIRE( C.isApprox(Cx) );
    }

    //-------------------------------------------------------------------------
    // Testing EquilibriumProblem::objective
    //-------------------------------------------------------------------------
    ChemicalProps props0(system);

    auto obj = problem.objective(props0);

    VectorXd x = VectorXd::Ones(Nx);

    VectorXd g(Nx);
    MatrixXd H(Nx, Nx);

    REQUIRE_NOTHROW( obj.f(x) );
    REQUIRE_NOTHROW( obj.g(x, g) );
    REQUIRE_NOTHROW( obj.H(x, H) );
}
