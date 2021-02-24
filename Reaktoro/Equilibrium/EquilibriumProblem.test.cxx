// // Reaktoro is a unified framework for modeling chemically reactive systems.
// //
// // Copyright (C) 2014-2020 Allan Leal
// //
// // This library is free software; you can redistribute it and/or
// // modify it under the terms of the GNU Lesser General Public
// // License as published by the Free Software Foundation; either
// // version 2.1 of the License, or (at your option) any later version.
// //
// // This library is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// // Lesser General Public License for more details.
// //
// // You should have received a copy of the GNU Lesser General Public License
// // along with this library. If not, see <http://www.gnu.org/licenses/>.

// // Catch includes
// #include <catch2/catch.hpp>

// // Reaktoro includes
// #include <Reaktoro/Core/ChemicalState.hpp>
// #include <Reaktoro/Core/ChemicalSystem.hpp>
// #include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
// #include <Reaktoro/Equilibrium/EquilibriumDims.hpp>
// #include <Reaktoro/Equilibrium/EquilibriumProblem.hpp>
// using namespace Reaktoro;

// namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

// TEST_CASE("Testing EquilibriumProblem", "[EquilibriumProblem]")
// {
//     ChemicalSystem system = test::createChemicalSystem();

//     EquilibriumConditions conditions(system);

//     // conditions.temperature();
//     // conditions.pressure();

//     // REQUIRE_THROWS( EquilibriumProblem(conditions) ); // two functional constraints are needed for the two control variables

//     conditions.volume(1.0, "m3");
//     conditions.internalEnergy(1.0, "kJ");

//     conditions.pH(3.0);
//     conditions.fugacity("O2(g)", 1.0, "bar");

//     conditions.cannotReact("H2O(aq) = H+(aq) + OH-(aq)");
//     conditions.cannotReact("CO2(g) = CO2(aq)");
//     conditions.cannotReact("CaCO3(s)");
//     conditions.cannotIncrease("NaCl(s)");
//     conditions.cannotDecrease("SiO2(s)");

//     EquilibriumProblem problem(conditions);

//     //-------------------------------------------------------------------------
//     // Testing EquilibriumProblem::conservationMatrix
//     //-------------------------------------------------------------------------
//     const auto dims = problem.dims();

//     const auto C = problem.conservationMatrix();

//     REQUIRE( C.rows() == dims.Nb ); // the number of components
//     REQUIRE( C.cols() == dims.Nx ); // the number of variables

//     const auto Cx = MatrixXd // the expected conservation matrix!
//     {{
//         { 2,  1,  1,  2,  0,  0,  0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2,  2,  4,  0,  0,  0,  0,  0,  0,  1,  0 },
//         { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  0,  0,  1,  0,  0,  0,  1,  1,  0,  1,  0,  0,  0,  0,  0 },
//         { 1,  0,  1,  0,  2,  0,  0,  0,  0,  1,  0,  0,  2,  3,  3,  0,  0,  2,  2,  2,  0,  1,  0,  1,  0,  3,  2,  0,  0,  0,  2 },
//         { 0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
//         { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
//         { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0 },
//         { 0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 },
//         { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0 },
//         { 0,  1, -1,  0,  0,  1, -1,  0,  0,  0,  2,  2,  0, -1, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0 },
//         {-1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
//         { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
//     }};

//     {
//         INFO("C(assembled) = (transpose)\n" << C.transpose());
//         INFO("C(expected) = (transpose)\n" << Cx.transpose());
//         REQUIRE( C.isApprox(Cx) );
//     }

//     //-------------------------------------------------------------------------
//     // Testing EquilibriumProblem::objective
//     //-------------------------------------------------------------------------
//     ChemicalState state0(system);
//     state0.setTemperature(298.15);
//     state0.setPressure(1.0e5);
//     state0.setSpeciesAmounts(1.0);

//     auto obj = problem.objective(state0);

//     VectorXd x = VectorXd::Ones(dims.Nx);

//     VectorXd g(dims.Nx);
//     MatrixXd H(dims.Nx, dims.Nx);

//     REQUIRE_NOTHROW( obj.f(x) );
//     REQUIRE_NOTHROW( obj.g(x, g) );
//     REQUIRE_NOTHROW( obj.H(x, H) );

//     //-------------------------------------------------------------------------
//     // Testing EquilibriumProblem::xlower
//     //-------------------------------------------------------------------------
//     ArrayXd xlower(dims.Nx);

//     problem.xlower(state0, xlower);

//     REQUIRE( xlower[system.species().index("CaCO3(s)")] == state0.speciesAmount("CaCO3(s)") );
//     REQUIRE( xlower[system.species().index("SiO2(s)")]  == state0.speciesAmount("SiO2(s)")  );

//     //-------------------------------------------------------------------------
//     // Testing EquilibriumProblem::xupper
//     //-------------------------------------------------------------------------
//     ArrayXd xupper(dims.Nx);

//     problem.xupper(state0, xupper);

//     REQUIRE( xupper[system.species().index("CaCO3(s)")] == state0.speciesAmount("CaCO3(s)") );
//     REQUIRE( xupper[system.species().index("NaCl(s)")]  == state0.speciesAmount("NaCl(s)")  );
// }
