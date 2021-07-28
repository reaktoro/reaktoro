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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumJacobian.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumJacobian", "[EquilibriumJacobian]")
{
    ChemicalSystem system = test::createChemicalSystem();

    ChemicalState state(system);
    state.temperature(3.0);
    state.pressure(5.0);
    for(auto i = 0; i < Nn; ++i)
        state.setSpeciesAmount(i, 0.1 * i);
    state.set("H2O(aq)", 1.0, "kg");

    const auto T = state.temperature();
    const auto P = state.pressure();
    const auto n = state.speciesAmounts();

    EquilibriumJacobian jac(system);

    auto dudn_exact_expected_fn = [&](VectorXr n)
    {
        ChemicalProps props(system);
        auto u = [&](VectorXrConstRef n)
        {
            props.update(T, P, n);
            return props.chemicalPotentials();
        };
        return jacobian(u, wrt(n), at(n));
    };

    auto dudn_approx_expected_fn = [&](VectorXr n)
    {
        ChemicalProps props(system);
        auto u = [&](VectorXrConstRef n)
        {
            props.updateIdeal(T, P, n);
            return props.chemicalPotentials();
        };
        return jacobian(u, wrt(n), at(n));
    };

    MatrixXd dudn_exact = jac.dudnExact(T, P, n);
    MatrixXd dudn_exact_expected = dudn_exact_expected_fn(n);

    INFO("dudn_exact = \n" << dudn_exact);
    INFO("dudn_exact(expected) = \n" << dudn_exact_expected);
    CHECK( dudn_exact.isApprox(dudn_exact_expected) );

    MatrixXd dudn_approx = jac.dudnApproximate(n);
    MatrixXd dudn_approx_expected = dudn_approx_expected_fn(n);

    INFO("dudn_approx = \n" << dudn_approx);
    INFO("dudn_approx(expected) = \n" << dudn_approx_expected);
    CHECK( dudn_approx.isApprox(dudn_approx_expected) );

    MatrixXd dudn_diag = jac.dudnDiagonal(n);
    MatrixXd dudn_diag_expected = dudn_approx_expected.diagonal().asDiagonal();

    INFO("dudn_diag = \n" << dudn_diag);
    INFO("dudn_diag(expected) = \n" << dudn_diag_expected);
    CHECK( dudn_diag.isApprox(dudn_diag_expected) );

    VectorXl idxs = VectorXl{{0, 2, 4, 7}};
    MatrixXd dudn_partially_exact = jac.dudnPartiallyExact(T, P, n, idxs);
    MatrixXd dudn_partially_exact_expected = dudn_approx_expected;
    dudn_partially_exact_expected(Eigen::all, idxs) = dudn_exact_expected(Eigen::all, idxs);

    INFO("dudn_partially_exact = \n" << dudn_partially_exact);
    INFO("dudn_partially_exact(expected) = \n" << dudn_partially_exact_expected);
    CHECK( dudn_partially_exact.isApprox(dudn_partially_exact_expected) );
}
