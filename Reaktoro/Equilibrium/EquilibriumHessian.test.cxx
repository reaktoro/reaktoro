// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
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
#include <Reaktoro/Equilibrium/EquilibriumHessian.hpp>
using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing EquilibriumHessian", "[EquilibriumHessian]")
{
    ChemicalSystem system = test::createChemicalSystem();

    const auto Nn = system.species().size();

    ChemicalState state(system);
    state.temperature(3.0);
    state.pressure(5.0);
    for(auto i = 0; i < Nn; ++i)
        state.setSpeciesAmount(i, 0.1 * i, "mol");
    state.set("H2O(aq)", 1.0, "kg");

    const auto T = state.temperature();
    const auto P = state.pressure();
    const auto n = state.speciesAmounts();
    const auto RT = universalGasConstant * T;

    EquilibriumHessian hessian(system);

    auto dudn_exact_expected_fn = [&](VectorXr n) -> MatrixXd
    {
        ChemicalProps props(system);
        auto u = [&](VectorXrConstRef n) -> VectorXr
        {
            props.update(T, P, n);
            return props.speciesChemicalPotentials()/RT;
        };
        return jacobian(u, wrt(n), at(n));
    };

    auto dudn_approx_expected_fn = [&](VectorXr n) -> MatrixXd
    {
        ChemicalProps props(system);
        auto u = [&](VectorXrConstRef n) -> VectorXr
        {
            props.updateIdeal(T, P, n);
            return props.speciesChemicalPotentials()/RT;
        };
        return jacobian(u, wrt(n), at(n));
    };

    MatrixXd dudn_exact = hessian.exact(T, P, n);
    MatrixXd dudn_exact_expected = dudn_exact_expected_fn(n);

    MatrixXd dudn_approx = hessian.approximate(n);
    MatrixXd dudn_approx_expected = dudn_approx_expected_fn(n);

    MatrixXd dudn_diag = hessian.diagonal(n);
    MatrixXd dudn_diag_expected = dudn_approx_expected.diagonal().asDiagonal();

    VectorXl idxs = VectorXl{{0, 2, 4, 7}};
    MatrixXd dudn_partially_exact = hessian.partiallyExact(T, P, n, idxs);
    MatrixXd dudn_partially_exact_expected = dudn_approx_expected;
    dudn_partially_exact_expected(Eigen::all, idxs) = dudn_exact_expected(Eigen::all, idxs);

    SECTION("testing EquilibriumHessian::dudnExact")
    {
        INFO("dudn_exact = \n" << dudn_exact);
        INFO("dudn_exact(expected) = \n" << dudn_exact_expected);
        CHECK( dudn_exact.isApprox(dudn_exact_expected) );
    }

    SECTION("testing EquilibriumHessian::dudnApproximate")
    {

        INFO("dudn_approx = \n" << dudn_approx);
        INFO("dudn_approx(expected) = \n" << dudn_approx_expected);
        CHECK( dudn_approx.isApprox(dudn_approx_expected) );
    }

    SECTION("testing EquilibriumHessian::dudnDiagonal")
    {
        INFO("dudn_diag = \n" << dudn_diag);
        INFO("dudn_diag(expected) = \n" << dudn_diag_expected);
        CHECK( dudn_diag.isApprox(dudn_diag_expected) );
    }

    SECTION("testing EquilibriumHessian::dudnPartiallyExact")
    {
        INFO("dudn_partially_exact = \n" << dudn_partially_exact);
        INFO("dudn_partially_exact(expected) = \n" << dudn_partially_exact_expected);
        CHECK( dudn_partially_exact.isApprox(dudn_partially_exact_expected) );
    }
}
