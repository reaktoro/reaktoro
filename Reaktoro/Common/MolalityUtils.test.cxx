// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2022 Allan Leal
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
#include <Reaktoro/Common/MolalityUtils.hpp>
using namespace Reaktoro;

TEST_CASE("Testing MolalityUtils module", "[MolalityUtils]")
{
    // The possibilities for an array of species amounts
    auto n           = ArrayXr{{1.0, 2.0, 3.0, 4.0}};
    auto nzero       = ArrayXr{{0.0, 0.0, 0.0, 0.0}};
    auto nsingle     = ArrayXr{{2.0}};
    auto nsinglezero = ArrayXr{{0.0}};

    // The index of water species within a vector of species amounts
    auto iH2O = 0;

    //-------------------------------------------------------------------------
    // TESTING METHOD: molalities
    //-------------------------------------------------------------------------
    ArrayXr m;

    iH2O = 0; m = molalities(n, iH2O);

    CHECK( m[0] == Approx(1.0 / 10.0) );
    CHECK( m[1] == Approx(2.0 / (1.0 * waterMolarMass)) );
    CHECK( m[2] == Approx(3.0 / (1.0 * waterMolarMass)) );
    CHECK( m[3] == Approx(4.0 / (1.0 * waterMolarMass)) );

    iH2O = 2; m = molalities(n, iH2O);

    CHECK( m[0] == Approx(1.0 / (3.0 * waterMolarMass)) );
    CHECK( m[1] == Approx(2.0 / (3.0 * waterMolarMass)) );
    CHECK( m[2] == Approx(3.0 / 10.0) );
    CHECK( m[3] == Approx(4.0 / (3.0 * waterMolarMass)) );

    iH2O = 3; m = molalities(nzero, iH2O);

    CHECK( m[0] == Approx(0.0) );
    CHECK( m[1] == Approx(0.0) );
    CHECK( m[2] == Approx(0.0) );
    CHECK( m[3] == Approx(1.0) );

    iH2O = 0; m = molalities(nsingle, iH2O);

    CHECK( m[0] == Approx(1.0) );

    iH2O = 0; m = molalities(nsinglezero, iH2O);

    CHECK( m[0] == Approx(1.0) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: molalitiesJacobian
    //-------------------------------------------------------------------------
    MatrixXd dmdn(4, 4);
    MatrixXd dmdn_expected(4, 4);

    iH2O = 2; dmdn = molalitiesJacobian(n, iH2O);

    const auto nsum = n.sum();
    const auto xH2O = n[iH2O]/nsum;

    dmdn_expected.fill(0.0);
    dmdn_expected.diagonal().array() = 1.0 / (n[iH2O]*waterMolarMass);
    dmdn_expected.col(iH2O) = -n/(n[iH2O]*n[iH2O]*waterMolarMass);
    dmdn_expected.row(iH2O).array() = -xH2O/nsum;
    dmdn_expected(iH2O, iH2O) += double(xH2O/n[iH2O]);

    INFO("dmdn = \n" << dmdn);
    INFO("dmdn(expected) = \n" << dmdn_expected);
    CHECK( dmdn.isApprox(dmdn_expected) );

    iH2O = 0; dmdn = molalitiesJacobian(nzero, iH2O);

    INFO("dmdn = \n" << dmdn);
    CHECK( dmdn.isApprox(MatrixXd::Zero(4, 4)) );

    iH2O = 0; dmdn = molalitiesJacobian(nsingle, iH2O);

    INFO("dmdn = \n" << dmdn);
    CHECK( dmdn.isApprox(MatrixXd::Zero(1, 1)) );

    iH2O = 0; dmdn = molalitiesJacobian(nsinglezero, iH2O);

    INFO("dmdn = \n" << dmdn);
    CHECK( dmdn.isApprox(MatrixXd::Zero(1, 1)) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: lnMolalitiesJacobian
    //-------------------------------------------------------------------------
    MatrixXd dlnmdn(4, 4);
    MatrixXd dlnmdn_expected(4, 4);

    iH2O = 1; dlnmdn = lnMolalitiesJacobian(n, iH2O);

    dlnmdn_expected.fill(0.0);
    dlnmdn_expected.diagonal().array() = 1.0 / n.array();
    dlnmdn_expected.col(iH2O).array() = -1.0/n[iH2O];
    dlnmdn_expected.row(iH2O).array() = -1.0/nsum;
    dlnmdn_expected(iH2O, iH2O) += double(1.0/n[iH2O]);

    INFO("dlnmdn = \n" << dlnmdn);
    INFO("dlnmdn(expected) = \n" << dlnmdn_expected);
    CHECK( dlnmdn.isApprox(dlnmdn_expected));

    iH2O = 0; dlnmdn = lnMolalitiesJacobian(nzero, iH2O);

    INFO("dlnmdn = \n" << dlnmdn);
    CHECK( dlnmdn.isApprox(MatrixXd::Zero(4, 4)) );

    iH2O = 0; dlnmdn = lnMolalitiesJacobian(nsingle, iH2O);

    INFO("dlnmdn = \n" << dlnmdn);
    CHECK( dlnmdn.isApprox(MatrixXd::Zero(1, 1)) );

    iH2O = 0; dlnmdn = lnMolalitiesJacobian(nsinglezero, iH2O);

    INFO("dlnmdn = \n" << dlnmdn);
    CHECK( dlnmdn.isApprox(MatrixXd::Zero(1, 1)) );

    //-------------------------------------------------------------------------
    // TESTING METHOD: lnMolalitiesJacobianDiagonal
    //-------------------------------------------------------------------------
    ArrayXd dlnmdn_diag(4);
    ArrayXd dlnmdn_diag_expected(4);

    iH2O = 1; dlnmdn_diag = lnMolalitiesJacobianDiagonal(n, iH2O);

    dlnmdn_diag_expected = 1/n;
    dlnmdn_diag_expected[iH2O] = 1.0/n[iH2O] - 1.0/nsum;

    INFO("dlnmdn_diag = " << dlnmdn_diag.transpose());
    INFO("dlnmdn_diag(expected) = " << dlnmdn_diag_expected.transpose());
    CHECK( dlnmdn_diag.isApprox(dlnmdn_diag_expected) );
}
