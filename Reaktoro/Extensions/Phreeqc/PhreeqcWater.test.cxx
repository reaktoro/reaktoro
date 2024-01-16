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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Extensions/Phreeqc/PhreeqcWater.hpp>
using namespace Reaktoro;

TEST_CASE("Testing PhreeqcWater module", "[PhreeqcWater]")
{
    const auto T = 298.15;
    const auto P = 1.0e5;

    const auto wtp = PhreeqcUtils::waterThermoProps(T, P);
    const auto wep = PhreeqcUtils::waterElectroProps(T, P, wtp);

    CHECK( wtp.rho_0   == Approx(0.997042)    );
    CHECK( wtp.kappa_0 == Approx(4.51729e-05) );

    CHECK( wep.eps_r == Approx(78.3844)       );
    CHECK( wep.DH_A  == Approx(0.510025)      );
    CHECK( wep.DH_B  == Approx(0.328491)      );
    CHECK( wep.DH_Av == Approx(1.88731)       );
    CHECK( wep.ZBrn  == Approx(41.3063)       );
    CHECK( wep.QBrn  == Approx(2.5234e-05)    );

    const auto wprops = PhreeqcUtils::waterProps(T, P);

    CHECK( wprops.wtp.rho_0   == wtp.rho_0   );
    CHECK( wprops.wtp.kappa_0 == wtp.kappa_0 );

    CHECK( wprops.wep.eps_r == wep.eps_r );
    CHECK( wprops.wep.DH_A  == wep.DH_A  );
    CHECK( wprops.wep.DH_B  == wep.DH_B  );
    CHECK( wprops.wep.DH_Av == wep.DH_Av );
    CHECK( wprops.wep.ZBrn  == wep.ZBrn  );
    CHECK( wprops.wep.QBrn  == wep.QBrn  );
}
