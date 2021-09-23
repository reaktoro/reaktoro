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
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>
#include <Reaktoro/Thermodynamics/Surface/IonExchangeSurface.hpp>

using namespace Reaktoro;

TEST_CASE("Testing IonExchangeSurface", "[IonExchangeSurface]")
{
    // Load phreeqc database
    PhreeqcDatabase db("phreeqc.dat");

    // Define ion exchange species list
    // Expected species: X- AlOHX2 AlX3 BaX2 CaX2 CdX2 CuX2 FeX2 KX LiX MgX2 MnX2 NH4X NaX PbX2 SrX2 ZnX2
    SpeciesList species = db.species().withAggregateState(AggregateState::IonExchange);

    SECTION("Checking the charges of the species in IonExchangeSurface")
    {
        // Create the aqueous mixture
        IonExchangeSurface surface(species);

        // The numbers of exchanger's equivalents for exchange species
        ArrayXd ze = surface.ze();

        CHECK( ze[0]  == 0 ); // X-
        CHECK( ze[1]  == 2 ); // AlOHX2
        CHECK( ze[2]  == 3 ); // AlX3
        CHECK( ze[3]  == 2 ); // BaX2
        CHECK( ze[4]  == 2 ); // CaX2
        CHECK( ze[5]  == 2 ); // CdX2
        CHECK( ze[6]  == 2 ); // CuX2
        CHECK( ze[7]  == 2 ); // FeX2
        CHECK( ze[8]  == 1 ); // KX
        CHECK( ze[9]  == 1 ); // LiX
        CHECK( ze[10] == 2 ); // MgX2
        CHECK( ze[11] == 2 ); // MnX2
        CHECK( ze[12] == 1 ); // NH4X
        CHECK( ze[13] == 1 ); // NaX
        CHECK( ze[14] == 2 ); // PbX2
        CHECK( ze[15] == 2 ); // SrX2
        CHECK( ze[16] == 2 ); // ZnX2
    }
}
