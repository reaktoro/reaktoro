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
#include <Reaktoro/Extensions/Geochemistry/AqueousMixture.hpp>
#include <Reaktoro/Extensions/Water/WaterConstants.hpp>
using namespace Reaktoro;

auto moleFractions(Index size) -> ArrayXr
{
    const auto n = ArrayXr::Random(size);
    return n / n.sum();
}

TEST_CASE("Testing AqueousMixture class", "[AqueousMixture]")
{
    WHEN("When the aqueous mixture is setup correctly")
    {
        SpeciesList species("H2O H+ OH- Na+ Cl- Ca++ Mg++ HCO3- CO3-- K+ CO2 HCl NaCl NaOH CaCl2 MgCl2 CaCO3 MgCO3");

        AqueousMixture mixture(species);

        REQUIRE( mixture.species().size() == species.size() );

        // Test H2O is not classified as neutral solute species
        REQUIRE_THROWS( mixture.neutral().index("H2O") );

        for(auto x : species)
        {
            if(x.name() == "H2O") continue;

            // Test AqueousMixture::neutral|charged|cations|anions methods
            if(x.charge() == 0.0) REQUIRE_NOTHROW( mixture.neutral().index(x.name()) );
            if(x.charge() != 0.0) REQUIRE_NOTHROW( mixture.charged().index(x.name()) );
            if(x.charge()  > 0.0) REQUIRE_NOTHROW( mixture.cations().index(x.name()) );
            if(x.charge()  < 0.0) REQUIRE_NOTHROW( mixture.anions().index(x.name()) );

            // Test AqueousMixture::indicesXYZ methods
            if(x.charge() == 0.0) REQUIRE( contains(mixture.indicesNeutral(), species.index(x.name())) );
            if(x.charge() != 0.0) REQUIRE( contains(mixture.indicesCharged(), species.index(x.name())) );
            if(x.charge()  > 0.0) REQUIRE( contains(mixture.indicesCations(), species.index(x.name())) );
            if(x.charge()  < 0.0) REQUIRE( contains(mixture.indicesAnions(), species.index(x.name())) );
        }

        // Test AqueousMixture::indexWater method
        REQUIRE( mixture.indexWater() == species.index("H2O") );

        // Test AqueousMixture::charges method
        REQUIRE( mixture.charges().size() == species.size() );
        for(auto i = 0; i < species.size(); ++i)
            REQUIRE( mixture.charges()[i] == species[i].charge() );

        // Test AqueousMixture::dissociationMatrix method
        REQUIRE( mixture.dissociationMatrix().rows() == mixture.neutral().size() );
        REQUIRE( mixture.dissociationMatrix().cols() == mixture.charged().size() );

        // Assemble the expected dissociation matrix for the constructed aqueous mixture
        const MatrixXd M
        { //  H+     OH-    Na+    Cl-    Ca++   Mg++   HCO3-  CO3--  K+
            {  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0 }, // CO2
            {  1.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0 }, // HCl
            {  0.0,   0.0,   1.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0 }, // NaCl
            {  0.0,   1.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0 }, // NaOH
            {  0.0,   0.0,   0.0,   2.0,   1.0,   0.0,   0.0,   0.0,   0.0 }, // CaCl2
            {  0.0,   0.0,   0.0,   2.0,   0.0,   1.0,   0.0,   0.0,   0.0 }, // MgCl2
            {  0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   1.0,   0.0 }, // CaCO3
            {  0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   1.0,   0.0 }, // MgCO3
        };

        INFO( "dissociation matrix is\n" << mixture.dissociationMatrix() << "\nbut expected is\n" << M );
        REQUIRE( mixture.dissociationMatrix().isApprox(M) );

        // The temperature (in K), pressure (in Pa) and mole fractions of the aqueous species
        const auto T = 345.67;
        const auto P = 123.4e+5;
        const auto x = moleFractions(species.size());

        const auto state = mixture.state(T, P, x);

        const ArrayXd z  = mixture.charges();                            // charges of the species
        const ArrayXd zc = mixture.charges()(mixture.indicesCharged());  // charges of charged species

        const auto m  = x/(waterMolarMass * x[mixture.indexWater()]); // molalities of spcies
        const auto mc = m(mixture.indicesCharged()).matrix();         // molalities of charged solutes
        const auto mn = m(mixture.indicesNeutral()).matrix();         // molalities of neutral solutes
        const auto ms = (mc + M.transpose() * mn).array();            // stoichiometric molalities (as array)
        const auto Ie = 0.5 * (m * z * z).sum();                      // effective ionic strength
        const auto Is = 0.5 * (ms * zc * zc).sum();                   // stoichiometric ionic strength

        REQUIRE( state.T       == T                      );
        REQUIRE( state.P       == P                      );
        REQUIRE( state.Ie      == Approx(Ie)             );
        REQUIRE( state.Is      == Approx(Is)             );
        REQUIRE( state.rho     == Approx(997.0470390177) );
        REQUIRE( state.epsilon == Approx(78.2451448082)  );

        REQUIRE( state.m.isApprox(m)   );
        REQUIRE( state.ms.isApprox(ms) );
    }
}
