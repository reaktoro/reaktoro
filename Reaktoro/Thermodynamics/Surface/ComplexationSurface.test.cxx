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
#include <Reaktoro/Extensions/Phreeqc/PhreeqcDatabase.hpp>
#include <Reaktoro/Thermodynamics/Surface/ComplexationSurface.hpp>
#include <Reaktoro/Common/Constants.hpp>

using namespace Reaktoro;

namespace test {

extern auto createDatabasePhases() -> Database;

auto getPhreeqcDatabase(const String& name) -> PhreeqcDatabase;

} // namespace test

TEST_CASE("Testing ComplexationSurface", "[ComplexationSurface]")
{
    // Initialize the database corresponding to the string `phreeqc.dat` has been already initialized
    auto dbphreeqc = test::getPhreeqcDatabase("phreeqc.dat");

    // Create complexation surface
    ComplexationSurface surface_Hfo("Hfo");
    surface_Hfo.setSpecificSurfaceArea(60, "m2/g")
        .setMass(4.45, "g");

    // Defined strong site of the complexation surface
    surface_Hfo.addSite("Hfo_w", "_w")
        .setAmount(1e-3, "mol");

    // Defined weak site of the complexation surface
    ComplexationSurfaceSite site_Hfo_s;
    site_Hfo_s.setName("Hfo_s")
        .setAmount(0.025e-3, "mol");
    surface_Hfo.addSite(site_Hfo_s);

    // Fetch all the species with Adsorbed aggregate state
    SpeciesList adsorbed_species = dbphreeqc.species().withAggregateState(AggregateState::Adsorbed);

    // Defined and add surface species
    String selected_absorbed_species = "Hfo_sOH Hfo_sOHCa+2 Hfo_sOH2+ Hfo_sO- Hfo_sOHSr+2 "
                                       "Hfo_wOH Hfo_wOH2+ Hfo_wO- Hfo_wOCa+ Hfo_wOSr+ Hfo_wOSrOH";
    // Select the names of considered absorbed species
    SpeciesList species_list = adsorbed_species.withNames(selected_absorbed_species);

    surface_Hfo.addSurfaceSpecies(species_list);

    SECTION("Checking the species in ComplexationSurface")
    {
        CHECK(surface_Hfo.species()[0].name()  == "Hfo_sOH"     ); // Hfo_sOH
        CHECK(surface_Hfo.species()[1].name()  == "Hfo_sOHCa+2" ); // Hfo_sOHCa+2
        CHECK(surface_Hfo.species()[2].name()  == "Hfo_sOH2+"   ); // Hfo_sOH2+
        CHECK(surface_Hfo.species()[3].name()  == "Hfo_sO-"     ); // Hfo_sO-
        CHECK(surface_Hfo.species()[4].name()  == "Hfo_sOHSr+2" ); // Hfo_sOHSr+2
        CHECK(surface_Hfo.species()[5].name()  == "Hfo_wOH"     ); // Hfo_wOH
        CHECK(surface_Hfo.species()[6].name()  == "Hfo_wOH2+"   ); // Hfo_wOH2+
        CHECK(surface_Hfo.species()[7].name()  == "Hfo_wO-"     ); // Hfo_wO-
        CHECK(surface_Hfo.species()[8].name()  == "Hfo_wOCa+"   ); // Hfo_wOCa+
        CHECK(surface_Hfo.species()[9].name()  == "Hfo_wOSr+"   ); // Hfo_wOSr+
        CHECK(surface_Hfo.species()[10].name() == "Hfo_wOSrOH"  ); // Hfo_wOSrOH
    }

    SECTION("Checking the charges of the species in ComplexationSurface")
    {
        // The numbers of exchanger's equivalents for exchange species
        ArrayXd z = surface_Hfo.charges();

        CHECK( z[0]  ==  0.0 ); // Hfo_sOH
        CHECK( z[1]  ==  2.0 ); // Hfo_sOHCa+2
        CHECK( z[2]  ==  1.0 ); // Hfo_sOH2+
        CHECK( z[3]  == -1.0 ); // Hfo_sO-
        CHECK( z[4]  ==  2.0 ); // Hfo_sOHSr+2
        CHECK( z[5]  ==  0.0 ); // Hfo_wOH
        CHECK( z[6]  ==  1.0 ); // Hfo_wOH2+
        CHECK( z[7]  == -1.0 ); // Hfo_wO-
        CHECK( z[8]  ==  1.0 ); // Hfo_wOCa+
        CHECK( z[9]  ==  1.0 ); // Hfo_wOSr+
        CHECK( z[10] ==  0.0 ); // Hfo_wOSrOH
    }

    // Create surface complexation state
    ComplexationSurfaceState state;
    state.T = 273.15 + 25;

    // Data for this test is taken from Appelo etal, 2005, Example 6.12, p. 292
    state.sigma = faradayConstant * 1.07e-6;
    real I = 0.1;
    state.updatePotential(I);

    // y = 2.78084626624
    // arcsinhy = 1.74676718758
    // psi = 0.08975790202

    SECTION("Checking the potential of ComplexationSurface for surface charge 1.07e-6 eq/m2 "
            "and an ionic strength 0.1 molal")
    {
        CHECK(state.psi == Approx(0.08975790202) );
    }

    I = 0.025;
    state.updatePotential(I);

    // y = 5.56169253249
    // arcsinhy = 2.41703553095
    // psi = 0.08975790202

    SECTION("Checking the potential of ComplexationSurface for surface charge 1.07e-6 eq/m2 "
            "and an ionic strength 0.025 molal")
    {
        CHECK(state.psi == Approx(0.12419974449) );
    }

    // Data for this test is taken from Appelo etal, 2005, Example 7.6, p. 337
    state.sigma = faradayConstant * 2.6e-7;
    I = 0.1;
    state.updatePotential(I);

    // y = 0.67571965348
    // arcsinhy = 0.63266191582
    // psi = 0.03250943037

    SECTION("Checking the potential of ComplexationSurface for surface charge 2.6e-7 eq/m2 "
            "and an ionic strength 0.1 molal")
    {
        CHECK(state.psi == Approx(0.03250943037) );
    }
}
