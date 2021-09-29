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
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Thermodynamics/Aqueous/AqueousProps.hpp>

using namespace Reaktoro;

namespace test { extern auto createChemicalSystem() -> ChemicalSystem; }

TEST_CASE("Testing AqueousProps class", "[AqueousProps]")
{
    ChemicalSystem system = test::createChemicalSystem();

    EquilibriumSolver solver(system);

    AqueousProps aqprops(system);

    auto phase = aqprops.phase();

    const auto& species = system.species();
    const auto& elements = system.elements();
    const auto& aqspecies = phase.species();
    const auto& aqelements = phase.elements();

    auto num_species = species.size();
    auto num_elements = elements.size();

    SECTION("Testing correct initialization of the `AqueousProps` instance")
    {
        CHECK( phase.species().size()   == 19           );
        CHECK( phase.elements().size()  == num_elements );
        CHECK( phase.species(0).name()  == "H2O(aq)"    );
        CHECK( phase.species(1).name()  == "H+(aq)"     );
        CHECK( phase.species(2).name()  == "OH-(aq)"    );
        CHECK( phase.species(3).name()  == "H2(aq)"     );
        CHECK( phase.species(4).name()  == "O2(aq)"     );
        CHECK( phase.species(5).name()  == "Na+(aq)"    );
        CHECK( phase.species(6).name()  == "Cl-(aq)"    );
        CHECK( phase.species(7).name()  == "NaCl(aq)"   );
        CHECK( phase.species(8).name()  == "HCl(aq)"    );
        CHECK( phase.species(9).name()  == "NaOH(aq)"   );
        CHECK( phase.species(10).name() == "Ca++(aq)"   );
        CHECK( phase.species(11).name() == "Mg++(aq)"   );
        CHECK( phase.species(12).name() == "CO2(aq)"    );
        CHECK( phase.species(13).name() == "HCO3-(aq)"  );
        CHECK( phase.species(14).name() == "CO3--(aq)"  );
        CHECK( phase.species(15).name() == "CaCl2(aq)"  );
        CHECK( phase.species(16).name() == "MgCl2(aq)"  );
        CHECK( phase.species(17).name() == "SiO2(aq)"   );
        CHECK( phase.species(18).name() == "e-(aq)"     );
    }

    const real T = 25.0; // celsius
    const real P = 1.0;  // bar

    SECTION("Testing when species have zero amounts")
    {
        const ArrayXr n = ArrayXr::Zero(num_species);

        ChemicalState state(system);
        state.setTemperature(T);
        state.setPressure(P);
        state.setSpeciesAmounts(n);

        CHECK_THROWS(aqprops.update(state));
    }

    SECTION("Testing when species have nonzero amounts")
    {
        ArrayXd n = ArrayXd::Ones(num_species);

        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesAmounts(n);

        aqprops.update(state);

        // Check values of aqueous properties
        CHECK( aqprops.pressure()                    == Approx(P*1e5)        );
        CHECK( aqprops.temperature()                 == Approx(T + 273.15)   );
        CHECK( aqprops.Eh()                          == Approx(-0.000671794) );
        CHECK( aqprops.pH()                          == Approx(-0.0205718)   );
        CHECK( aqprops.ionicStrength()               == Approx(499.576)      );
        CHECK( aqprops.ionicStrengthEffective()      == Approx(499.576)      );
        CHECK( aqprops.ionicStrengthStoichiometric() == Approx(999.152)      );

        // Check molalities of all species
        for (Index i = 0; i < aqspecies.size() ; i++)
            CHECK( aqprops.speciesMolalities()[i] == Approx(55.5085) );

        // Check molalities of all elements
        CHECK( aqprops.elementMolality("H")  == Approx(499.576) );
        CHECK( aqprops.elementMolality("C")  == Approx(166.525) );
        CHECK( aqprops.elementMolality("O")  == Approx(832.627) );
        CHECK( aqprops.elementMolality("Na") == Approx(166.525) );
        CHECK( aqprops.elementMolality("Mg") == Approx(111.017) );
        CHECK( aqprops.elementMolality("Si") == Approx(55.5085) );
        CHECK( aqprops.elementMolality("Cl") == Approx(388.559) );
        CHECK( aqprops.elementMolality("Ca") == Approx(111.017) );
    }

    SECTION("Testing when state is a brine")
    {
        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.setSpeciesMass("H2O(aq)" , 1.00, "kg");
        state.setSpeciesMass("Na+(aq)" , 2.05, "mg");
        state.setSpeciesMass("Ca++(aq)", 1.42, "mg");
        state.setSpeciesMass("Mg++(aq)", 0.39, "mg");
        state.setSpeciesMass("Cl-(aq)" , 3.47, "mg");

        EquilibriumResult result = solver.solve(state);

        aqprops.update(state);

        // Check values of aqueous properties
        CHECK( aqprops.pressure()                    == Approx(P*1e5)        );
        CHECK( aqprops.temperature()                 == Approx(T + 273.15)   );
        CHECK( aqprops.ionicStrength()               == Approx(27.2725)      );
        CHECK( aqprops.Eh()                          == Approx(-0.000754246) );
        CHECK( aqprops.pH()                          == Approx(-0.0605375)   );
        CHECK( aqprops.ionicStrengthEffective()      == Approx(27.2725)      );
        CHECK( aqprops.ionicStrengthStoichiometric() == Approx(27.273)       );

        // Check molalities of all species
        CHECK( aqprops.speciesMolality("H2O(aq)"  ) == Approx(55.5085)     );
        CHECK( aqprops.speciesMolality("H+(aq)"   ) == Approx(27.2723)     );
        CHECK( aqprops.speciesMolality("OH-(aq)"  ) == Approx(27.2723)     );
        CHECK( aqprops.speciesMolality("H2(aq)"   ) == Approx(44.0212)     );
        CHECK( aqprops.speciesMolality("O2(aq)"   ) == Approx(22.0106)     );
        CHECK( aqprops.speciesMolality("Na+(aq)"  ) == Approx(2.28438e-16) );
        CHECK( aqprops.speciesMolality("Cl-(aq)"  ) == Approx(2.28438e-16) );
        CHECK( aqprops.speciesMolality("NaCl(aq)" ) == Approx(2.28438e-16) );
        CHECK( aqprops.speciesMolality("HCl(aq)"  ) == Approx(0.000223584) );
        CHECK( aqprops.speciesMolality("NaOH(aq)" ) == Approx(0.000203703) );
        CHECK( aqprops.speciesMolality("Ca++(aq)" ) == Approx(8.09398e-05) );
        CHECK( aqprops.speciesMolality("Mg++(aq)" ) == Approx(3.6657e-05)  );
        CHECK( aqprops.speciesMolality("CO2(aq)"  ) == Approx(2.28438e-16) );
        CHECK( aqprops.speciesMolality("HCO3-(aq)") == Approx(2.28438e-16) );
        CHECK( aqprops.speciesMolality("CO3--(aq)") == Approx(2.28438e-16) );
        CHECK( aqprops.speciesMolality("CaCl2(aq)") == Approx(2.28438e-16) );
        CHECK( aqprops.speciesMolality("MgCl2(aq)") == Approx(2.28438e-16) );
        CHECK( aqprops.speciesMolality("SiO2(aq)" ) == Approx(2.28438e-16) );
        CHECK( aqprops.speciesMolality("e-(aq)"   ) == Approx(2.28438e-16) );

        // Check molalities of all elements
        CHECK( aqprops.elementMolality("H")  == Approx(253.604)     );
        CHECK( aqprops.elementMolality("C")  == Approx(6.85313e-16) );
        CHECK( aqprops.elementMolality("O")  == Approx(126.802)     );
        CHECK( aqprops.elementMolality("Na") == Approx(0.000203703) );
        CHECK( aqprops.elementMolality("Mg") == Approx(3.6657e-05)  );
        CHECK( aqprops.elementMolality("Si") == Approx(2.28438e-16) );
        CHECK( aqprops.elementMolality("Cl") == Approx(0.000223584) );
        CHECK( aqprops.elementMolality("Ca") == Approx(8.09398e-05) );

        // Test convenience methods species and element molalities
        for(const auto& s : aqspecies)
        {
            const auto name = s.name();
            const auto idx = aqspecies.index(name);
            CHECK( aqprops.speciesMolality(name)     == Approx(aqprops.speciesMolalities()[idx]) );
        }
        for(const auto& e : aqelements)
        {
            const auto symbol = e.symbol();
            const auto idx = aqelements.index(symbol);
            CHECK( aqprops.elementMolality(symbol)     == Approx(aqprops.elementMolalities()[idx]) );
        }
    }
}
