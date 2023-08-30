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
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Warnings.hpp>
#include <Reaktoro/Core/ChemicalState.hpp>
#include <Reaktoro/Core/ChemicalSystem.hpp>
#include <Reaktoro/Equilibrium/EquilibriumResult.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelCubicEOS.hpp>
#include <Reaktoro/Utils/AqueousProps.hpp>

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

    const auto iH2O = aqspecies.indexWithFormula("H2O");

    auto num_species = species.size();
    auto num_elements = elements.size();

    SECTION("Testing correct initialization of the `AqueousProps` instance")
    {
        CHECK( aqprops.phase().species().size()   == 19           );
        CHECK( aqprops.phase().elements().size()  == num_elements );
        CHECK( aqprops.phase().species(0).name()  == "H2O(aq)"    );
        CHECK( aqprops.phase().species(1).name()  == "H+(aq)"     );
        CHECK( aqprops.phase().species(2).name()  == "OH-(aq)"    );
        CHECK( aqprops.phase().species(3).name()  == "H2(aq)"     );
        CHECK( aqprops.phase().species(4).name()  == "O2(aq)"     );
        CHECK( aqprops.phase().species(5).name()  == "Na+(aq)"    );
        CHECK( aqprops.phase().species(6).name()  == "Cl-(aq)"    );
        CHECK( aqprops.phase().species(7).name()  == "NaCl(aq)"   );
        CHECK( aqprops.phase().species(8).name()  == "HCl(aq)"    );
        CHECK( aqprops.phase().species(9).name()  == "NaOH(aq)"   );
        CHECK( aqprops.phase().species(10).name() == "Ca++(aq)"   );
        CHECK( aqprops.phase().species(11).name() == "Mg++(aq)"   );
        CHECK( aqprops.phase().species(12).name() == "CO2(aq)"    );
        CHECK( aqprops.phase().species(13).name() == "HCO3-(aq)"  );
        CHECK( aqprops.phase().species(14).name() == "CO3--(aq)"  );
        CHECK( aqprops.phase().species(15).name() == "CaCl2(aq)"  );
        CHECK( aqprops.phase().species(16).name() == "MgCl2(aq)"  );
        CHECK( aqprops.phase().species(17).name() == "SiO2(aq)"   );
        CHECK( aqprops.phase().species(18).name() == "e-(aq)"     );

        const auto saturation_species = aqprops.saturationSpecies();

        CHECK( saturation_species.size() == 11 );
        CHECK( saturation_species[0].name()  == "CO2(g)"        );
        CHECK( saturation_species[1].name()  == "O2(g)"         );
        CHECK( saturation_species[2].name()  == "H2(g)"         );
        CHECK( saturation_species[3].name()  == "H2O(g)"        );
        CHECK( saturation_species[4].name()  == "CH4(g)"        );
        CHECK( saturation_species[5].name()  == "CO(g)"         );
        CHECK( saturation_species[6].name()  == "NaCl(s)"       );
        CHECK( saturation_species[7].name()  == "CaCO3(s)"      );
        CHECK( saturation_species[8].name()  == "MgCO3(s)"      );
        CHECK( saturation_species[9].name()  == "CaMg(CO3)2(s)" );
        CHECK( saturation_species[10].name() == "SiO2(s)"       );
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
        CHECK( aqprops.Eh()                          == Approx(0.000671794)  );
        CHECK( aqprops.pH()                          == Approx(-0.0205718)   );
        CHECK( aqprops.ionicStrength()               == Approx(499.576)      );
        CHECK( aqprops.ionicStrengthEffective()      == Approx(499.576)      );
        CHECK( aqprops.ionicStrengthStoichiometric() == Approx(999.152)      );

        // Check electric charge properties
        CHECK( aqprops.charge() == Approx(state.charge()) );
        CHECK( aqprops.chargeMolality() == Approx(state.charge() / state.speciesMass(iH2O)) );

        // Check molalities of all species
        for(auto i = 0; i < aqspecies.size(); i++)
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

        // Check saturation ratios of the non-aqueous species
        auto lnOmega = aqprops.saturationRatiosLn();

        CHECK( lnOmega[0]  == Approx(-1.474090) );
        CHECK( lnOmega[1]  == Approx(-1.447940) );
        CHECK( lnOmega[2]  == Approx(-1.447940) );
        CHECK( lnOmega[3]  == Approx(-1.421790) );
        CHECK( lnOmega[4]  == Approx(-1.421790) );
        CHECK( lnOmega[5]  == Approx(-1.500230) );
        CHECK( lnOmega[6]  == Approx(-9.050290) );
        CHECK( lnOmega[7]  == Approx(-9.050290) );
        CHECK( lnOmega[8]  == Approx(-9.050290) );
        CHECK( lnOmega[9]  == Approx( 0.102009) );
        CHECK( lnOmega[10] == Approx(-9.050290) );

        // Set activity model of gases to that of Peng-Robinson and for solids,
        // ideal model. Note: the reason the saturation indices below for
        // solids differ from those above is because the mock chemical system
        // does consider unit activities for the solid species!
        aqprops.setActivityModel("CO2(g)"       , ActivityModelPengRobinson());
        aqprops.setActivityModel("O2(g)"        , ActivityModelPengRobinson());
        aqprops.setActivityModel("H2(g)"        , ActivityModelPengRobinson());
        aqprops.setActivityModel("H2O(g)"       , ActivityModelPengRobinson());
        aqprops.setActivityModel("CH4(g)"       , ActivityModelPengRobinson());
        aqprops.setActivityModel("CO(g)"        , ActivityModelPengRobinson());
        aqprops.setActivityModel("NaCl(s)"      , ActivityModelIdealSolution(StateOfMatter::Solid));
        aqprops.setActivityModel("CaCO3(s)"     , ActivityModelIdealSolution(StateOfMatter::Solid));
        aqprops.setActivityModel("MgCO3(s)"     , ActivityModelIdealSolution(StateOfMatter::Solid));
        aqprops.setActivityModel("CaMg(CO3)2(s)", ActivityModelIdealSolution(StateOfMatter::Solid));
        aqprops.setActivityModel("SiO2(s)"      , ActivityModelIdealSolution(StateOfMatter::Solid));

        lnOmega = aqprops.saturationRatiosLn();

        CHECK( lnOmega[0]  == Approx(0.031383400) );
        CHECK( lnOmega[1]  == Approx(0.052984200) );
        CHECK( lnOmega[2]  == Approx(0.051783300) );
        CHECK( lnOmega[3]  == Approx(3.703000000) );
        CHECK( lnOmega[4]  == Approx(0.080426100) );
        CHECK( lnOmega[5]  == Approx(0.000339846) );
        CHECK( lnOmega[6]  == Approx(0.049714300) );
        CHECK( lnOmega[7]  == Approx(0.049714300) );
        CHECK( lnOmega[8]  == Approx(0.049714300) );
        CHECK( lnOmega[9]  == Approx(0.102009000) );
        CHECK( lnOmega[10] == Approx(0.049714300) );

        CHECK( aqprops.saturationRatio(5)       == Approx(exp(0.000339846)) );
        CHECK( aqprops.saturationRatio("CO(g)") == Approx(exp(0.000339846)) );

        CHECK( aqprops.saturationIndex(5)       == Approx(0.000339846/ln10) );
        CHECK( aqprops.saturationIndex("CO(g)") == Approx(0.000339846/ln10) );
    }

    SECTION("Testing when state is a brine")
    {
        ChemicalState state(system);
        state.setTemperature(T, "celsius");
        state.setPressure(P, "bar");
        state.set("H2O(aq)" , 1.00, "kg");
        state.set("Na+(aq)" , 2.05, "mg");
        state.set("Ca++(aq)", 1.42, "mg");
        state.set("Mg++(aq)", 0.39, "mg");
        state.set("Cl-(aq)" , 3.47, "mg");

        EquilibriumResult result = solver.solve(state);

        aqprops.update(state);

        // Check values of aqueous properties
        CHECK( aqprops.pressure().val()                    == Approx(P*1e5)         );
        CHECK( aqprops.temperature().val()                 == Approx(T + 273.15)    );
        CHECK( aqprops.ionicStrength().val()               == Approx(27.2721122969) );
        CHECK( aqprops.Eh().val()                          == Approx(0.000754321)   );
        CHECK( aqprops.pH().val()                          == Approx(-0.0605375)    );
        CHECK( aqprops.ionicStrengthEffective().val()      == Approx(27.2721122969) );
        CHECK( aqprops.ionicStrengthStoichiometric().val() == Approx(27.2725395818) );

        // Check water properties
        CHECK( aqprops.waterAmount().val() == Approx(state.speciesAmount(iH2O))  );
        CHECK( aqprops.waterMass().val()   == Approx(state.speciesMass(iH2O))    );

        // Check electric charge properties
        CHECK( aqprops.charge().val() == Approx(state.charge()) );
        CHECK( aqprops.chargeMolality().val() == Approx(state.charge() / state.speciesMass(iH2O)) );

        // Check molalities of all species
        CHECK( aqprops.speciesMolality("H2O(aq)"  ).val() == Approx(55.5085)       );
        CHECK( aqprops.speciesMolality("H+(aq)"   ).val() == Approx(27.2718671642) );
        CHECK( aqprops.speciesMolality("OH-(aq)"  ).val() == Approx(27.2718671642) );
        CHECK( aqprops.speciesMolality("H2(aq)"   ).val() == Approx(44.0212)       );
        CHECK( aqprops.speciesMolality("O2(aq)"   ).val() == Approx(22.0106)       );
        CHECK( aqprops.speciesMolality("Na+(aq)"  ).val() == Approx(2.28438e-16)   );
        CHECK( aqprops.speciesMolality("Cl-(aq)"  ).val() == Approx(2.28438e-16)   );
        CHECK( aqprops.speciesMolality("NaCl(aq)" ).val() == Approx(2.28438e-16)   );
        CHECK( aqprops.speciesMolality("HCl(aq)"  ).val() == Approx(0.000223584)   );
        CHECK( aqprops.speciesMolality("NaOH(aq)" ).val() == Approx(0.000203703)   );
        CHECK( aqprops.speciesMolality("Ca++(aq)" ).val() == Approx(8.09398e-05)   );
        CHECK( aqprops.speciesMolality("Mg++(aq)" ).val() == Approx(3.66570e-05)   );
        CHECK( aqprops.speciesMolality("CO2(aq)"  ).val() == Approx(2.28438e-16)   );
        CHECK( aqprops.speciesMolality("HCO3-(aq)").val() == Approx(2.28438e-16)   );
        CHECK( aqprops.speciesMolality("CO3--(aq)").val() == Approx(2.28438e-16)   );
        CHECK( aqprops.speciesMolality("CaCl2(aq)").val() == Approx(2.28438e-16)   );
        CHECK( aqprops.speciesMolality("MgCl2(aq)").val() == Approx(2.28438e-16)   );
        CHECK( aqprops.speciesMolality("SiO2(aq)" ).val() == Approx(2.28438e-16)   );
        CHECK( aqprops.speciesMolality("e-(aq)"   ).val() == Approx(2.28438e-16)   );

        // Check molalities of all elements
        CHECK( aqprops.elementMolality("H").val()  == Approx(253.604)     );
        CHECK( aqprops.elementMolality("C").val()  == Approx(6.85313e-16) );
        CHECK( aqprops.elementMolality("O").val()  == Approx(126.802)     );
        CHECK( aqprops.elementMolality("Na").val() == Approx(0.000203703) );
        CHECK( aqprops.elementMolality("Mg").val() == Approx(3.6657e-05)  );
        CHECK( aqprops.elementMolality("Si").val() == Approx(2.28438e-16) );
        CHECK( aqprops.elementMolality("Cl").val() == Approx(0.000223584) );
        CHECK( aqprops.elementMolality("Ca").val() == Approx(8.09398e-05) );

        // Test convenience methods species and element molalities
        for(const auto& s : aqspecies)
        {
            const auto name = s.name();
            const auto idx = aqspecies.index(name);
            CHECK( aqprops.speciesMolality(name).val() == Approx(aqprops.speciesMolalities()[idx]) );
            CHECK( aqprops.speciesMolality(idx).val()  == Approx(aqprops.speciesMolalities()[idx]) );
        }

        for(const auto& e : aqelements)
        {
            const auto symbol = e.symbol();
            const auto idx = aqelements.index(symbol);
            CHECK( aqprops.elementMolality(symbol).val() == Approx(aqprops.elementMolalities()[idx]) );
            CHECK( aqprops.elementMolality(idx).val()    == Approx(aqprops.elementMolalities()[idx]) );
        }
    }

    SECTION("Testing static method AqueousProps::compute")
    {
        ChemicalState state(system);
        state.temperature(300);
        state.setPressure(3e5);
        state.setSpeciesAmounts(1.0);

        ChemicalProps props(system);
        props.update(state); // stateid is 1

        aqprops = AqueousProps::compute(props);
        CHECK( aqprops.props().stateid() == 1 );
        CHECK( aqprops.props().speciesAmount(0) == 1.0 );

        state.setSpeciesAmount(0, 99.0);
        props.update(state); // stateid is now 2
        aqprops = AqueousProps::compute(props);
        CHECK( aqprops.props().stateid() == 2 );
        CHECK( aqprops.props().speciesAmount(0) == 99.0 );

        // Force internal change of amount for species at index 0, which will
        // not cause props.stateid() to change!
        const_cast<real*>(props.speciesAmounts().data())[0] = 123.0;

        aqprops = AqueousProps::compute(props);
        CHECK( aqprops.props().stateid() == 2 );
        CHECK( aqprops.props().speciesAmount(0) == 99.0 ); // this shows that memoization is working; the AqueousProps object has not been updated with the given ChemicalProps object

        ChemicalProps clone_props(props);

        aqprops = AqueousProps::compute(clone_props);
        CHECK( aqprops.props().stateid() == 2 ); // same stateid, because of clone
        CHECK( aqprops.props().speciesAmount(0) == 123.0 ); // this shows that memoization did not work with the cloned props because it is a different ChemicalProps object with different pointer
    }
}
